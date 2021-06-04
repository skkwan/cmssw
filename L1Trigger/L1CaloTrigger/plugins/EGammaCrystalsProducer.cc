/* 
 * Description: Read crystal-level information
 */

// system include files
#include <ap_int.h>
#include <array>
#include <cmath>
// #include <cstdint>
#include <iostream>
#include <memory>

// user include files
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

// ECAL TPs
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// HCAL TPs
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"

// Output tower collection
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloTower.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"

#include "L1Trigger/L1CaloTrigger/interface/ParametricCalibration.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

static constexpr int n_towers_Eta = 34;
static constexpr int n_towers_Phi = 72;
static constexpr int n_towers_halfPhi = 36;
static constexpr int n_towers_cardEta = 17;   // new: equivalent to n_towers_per_link
static constexpr int n_towers_cardPhi = 4;    
static constexpr int n_crystals_cardEta = (n_towers_Eta * n_towers_cardEta);
static constexpr int n_crystals_cardPhi = (n_towers_Phi * n_towers_cardPhi);

// outputs
static constexpr int n_links_card = 4;
static constexpr int n_clusters_link = 3;
static constexpr int n_towers_per_link = 17;
static constexpr int n_clusters_4link = 4 * 3;

static constexpr int CRYSTALS_IN_TOWER_ETA = 5;
static constexpr int CRYSTALS_IN_TOWER_PHI = 5;

static constexpr int TOWER_IN_ETA = 3;      // number of towers in eta, in one 3x4 region (barrel)
static constexpr int TOWER_IN_PHI = 4;      // number of towers in phi, in one 3x4 region (barrel)

//static constexpr int TOWER_IN_ETA_OVERLAP = 2; // number of towers in eta, in one 2x4 region (overlap)
//static constexpr int TOWER_IN_PHI_OVERLAP = 4; // number of towers in phi, in one 2x4 region (overlap)

static constexpr int CRYSTAL_IN_ETA = 15;   // number of crystals in eta, in one 3x4 region (barrel)
static constexpr int CRYSTAL_IN_PHI = 20;   // number of crystals in phi, in one 3x4 region (barrel)

static constexpr float ECAL_eta_range = 1.4841;
static constexpr float cut_500_MeV = 0.5;

static constexpr int N_CLUSTERS_PER_REGION = 4;       // number of clusters per ECAL region
static constexpr int N_REGIONS_PER_CARD = 6;          // number of ECAL regions per card


// Assert that the card index is within bounds. (Valid cc: 0 to 35, since there are 36 RCT cards)
bool isValidCard(int cc) {
  return ((cc > -1) && (cc < 36));
}

// 

// Get crystal's iEta from real eta. (identical to getCrystal_etaID in L1EGammaCrystalsEmulatorProducer.cc)
// This "global" iEta ranges from 0 to (33*5) since there are 34 towers in eta in the full detector, 
// each with five crystals in eta.
int getCrystal_iEta(float eta) {
  float size_cell = 2 * ECAL_eta_range / (CRYSTALS_IN_TOWER_ETA * n_towers_Eta);
  int iEta = int((eta + ECAL_eta_range) / size_cell);
  return iEta;
}

// Get crystal's iPhi from real phi. (identical to getCrystal_phiID in L1EGammaCrystalsEmulatorProducer.cc)
// This "global" iPhi ranges from 0 to (71*5) since there are 72 towers in phi in the full detector,
// each with five crystals in eta.
int getCrystal_iPhi(float phi) {
  float size_cell = 2 * M_PI / (CRYSTALS_IN_TOWER_PHI * n_towers_Phi);
  int iPhi = int((phi + M_PI) / size_cell);
  return iPhi;
}

// For a card (ranging from 0 to 35, since there are 36 cards), return the iEta of the crystal with max iEta.
// This represents the card boundaries in eta (identical to getEtaMax_card in the original emulator)
int getCard_iEtaMax(int cc) {
  assert(isValidCard(cc));
  
  int etamax = 0;
  if (cc % 2 == 0)   // Even card: negative eta
    etamax = (n_towers_cardEta * CRYSTALS_IN_TOWER_ETA - 1);  // First eta half. 5 crystals in eta in 1 tower.
  else               // Odd card: positive eta
    etamax = (n_towers_Eta * CRYSTALS_IN_TOWER_ETA - 1);
  return etamax;
}

// Same as above but for minimum iEta.
int getCard_iEtaMin(int cc) {
  int etamin = 0;
  if (cc % 2 == 0)  // Even card: negative eta 
    etamin = (0);  
  else                // Odd card: positive eta
    etamin = (n_towers_cardEta * CRYSTALS_IN_TOWER_ETA);
  return etamin;
}

// Same as above but for maximum iPhi. 
int getCard_iPhiMax(int cc) {
  int phimax = ((cc / 2) + 1) * 4 * CRYSTALS_IN_TOWER_PHI - 1;
  return phimax;
}

// Same as above but for minimum iPhi.
int getCard_iPhiMin(int cc) {
  int phimin = (cc / 2) * 4 * CRYSTALS_IN_TOWER_PHI;
  return phimin;
}

// For a real eta, get the tower absolute Eta index (possible values are 0-33, since there
// are 34 towers in eta. (Adapted from getTower_absoluteEtaID)
int getTower_absEtaID(float eta) {
  float size_cell = 2 * ECAL_eta_range / n_towers_Eta;
  int etaID = int((eta + ECAL_eta_range) / size_cell);
  return etaID;
}

// Same as above, but for phi.
// Possible values range from 0-71 (Adapted from getTower_absolutePhiID)
int getTower_absPhiID(float phi) {
  float size_cell = 2 * M_PI / n_towers_Phi;
  int phiID = int((phi + M_PI) / size_cell);
  return phiID;
}

// Given the RCT card number (0-35), get the TOWER iPhi of the "bottom left" corner (0-71)
int getCard_refTower_iPhi(int cc) {
  assert(isValidCard(cc));
      
  if ((cc % 2) == 1) { // if cc is odd
    return (int(cc / 2) * 4);
  }
  else {  // if cc is even
    return ((int(cc / 2) * 4) + 3);
  }
  
}
 
// Given the RCT card number (0-35), get the TOWER iEta of the "bottom left" corner (0-33)
 int getCard_refTower_iEta(int cc) {
   assert(isValidCard(cc));
   
   if ((cc % 2) == 1) { // if cc is odd
     return 17;
   }
   else {  // if cc is even 
     return 16;
   }
   
 }

// Given the RCT card number (0-35), get the crystal iPhi of the "bottom left" corner (0- 71*5)
int getCard_ref_iPhi(int cc) {
  assert(isValidCard(cc));

  int ref_iPhi;
  if ((cc % 2) == 1) { // if cc is odd                   
    ref_iPhi = int(cc / 2) * TOWER_IN_PHI * CRYSTALS_IN_TOWER_PHI;
  }
  else {  // if cc is even, the bottom left corner is further in iPhi, hence the +4                
    ref_iPhi = ((int(cc / 2) * TOWER_IN_PHI) + 4) * CRYSTALS_IN_TOWER_PHI;
  }

  return ref_iPhi;
}

// Given the RCT card number (0-35), get the crystal iEta of the "bottom left" corner
int getCard_ref_iEta(int cc) {

  assert(isValidCard(cc));
  if ((cc % 2) == 1) {  // if cc is odd (positive eta)
    return (17 * CRYSTALS_IN_TOWER_ETA);
  }
  else {   // if cc is even (negative eta)
    return ((16 * CRYSTALS_IN_TOWER_ETA) + 4);
  }
}

// For a crystal with real (eta, phi) and falling in card cc, get its local iEta 
// relative to the bottom left corner of the card (possible local iEta ranges from 0 to 17 * 5,
// since in one card, there are 17 towers in eta, each with 5 crystals in eta.
int getCrystal_local_iEta(float hitEta, int cc) {
  assert(isValidCard(cc));

  // Functionally the same thing as an absolute value: 
  int diff = (getCard_ref_iEta(cc) - getCrystal_iEta(hitEta));

  return abs(diff);
}

// Same as above, but for iPhi (possible local iPhi ranges from 0 to (3*5), since in one card,
// there are 4 towers in phi, each with 5 crystals in phi.
int getCrystal_local_iPhi(float hitPhi, int cc) {
  assert(isValidCard(cc));

  int diff = (getCard_ref_iPhi(cc) - getCrystal_iPhi(hitPhi));
  return abs(diff);
}



// Get the RCT card region that a crystal is in, given the "local" iEta of the crystal 
// 0 is region closest to eta = 0. Regions 0, 1, 2, 3, 4 are in the barrel, Region 5 is in overlap
int getRegionNumber(int local_iEta) {
  int no = int(local_iEta / (TOWER_IN_ETA * CRYSTALS_IN_TOWER_ETA));
  assert(no < 6); 
  return no;
}

// Helper functions from the emulator
int getEtaMin_card_emulator(int card) {
  int etamin = 0;
  if (card % 2 == 0)
    etamin = 0 * CRYSTALS_IN_TOWER_ETA;  // First eta half. 5 crystals in eta in 1 tower.
  else
    etamin = n_towers_cardEta * CRYSTALS_IN_TOWER_ETA;
  return etamin;
}

int getPhiMin_card_emulator(int card) {
  int phimin = (card / 2) * 4 * CRYSTALS_IN_TOWER_PHI;
  return phimin;
}

// From the etaID (0- 33*5) and phiID (0 - 71*5) of a crystal, get the crystal ID (0-24) within the
// tower it's in, where 0-24 is "unwrapping" the 5x5 array, using rows in eta.
int getCrystalIDInTower_emulator(int etaID, int phiID) {
  return int(CRYSTALS_IN_TOWER_PHI * (phiID % CRYSTALS_IN_TOWER_PHI) + (etaID % CRYSTALS_IN_TOWER_ETA));
}

//////////////////////////////////////////////////////////////////////////  

// Declare the EGammaCrystalsProducer class and its methods

class EGammaCrystalsProducer : public edm::stream::EDProducer<> {
public:
  explicit EGammaCrystalsProducer(const edm::ParameterSet&);
  ~EGammaCrystalsProducer() override;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<EcalEBTrigPrimDigiCollection> ecalTPEBToken_;
  edm::EDGetTokenT<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hcalTPToken_;
  edm::ESHandle<CaloTPGTranscoder> decoder_;

  l1tp2::ParametricCalibration calib_;

  edm::ESHandle<CaloGeometry> caloGeometry_;
  const CaloSubdetectorGeometry* ebGeometry;
  const CaloSubdetectorGeometry* hbGeometry;
  edm::ESHandle<HcalTopology> hbTopology;
  const HcalTopology* hcTopology_;
};

//////////////////////////////////////////////////////////////////////////

/*
 * Declare the SimpleCaloHit class (taken from previous emulator)
 */

class SimpleCaloHit {
private:
  float pt_ = 0;
  float energy_ = 0.;
  ap_uint<10> et_uint_;
  GlobalVector position_;  // As opposed to GlobalPoint, so we can add them (for weighted average)
  HcalDetId id_hcal_;
  EBDetId id_;
  
public:
  // tool functions
  inline void setPt() { pt_ = (position_.mag2() > 0) ? energy_ * sin(position_.theta()) : 0; };
  inline void setEnergy(float et) { energy_ = et / sin(position_.theta()); };
  inline void setEt_uint(ap_uint<10> et_uint) { et_uint_ = et_uint; }
  inline void setPosition(const GlobalVector& pos) { position_ = pos; };
  inline void setIdHcal(const HcalDetId& idhcal) { id_hcal_ = idhcal; };
  inline void setId(const EBDetId& id) { id_ = id; };

  inline float pt() const { return pt_; };
  inline float energy() const { return energy_; };
  inline ap_uint<10> et_uint() const { return et_uint_; };
  inline const GlobalVector& position() const { return position_; };
  inline const EBDetId& id() const { return id_; };
};

/*******************************************************************/

/*
 * linkECAL class: represents one ECAL link (one tower: 5x5 crystals)
 */

class linkECAL {
private:
  ap_uint<10> crystalE[CRYSTALS_IN_TOWER_ETA][CRYSTALS_IN_TOWER_PHI];

public:
  // constructor                                                                                               
  linkECAL() { crystalE[CRYSTALS_IN_TOWER_ETA][CRYSTALS_IN_TOWER_PHI] = {};  }

  // Set members
  inline void zeroOut() {  // zero out the crystalE array
    for (int i = 0; i < CRYSTALS_IN_TOWER_ETA; i++) {
      for (int j = 0; j < CRYSTALS_IN_TOWER_PHI; j++) { 
	crystalE[i][j] = 0;
      }}};
  inline void setCrystalE(int iEta, int iPhi, ap_uint<10> energy) {
    assert(iEta < CRYSTALS_IN_TOWER_ETA); assert(iPhi < CRYSTALS_IN_TOWER_PHI); crystalE[iEta][iPhi] = energy; };
  inline void addCrystalE(int iEta, int iPhi, ap_uint<10> energy) { 
    assert(iEta < CRYSTALS_IN_TOWER_ETA); assert(iPhi < CRYSTALS_IN_TOWER_PHI); crystalE[iEta][iPhi] += energy; };

  // Access members
  inline ap_uint<10> getCrystalE(int iEta, int iPhi) { assert(iEta < 5); assert(iPhi < 5); return crystalE[iEta][iPhi]; };

};

/*******************************************************************/

/*
 * region3x4 class: represents one 3x4 ECAL region. The region stores no
 *                  information about which card it is located in.
 *                  idx: 0-4. Region 0 is the one closest to eta = 0, counting outwards in eta  
 */

class region3x4 {
private:
  int idx_ = -1; 
  linkECAL linksECAL[TOWER_IN_ETA][TOWER_IN_PHI];  // 3x4 in towers

public:
  // constructor                                                                            
  region3x4() { idx_ = -1; }

  // copy constructor
  region3x4(const region3x4& other) {
    idx_ = other.idx_;
    for (int i = 0; i < TOWER_IN_ETA; i++) {
      for (int j = 0; j < TOWER_IN_PHI; j++ ) {
        linksECAL[i][j] = other.linksECAL[i][j];
      }
    }
  }
  
  // overload operator= to use copy constructor                                                            
  region3x4 operator=(const region3x4& other) {
    region3x4 newRegion(other);
    return newRegion;
  };

  // set members
  inline void zeroOut() { 
    for (int i = 0; i < TOWER_IN_ETA; i++) { 
      for (int j = 0; j < TOWER_IN_PHI; j++) { 
	linksECAL[i][j].zeroOut();
      }}
  };
  inline void setIdx(int idx) { idx_ = idx; };

  // get members
  inline float getIdx() const { return idx_; };
  inline linkECAL& getLinkECAL (int iEta, int iPhi) { return linksECAL[iEta][iPhi]; };
};

/*******************************************************************/

/*
 * towerHCAL class: represents one HCAL tower
 */

class towerHCAL {
private:
  ap_uint<10> et;
  ap_uint<6> fb;

public:
  // constructor
  towerHCAL() { et = 0;  fb = 0; };

  // copy constructor
  towerHCAL(const towerHCAL &other) {
    et = other.et; fb = other.fb;   };

  // set members
  inline void zeroOut() { et = 0; fb = 0; };
  inline void addEt(ap_uint<6> newEt) { et += newEt; };

  // get members
  inline ap_uint<10> getEt() { return et; };

};

/*******************************************************************/

/*
 * towers3x4 class: represents 3x4 array of HCAL towers. idx = 0, 1, ... 4 are the barrel gion
 */

class towers3x4 {
private:
  int idx_ = -1; 
  towerHCAL towersHCAL[TOWER_IN_ETA][TOWER_IN_PHI];  // 3x4 in towers

public:
  // constructor                                                                            
  towers3x4() { idx_ = -1; };

  // copy constructor
  towers3x4(const towers3x4& other) {
    idx_ = other.idx_;
    for (int i = 0; i < TOWER_IN_ETA; i++) {
      for (int j = 0; j < TOWER_IN_PHI; j++ ) {
	towersHCAL[i][j] = other.towersHCAL[i][j]; }}; };

  // set members
  inline void zeroOut() { 
    for (int i = 0; i < TOWER_IN_ETA; i++) { 
      for (int j = 0; j < TOWER_IN_PHI; j++) { 
	towersHCAL[i][j].zeroOut();      }}   };
  inline void setIdx(int idx) { idx_ = idx; };

  // get members
  inline float getIdx() const { return idx_; };
  inline towerHCAL& getTowerHCAL ( int iEta, int iPhi ){ return towersHCAL[iEta][iPhi]; }; 
  
};

/*******************************************************************/

/* 
 * card class: represents one RCT card. Each card has five 3x4 regions and one 2x4 region,
 *             which is represented by a 3x4 region with its third row zero'd out.
 *             idx 0-35: odd values of cardIdx span eta = 0 to eta = 1.41 
 *                       even values of cardIdx span eta = -1.41 to eta = 0
 *             The realEta and realPhi arrays store the (eta, phi) of the center of the towers.
 */

class card {
private:
  int idx_ = -1 ; 
  region3x4 card3x4Regions[N_REGIONS_PER_CARD];
  towers3x4 card3x4Towers[N_REGIONS_PER_CARD];
  
public:
  // constructor
  card() {  
    idx_ = -1; 
    for (int i = 0; i < N_REGIONS_PER_CARD; i++) {
      card3x4Regions[i].setIdx(i);  
      card3x4Regions[i].zeroOut();
      card3x4Towers[i].setIdx(i);
      card3x4Towers[i].zeroOut();
    }};
  
  // copy constructor
  card(const card& other) {
    idx_ = other.idx_;
    for (int i = 0; i < N_REGIONS_PER_CARD; i++) {
      card3x4Regions[i] = other.card3x4Regions[i]; 
      card3x4Towers[i]  = other.card3x4Towers[i];  }};

  // overload operator= to use copy constructor
  card operator=(const card& other) {
    card newCard(other);
    return newCard;  };
  
  // set members
  inline void setIdx(int idx) { idx_ = idx; };
  inline void zeroOut() { 
    for (int i = 0; i < N_REGIONS_PER_CARD; i++) { card3x4Regions[i].zeroOut(); card3x4Towers[i].zeroOut();    }; };

  // get members
  inline float getIdx() const { return idx_; };
  inline region3x4& getRegion3x4(int idx) { assert(idx < N_REGIONS_PER_CARD); return card3x4Regions[idx]; }
  inline towers3x4& getTowers3x4(int idx) { assert(idx < N_REGIONS_PER_CARD); return card3x4Towers[idx]; }

};

/*******************************************************************/

/*
 *  crystal class: 
 */

class crystal{
public:
  ap_uint<10> energy;   // formerly ap_uint<10>
  //  uint8_t  timing; // formerly ap_uint<4>

  crystal(){
    energy = 0;
    //    timing = 0;
  }

  crystal(ap_uint<10> energy){  // To-do: add timing information
    this->energy = energy;
    //    this->timing = 0; 
  }

  crystal& operator=(const crystal& rhs){
    energy = rhs.energy;
    //    timing = rhs.timing;
    return *this;
  }
};

/*
 * crystalMax class: adapted from algo_top.h
 */ 
class crystalMax{
public:
  ap_uint<10> energy;
  uint8_t phiMax;
  uint8_t etaMax;

  crystalMax(){
    energy = 0;
    phiMax = 0;
    etaMax = 0;
  }

  crystalMax& operator=(const crystalMax& rhs){
    energy = rhs.energy;
    phiMax = rhs.phiMax;
    etaMax = rhs.etaMax;
    return *this;
  }
};
/*******************************************************************/

// Copied from algo_top.h, with the exception of changing out the ap uints to ap_uint<10>

typedef struct {
  ap_uint<10> energy;
  ap_uint<5> phi;
  ap_uint<5> eta;
} ecaltp_t ;

typedef struct {
  ecaltp_t cr0;
  ecaltp_t cr1 ;
  ecaltp_t cr2 ;
  ecaltp_t cr3 ;
  ecaltp_t cr4 ;
  ecaltp_t cr5 ;
  ecaltp_t cr6 ;
  ecaltp_t cr7 ;
  ecaltp_t cr8 ;
  ecaltp_t cr9 ;
  ecaltp_t cr10 ;
  ecaltp_t cr11 ;
  ecaltp_t cr12 ;
  ecaltp_t cr13 ;
  ecaltp_t cr14 ;
  ecaltp_t cr15 ;
  ecaltp_t cr16 ;
  ecaltp_t cr17 ;
  ecaltp_t cr18 ;
  ecaltp_t cr19 ;
} etaStrip_t ;

typedef struct {
  etaStrip_t etaStrip0 ;
  etaStrip_t etaStrip1 ;
  etaStrip_t etaStrip2 ;
  etaStrip_t etaStrip3 ;
  etaStrip_t etaStrip4 ;
  etaStrip_t etaStrip5 ;
  etaStrip_t etaStrip6 ;
  etaStrip_t etaStrip7 ;
  etaStrip_t etaStrip8 ;
  etaStrip_t etaStrip9 ;
  etaStrip_t etaStrip10 ;
  etaStrip_t etaStrip11 ;
  etaStrip_t etaStrip12 ;
  etaStrip_t etaStrip13 ;
  etaStrip_t etaStrip14 ;
} ecalRegion_t ;


typedef struct {
  ecaltp_t pk0;
  ecaltp_t pk1;
  ecaltp_t pk2;
  ecaltp_t pk3;
  ecaltp_t pk4;
  ecaltp_t pk5;
  ecaltp_t pk6;
  ecaltp_t pk7;
  ecaltp_t pk8;
  ecaltp_t pk9;
  ecaltp_t pk10;
  ecaltp_t pk11;
  ecaltp_t pk12;
  ecaltp_t pk13;
  ecaltp_t pk14;
} etaStripPeak_t ;

class tower_t {
public:
  ap_uint<16> data;

  tower_t(){
    data = 0;
  }
  tower_t& operator=(const tower_t& rhs){
    data = rhs.data;
    return *this;
  }

  tower_t(ap_uint<12> et, ap_uint<3> hoe, ap_uint<1> satur){
    data = (et) | 
      (((ap_uint<16>) hoe)  << 12) | 
      (((ap_uint<16>) satur)  << 15);
  }

  ap_uint<12> et() {return (data & 0xFFF);}
  ap_uint<3> hoe() {return ((data >> 12) & 0x7);}
  ap_uint<1> satur() {return ((data >> 15) & 0x2);}
  
  operator uint16_t() {return (uint16_t) data;}

};

/*******************************************************************/

//////////////////////////////////////////////////////////////////////////  

// Helper functions

//--------------------------------------------------------// 

ecalRegion_t initStructure(crystal temporary[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI]){

  ap_uint<5> Eta = 0x0 ;
  ap_uint<5> Phi = 0x0 ;

  ecalRegion_t out;

  out.etaStrip0.cr0.energy=temporary[Eta+0][Phi+0].energy; out.etaStrip0.cr0.eta=0; out.etaStrip0.cr0.phi=0;
  out.etaStrip0.cr1.energy=temporary[Eta+0][Phi+1].energy; out.etaStrip0.cr1.eta=0; out.etaStrip0.cr1.phi=1;
  out.etaStrip0.cr2.energy=temporary[Eta+0][Phi+2].energy; out.etaStrip0.cr2.eta=0; out.etaStrip0.cr2.phi=2;
  out.etaStrip0.cr3.energy=temporary[Eta+0][Phi+3].energy; out.etaStrip0.cr3.eta=0; out.etaStrip0.cr3.phi=3;
  out.etaStrip0.cr4.energy=temporary[Eta+0][Phi+4].energy; out.etaStrip0.cr4.eta=0; out.etaStrip0.cr4.phi=4;
  out.etaStrip0.cr5.energy=temporary[Eta+0][Phi+5].energy; out.etaStrip0.cr5.eta=0; out.etaStrip0.cr5.phi=5;
  out.etaStrip0.cr6.energy=temporary[Eta+0][Phi+6].energy; out.etaStrip0.cr6.eta=0; out.etaStrip0.cr6.phi=6;
  out.etaStrip0.cr7.energy=temporary[Eta+0][Phi+7].energy; out.etaStrip0.cr7.eta=0; out.etaStrip0.cr7.phi=7;
  out.etaStrip0.cr8.energy=temporary[Eta+0][Phi+8].energy; out.etaStrip0.cr8.eta=0; out.etaStrip0.cr8.phi=8;
  out.etaStrip0.cr9.energy=temporary[Eta+0][Phi+9].energy; out.etaStrip0.cr9.eta=0; out.etaStrip0.cr9.phi=9;
  out.etaStrip0.cr10.energy=temporary[Eta+0][Phi+10].energy; out.etaStrip0.cr10.eta=0; out.etaStrip0.cr10.phi=10;
  out.etaStrip0.cr11.energy=temporary[Eta+0][Phi+11].energy; out.etaStrip0.cr11.eta=0; out.etaStrip0.cr11.phi=11;
  out.etaStrip0.cr12.energy=temporary[Eta+0][Phi+12].energy; out.etaStrip0.cr12.eta=0; out.etaStrip0.cr12.phi=12;
  out.etaStrip0.cr13.energy=temporary[Eta+0][Phi+13].energy; out.etaStrip0.cr13.eta=0; out.etaStrip0.cr13.phi=13;
  out.etaStrip0.cr14.energy=temporary[Eta+0][Phi+14].energy; out.etaStrip0.cr14.eta=0; out.etaStrip0.cr14.phi=14;
  out.etaStrip0.cr15.energy=temporary[Eta+0][Phi+15].energy; out.etaStrip0.cr15.eta=0; out.etaStrip0.cr15.phi=15;
  out.etaStrip0.cr16.energy=temporary[Eta+0][Phi+16].energy; out.etaStrip0.cr16.eta=0; out.etaStrip0.cr16.phi=16;
  out.etaStrip0.cr17.energy=temporary[Eta+0][Phi+17].energy; out.etaStrip0.cr17.eta=0; out.etaStrip0.cr17.phi=17;
  out.etaStrip0.cr18.energy=temporary[Eta+0][Phi+18].energy; out.etaStrip0.cr18.eta=0; out.etaStrip0.cr18.phi=18;
  out.etaStrip0.cr19.energy=temporary[Eta+0][Phi+19].energy; out.etaStrip0.cr19.eta=0; out.etaStrip0.cr19.phi=19;

  out.etaStrip1.cr0.energy=temporary[Eta+1][Phi+0].energy; out.etaStrip1.cr0.eta=1; out.etaStrip1.cr0.phi=0;
  out.etaStrip1.cr1.energy=temporary[Eta+1][Phi+1].energy; out.etaStrip1.cr1.eta=1; out.etaStrip1.cr1.phi=1;
  out.etaStrip1.cr2.energy=temporary[Eta+1][Phi+2].energy; out.etaStrip1.cr2.eta=1; out.etaStrip1.cr2.phi=2;
  out.etaStrip1.cr3.energy=temporary[Eta+1][Phi+3].energy; out.etaStrip1.cr3.eta=1; out.etaStrip1.cr3.phi=3;
  out.etaStrip1.cr4.energy=temporary[Eta+1][Phi+4].energy; out.etaStrip1.cr4.eta=1; out.etaStrip1.cr4.phi=4;
  out.etaStrip1.cr5.energy=temporary[Eta+1][Phi+5].energy; out.etaStrip1.cr5.eta=1; out.etaStrip1.cr5.phi=5;
  out.etaStrip1.cr6.energy=temporary[Eta+1][Phi+6].energy; out.etaStrip1.cr6.eta=1; out.etaStrip1.cr6.phi=6;
  out.etaStrip1.cr7.energy=temporary[Eta+1][Phi+7].energy; out.etaStrip1.cr7.eta=1; out.etaStrip1.cr7.phi=7;
  out.etaStrip1.cr8.energy=temporary[Eta+1][Phi+8].energy; out.etaStrip1.cr8.eta=1; out.etaStrip1.cr8.phi=8;
  out.etaStrip1.cr9.energy=temporary[Eta+1][Phi+9].energy; out.etaStrip1.cr9.eta=1; out.etaStrip1.cr9.phi=9;
  out.etaStrip1.cr10.energy=temporary[Eta+1][Phi+10].energy; out.etaStrip1.cr10.eta=1; out.etaStrip1.cr10.phi=10;
  out.etaStrip1.cr11.energy=temporary[Eta+1][Phi+11].energy; out.etaStrip1.cr11.eta=1; out.etaStrip1.cr11.phi=11;
  out.etaStrip1.cr12.energy=temporary[Eta+1][Phi+12].energy; out.etaStrip1.cr12.eta=1; out.etaStrip1.cr12.phi=12;
  out.etaStrip1.cr13.energy=temporary[Eta+1][Phi+13].energy; out.etaStrip1.cr13.eta=1; out.etaStrip1.cr13.phi=13;
  out.etaStrip1.cr14.energy=temporary[Eta+1][Phi+14].energy; out.etaStrip1.cr14.eta=1; out.etaStrip1.cr14.phi=14;
  out.etaStrip1.cr15.energy=temporary[Eta+1][Phi+15].energy; out.etaStrip1.cr15.eta=1; out.etaStrip1.cr15.phi=15;
  out.etaStrip1.cr16.energy=temporary[Eta+1][Phi+16].energy; out.etaStrip1.cr16.eta=1; out.etaStrip1.cr16.phi=16;
  out.etaStrip1.cr17.energy=temporary[Eta+1][Phi+17].energy; out.etaStrip1.cr17.eta=1; out.etaStrip1.cr17.phi=17;
  out.etaStrip1.cr18.energy=temporary[Eta+1][Phi+18].energy; out.etaStrip1.cr18.eta=1; out.etaStrip1.cr18.phi=18;
  out.etaStrip1.cr19.energy=temporary[Eta+1][Phi+19].energy; out.etaStrip1.cr19.eta=1; out.etaStrip1.cr19.phi=19;

  out.etaStrip2.cr0.energy=temporary[Eta+2][Phi+0].energy; out.etaStrip2.cr0.eta=2; out.etaStrip2.cr0.phi=0;
  out.etaStrip2.cr1.energy=temporary[Eta+2][Phi+1].energy; out.etaStrip2.cr1.eta=2; out.etaStrip2.cr1.phi=1;
  out.etaStrip2.cr2.energy=temporary[Eta+2][Phi+2].energy; out.etaStrip2.cr2.eta=2; out.etaStrip2.cr2.phi=2;
  out.etaStrip2.cr3.energy=temporary[Eta+2][Phi+3].energy; out.etaStrip2.cr3.eta=2; out.etaStrip2.cr3.phi=3;
  out.etaStrip2.cr4.energy=temporary[Eta+2][Phi+4].energy; out.etaStrip2.cr4.eta=2; out.etaStrip2.cr4.phi=4;
  out.etaStrip2.cr5.energy=temporary[Eta+2][Phi+5].energy; out.etaStrip2.cr5.eta=2; out.etaStrip2.cr5.phi=5;
  out.etaStrip2.cr6.energy=temporary[Eta+2][Phi+6].energy; out.etaStrip2.cr6.eta=2; out.etaStrip2.cr6.phi=6;
  out.etaStrip2.cr7.energy=temporary[Eta+2][Phi+7].energy; out.etaStrip2.cr7.eta=2; out.etaStrip2.cr7.phi=7;
  out.etaStrip2.cr8.energy=temporary[Eta+2][Phi+8].energy; out.etaStrip2.cr8.eta=2; out.etaStrip2.cr8.phi=8;
  out.etaStrip2.cr9.energy=temporary[Eta+2][Phi+9].energy; out.etaStrip2.cr9.eta=2; out.etaStrip2.cr9.phi=9;
  out.etaStrip2.cr10.energy=temporary[Eta+2][Phi+10].energy; out.etaStrip2.cr10.eta=2; out.etaStrip2.cr10.phi=10;
  out.etaStrip2.cr11.energy=temporary[Eta+2][Phi+11].energy; out.etaStrip2.cr11.eta=2; out.etaStrip2.cr11.phi=11;
  out.etaStrip2.cr12.energy=temporary[Eta+2][Phi+12].energy; out.etaStrip2.cr12.eta=2; out.etaStrip2.cr12.phi=12;
  out.etaStrip2.cr13.energy=temporary[Eta+2][Phi+13].energy; out.etaStrip2.cr13.eta=2; out.etaStrip2.cr13.phi=13;
  out.etaStrip2.cr14.energy=temporary[Eta+2][Phi+14].energy; out.etaStrip2.cr14.eta=2; out.etaStrip2.cr14.phi=14;
  out.etaStrip2.cr15.energy=temporary[Eta+2][Phi+15].energy; out.etaStrip2.cr15.eta=2; out.etaStrip2.cr15.phi=15;
  out.etaStrip2.cr16.energy=temporary[Eta+2][Phi+16].energy; out.etaStrip2.cr16.eta=2; out.etaStrip2.cr16.phi=16;
  out.etaStrip2.cr17.energy=temporary[Eta+2][Phi+17].energy; out.etaStrip2.cr17.eta=2; out.etaStrip2.cr17.phi=17;
  out.etaStrip2.cr18.energy=temporary[Eta+2][Phi+18].energy; out.etaStrip2.cr18.eta=2; out.etaStrip2.cr18.phi=18;
  out.etaStrip2.cr19.energy=temporary[Eta+2][Phi+19].energy; out.etaStrip2.cr19.eta=2; out.etaStrip2.cr19.phi=19;

  out.etaStrip3.cr0.energy=temporary[Eta+3][Phi+0].energy; out.etaStrip3.cr0.eta=3; out.etaStrip3.cr0.phi=0;
  out.etaStrip3.cr1.energy=temporary[Eta+3][Phi+1].energy; out.etaStrip3.cr1.eta=3; out.etaStrip3.cr1.phi=1;
  out.etaStrip3.cr2.energy=temporary[Eta+3][Phi+2].energy; out.etaStrip3.cr2.eta=3; out.etaStrip3.cr2.phi=2;
  out.etaStrip3.cr3.energy=temporary[Eta+3][Phi+3].energy; out.etaStrip3.cr3.eta=3; out.etaStrip3.cr3.phi=3;
  out.etaStrip3.cr4.energy=temporary[Eta+3][Phi+4].energy; out.etaStrip3.cr4.eta=3; out.etaStrip3.cr4.phi=4;
  out.etaStrip3.cr5.energy=temporary[Eta+3][Phi+5].energy; out.etaStrip3.cr5.eta=3; out.etaStrip3.cr5.phi=5;
  out.etaStrip3.cr6.energy=temporary[Eta+3][Phi+6].energy; out.etaStrip3.cr6.eta=3; out.etaStrip3.cr6.phi=6;
  out.etaStrip3.cr7.energy=temporary[Eta+3][Phi+7].energy; out.etaStrip3.cr7.eta=3; out.etaStrip3.cr7.phi=7;
  out.etaStrip3.cr8.energy=temporary[Eta+3][Phi+8].energy; out.etaStrip3.cr8.eta=3; out.etaStrip3.cr8.phi=8;
  out.etaStrip3.cr9.energy=temporary[Eta+3][Phi+9].energy; out.etaStrip3.cr9.eta=3; out.etaStrip3.cr9.phi=9;
  out.etaStrip3.cr10.energy=temporary[Eta+3][Phi+10].energy; out.etaStrip3.cr10.eta=3; out.etaStrip3.cr10.phi=10;
  out.etaStrip3.cr11.energy=temporary[Eta+3][Phi+11].energy; out.etaStrip3.cr11.eta=3; out.etaStrip3.cr11.phi=11;
  out.etaStrip3.cr12.energy=temporary[Eta+3][Phi+12].energy; out.etaStrip3.cr12.eta=3; out.etaStrip3.cr12.phi=12;
  out.etaStrip3.cr13.energy=temporary[Eta+3][Phi+13].energy; out.etaStrip3.cr13.eta=3; out.etaStrip3.cr13.phi=13;
  out.etaStrip3.cr14.energy=temporary[Eta+3][Phi+14].energy; out.etaStrip3.cr14.eta=3; out.etaStrip3.cr14.phi=14;
  out.etaStrip3.cr15.energy=temporary[Eta+3][Phi+15].energy; out.etaStrip3.cr15.eta=3; out.etaStrip3.cr15.phi=15;
  out.etaStrip3.cr16.energy=temporary[Eta+3][Phi+16].energy; out.etaStrip3.cr16.eta=3; out.etaStrip3.cr16.phi=16;
  out.etaStrip3.cr17.energy=temporary[Eta+3][Phi+17].energy; out.etaStrip3.cr17.eta=3; out.etaStrip3.cr17.phi=17;
  out.etaStrip3.cr18.energy=temporary[Eta+3][Phi+18].energy; out.etaStrip3.cr18.eta=3; out.etaStrip3.cr18.phi=18;
  out.etaStrip3.cr19.energy=temporary[Eta+3][Phi+19].energy; out.etaStrip3.cr19.eta=3; out.etaStrip3.cr19.phi=19;

  out.etaStrip4.cr0.energy=temporary[Eta+4][Phi+0].energy; out.etaStrip4.cr0.eta=4; out.etaStrip4.cr0.phi=0;
  out.etaStrip4.cr1.energy=temporary[Eta+4][Phi+1].energy; out.etaStrip4.cr1.eta=4; out.etaStrip4.cr1.phi=1;
  out.etaStrip4.cr2.energy=temporary[Eta+4][Phi+2].energy; out.etaStrip4.cr2.eta=4; out.etaStrip4.cr2.phi=2;
  out.etaStrip4.cr3.energy=temporary[Eta+4][Phi+3].energy; out.etaStrip4.cr3.eta=4; out.etaStrip4.cr3.phi=3;
  out.etaStrip4.cr4.energy=temporary[Eta+4][Phi+4].energy; out.etaStrip4.cr4.eta=4; out.etaStrip4.cr4.phi=4;
  out.etaStrip4.cr5.energy=temporary[Eta+4][Phi+5].energy; out.etaStrip4.cr5.eta=4; out.etaStrip4.cr5.phi=5;
  out.etaStrip4.cr6.energy=temporary[Eta+4][Phi+6].energy; out.etaStrip4.cr6.eta=4; out.etaStrip4.cr6.phi=6;
  out.etaStrip4.cr7.energy=temporary[Eta+4][Phi+7].energy; out.etaStrip4.cr7.eta=4; out.etaStrip4.cr7.phi=7;
  out.etaStrip4.cr8.energy=temporary[Eta+4][Phi+8].energy; out.etaStrip4.cr8.eta=4; out.etaStrip4.cr8.phi=8;
  out.etaStrip4.cr9.energy=temporary[Eta+4][Phi+9].energy; out.etaStrip4.cr9.eta=4; out.etaStrip4.cr9.phi=9;
  out.etaStrip4.cr10.energy=temporary[Eta+4][Phi+10].energy; out.etaStrip4.cr10.eta=4; out.etaStrip4.cr10.phi=10;
  out.etaStrip4.cr11.energy=temporary[Eta+4][Phi+11].energy; out.etaStrip4.cr11.eta=4; out.etaStrip4.cr11.phi=11;
  out.etaStrip4.cr12.energy=temporary[Eta+4][Phi+12].energy; out.etaStrip4.cr12.eta=4; out.etaStrip4.cr12.phi=12;
  out.etaStrip4.cr13.energy=temporary[Eta+4][Phi+13].energy; out.etaStrip4.cr13.eta=4; out.etaStrip4.cr13.phi=13;
  out.etaStrip4.cr14.energy=temporary[Eta+4][Phi+14].energy; out.etaStrip4.cr14.eta=4; out.etaStrip4.cr14.phi=14;
  out.etaStrip4.cr15.energy=temporary[Eta+4][Phi+15].energy; out.etaStrip4.cr15.eta=4; out.etaStrip4.cr15.phi=15;
  out.etaStrip4.cr16.energy=temporary[Eta+4][Phi+16].energy; out.etaStrip4.cr16.eta=4; out.etaStrip4.cr16.phi=16;
  out.etaStrip4.cr17.energy=temporary[Eta+4][Phi+17].energy; out.etaStrip4.cr17.eta=4; out.etaStrip4.cr17.phi=17;
  out.etaStrip4.cr18.energy=temporary[Eta+4][Phi+18].energy; out.etaStrip4.cr18.eta=4; out.etaStrip4.cr18.phi=18;
  out.etaStrip4.cr19.energy=temporary[Eta+4][Phi+19].energy; out.etaStrip4.cr19.eta=4; out.etaStrip4.cr19.phi=19;

  out.etaStrip5.cr0.energy=temporary[Eta+5][Phi+0].energy; out.etaStrip5.cr0.eta=5; out.etaStrip5.cr0.phi=0;
  out.etaStrip5.cr1.energy=temporary[Eta+5][Phi+1].energy; out.etaStrip5.cr1.eta=5; out.etaStrip5.cr1.phi=1;
  out.etaStrip5.cr2.energy=temporary[Eta+5][Phi+2].energy; out.etaStrip5.cr2.eta=5; out.etaStrip5.cr2.phi=2;
  out.etaStrip5.cr3.energy=temporary[Eta+5][Phi+3].energy; out.etaStrip5.cr3.eta=5; out.etaStrip5.cr3.phi=3;
  out.etaStrip5.cr4.energy=temporary[Eta+5][Phi+4].energy; out.etaStrip5.cr4.eta=5; out.etaStrip5.cr4.phi=4;
  out.etaStrip5.cr5.energy=temporary[Eta+5][Phi+5].energy; out.etaStrip5.cr5.eta=5; out.etaStrip5.cr5.phi=5;
  out.etaStrip5.cr6.energy=temporary[Eta+5][Phi+6].energy; out.etaStrip5.cr6.eta=5; out.etaStrip5.cr6.phi=6;
  out.etaStrip5.cr7.energy=temporary[Eta+5][Phi+7].energy; out.etaStrip5.cr7.eta=5; out.etaStrip5.cr7.phi=7;
  out.etaStrip5.cr8.energy=temporary[Eta+5][Phi+8].energy; out.etaStrip5.cr8.eta=5; out.etaStrip5.cr8.phi=8;
  out.etaStrip5.cr9.energy=temporary[Eta+5][Phi+9].energy; out.etaStrip5.cr9.eta=5; out.etaStrip5.cr9.phi=9;
  out.etaStrip5.cr10.energy=temporary[Eta+5][Phi+10].energy; out.etaStrip5.cr10.eta=5; out.etaStrip5.cr10.phi=10;
  out.etaStrip5.cr11.energy=temporary[Eta+5][Phi+11].energy; out.etaStrip5.cr11.eta=5; out.etaStrip5.cr11.phi=11;
  out.etaStrip5.cr12.energy=temporary[Eta+5][Phi+12].energy; out.etaStrip5.cr12.eta=5; out.etaStrip5.cr12.phi=12;
  out.etaStrip5.cr13.energy=temporary[Eta+5][Phi+13].energy; out.etaStrip5.cr13.eta=5; out.etaStrip5.cr13.phi=13;
  out.etaStrip5.cr14.energy=temporary[Eta+5][Phi+14].energy; out.etaStrip5.cr14.eta=5; out.etaStrip5.cr14.phi=14;
  out.etaStrip5.cr15.energy=temporary[Eta+5][Phi+15].energy; out.etaStrip5.cr15.eta=5; out.etaStrip5.cr15.phi=15;
  out.etaStrip5.cr16.energy=temporary[Eta+5][Phi+16].energy; out.etaStrip5.cr16.eta=5; out.etaStrip5.cr16.phi=16;
  out.etaStrip5.cr17.energy=temporary[Eta+5][Phi+17].energy; out.etaStrip5.cr17.eta=5; out.etaStrip5.cr17.phi=17;
  out.etaStrip5.cr18.energy=temporary[Eta+5][Phi+18].energy; out.etaStrip5.cr18.eta=5; out.etaStrip5.cr18.phi=18;
  out.etaStrip5.cr19.energy=temporary[Eta+5][Phi+19].energy; out.etaStrip5.cr19.eta=5; out.etaStrip5.cr19.phi=19;

  out.etaStrip6.cr0.energy=temporary[Eta+6][Phi+0].energy; out.etaStrip6.cr0.eta=6; out.etaStrip6.cr0.phi=0;
  out.etaStrip6.cr1.energy=temporary[Eta+6][Phi+1].energy; out.etaStrip6.cr1.eta=6; out.etaStrip6.cr1.phi=1;
  out.etaStrip6.cr2.energy=temporary[Eta+6][Phi+2].energy; out.etaStrip6.cr2.eta=6; out.etaStrip6.cr2.phi=2;
  out.etaStrip6.cr2.energy=temporary[Eta+6][Phi+2].energy; out.etaStrip6.cr2.eta=6; out.etaStrip6.cr2.phi=2;
  out.etaStrip6.cr3.energy=temporary[Eta+6][Phi+3].energy; out.etaStrip6.cr3.eta=6; out.etaStrip6.cr3.phi=3;
  out.etaStrip6.cr4.energy=temporary[Eta+6][Phi+4].energy; out.etaStrip6.cr4.eta=6; out.etaStrip6.cr4.phi=4;
  out.etaStrip6.cr5.energy=temporary[Eta+6][Phi+5].energy; out.etaStrip6.cr5.eta=6; out.etaStrip6.cr5.phi=5;
  out.etaStrip6.cr6.energy=temporary[Eta+6][Phi+6].energy; out.etaStrip6.cr6.eta=6; out.etaStrip6.cr6.phi=6;
  out.etaStrip6.cr7.energy=temporary[Eta+6][Phi+7].energy; out.etaStrip6.cr7.eta=6; out.etaStrip6.cr7.phi=7;
  out.etaStrip6.cr8.energy=temporary[Eta+6][Phi+8].energy; out.etaStrip6.cr8.eta=6; out.etaStrip6.cr8.phi=8;
  out.etaStrip6.cr9.energy=temporary[Eta+6][Phi+9].energy; out.etaStrip6.cr9.eta=6; out.etaStrip6.cr9.phi=9;
  out.etaStrip6.cr10.energy=temporary[Eta+6][Phi+10].energy; out.etaStrip6.cr10.eta=6; out.etaStrip6.cr10.phi=10;
  out.etaStrip6.cr11.energy=temporary[Eta+6][Phi+11].energy; out.etaStrip6.cr11.eta=6; out.etaStrip6.cr11.phi=11;
  out.etaStrip6.cr12.energy=temporary[Eta+6][Phi+12].energy; out.etaStrip6.cr12.eta=6; out.etaStrip6.cr12.phi=12;
  out.etaStrip6.cr13.energy=temporary[Eta+6][Phi+13].energy; out.etaStrip6.cr13.eta=6; out.etaStrip6.cr13.phi=13;
  out.etaStrip6.cr14.energy=temporary[Eta+6][Phi+14].energy; out.etaStrip6.cr14.eta=6; out.etaStrip6.cr14.phi=14;
  out.etaStrip6.cr15.energy=temporary[Eta+6][Phi+15].energy; out.etaStrip6.cr15.eta=6; out.etaStrip6.cr15.phi=15;
  out.etaStrip6.cr16.energy=temporary[Eta+6][Phi+16].energy; out.etaStrip6.cr16.eta=6; out.etaStrip6.cr16.phi=16;
  out.etaStrip6.cr17.energy=temporary[Eta+6][Phi+17].energy; out.etaStrip6.cr17.eta=6; out.etaStrip6.cr17.phi=17;
  out.etaStrip6.cr18.energy=temporary[Eta+6][Phi+18].energy; out.etaStrip6.cr18.eta=6; out.etaStrip6.cr18.phi=18;
  out.etaStrip6.cr19.energy=temporary[Eta+6][Phi+19].energy; out.etaStrip6.cr19.eta=6; out.etaStrip6.cr19.phi=19;

  out.etaStrip7.cr0.energy=temporary[Eta+7][Phi+0].energy; out.etaStrip7.cr0.eta=7; out.etaStrip7.cr0.phi=0;
  out.etaStrip7.cr1.energy=temporary[Eta+7][Phi+1].energy; out.etaStrip7.cr1.eta=7; out.etaStrip7.cr1.phi=1;
  out.etaStrip7.cr2.energy=temporary[Eta+7][Phi+2].energy; out.etaStrip7.cr2.eta=7; out.etaStrip7.cr2.phi=2;
  out.etaStrip7.cr3.energy=temporary[Eta+7][Phi+3].energy; out.etaStrip7.cr3.eta=7; out.etaStrip7.cr3.phi=3;
  out.etaStrip7.cr4.energy=temporary[Eta+7][Phi+4].energy; out.etaStrip7.cr4.eta=7; out.etaStrip7.cr4.phi=4;
  out.etaStrip7.cr5.energy=temporary[Eta+7][Phi+5].energy; out.etaStrip7.cr5.eta=7; out.etaStrip7.cr5.phi=5;
  out.etaStrip7.cr6.energy=temporary[Eta+7][Phi+6].energy; out.etaStrip7.cr6.eta=7; out.etaStrip7.cr6.phi=6;
  out.etaStrip7.cr7.energy=temporary[Eta+7][Phi+7].energy; out.etaStrip7.cr7.eta=7; out.etaStrip7.cr7.phi=7;
  out.etaStrip7.cr8.energy=temporary[Eta+7][Phi+8].energy; out.etaStrip7.cr8.eta=7; out.etaStrip7.cr8.phi=8;
  out.etaStrip7.cr9.energy=temporary[Eta+7][Phi+9].energy; out.etaStrip7.cr9.eta=7; out.etaStrip7.cr9.phi=9;
  out.etaStrip7.cr10.energy=temporary[Eta+7][Phi+10].energy; out.etaStrip7.cr10.eta=7; out.etaStrip7.cr10.phi=10;
  out.etaStrip7.cr11.energy=temporary[Eta+7][Phi+11].energy; out.etaStrip7.cr11.eta=7; out.etaStrip7.cr11.phi=11;
  out.etaStrip7.cr12.energy=temporary[Eta+7][Phi+12].energy; out.etaStrip7.cr12.eta=7; out.etaStrip7.cr12.phi=12;
  out.etaStrip7.cr13.energy=temporary[Eta+7][Phi+13].energy; out.etaStrip7.cr13.eta=7; out.etaStrip7.cr13.phi=13;
  out.etaStrip7.cr14.energy=temporary[Eta+7][Phi+14].energy; out.etaStrip7.cr14.eta=7; out.etaStrip7.cr14.phi=14;
  out.etaStrip7.cr15.energy=temporary[Eta+7][Phi+15].energy; out.etaStrip7.cr15.eta=7; out.etaStrip7.cr15.phi=15;
  out.etaStrip7.cr16.energy=temporary[Eta+7][Phi+16].energy; out.etaStrip7.cr16.eta=7; out.etaStrip7.cr16.phi=16;
  out.etaStrip7.cr17.energy=temporary[Eta+7][Phi+17].energy; out.etaStrip7.cr17.eta=7; out.etaStrip7.cr17.phi=17;
  out.etaStrip7.cr18.energy=temporary[Eta+7][Phi+18].energy; out.etaStrip7.cr18.eta=7; out.etaStrip7.cr18.phi=18;
  out.etaStrip7.cr19.energy=temporary[Eta+7][Phi+19].energy; out.etaStrip7.cr19.eta=7; out.etaStrip7.cr19.phi=19;

  out.etaStrip8.cr0.energy=temporary[Eta+8][Phi+0].energy; out.etaStrip8.cr0.eta=8; out.etaStrip8.cr0.phi=0;
  out.etaStrip8.cr1.energy=temporary[Eta+8][Phi+1].energy; out.etaStrip8.cr1.eta=8; out.etaStrip8.cr1.phi=1;
  out.etaStrip8.cr2.energy=temporary[Eta+8][Phi+2].energy; out.etaStrip8.cr2.eta=8; out.etaStrip8.cr2.phi=2;
  out.etaStrip8.cr3.energy=temporary[Eta+8][Phi+3].energy; out.etaStrip8.cr3.eta=8; out.etaStrip8.cr3.phi=3;
  out.etaStrip8.cr4.energy=temporary[Eta+8][Phi+4].energy; out.etaStrip8.cr4.eta=8; out.etaStrip8.cr4.phi=4;
  out.etaStrip8.cr5.energy=temporary[Eta+8][Phi+5].energy; out.etaStrip8.cr5.eta=8; out.etaStrip8.cr5.phi=5;
  out.etaStrip8.cr6.energy=temporary[Eta+8][Phi+6].energy; out.etaStrip8.cr6.eta=8; out.etaStrip8.cr6.phi=6;
  out.etaStrip8.cr7.energy=temporary[Eta+8][Phi+7].energy; out.etaStrip8.cr7.eta=8; out.etaStrip8.cr7.phi=7;
  out.etaStrip8.cr8.energy=temporary[Eta+8][Phi+8].energy; out.etaStrip8.cr8.eta=8; out.etaStrip8.cr8.phi=8;
  out.etaStrip8.cr9.energy=temporary[Eta+8][Phi+9].energy; out.etaStrip8.cr9.eta=8; out.etaStrip8.cr9.phi=9;
  out.etaStrip8.cr10.energy=temporary[Eta+8][Phi+10].energy; out.etaStrip8.cr10.eta=8; out.etaStrip8.cr10.phi=10;
  out.etaStrip8.cr11.energy=temporary[Eta+8][Phi+11].energy; out.etaStrip8.cr11.eta=8; out.etaStrip8.cr11.phi=11;
  out.etaStrip8.cr12.energy=temporary[Eta+8][Phi+12].energy; out.etaStrip8.cr12.eta=8; out.etaStrip8.cr12.phi=12;
  out.etaStrip8.cr13.energy=temporary[Eta+8][Phi+13].energy; out.etaStrip8.cr13.eta=8; out.etaStrip8.cr13.phi=13;
  out.etaStrip8.cr14.energy=temporary[Eta+8][Phi+14].energy; out.etaStrip8.cr14.eta=8; out.etaStrip8.cr14.phi=14;
  out.etaStrip8.cr15.energy=temporary[Eta+8][Phi+15].energy; out.etaStrip8.cr15.eta=8; out.etaStrip8.cr15.phi=15;
  out.etaStrip8.cr16.energy=temporary[Eta+8][Phi+16].energy; out.etaStrip8.cr16.eta=8; out.etaStrip8.cr16.phi=16;
  out.etaStrip8.cr17.energy=temporary[Eta+8][Phi+17].energy; out.etaStrip8.cr17.eta=8; out.etaStrip8.cr17.phi=17;
  out.etaStrip8.cr18.energy=temporary[Eta+8][Phi+18].energy; out.etaStrip8.cr18.eta=8; out.etaStrip8.cr18.phi=18;
  out.etaStrip8.cr19.energy=temporary[Eta+8][Phi+19].energy; out.etaStrip8.cr19.eta=8; out.etaStrip8.cr19.phi=19;

  out.etaStrip9.cr0.energy=temporary[Eta+9][Phi+0].energy; out.etaStrip9.cr0.eta=9; out.etaStrip9.cr0.phi=0;
  out.etaStrip9.cr1.energy=temporary[Eta+9][Phi+1].energy; out.etaStrip9.cr1.eta=9; out.etaStrip9.cr1.phi=1;
  out.etaStrip9.cr2.energy=temporary[Eta+9][Phi+2].energy; out.etaStrip9.cr2.eta=9; out.etaStrip9.cr2.phi=2;
  out.etaStrip9.cr3.energy=temporary[Eta+9][Phi+3].energy; out.etaStrip9.cr3.eta=9; out.etaStrip9.cr3.phi=3;
  out.etaStrip9.cr4.energy=temporary[Eta+9][Phi+4].energy; out.etaStrip9.cr4.eta=9; out.etaStrip9.cr4.phi=4;
  out.etaStrip9.cr5.energy=temporary[Eta+9][Phi+5].energy; out.etaStrip9.cr5.eta=9; out.etaStrip9.cr5.phi=5;
  out.etaStrip9.cr6.energy=temporary[Eta+9][Phi+6].energy; out.etaStrip9.cr6.eta=9; out.etaStrip9.cr6.phi=6;
  out.etaStrip9.cr7.energy=temporary[Eta+9][Phi+7].energy; out.etaStrip9.cr7.eta=9; out.etaStrip9.cr7.phi=7;
  out.etaStrip9.cr8.energy=temporary[Eta+9][Phi+8].energy; out.etaStrip9.cr8.eta=9; out.etaStrip9.cr8.phi=8;
  out.etaStrip9.cr9.energy=temporary[Eta+9][Phi+9].energy; out.etaStrip9.cr9.eta=9; out.etaStrip9.cr9.phi=9;
  out.etaStrip9.cr10.energy=temporary[Eta+9][Phi+10].energy; out.etaStrip9.cr10.eta=9; out.etaStrip9.cr10.phi=10;
  out.etaStrip9.cr11.energy=temporary[Eta+9][Phi+11].energy; out.etaStrip9.cr11.eta=9; out.etaStrip9.cr11.phi=11;
  out.etaStrip9.cr12.energy=temporary[Eta+9][Phi+12].energy; out.etaStrip9.cr12.eta=9; out.etaStrip9.cr12.phi=12;
  out.etaStrip9.cr13.energy=temporary[Eta+9][Phi+13].energy; out.etaStrip9.cr13.eta=9; out.etaStrip9.cr13.phi=13;
  out.etaStrip9.cr14.energy=temporary[Eta+9][Phi+14].energy; out.etaStrip9.cr14.eta=9; out.etaStrip9.cr14.phi=14;
  out.etaStrip9.cr15.energy=temporary[Eta+9][Phi+15].energy; out.etaStrip9.cr15.eta=9; out.etaStrip9.cr15.phi=15;
  out.etaStrip9.cr16.energy=temporary[Eta+9][Phi+16].energy; out.etaStrip9.cr16.eta=9; out.etaStrip9.cr16.phi=16;
  out.etaStrip9.cr17.energy=temporary[Eta+9][Phi+17].energy; out.etaStrip9.cr17.eta=9; out.etaStrip9.cr17.phi=17;
  out.etaStrip9.cr18.energy=temporary[Eta+9][Phi+18].energy; out.etaStrip9.cr18.eta=9; out.etaStrip9.cr18.phi=18;
  out.etaStrip9.cr19.energy=temporary[Eta+9][Phi+19].energy; out.etaStrip9.cr19.eta=9; out.etaStrip9.cr19.phi=19;

  out.etaStrip10.cr0.energy=temporary[Eta+10][Phi+0].energy; out.etaStrip10.cr0.eta=10; out.etaStrip10.cr0.phi=0;
  out.etaStrip10.cr1.energy=temporary[Eta+10][Phi+1].energy; out.etaStrip10.cr1.eta=10; out.etaStrip10.cr1.phi=1;
  out.etaStrip10.cr2.energy=temporary[Eta+10][Phi+2].energy; out.etaStrip10.cr2.eta=10; out.etaStrip10.cr2.phi=2;
  out.etaStrip10.cr3.energy=temporary[Eta+10][Phi+3].energy; out.etaStrip10.cr3.eta=10; out.etaStrip10.cr3.phi=3;
  out.etaStrip10.cr4.energy=temporary[Eta+10][Phi+4].energy; out.etaStrip10.cr4.eta=10; out.etaStrip10.cr4.phi=4;
  out.etaStrip10.cr5.energy=temporary[Eta+10][Phi+5].energy; out.etaStrip10.cr5.eta=10; out.etaStrip10.cr5.phi=5;
  out.etaStrip10.cr6.energy=temporary[Eta+10][Phi+6].energy; out.etaStrip10.cr6.eta=10; out.etaStrip10.cr6.phi=6;
  out.etaStrip10.cr7.energy=temporary[Eta+10][Phi+7].energy; out.etaStrip10.cr7.eta=10; out.etaStrip10.cr7.phi=7;
  out.etaStrip10.cr8.energy=temporary[Eta+10][Phi+8].energy; out.etaStrip10.cr8.eta=10; out.etaStrip10.cr8.phi=8;
  out.etaStrip10.cr9.energy=temporary[Eta+10][Phi+9].energy; out.etaStrip10.cr9.eta=10; out.etaStrip10.cr9.phi=9;
  out.etaStrip10.cr10.energy=temporary[Eta+10][Phi+10].energy; out.etaStrip10.cr10.eta=10; out.etaStrip10.cr10.phi=10;
  out.etaStrip10.cr11.energy=temporary[Eta+10][Phi+11].energy; out.etaStrip10.cr11.eta=10; out.etaStrip10.cr11.phi=11;
  out.etaStrip10.cr12.energy=temporary[Eta+10][Phi+12].energy; out.etaStrip10.cr12.eta=10; out.etaStrip10.cr12.phi=12;
  out.etaStrip10.cr13.energy=temporary[Eta+10][Phi+13].energy; out.etaStrip10.cr13.eta=10; out.etaStrip10.cr13.phi=13;
  out.etaStrip10.cr14.energy=temporary[Eta+10][Phi+14].energy; out.etaStrip10.cr14.eta=10; out.etaStrip10.cr14.phi=14;
  out.etaStrip10.cr15.energy=temporary[Eta+10][Phi+15].energy; out.etaStrip10.cr15.eta=10; out.etaStrip10.cr15.phi=15;
  out.etaStrip10.cr16.energy=temporary[Eta+10][Phi+16].energy; out.etaStrip10.cr16.eta=10; out.etaStrip10.cr16.phi=16;
  out.etaStrip10.cr17.energy=temporary[Eta+10][Phi+17].energy; out.etaStrip10.cr17.eta=10; out.etaStrip10.cr17.phi=17;
  out.etaStrip10.cr18.energy=temporary[Eta+10][Phi+18].energy; out.etaStrip10.cr18.eta=10; out.etaStrip10.cr18.phi=18;
  out.etaStrip10.cr19.energy=temporary[Eta+10][Phi+19].energy; out.etaStrip10.cr19.eta=10; out.etaStrip10.cr19.phi=19;

  out.etaStrip11.cr0.energy=temporary[Eta+11][Phi+0].energy; out.etaStrip11.cr0.eta=11; out.etaStrip11.cr0.phi=0;
  out.etaStrip11.cr1.energy=temporary[Eta+11][Phi+1].energy; out.etaStrip11.cr1.eta=11; out.etaStrip11.cr1.phi=1;
  out.etaStrip11.cr2.energy=temporary[Eta+11][Phi+2].energy; out.etaStrip11.cr2.eta=11; out.etaStrip11.cr2.phi=2;
  out.etaStrip11.cr3.energy=temporary[Eta+11][Phi+3].energy; out.etaStrip11.cr3.eta=11; out.etaStrip11.cr3.phi=3;
  out.etaStrip11.cr4.energy=temporary[Eta+11][Phi+4].energy; out.etaStrip11.cr4.eta=11; out.etaStrip11.cr4.phi=4;
  out.etaStrip11.cr5.energy=temporary[Eta+11][Phi+5].energy; out.etaStrip11.cr5.eta=11; out.etaStrip11.cr5.phi=5;
  out.etaStrip11.cr6.energy=temporary[Eta+11][Phi+6].energy; out.etaStrip11.cr6.eta=11; out.etaStrip11.cr6.phi=6;
  out.etaStrip11.cr7.energy=temporary[Eta+11][Phi+7].energy; out.etaStrip11.cr7.eta=11; out.etaStrip11.cr7.phi=7;
  out.etaStrip11.cr8.energy=temporary[Eta+11][Phi+8].energy; out.etaStrip11.cr8.eta=11; out.etaStrip11.cr8.phi=8;
  out.etaStrip11.cr9.energy=temporary[Eta+11][Phi+9].energy; out.etaStrip11.cr9.eta=11; out.etaStrip11.cr9.phi=9;
  out.etaStrip11.cr10.energy=temporary[Eta+11][Phi+10].energy; out.etaStrip11.cr10.eta=11; out.etaStrip11.cr10.phi=10;
  out.etaStrip11.cr11.energy=temporary[Eta+11][Phi+11].energy; out.etaStrip11.cr11.eta=11; out.etaStrip11.cr11.phi=11;
  out.etaStrip11.cr12.energy=temporary[Eta+11][Phi+12].energy; out.etaStrip11.cr12.eta=11; out.etaStrip11.cr12.phi=12;
  out.etaStrip11.cr13.energy=temporary[Eta+11][Phi+13].energy; out.etaStrip11.cr13.eta=11; out.etaStrip11.cr13.phi=13;
  out.etaStrip11.cr14.energy=temporary[Eta+11][Phi+14].energy; out.etaStrip11.cr14.eta=11; out.etaStrip11.cr14.phi=14;
  out.etaStrip11.cr15.energy=temporary[Eta+11][Phi+15].energy; out.etaStrip11.cr15.eta=11; out.etaStrip11.cr15.phi=15;
  out.etaStrip11.cr16.energy=temporary[Eta+11][Phi+16].energy; out.etaStrip11.cr16.eta=11; out.etaStrip11.cr16.phi=16;
  out.etaStrip11.cr17.energy=temporary[Eta+11][Phi+17].energy; out.etaStrip11.cr17.eta=11; out.etaStrip11.cr17.phi=17;
  out.etaStrip11.cr18.energy=temporary[Eta+11][Phi+18].energy; out.etaStrip11.cr18.eta=11; out.etaStrip11.cr18.phi=18;
  out.etaStrip11.cr19.energy=temporary[Eta+11][Phi+19].energy; out.etaStrip11.cr19.eta=11; out.etaStrip11.cr19.phi=19;

  out.etaStrip12.cr0.energy=temporary[Eta+12][Phi+0].energy; out.etaStrip12.cr0.eta=12; out.etaStrip12.cr0.phi=0;
  out.etaStrip12.cr1.energy=temporary[Eta+12][Phi+1].energy; out.etaStrip12.cr1.eta=12; out.etaStrip12.cr1.phi=1;
  out.etaStrip12.cr2.energy=temporary[Eta+12][Phi+2].energy; out.etaStrip12.cr2.eta=12; out.etaStrip12.cr2.phi=2;
  out.etaStrip12.cr3.energy=temporary[Eta+12][Phi+3].energy; out.etaStrip12.cr3.eta=12; out.etaStrip12.cr3.phi=3;
  out.etaStrip12.cr4.energy=temporary[Eta+12][Phi+4].energy; out.etaStrip12.cr4.eta=12; out.etaStrip12.cr4.phi=4;
  out.etaStrip12.cr5.energy=temporary[Eta+12][Phi+5].energy; out.etaStrip12.cr5.eta=12; out.etaStrip12.cr5.phi=5;
  out.etaStrip12.cr6.energy=temporary[Eta+12][Phi+6].energy; out.etaStrip12.cr6.eta=12; out.etaStrip12.cr6.phi=6;
  out.etaStrip12.cr7.energy=temporary[Eta+12][Phi+7].energy; out.etaStrip12.cr7.eta=12; out.etaStrip12.cr7.phi=7;
  out.etaStrip12.cr8.energy=temporary[Eta+12][Phi+8].energy; out.etaStrip12.cr8.eta=12; out.etaStrip12.cr8.phi=8;
  out.etaStrip12.cr9.energy=temporary[Eta+12][Phi+9].energy; out.etaStrip12.cr9.eta=12; out.etaStrip12.cr9.phi=9;
  out.etaStrip12.cr10.energy=temporary[Eta+12][Phi+10].energy; out.etaStrip12.cr10.eta=12; out.etaStrip12.cr10.phi=10;
  out.etaStrip12.cr11.energy=temporary[Eta+12][Phi+11].energy; out.etaStrip12.cr11.eta=12; out.etaStrip12.cr11.phi=11;
  out.etaStrip12.cr12.energy=temporary[Eta+12][Phi+12].energy; out.etaStrip12.cr12.eta=12; out.etaStrip12.cr12.phi=12;
  out.etaStrip12.cr13.energy=temporary[Eta+12][Phi+13].energy; out.etaStrip12.cr13.eta=12; out.etaStrip12.cr13.phi=13;
  out.etaStrip12.cr14.energy=temporary[Eta+12][Phi+14].energy; out.etaStrip12.cr14.eta=12; out.etaStrip12.cr14.phi=14;
  out.etaStrip12.cr15.energy=temporary[Eta+12][Phi+15].energy; out.etaStrip12.cr15.eta=12; out.etaStrip12.cr15.phi=15;
  out.etaStrip12.cr16.energy=temporary[Eta+12][Phi+16].energy; out.etaStrip12.cr16.eta=12; out.etaStrip12.cr16.phi=16;
  out.etaStrip12.cr17.energy=temporary[Eta+12][Phi+17].energy; out.etaStrip12.cr17.eta=12; out.etaStrip12.cr17.phi=17;
  out.etaStrip12.cr18.energy=temporary[Eta+12][Phi+18].energy; out.etaStrip12.cr18.eta=12; out.etaStrip12.cr18.phi=18;
  out.etaStrip12.cr19.energy=temporary[Eta+12][Phi+19].energy; out.etaStrip12.cr19.eta=12; out.etaStrip12.cr19.phi=19;

  out.etaStrip13.cr0.energy=temporary[Eta+13][Phi+0].energy; out.etaStrip13.cr0.eta=13; out.etaStrip13.cr0.phi=0;
  out.etaStrip13.cr1.energy=temporary[Eta+13][Phi+1].energy; out.etaStrip13.cr1.eta=13; out.etaStrip13.cr1.phi=1;
  out.etaStrip13.cr2.energy=temporary[Eta+13][Phi+2].energy; out.etaStrip13.cr2.eta=13; out.etaStrip13.cr2.phi=2;
  out.etaStrip13.cr3.energy=temporary[Eta+13][Phi+3].energy; out.etaStrip13.cr3.eta=13; out.etaStrip13.cr3.phi=3;
  out.etaStrip13.cr4.energy=temporary[Eta+13][Phi+4].energy; out.etaStrip13.cr4.eta=13; out.etaStrip13.cr4.phi=4;
  out.etaStrip13.cr5.energy=temporary[Eta+13][Phi+5].energy; out.etaStrip13.cr5.eta=13; out.etaStrip13.cr5.phi=5;
  out.etaStrip13.cr6.energy=temporary[Eta+13][Phi+6].energy; out.etaStrip13.cr6.eta=13; out.etaStrip13.cr6.phi=6;
  out.etaStrip13.cr7.energy=temporary[Eta+13][Phi+7].energy; out.etaStrip13.cr7.eta=13; out.etaStrip13.cr7.phi=7;
  out.etaStrip13.cr8.energy=temporary[Eta+13][Phi+8].energy; out.etaStrip13.cr8.eta=13; out.etaStrip13.cr8.phi=8;
  out.etaStrip13.cr9.energy=temporary[Eta+13][Phi+9].energy; out.etaStrip13.cr9.eta=13; out.etaStrip13.cr9.phi=9;
  out.etaStrip13.cr10.energy=temporary[Eta+13][Phi+10].energy; out.etaStrip13.cr10.eta=13; out.etaStrip13.cr10.phi=10;
  out.etaStrip13.cr11.energy=temporary[Eta+13][Phi+11].energy; out.etaStrip13.cr11.eta=13; out.etaStrip13.cr11.phi=11;
  out.etaStrip13.cr12.energy=temporary[Eta+13][Phi+12].energy; out.etaStrip13.cr12.eta=13; out.etaStrip13.cr12.phi=12;
  out.etaStrip13.cr13.energy=temporary[Eta+13][Phi+13].energy; out.etaStrip13.cr13.eta=13; out.etaStrip13.cr13.phi=13;
  out.etaStrip13.cr14.energy=temporary[Eta+13][Phi+14].energy; out.etaStrip13.cr14.eta=13; out.etaStrip13.cr14.phi=14;
  out.etaStrip13.cr15.energy=temporary[Eta+13][Phi+15].energy; out.etaStrip13.cr15.eta=13; out.etaStrip13.cr15.phi=15;
  out.etaStrip13.cr16.energy=temporary[Eta+13][Phi+16].energy; out.etaStrip13.cr16.eta=13; out.etaStrip13.cr16.phi=16;
  out.etaStrip13.cr17.energy=temporary[Eta+13][Phi+17].energy; out.etaStrip13.cr17.eta=13; out.etaStrip13.cr17.phi=17;
  out.etaStrip13.cr18.energy=temporary[Eta+13][Phi+18].energy; out.etaStrip13.cr18.eta=13; out.etaStrip13.cr18.phi=18;
  out.etaStrip13.cr19.energy=temporary[Eta+13][Phi+19].energy; out.etaStrip13.cr19.eta=13; out.etaStrip13.cr19.phi=19;
  
  out.etaStrip14.cr0.energy=temporary[Eta+14][Phi+0].energy; out.etaStrip14.cr0.eta=14; out.etaStrip14.cr0.phi=0;
  out.etaStrip14.cr1.energy=temporary[Eta+14][Phi+1].energy; out.etaStrip14.cr1.eta=14; out.etaStrip14.cr1.phi=1;
  out.etaStrip14.cr2.energy=temporary[Eta+14][Phi+2].energy; out.etaStrip14.cr2.eta=14; out.etaStrip14.cr2.phi=2;
  out.etaStrip14.cr3.energy=temporary[Eta+14][Phi+3].energy; out.etaStrip14.cr3.eta=14; out.etaStrip14.cr3.phi=3;
  out.etaStrip14.cr4.energy=temporary[Eta+14][Phi+4].energy; out.etaStrip14.cr4.eta=14; out.etaStrip14.cr4.phi=4;
  out.etaStrip14.cr5.energy=temporary[Eta+14][Phi+5].energy; out.etaStrip14.cr5.eta=14; out.etaStrip14.cr5.phi=5;
  out.etaStrip14.cr6.energy=temporary[Eta+14][Phi+6].energy; out.etaStrip14.cr6.eta=14; out.etaStrip14.cr6.phi=6;
  out.etaStrip14.cr7.energy=temporary[Eta+14][Phi+7].energy; out.etaStrip14.cr7.eta=14; out.etaStrip14.cr7.phi=7;
  out.etaStrip14.cr8.energy=temporary[Eta+14][Phi+8].energy; out.etaStrip14.cr8.eta=14; out.etaStrip14.cr8.phi=8;
  out.etaStrip14.cr9.energy=temporary[Eta+14][Phi+9].energy; out.etaStrip14.cr9.eta=14; out.etaStrip14.cr9.phi=9;
  out.etaStrip14.cr10.energy=temporary[Eta+14][Phi+10].energy; out.etaStrip14.cr10.eta=14; out.etaStrip14.cr10.phi=10;
  out.etaStrip14.cr11.energy=temporary[Eta+14][Phi+11].energy; out.etaStrip14.cr11.eta=14; out.etaStrip14.cr11.phi=11;
  out.etaStrip14.cr12.energy=temporary[Eta+14][Phi+12].energy; out.etaStrip14.cr12.eta=14; out.etaStrip14.cr12.phi=12;
  out.etaStrip14.cr13.energy=temporary[Eta+14][Phi+13].energy; out.etaStrip14.cr13.eta=14; out.etaStrip14.cr13.phi=13;
  out.etaStrip14.cr14.energy=temporary[Eta+14][Phi+14].energy; out.etaStrip14.cr14.eta=14; out.etaStrip14.cr14.phi=14;
  out.etaStrip14.cr15.energy=temporary[Eta+14][Phi+15].energy; out.etaStrip14.cr15.eta=14; out.etaStrip14.cr15.phi=15;
  out.etaStrip14.cr16.energy=temporary[Eta+14][Phi+16].energy; out.etaStrip14.cr16.eta=14; out.etaStrip14.cr16.phi=16;
  out.etaStrip14.cr17.energy=temporary[Eta+14][Phi+17].energy; out.etaStrip14.cr17.eta=14; out.etaStrip14.cr17.phi=17;
  out.etaStrip14.cr18.energy=temporary[Eta+14][Phi+18].energy; out.etaStrip14.cr18.eta=14; out.etaStrip14.cr18.phi=18;
  out.etaStrip14.cr19.energy=temporary[Eta+14][Phi+19].energy; out.etaStrip14.cr19.eta=14; out.etaStrip14.cr19.phi=19;
  
  return out;
}

//--------------------------------------------------------// 

// Compare two ecaltp_t and return the one with the larger pT.
ecaltp_t bestOf2(const ecaltp_t& ecaltp0, const ecaltp_t& ecaltp1) {
  ecaltp_t x;
  x = (ecaltp0.energy > ecaltp1.energy)?ecaltp0:ecaltp1;

  return x;
}

//--------------------------------------------------------// 

// For a given etaStrip_t, find the ecaltp_t (out of 20 of them) with the largest pT, using pairwise comparison
ecaltp_t getPeakBin20N(const etaStrip_t& etaStrip){

  ecaltp_t best01 = bestOf2(etaStrip.cr0,etaStrip.cr1) ;
  ecaltp_t best23 = bestOf2(etaStrip.cr2,etaStrip.cr3) ;
  ecaltp_t best45 = bestOf2(etaStrip.cr4,etaStrip.cr5) ;
  ecaltp_t best67 = bestOf2(etaStrip.cr6,etaStrip.cr7) ;
  ecaltp_t best89 = bestOf2(etaStrip.cr8,etaStrip.cr9) ;
  ecaltp_t best1011 = bestOf2(etaStrip.cr10,etaStrip.cr11) ;
  ecaltp_t best1213 = bestOf2(etaStrip.cr12,etaStrip.cr13) ;
  ecaltp_t best1415 = bestOf2(etaStrip.cr14,etaStrip.cr15) ;
  ecaltp_t best1617 = bestOf2(etaStrip.cr16,etaStrip.cr17) ;
  ecaltp_t best1819 = bestOf2(etaStrip.cr18,etaStrip.cr19) ;

  ecaltp_t best0123 = bestOf2(best01,best23) ;
  ecaltp_t best4567 = bestOf2(best45,best67) ;
  ecaltp_t best891011 = bestOf2(best89,best1011) ;
  ecaltp_t best12131415 = bestOf2(best1213,best1415) ;
  ecaltp_t best16171819 = bestOf2(best1617,best1819) ;

  ecaltp_t best01234567 = bestOf2(best0123,best4567) ;
  ecaltp_t best89101112131415 = bestOf2(best891011,best12131415) ;

  ecaltp_t best0to15 = bestOf2(best01234567,best89101112131415) ;
  ecaltp_t bestOf20 = bestOf2(best0to15,best16171819) ;

  return bestOf20 ;
}

//--------------------------------------------------------//

// For a given etaStripPeak_t (representing the 15 crystals, one per row in eta, not necessarily with the same phi),
// return the crystal with the highest pT).
// This is used in getClusterPosition to get the cluster seed.

crystalMax getPeakBin15N(const etaStripPeak_t& etaStrip){
  crystalMax x;

  ecaltp_t best01 = bestOf2(etaStrip.pk0,etaStrip.pk1) ;
  ecaltp_t best23 = bestOf2(etaStrip.pk2,etaStrip.pk3) ;
  ecaltp_t best45 = bestOf2(etaStrip.pk4,etaStrip.pk5) ;
  ecaltp_t best67 = bestOf2(etaStrip.pk6,etaStrip.pk7) ;
  ecaltp_t best89 = bestOf2(etaStrip.pk8,etaStrip.pk9) ;
  ecaltp_t best1011 = bestOf2(etaStrip.pk10,etaStrip.pk11) ;
  ecaltp_t best1213 = bestOf2(etaStrip.pk12,etaStrip.pk13) ;

  ecaltp_t best0123 = bestOf2(best01,best23) ;
  ecaltp_t best4567 = bestOf2(best45,best67) ;
  ecaltp_t best891011 = bestOf2(best89,best1011) ;
  ecaltp_t best121314 = bestOf2(best1213,etaStrip.pk14) ;

  ecaltp_t best01234567 = bestOf2(best0123,best4567);
  ecaltp_t best891011121314 = bestOf2(best891011,best121314) ;

  ecaltp_t bestOf15 = bestOf2(best01234567,best891011121314) ;

  x.energy = bestOf15.energy ;
  x.etaMax = bestOf15.eta ;
  x.phiMax = bestOf15.phi ;

  return x ;
}

//--------------------------------------------------------// 

// Take a 3x4 ECAL region (i.e. 15x20 in crystals, add crystal energies in squares of 5x5, giving
// 3x4 = 12 ECAL tower sums.) Store these 12 values in towerEt.

void getECALTowersEt(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<12> towerEt[12]){

    ap_uint<10> temp[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI] ;
    ap_uint<12> towerEtN[3][4][5] ;
    for (int i = 0; i < CRYSTAL_IN_ETA; i++){
      for(int k = 0; k < CRYSTAL_IN_PHI; k++){
	temp[i][k] = tempX[i][k].energy ;
      }}
    
    for(int i=0; i<CRYSTAL_IN_ETA; i=i+5){
      for(int k=0; k<CRYSTAL_IN_PHI; k=k+5){
	towerEtN[i/5][k/5][0] = temp[i][k] + temp[i][k+1] + temp[i][k+2] + temp[i][k+3] + temp[i][k+4] ;
	towerEtN[i/5][k/5][1] = temp[i+1][k] + temp[i+1][k+1] + temp[i+1][k+2] + temp[i+1][k+3] + temp[i+1][k+4] ;
	towerEtN[i/5][k/5][2] = temp[i+2][k] + temp[i+2][k+1] + temp[i+2][k+2] + temp[i+2][k+3] + temp[i+2][k+4] ;
	towerEtN[i/5][k/5][3] = temp[i+3][k] + temp[i+3][k+1] + temp[i+3][k+2] + temp[i+3][k+3] + temp[i+3][k+4] ;
	towerEtN[i/5][k/5][4] = temp[i+4][k] + temp[i+4][k+1] + temp[i+4][k+2] + temp[i+4][k+3] + temp[i+4][k+4] ;
      }}

    towerEt[0]= towerEtN[0][0][0] + towerEtN[0][0][1] + towerEtN[0][0][2] + towerEtN[0][0][3] + towerEtN[0][0][4] ;
    towerEt[1]= towerEtN[0][1][0] + towerEtN[0][1][1] + towerEtN[0][1][2] + towerEtN[0][1][3] + towerEtN[0][1][4] ;
    towerEt[2]= towerEtN[0][2][0] + towerEtN[0][2][1] + towerEtN[0][2][2] + towerEtN[0][2][3] + towerEtN[0][2][4] ;
    towerEt[3]= towerEtN[0][3][0] + towerEtN[0][3][1] + towerEtN[0][3][2] + towerEtN[0][3][3] + towerEtN[0][3][4] ;
    towerEt[4]= towerEtN[1][0][0] + towerEtN[1][0][1] + towerEtN[1][0][2] + towerEtN[1][0][3] + towerEtN[1][0][4] ;
    towerEt[5]= towerEtN[1][1][0] + towerEtN[1][1][1] + towerEtN[1][1][2] + towerEtN[1][1][3] + towerEtN[1][1][4] ;
    towerEt[6]= towerEtN[1][2][0] + towerEtN[1][2][1] + towerEtN[1][2][2] + towerEtN[1][2][3] + towerEtN[1][2][4] ;
    towerEt[7]= towerEtN[1][3][0] + towerEtN[1][3][1] + towerEtN[1][3][2] + towerEtN[1][3][3] + towerEtN[1][3][4] ;
    towerEt[8]= towerEtN[2][0][0] + towerEtN[2][0][1] + towerEtN[2][0][2] + towerEtN[2][0][3] + towerEtN[2][0][4] ;
    towerEt[9]= towerEtN[2][1][0] + towerEtN[2][1][1] + towerEtN[2][1][2] + towerEtN[2][1][3] + towerEtN[2][1][4] ;
    towerEt[10]= towerEtN[2][2][0] + towerEtN[2][2][1] + towerEtN[2][2][2] + towerEtN[2][2][3] + towerEtN[2][2][4] ;
    towerEt[11]= towerEtN[2][3][0] + towerEtN[2][3][1] + towerEtN[2][3][2] + towerEtN[2][3][3] + towerEtN[2][3][4] ;

    std::cout << "getECALTowersEt (after subtracting clusters): ";
    for (int i = 0; i < 12; i++) {
      std::cout << towerEt[i] << " ";
    }
    std::cout << std::endl;
    
}

//--------------------------------------------------------//  

class clusterInfo{
public:
  ap_uint<10> seedEnergy;
  ap_uint<15> energy;
  ap_uint<15> et5x5;
  ap_uint<15> et2x5;
  ap_uint<5> phiMax;
  ap_uint<5> etaMax;
  ap_uint<2> brems;

  clusterInfo(){
    seedEnergy = 0;
    energy = 0;
    et5x5 = 0;
    et2x5 = 0;
    phiMax = 0;
    etaMax = 0;
    brems = 0;
  }

  clusterInfo& operator=(const clusterInfo& rhs){
    seedEnergy = rhs.seedEnergy;
    energy = rhs.energy;
    et5x5 = rhs.et5x5;
    et2x5 = rhs.et2x5;
    phiMax = rhs.phiMax;
    etaMax = rhs.etaMax;
    brems = rhs.brems;
    return *this;
  }
};

//--------------------------------------------------------//

class Cluster{
public:
  ap_uint<28> data;

  Cluster(){
    data = 0;
  }

  Cluster& operator=(const Cluster& rhs){
    data = rhs.data;
    return *this;
  }

  Cluster(ap_uint<12> clusterEnergy, ap_uint<5> towerEta, ap_uint<2> towerPhi, ap_uint<3> clusterEta, ap_uint<3> clusterPhi, ap_uint<3> satur){
    data = (clusterEnergy) | 
      (((ap_uint<32>) towerEta)  << 12) | 
      (((ap_uint<32>) towerPhi)  << 17) | 
      (((ap_uint<32>) clusterEta)  << 19) | 
      (((ap_uint<32>) clusterPhi) << 22) | 
      (((ap_uint<32>) satur)       << 25);
  }

  ap_uint<12> clusterEnergy() {return (data & 0xFFF);}
  ap_uint<5> towerEta() {return ((data >> 12) & 0x37);}
  ap_uint<2> towerPhi() {return ((data >> 17) & 0x3);}
  ap_uint<3> clusterEta() {return ((data >> 19) & 0x7);}
  ap_uint<3> clusterPhi() {return ((data >> 22) & 0x7);}
  ap_uint<3> satur() {return ((data >> 25) & 0x7);}
  
  operator uint32_t() {return (ap_uint<28>) data;}

};

//--------------------------------------------------------// 

// Compare the ET of two clusters (pass this to std::sort to get clusters sorted in decreasing ET).

bool compareClusterET(Cluster& lhs, Cluster& rhs) {
  return ( lhs.clusterEnergy() > rhs.clusterEnergy() );
}

//--------------------------------------------------------//

clusterInfo getClusterPosition(const ecalRegion_t& ecalRegion){
  etaStripPeak_t etaStripPeak;
  clusterInfo cluster;

  etaStripPeak.pk0 = getPeakBin20N(ecalRegion.etaStrip0);
  etaStripPeak.pk1 = getPeakBin20N(ecalRegion.etaStrip1);
  etaStripPeak.pk2 = getPeakBin20N(ecalRegion.etaStrip2);
  etaStripPeak.pk3 = getPeakBin20N(ecalRegion.etaStrip3);
  etaStripPeak.pk4 = getPeakBin20N(ecalRegion.etaStrip4);
  etaStripPeak.pk5 = getPeakBin20N(ecalRegion.etaStrip5);
  etaStripPeak.pk6 = getPeakBin20N(ecalRegion.etaStrip6);
  etaStripPeak.pk7 = getPeakBin20N(ecalRegion.etaStrip7);
  etaStripPeak.pk8 = getPeakBin20N(ecalRegion.etaStrip8);
  etaStripPeak.pk9 = getPeakBin20N(ecalRegion.etaStrip9);
  etaStripPeak.pk10 = getPeakBin20N(ecalRegion.etaStrip10);
  etaStripPeak.pk11 = getPeakBin20N(ecalRegion.etaStrip11);
  etaStripPeak.pk12 = getPeakBin20N(ecalRegion.etaStrip12);
  etaStripPeak.pk13 = getPeakBin20N(ecalRegion.etaStrip13);
  etaStripPeak.pk14 = getPeakBin20N(ecalRegion.etaStrip14);

  crystalMax peakIn15;
  peakIn15 = getPeakBin15N(etaStripPeak);
  
  cluster.seedEnergy = peakIn15.energy;
  cluster.energy = 0;
  cluster.etaMax = peakIn15.etaMax;
  cluster.phiMax = peakIn15.phiMax;
  cluster.brems = 0;
  cluster.et5x5 = 0;
  cluster.et2x5 = 0;

  return cluster;
}

//--------------------------------------------------------//

Cluster packCluster(ap_uint<15>& clusterEt, ap_uint<5>& etaMax_t, ap_uint<5>& phiMax_t){
  

  //  std::cout << "[-->] packCluster: clusterEt " << clusterEt << ", eta and phi: " << etaMax_t << ", " << phiMax_t << std::endl;
  
  ap_uint<12> peggedEt;
  Cluster pack;

  ap_uint<5> towerEta = (etaMax_t)/5;
  ap_uint<2> towerPhi = (phiMax_t)/5;
  ap_uint<3> clusterEta = etaMax_t - 5*towerEta;
  ap_uint<3> clusterPhi = phiMax_t - 5*towerPhi;

  peggedEt = (clusterEt > 0xFFF)? (ap_uint<12>)0xFFF : (ap_uint<12>) clusterEt;

  
  pack = Cluster(peggedEt, towerEta, towerPhi, clusterEta, clusterPhi, 0);
  std::cout << "[-->] packCluster: clusterEt " << clusterEt << ", clusterEta and Phi: " << clusterEta << ", "
	    << clusterPhi << std::endl;

  return pack;
}


//--------------------------------------------------------// 

// Given the cluster seed_eta, seed_phi, and brems, remove the cluster energy
// from the given crystal array temp. Functionally identical to "RemoveTmp".

void removeClusterFromCrystal(crystal temp[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi, ap_uint<2> brems) {

  // Zero out the crystal energies in a 3 (eta) by 5 (phi) window (the clusters are 3x5 in crystals)
  for (int i = 0; i < CRYSTAL_IN_ETA; i++){
    for (int k = 0; k < CRYSTAL_IN_PHI; k++){
      if ((i >= (seed_eta-1)) && (i <= (seed_eta+1)) && (k >= (seed_phi-2)) && (k <= (seed_phi+2)))  temp[i][k].energy = 0;
    }
  }

  // If brems flag is 1, *also* zero the energies in the 3x5 window to the "left" of the cluster 
  // N.B. in the positive eta cards, "left" in the region = towards negative phi, 
  // but for negative eta cards, everything is flipped, so "left" in the region" = towards positive phi
  if (brems == 1) {
    for (int i = 0; i < CRYSTAL_IN_ETA; i++){
      for (int k = 0; k < CRYSTAL_IN_PHI; k++){
	if ((i >= (seed_eta-1)) && (i <= (seed_eta+1)) && (k >= (seed_phi-2-5)) && (k <= (seed_phi+2-5)))  temp[i][k].energy = 0;
      }
    }
  }
  
  // If brems flag is 2, *also* zero the energies in the 3x5 window to the "right" of the cluster
  // N.B. in the positive eta cards, "right" in the region = towards POSITIVE phi,
  // but for negative eta cards, everything is flipped, so "right" in the region = towards NEGATIVE phi
  if (brems == 2) {
    for (int i = 0; i < CRYSTAL_IN_ETA; i++){
      for (int k = 0; k < CRYSTAL_IN_PHI; k++){
	if ((i >= (seed_eta-1)) && (i <= (seed_eta+1)) && (k >= (seed_phi-2+5)) && (k <= (seed_phi+2+5)))  temp[i][k].energy = 0;
      }
    }
  }
  
}

//--------------------------------------------------------//

// Given a 15x20 crystal tempX, and a seed with seed_eta and seed_phi, return a clusterInfo containing
// the cluster energy for a positive bremmstrahulung shift 

clusterInfo getBremsValuesPos(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi ){

  ap_uint<12> temp[CRYSTAL_IN_ETA+2][CRYSTAL_IN_PHI+4];
  ap_uint<12> phi0eta[3], phi1eta[3], phi2eta[3], phi3eta[3], phi4eta[3];
  ap_uint<12> eta_slice[3];
  clusterInfo cluster_tmp;
  
  // Set all entries in a new ((15+2)x(20+4)) array to be zero.
  for (int i = 0; i < (CRYSTAL_IN_ETA + 2); i++) {
    for (int j = 0; j < (CRYSTAL_IN_PHI + 4); j++) {
      temp[i][j] = 0;
    }}
  
  // Read the energies of the input crystal tempX into the slightly larger array temp, with an offset so temp is tempX
  // except shifted +1 in eta, and -3 in phi. 
  for (int i = 0; i < (CRYSTAL_IN_ETA); i++) {
    for (int j = 0; j < (CRYSTAL_IN_PHI); j++) {
      temp[i+1][j-3] = tempX[i][j].energy;   
    }}
  
  ap_uint<6> seed_eta1, seed_phi1;
  seed_eta1 = seed_eta; //to start from corner
  seed_phi1 = seed_phi; //to start from corner

  // now we are in the left bottom corner 
  // Loop over the shifted array, and at the original location of the seed (seed_eta1/seed_phi1), 
  // read a 3 (eta) x 5 (phi) rectangle of crystals where the original location of the seed is in the bottom left corner
  for (int j = 0; j < CRYSTAL_IN_ETA; j++) {
    if (j == seed_eta1) {
      for (int k = 0; k < CRYSTAL_IN_PHI; k++) {
	if (k == seed_phi1) {
	  // Same eta as the seed, read next five crystals in phi
	  phi0eta[0] = temp[j][k];
	  phi1eta[0] = temp[j][k+1];
	  phi2eta[0] = temp[j][k+2];
	  phi3eta[0] = temp[j][k+3];
	  phi4eta[0] = temp[j][k+4];
	  
	  // +1 eta from the seed, read next five crystals in phi
	  phi0eta[1] = temp[j+1][k];
	  phi1eta[1] = temp[j+1][k+1];
	  phi2eta[1] = temp[j+1][k+2];
	  phi3eta[1] = temp[j+1][k+3];
	  phi4eta[1] = temp[j+1][k+4];
                        
	  // +2 eta from the seed, read next five crystals in phi
	  phi0eta[2] = temp[j+2][k];
	  phi1eta[2] = temp[j+2][k+1];
	  phi2eta[2] = temp[j+2][k+2];
	  phi3eta[2] = temp[j+2][k+3];
	  phi4eta[2] = temp[j+2][k+4];
	  
	  continue;
	}}
    }}

  // Add up the energies in this 3x5 of crystals, initialize a cluster_tmp, and return it
  for (int i = 0; i < 3; i++) { 
    eta_slice[i] = phi0eta[i] + phi1eta[i] + phi2eta[i] + phi3eta[i] + phi4eta[i];
  }
  cluster_tmp.energy = (eta_slice[0] + eta_slice[1] + eta_slice[2]);

  //  std::cout << "getBremsValuesPos: energy, seed eta/phi = " << cluster_tmp.energy << ", " << seed_eta << ", " << seed_phi << std::endl;
  return cluster_tmp;

}

//--------------------------------------------------------//  

// Given a 15x20 crystal tempX, and a seed with seed_eta and seed_phi, return a clusterInfo containing                                     
// the cluster energy for a *negative* bremmstrahlung shift 

clusterInfo getBremsValuesNeg(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi){

  ap_uint<12> temp[CRYSTAL_IN_ETA+2][CRYSTAL_IN_PHI+4];
  ap_uint<12> phi0eta[3], phi1eta[3], phi2eta[3], phi3eta[3], phi4eta[3];

  ap_uint<12> eta_slice[3];

  clusterInfo cluster_tmp;

  // Initialize all entries in a new ((15+2)x(20+4)) array to be zero.                                                                  
  for (int i = 0; i < (CRYSTAL_IN_ETA + 2); i++) {
    for (int j = 0; j < (CRYSTAL_IN_PHI + 4); j++) {
      temp[i][j] = 0;
    }}

  // Read the energies of the input crystal tempX into the slightly larger array temp, with an offset so temp is tempX
  // except shifted in +1 in eta and +7 in phi
  for (int i = 0; i < (CRYSTAL_IN_ETA); i++) {
    for (int j = 0; j < (CRYSTAL_IN_PHI - 5); j++) {
      temp[i+1][j+7] = tempX[i][j].energy;
    }}

  ap_uint<6> seed_eta1, seed_phi1;
  seed_eta1 = seed_eta; //to start from corner
  seed_phi1 = seed_phi; //to start from corner

  // Loop over the shifted array, and at the original location of the seed (seed_eta1/seed_phi1), 
  // read a 3 (eta) x 5 (phi) rectangle of crystals where the original location of the seed is in the bottom left corner
  for (int j = 0; j < CRYSTAL_IN_ETA; j++) {
    if (j == seed_eta1) {
      for (int k = 0; k < CRYSTAL_IN_PHI; k++) {
  	if (k == seed_phi1) {
  	  // Same eta as the seed, read next five crystals in phi
  	  phi0eta[0] = temp[j][k];
  	  phi1eta[0] = temp[j][k+1];
  	  phi2eta[0] = temp[j][k+2];
  	  phi3eta[0] = temp[j][k+3];
  	  phi4eta[0] = temp[j][k+4];
	    
  	  // +1 eta from the seed, read next five crystals in phi
  	  phi0eta[1] = temp[j+1][k];
  	  phi1eta[1] = temp[j+1][k+1];
  	  phi2eta[1] = temp[j+1][k+2];
  	  phi3eta[1] = temp[j+1][k+3];
  	  phi4eta[1] = temp[j+1][k+4];
                        
  	  // +2 eta from the seed, read next five crystals in phi
  	  phi0eta[2] = temp[j+2][k];
  	  phi1eta[2] = temp[j+2][k+1];
  	  phi2eta[2] = temp[j+2][k+2];
  	  phi3eta[2] = temp[j+2][k+3];
  	  phi4eta[2] = temp[j+2][k+4];
	    
  	  continue;
  	}}
    }}
  
  // Add up the energies in this 3x5 of crystals, initialize a cluster_tmp, and return it
  for (int i = 0; i < 3; i++) { 
    eta_slice[i] = phi0eta[i] + phi1eta[i] + phi2eta[i] + phi3eta[i] + phi4eta[i];
  }
  cluster_tmp.energy = (eta_slice[0] + eta_slice[1] + eta_slice[2]);
  
  //  std::cout << "getBremsValuesNeg: energy, seed eta/phi = " << cluster_tmp.energy << ", " << seed_eta << ", " << seed_phi << std::endl;
  
  return cluster_tmp;
  
}


//--------------------------------------------------------//                                                                           

// Given a 15x20 crystal tempX, and a seed with seed_eta and seed_phi, return a clusterInfo containing                              
// the cluster energy (central value)

clusterInfo getClusterValues(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi){

  ap_uint<12> temp[CRYSTAL_IN_ETA+4][CRYSTAL_IN_PHI+4];
  ap_uint<12> phi0eta[5], phi1eta[5], phi2eta[5], phi3eta[5], phi4eta[5];
  ap_uint<12> eta_slice[5];
  ap_uint<12> et2x5_1Tot, et2x5_2Tot, etSum2x5;
  ap_uint<12> et5x5Tot;

  clusterInfo cluster_tmp;
  // Initialize empty (15+4)x(20+4) array
  for (int i = 0; i < (CRYSTAL_IN_ETA + 4); i++){
    for (int k = 0; k< (CRYSTAL_IN_PHI + 4); k++){
      temp[i][k] = 0 ;
    }}

  // Copy input array energies into temp array with +2 eta and +2 phi offset.
  for (int i = 0; i < (CRYSTAL_IN_ETA); i++){
    for (int k = 0; k < (CRYSTAL_IN_PHI); k++){
      temp[i+2][k+2] = tempX[i][k].energy ;
    }}

  ap_uint<6> seed_eta1, seed_phi1;
  seed_eta1 = seed_eta; //to start from corner
  seed_phi1 = seed_phi; //to start from corner
  
  // now we are in the left bottom corner 
  // Loop over the shifted array, and at the original location of the seed (seed_eta1/seed_phi1), 
  // read a 5 (eta) x 5 (phi) rectangle of crystals where the original location of the seed is in the bottom left corner
  for (int j = 0; j < CRYSTAL_IN_ETA; j++) {
    if (j == seed_eta1) {
      for (int k = 0; k < CRYSTAL_IN_PHI; k++) {
	if (k == seed_phi1) {
	  // Same eta as the seed, read next five crystals in phi
	  phi0eta[0] = temp[j][k];
	  phi1eta[0] = temp[j][k+1];
	  phi2eta[0] = temp[j][k+2];
	  phi3eta[0] = temp[j][k+3];
	  phi4eta[0] = temp[j][k+4];
	    
	  // +1 eta from the seed, read next five crystals in phi
	  phi0eta[1] = temp[j+1][k];
	  phi1eta[1] = temp[j+1][k+1];
	  phi2eta[1] = temp[j+1][k+2];
	  phi3eta[1] = temp[j+1][k+3];
	  phi4eta[1] = temp[j+1][k+4];
                        
	  // +2 eta from the seed, read next five crystals in phi
	  phi0eta[2] = temp[j+2][k];
	  phi1eta[2] = temp[j+2][k+1];
	  phi2eta[2] = temp[j+2][k+2];
	  phi3eta[2] = temp[j+2][k+3];
	  phi4eta[2] = temp[j+2][k+4];

	  // +3 eta from the seed, read next five crystals in phi 
	  phi0eta[3] = temp[j+3][k];
	  phi1eta[3] = temp[j+3][k+1];
	  phi2eta[3] = temp[j+3][k+2];
	  phi3eta[3] = temp[j+3][k+3];
	  phi4eta[3] = temp[j+3][k+4];

	  // +4 eta from the seed, read next five crystals in phi
	  phi0eta[4] = temp[j+4][k];
	  phi1eta[4] = temp[j+4][k+1];
	  phi2eta[4] = temp[j+4][k+2];
	  phi3eta[4] = temp[j+4][k+3];
	  phi4eta[4] = temp[j+4][k+4];
	    
	  continue;
	}}
    }}
  
  // Add the first three eta strips into the cluster energy
  for (int i = 0; i < 5; i++) {
    eta_slice[i] = phi0eta[i] + phi1eta[i] + phi2eta[i] + phi3eta[i] + phi4eta[i];
  }
  cluster_tmp.energy = (eta_slice[1] + eta_slice[2] + eta_slice[3]);

  // Get the energy totals in the 5x5 and also in two 2x5 
  et5x5Tot   = (eta_slice[0] + eta_slice[1] + eta_slice[2] + eta_slice[3] + eta_slice[4]);
  et2x5_1Tot = (eta_slice[1] + eta_slice[2]);
  et2x5_2Tot = (eta_slice[2] + eta_slice[3]);

  if (et2x5_1Tot >= et2x5_2Tot) etSum2x5 = et2x5_1Tot;
  else                          etSum2x5 = et2x5_2Tot;

  cluster_tmp.et5x5 = et5x5Tot;
  cluster_tmp.et2x5 = etSum2x5;

  std::cout << "getClusterValues: energy, et5x5Tot, etSum2x5 = " << cluster_tmp.energy << ", " 
	    << cluster_tmp.et5x5 << ", " << cluster_tmp.et2x5 << std::endl;

  return cluster_tmp;
}

//--------------------------------------------------------// 

// Functionally identical to "getRegion3x4" from algo_top.cc (Renamed to avoid confusion with the card class method)
// In 15x20 crystal array temp, return the next cluster, and remove the cluster's energy
// from the crystal array.

Cluster getClusterFromRegion3x4(crystal temp[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI]){

  Cluster returnCluster;
  clusterInfo cluster_tmp;
  clusterInfo cluster_tmpCenter;
  clusterInfo cluster_tmpBneg;
  clusterInfo cluster_tmpBpos;

  ecalRegion_t ecalRegion;
  ecalRegion = initStructure(temp);

  cluster_tmp = getClusterPosition(ecalRegion);

  ap_uint<5> seed_phi = cluster_tmp.phiMax;
  ap_uint<5> seed_eta = cluster_tmp.etaMax;
  
  cluster_tmpCenter = getClusterValues(temp, seed_eta, seed_phi);
  cluster_tmpBneg   = getBremsValuesNeg(temp, seed_eta, seed_phi); 
  cluster_tmpBpos   = getBremsValuesPos(temp, seed_eta, seed_phi);

  cluster_tmp.energy = cluster_tmpCenter.energy;
  cluster_tmp.brems = 0;

  if ((cluster_tmpBneg.energy > cluster_tmpCenter.energy/8) && (cluster_tmpBneg.energy > cluster_tmpBpos.energy)) {    
    cluster_tmp.energy = (cluster_tmpCenter.energy + cluster_tmpBneg.energy);
    std::cout << "getClusterFromRegion3x4: Brems in negative phi direction: set cluster ET to " << cluster_tmp.energy << std::endl;
    cluster_tmp.brems = 1; }
  else if(cluster_tmpBpos.energy > cluster_tmpCenter.energy/8) {
    cluster_tmp.energy = (cluster_tmpCenter.energy + cluster_tmpBpos.energy);
    std::cout << "getClusterFromRegion3x4: Brems in positive phi direction: set cluster ET to " << cluster_tmp.energy << std::endl;
    cluster_tmp.brems = 2; }

  returnCluster = packCluster(cluster_tmp.energy, cluster_tmp.etaMax, cluster_tmp.phiMax);
  
  removeClusterFromCrystal(temp, seed_eta, seed_phi, cluster_tmp.brems);

  return returnCluster;
  
}

//--------------------------------------------------------// 

//////////////////////////////////////////////////////////////////////////  

// EGammaCrystalsProducer initializer, destructor, and produce methods

EGammaCrystalsProducer::EGammaCrystalsProducer(const edm::ParameterSet & iConfig)
  : ecalTPEBToken_(consumes<EcalEBTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("ecalTPEB"))),
    hcalTPToken_(
		 consumes<edm::SortedCollection<HcalTriggerPrimitiveDigi> >(iConfig.getParameter<edm::InputTag>("hcalTP"))),
  calib_(iConfig.getParameter<edm::ParameterSet>("calib")) {
}

EGammaCrystalsProducer::~EGammaCrystalsProducer() {}



void EGammaCrystalsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Detector geometry
  iSetup.get<CaloGeometryRecord>().get(caloGeometry_);
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  hbGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
  iSetup.get<HcalRecNumberingRecord>().get(hbTopology);
  hcTopology_ = hbTopology.product();
  HcalTrigTowerGeometry theTrigTowerGeometry(hcTopology_);
  iSetup.get<CaloTPGRecord>().get(decoder_);

  /////////////////////////////////////////////////////
  // Get all the ECAL hits
  /////////////////////////////////////////////////////
  edm::Handle<EcalEBTrigPrimDigiCollection> pcalohits;   // from l.354
  iEvent.getByToken(ecalTPEBToken_, pcalohits);

  std::vector<SimpleCaloHit> ecalhits;

  for (const auto& hit : *pcalohits.product()) {
    if (hit.encodedEt() > 0)  // hit.encodedEt() returns an int corresponding to 2x the crystal Et
      {
	// Et is 10 bit, by keeping the ADC saturation Et at 120 GeV it means that you have to divide by 8
	float et = hit.encodedEt() / 8.;
	if (et < cut_500_MeV)
	  continue;  // keep the 500 MeV ET Cut

	// std::cout << "ECAL hit et: " << et << std::endl;

	// Get cell coordinates and info
	auto cell = ebGeometry->getGeometry(hit.id());
	
	// std::cout << "Found ECAL cell/hit with coordinates " << cell->getPosition().x() << "," 
	// 	  << cell->getPosition().y() << "," 
	// 	  << cell->getPosition().z() << " and ET (GeV) " 
	// 	  << et << std::endl;

	SimpleCaloHit ehit;
	ehit.setId(hit.id());
	ehit.setPosition(GlobalVector(cell->getPosition().x(), cell->getPosition().y(), cell->getPosition().z()));
	ehit.setEnergy(et);
	ehit.setEt_uint((ap_uint<10>) hit.encodedEt());  // also save the 10-bit Et
	ehit.setPt();
	ecalhits.push_back(ehit);
	
      }
  }

  // Get all the HCAL hits
  std::vector<SimpleCaloHit> hcalhits;
  edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hbhecoll;
  iEvent.getByToken(hcalTPToken_, hbhecoll);
  
  for (const auto& hit : *hbhecoll.product()) {
    float et = decoder_->hcaletValue(hit.id(), hit.t0());
    ap_uint<10> encodedEt = hit.t0().compressedEt(); 
    // same thing as SOI_compressedEt() in HcalTriggerPrimitiveDigi.h///
    if (et <= 0)
      continue;
    
    if (!(hcTopology_->validHT(hit.id()))) {
      LogError("EGammaCrystalsProducer")
  	<< " -- Hcal hit DetID not present in HCAL Geom: " << hit.id() << std::endl;
      throw cms::Exception("EGammaCrystalsProducer");
      continue;
    }
    const std::vector<HcalDetId>& hcId = theTrigTowerGeometry.detIds(hit.id());
    if (hcId.empty()) {
      LogError("EGammaCrystalsProducer")
  	<< "Cannot find any HCalDetId corresponding to " << hit.id() << std::endl;
      throw cms::Exception("EGammaCrystalsProducer");
      continue;
    }
    if (hcId[0].subdetId() > 1)
      continue;
    GlobalVector hcal_tp_position = GlobalVector(0., 0., 0.);
    for (const auto& hcId_i : hcId) {
      if (hcId_i.subdetId() > 1)
        continue;
      // get the first HCAL TP/ cell
      auto cell = hbGeometry->getGeometry(hcId_i);
      if (cell == nullptr)
  	continue;
      GlobalVector tmpVector = GlobalVector(cell->getPosition().x(), cell->getPosition().y(), cell->getPosition().z());
      hcal_tp_position = tmpVector;
      
      // std::cout << "Found HCAL cell/TP with coordinates " << cell->getPosition().x() << ","
      //  		<< cell->getPosition().y() << ","
      //  		<< cell->getPosition().z() << " and ET (GeV) " << et
      // 		<< ", encoded Et " << encodedEt << std::endl;
      
      break;
    }
    SimpleCaloHit hhit;
    hhit.setId(hit.id());
    hhit.setIdHcal(hit.id());
    hhit.setPosition(hcal_tp_position);
    hhit.setEnergy(et);
    hhit.setPt();
    hhit.setEt_uint(encodedEt);
    hcalhits.push_back(hhit);
  }

  //*******************************************************************
  //*************** Declare Layer 1 outputs ***************************
  //*******************************************************************

  // Definition of L1 outputs
  // 36 L1 cards, each with 4x17 towers. All using CMSSW indexing convention, NOT firmware convention
  ap_uint<12> ECAL_tower_L1Card[n_links_card][n_towers_per_link][n_towers_halfPhi];
  ap_uint<12> HCAL_tower_L1Card[n_links_card][n_towers_per_link][n_towers_halfPhi];
  int iEta_tower_L1Card[n_links_card][n_towers_per_link][n_towers_halfPhi];  
  int iPhi_tower_L1Card[n_links_card][n_towers_per_link][n_towers_halfPhi];
  // 36 L1 cards send each 4 links with 3 clusters
  int crystalID_cluster_L1Card[n_links_card][n_clusters_link][n_towers_halfPhi];

  // Zero out the L1 outputs
  for (int ii = 0; ii < n_links_card; ++ii) {
    for (int jj = 0; jj < n_towers_per_link; ++jj) {
      for (int ll = 0; ll < n_towers_halfPhi; ++ll) {
        ECAL_tower_L1Card[ii][jj][ll] = 0;
        HCAL_tower_L1Card[ii][jj][ll] = 0;
        iPhi_tower_L1Card[ii][jj][ll] = -999;
        iEta_tower_L1Card[ii][jj][ll] = -999;
      }
    }
  }

  for (int ii = 0; ii < n_links_card; ++ii) {
    for (int jj = 0; jj < n_clusters_link; ++jj) {
      for (int ll = 0; ll < n_towers_halfPhi; ++ll) {
        crystalID_cluster_L1Card[ii][jj][ll] = 0;
      }
    }
  }
  //*******************************************************************
  //*************** Do RCT geometry (ECAL)  ***************************
  //*******************************************************************

  
  for (int cc = 0; cc < n_towers_halfPhi; ++cc) {  // Loop over 36 L1 cards
          
    // Initialize variables
    card rctCard;
    rctCard.setIdx(cc);

    for (const auto& hit : ecalhits) {
      // Check if the hit is in cards 0-35
      if ((getCrystal_iPhi(hit.position().phi()) <= getCard_iPhiMax(cc)) &&
  	  (getCrystal_iPhi(hit.position().phi()) >= getCard_iPhiMin(cc)) &&
  	  (getCrystal_iEta(hit.position().eta()) <= getCard_iEtaMax(cc)) &&
  	  (getCrystal_iEta(hit.position().eta()) >= getCard_iEtaMin(cc))) {
	
	// Get the crystal eta and phi, relative to the bottom left corner of the card 
	// (0 up to 17*5, 0 up to 4*5) 
	int local_iEta = getCrystal_local_iEta(hit.position().eta(), cc);
	int local_iPhi = getCrystal_local_iPhi(hit.position().phi(), cc);
	
	// Once we have the iEta and iPhi of the crystal relative to the bottom left corner,
	// everything else is determined.
	
	// Figure out what region (0-5) the hit falls into 
	// The region number (0-5) depends only on the local crystal iEta
	int regionNumber = getRegionNumber(local_iEta);
	
	// Get the tower eta and phi index inside the card (17x4)
	int inCard_tower_iEta = int(local_iEta / CRYSTALS_IN_TOWER_ETA); 
	int inCard_tower_iPhi = int(local_iPhi / CRYSTALS_IN_TOWER_PHI);
	
	// Get the tower eta and phi index inside the region (3x4)
	int inRegion_tower_iEta = inCard_tower_iEta % TOWER_IN_ETA;
	int inRegion_tower_iPhi = inCard_tower_iPhi % TOWER_IN_PHI;
	
	// Within the region, figure out which crystal and link the hit falls into                                         
        // Get the crystal eta and phi index inside the 3x4 region (15x20)
	int inRegion_crystal_iEta = local_iEta % (TOWER_IN_ETA * CRYSTALS_IN_TOWER_ETA);
	int inRegion_crystal_iPhi = local_iPhi;

	// Get the crystal eta and phi index inside the link (5x5). E.g. crystal (7, 9) is crystal number (2, 4) in its link
	int inLink_crystal_iEta = (inRegion_crystal_iEta % CRYSTALS_IN_TOWER_ETA);
	int inLink_crystal_iPhi = (inRegion_crystal_iPhi % CRYSTALS_IN_TOWER_PHI);

	// std::cout << "local_iEta, local_iPhi, regionNumber: "
	// 	  << local_iEta << ", " 
	// 	  << local_iPhi << ", "
	// 	  << regionNumber << std::endl;
	// std::cout << "inRegion_crystal_iEta, inRegion_crystal_iPhi (expecting 15x20), regionNumber (expecting 0-5): "
	// 	  << inRegion_crystal_iEta << ", "
	// 	  << inRegion_crystal_iPhi << ", "
	// 	  << regionNumber << std::endl;

	// Access the right region -> link -> crystal and increment the energy
	if (regionNumber < N_REGIONS_PER_CARD) {
	  region3x4& myRegion = rctCard.getRegion3x4(regionNumber);

	  //	  std::cout << "inRegion_tower_iEta, inRegion_tower_iPhi (expecting 3x4): " << inRegion_tower_iEta << ", "
	  //		    << inRegion_tower_iPhi << std::endl;
	  
	  // Get the right link
	  linkECAL& myLink = myRegion.getLinkECAL(inRegion_tower_iEta, inRegion_tower_iPhi);
	  
	  // Add the energy to the right crystal 
	  //	  std::cout << "inLink_crystal_iEta, inLink_crystal_iPhi (expecting 5x5): " 
	  //		    << inLink_crystal_iEta << ", " << inLink_crystal_iPhi << std::endl;
	  myLink.addCrystalE(inLink_crystal_iEta, inLink_crystal_iPhi, hit.et_uint());
	  
	}
	
	std::cout << "Card: " << cc << ", hit (Eta, phi, et): "
		  << hit.position().eta() << ", " << hit.position().phi() << ", " << hit.et_uint() << ", "
	  // << "local_iEta/iPhi: " << local_iEta << ", " << local_iPhi << ", " 
	  //	  << "inCard_tower_iEta/iPhi: " << inCard_tower_iEta     << ", " << inCard_tower_iPhi << ", "
		  << "region/inRegion_tower_iEta/iPhi: " << regionNumber << ", " << inRegion_tower_iEta << ", "
		  << inRegion_tower_iPhi << ", "
		  << "inLink crystal iEta/iPhi: " << inLink_crystal_iEta << ", " << inLink_crystal_iPhi << std::endl;
      }      
    } // end of loop over ECAL hits

    // Also initialize tower (iEta, iPhi) coordinates (code lifted from 
    // https://github.com/cms-l1t-offline/cmssw/blob/25a1610b718c4cf94c33afb6e23767b5d3a677d7/L1Trigger/L1CaloTrigger/plugins/L1EGammaCrystalsEmulatorProducer.cc#L717-L727), changing the const variable names and getTowers_absEtaID function name
    // n.b. iEta_tower_L1Card and iPhi_tower_L1Card (the CMSSW emulator Layer 1 outputs) use a different convention
    // so later on when we compute a 17x4 tower array, it is rotated for the negative eta cards.
    static constexpr float tower_width = 0.0873;
    for (int jj = 0; jj < n_links_card ; ++jj) {
      for (int ii = 0; ii < n_towers_per_link; ++ii) {
	float phi = getPhiMin_card_emulator(cc) * tower_width / CRYSTALS_IN_TOWER_PHI - M_PI + (jj + 0.5) * tower_width;
	float eta = getEtaMin_card_emulator(cc) * tower_width / CRYSTALS_IN_TOWER_ETA - n_towers_cardEta * tower_width +
	  (ii + 0.5) * tower_width;
	iEta_tower_L1Card[jj][ii][cc] = getTower_absEtaID(eta);
	iPhi_tower_L1Card[jj][ii][cc] = getTower_absPhiID(phi);
	std::cout << "Real (eta, phi): " << eta << ", " << phi << "; " 
		  << "iEta_tower_L1Card[" << jj << "][" << ii 
		  << "][" << cc << "] = " << getTower_absEtaID(eta) << ", "
		  << "iPhi_tower_L1Card[" << jj << "][" << ii
		  << "][" << cc << "] = " << getTower_absPhiID(phi) << std::endl;
      }
    }
    

    //*******************************************************************
    //************* Do RCT geometry (HCAL) ******************************
    //*******************************************************************

    // Same idea as the ECAL RCT geometry, except we only care about the ET in towers 

    // Loop over hcal hits to get the HCAL towers.
    for (const auto& hit : hcalhits) {
      if (getCrystal_iPhi(hit.position().phi()) <= getCard_iPhiMax(cc) &&
          getCrystal_iPhi(hit.position().phi()) >= getCard_iPhiMin(cc) &&
          getCrystal_iEta(hit.position().eta()) <= getCard_iEtaMax(cc) &&
          getCrystal_iEta(hit.position().eta()) >= getCard_iEtaMin(cc) && hit.pt() > 0) {

	// HCAL: Get the crystal eta and phi, relative to the bottom left corner of the card 
	// (0 up to 17*5, 0 up to 4*5) 
	int local_iEta = getCrystal_local_iEta(hit.position().eta(), cc);
	int local_iPhi = getCrystal_local_iPhi(hit.position().phi(), cc);
	
	// HCAL: Once we have the iEta and iPhi of the crystal relative to the bottom left corner,
	// everything else is determined.
	
	// HCAL: Figure out what region (0-5) the hit falls into 
	// The region number (0-5) depends only on the local crystal iEta
	int regionNumber = getRegionNumber(local_iEta);
	
	// HCAL: Get the tower eta and phi index inside the card (17x4)
	int inCard_tower_iEta = int(local_iEta / CRYSTALS_IN_TOWER_ETA); 
	int inCard_tower_iPhi = int(local_iPhi / CRYSTALS_IN_TOWER_PHI);
	
	// HCAL: Get the tower eta and phi index inside the region (3x4)
	int inRegion_tower_iEta = inCard_tower_iEta % TOWER_IN_ETA;
	int inRegion_tower_iPhi = inCard_tower_iPhi % TOWER_IN_PHI;

	std::cout << "HCAL hit eta/phi : "  << hit.position().eta() << ", " << hit.position().phi() << ", region: " << regionNumber << ", "
		  << "inRegion_tower_iEta and iPhi : " << inRegion_tower_iEta << ", " << inRegion_tower_iPhi 
		  << ", energy " << hit.et_uint() 
		  << std::endl;

	// Access the right HCAL region -> tower and increment the ET
	if (regionNumber < N_REGIONS_PER_CARD) {
	  towers3x4& myTowers3x4 = rctCard.getTowers3x4(regionNumber);
	  towerHCAL& myTower = myTowers3x4.getTowerHCAL(inRegion_tower_iEta, inRegion_tower_iPhi);
	  myTower.addEt(hit.et_uint());

	}
      }
    } // end of loop over hcal hits

    //*******************************************************************    
    //******* Within each ECAL region, read back the hits ***************
    //******************************************************************* 

    // Cluster sort_clusterIn[N_CLUSTERS_PER_REGION * N_REGIONS_PER_CARD];    // array of clusters per card, to be sorted by ET
    // Cluster sort_clusterOut[N_CLUSTERS_PER_REGION * N_REGIONS_PER_CARD];   // array of clusters per card, sorted by ET
    std::vector<Cluster> sort_clusterIn;                  // Vector of clusters in the card
    tower_t towerEt[n_towers_cardEta][n_towers_cardPhi];  // 17x4 array of tower_t structs, representing one card 
    tower_t towerEtHCAL[n_towers_cardEta][n_towers_cardPhi]; 
    tower_t towerEtECAL[n_towers_cardEta][n_towers_cardPhi]; 

    for (int idxRegion = 0; idxRegion < N_REGIONS_PER_CARD; idxRegion++) {
      
      crystal temporary[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI];       // ECAL crystal array (will be changed)
      ap_uint<12> towerEtHCAL[TOWER_IN_ETA * TOWER_IN_PHI];    // HCAL tower ET in the 3x4 region

      if (cc > -1) {
	region3x4& myRegion = rctCard.getRegion3x4(idxRegion);
	towers3x4& myTowers = rctCard.getTowers3x4(idxRegion);

	std::cout << std::endl << "[----] DOING CARD " << cc << " AND REGION IDX " << myRegion.getIdx() << std::endl;

	// In each 3x4 region, loop through the links (one link per tower)
	for (int iLinkEta = 0; iLinkEta < TOWER_IN_ETA; iLinkEta++) {
	  for (int iLinkPhi = 0; iLinkPhi < TOWER_IN_PHI; iLinkPhi++) {
	    
	    // Get the ECAL link (one link per tower) 
	    linkECAL& myLink = myRegion.getLinkECAL(iLinkEta, iLinkPhi);
	    
	    // Each link has a different bottom left corner in the (iEta, iPhi) of the ECAL region.
	    // This will help us fill the ECAL region 'temporary' array, going link by link 
	    int ref_iEta = (iLinkEta * CRYSTALS_IN_TOWER_ETA); 
	    int ref_iPhi = (iLinkPhi * CRYSTALS_IN_TOWER_PHI); 

	    // In the link, get the crystals (5x5 in each link) 
	    for (int iEta = 0; iEta < CRYSTALS_IN_TOWER_ETA; iEta++) {           	    
	      for (int iPhi = 0; iPhi < CRYSTALS_IN_TOWER_PHI; iPhi++) {    
		
		// Et as unsigned int
		ap_uint<10> uEnergy = myLink.getCrystalE(iEta, iPhi);

		if (uEnergy > 0) {
		  std::cout << "energy>0: Accessing temporary array: " << (ref_iEta + iEta) << ", " << (ref_iPhi + iPhi) 
			    << ", writing energy (uint:) " << uEnergy << std::endl;
		}

		// Fill the 'temporary' array with a crystal object 
		temporary[ref_iEta + iEta][ref_iPhi + iPhi] = crystal(uEnergy);
	      }
	    } // end of loop over crystals

	    ///// Get HCAL tower ET
	    towerHCAL& myTower = myTowers.getTowerHCAL(iLinkEta, iLinkPhi);
	    //	    std::cout << "towerEt element # " << (iLinkEta * TOWER_IN_PHI) + iLinkPhi << ": will fill with Et " << myTower.getEt() << std::endl;
	    towerEtHCAL[(iLinkEta * TOWER_IN_PHI) + iLinkPhi] = myTower.getEt();
	  }
	}
	
	// Iteratively find four clusters and remove them from 'temporary' as we go, and fill sort_clusterIn
	for (int c = 0; c < N_CLUSTERS_PER_REGION; c++) {
	  Cluster newCluster = getClusterFromRegion3x4(temporary); // iteratively remove energy from 'temporary' 
	  if (newCluster.clusterEnergy() > 0) {                    // only add clusters with >0 energy to the vector 
	    sort_clusterIn.push_back(newCluster);
	  }
	}

	// Create towers using remaining ECAL energy, and the HCAL towers were already calculated in towersEtHCAL[12] 
	ap_uint<12> towerEtECAL[12];
	getECALTowersEt(temporary, towerEtECAL);
	
	// Add the towerETHCAL and towerETECAL arrays, and fill the 17x4 array of towerEt 
	for (int i = 0; i < 12; i++) {
	  ap_uint<12> towerTotalEt = towerEtHCAL[i] + towerEtECAL[i];
	  // Get the tower's indices in a (17x4) card
	  int iEta = (idxRegion * TOWER_IN_ETA) + (i / TOWER_IN_PHI);   
	  int iPhi = (i % TOWER_IN_PHI);
	  //	  std::cout << "(" << iEta << "," << iPhi << ")";
	  towerEt[iEta][iPhi]     = tower_t(towerTotalEt, 0, 0);
	  towerEtHCAL[iEta][iPhi] = tower_t(towerEtHCAL[i], 0, 0);
	  towerEtECAL[iEta][iPhi] = tower_t(towerEtECAL[i], 0, 0);
	}

      } // end of "if"

    } // end of loop over regions

    //-------------------------------------------//
    // Sort the clusters                         //
    //-------------------------------------------//
    // std::cout << "card " << cc << ": unsorted ET: ";
    // for (auto & c : sort_clusterIn) {
    //   std::cout << c.clusterEnergy() << " (" << c.clusterEta() << ", " << c.clusterPhi() << ") ";
    // }
    // std::cout << std::endl;
    
    // Sort the clusters in decreasing ET
    std::sort(sort_clusterIn.begin(), sort_clusterIn.end(), compareClusterET);

    //-------------------------------------------//
    // Write the outputs to the L1 outputs       //
    //-------------------------------------------//

    // Print the sorted vector
    std::cout << "Sanity check: Card " << cc << ": SORTED ET: ";
    // for (auto & c : sort_clusterIn) {
    for (unsigned int jj = 0; jj < unsigned(sort_clusterIn.size()) && (jj < n_clusters_4link); ++jj) {
      Cluster c = sort_clusterIn[jj];

      std::cout << c.clusterEnergy() << " (" << c.clusterEta() << ", " << c.clusterPhi() << ") "
		<< ", setting crystal_ID_clusterL1Card[" << jj % n_links_card << "]["
		<< jj / n_links_card << "][" << cc << "] = " << getCrystalIDInTower_emulator(c.clusterEta(), c.clusterPhi()) 
		<< std::endl;
      
      // Distribute (up to 12) clusters across 4 links
      // crystalID_cluster_L1Card
      crystalID_cluster_L1Card[jj % n_links_card][jj / n_links_card][cc] = 
	getCrystalIDInTower_emulator(c.clusterEta(), c.clusterPhi());
      
      
    }
    
    std::cout << std::endl;

    // Do the towers. The firmware 17x4 array treats the "bottom left" corner of the card as (0, 0)
    // (rotating the negative eta cards so that the endcap region is pointing up)
    // while the emulator treats the 4x17 array as always starting in the top left corner if we look
    // at the usual RCT diagram.
    std::cout << "Sanity check: Card " << cc << ": towers (remaining ECAL, plus HCAL): ";
    for (int i = 0; i < n_towers_cardEta; i++) {
      for (int j = 0; j < n_towers_cardPhi; j++ ) {
        
      std::cout << towerEt[i][j].et() << " "; 
      
      // n.b. L1 output is 4*17*36, hence [j][i] instead of [i][j] on the L.H.S.                       
      if ((cc % 2) == 1) { // if cc is odd (positive eta)
	ECAL_tower_L1Card[j][i][cc] = towerEtECAL[i][j].et();
	HCAL_tower_L1Card[j][i][cc] = towerEtHCAL[i][j].et();
      }
      else {  // if cc is even (negative eta), we need to rotate the coordinates
	ECAL_tower_L1Card[j][i][cc] = towerEtECAL[16-i][3-j].et();
	HCAL_tower_L1Card[j][i][cc] = towerEtHCAL[16-i][3-j].et();
      }
      }
    }
    std::cout << std::endl;
    
  } // end of loop over cards
  
  std::cout << "I'm here!" << std::endl;
  
  
}

////////////////////////////////////////////////////////////////////////// 

//define this as a plug-in
DEFINE_FWK_MODULE(EGammaCrystalsProducer);

