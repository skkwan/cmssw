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

static constexpr int CRYSTALS_IN_TOWER_ETA = 5;
static constexpr int CRYSTALS_IN_TOWER_PHI = 5;

static constexpr int TOWER_IN_ETA = 3;      // number of towers in eta, in one 3x4 region
static constexpr int TOWER_IN_PHI = 4;      // number of towers in phi, in one 3x4 region
static constexpr int CRYSTAL_IN_ETA = 15;   // number of crystals in eta, in one 3x4 region
static constexpr int CRYSTAL_IN_PHI = 20;   // number of crystals in phi, in one 3x4 region 

static constexpr float ECAL_eta_range = 1.4841;
static constexpr float cut_500_MeV = 0.5;

static constexpr int NB_CLUSTERS = 5;       // number of clusters per region

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
  
  if ((cc % 2) == 1) { // if cc is odd                                                                                     
    return (int(cc / 2) * TOWER_IN_PHI * CRYSTALS_IN_TOWER_PHI);
  }
  else {  // if cc is even, the bottom left corner is further in iPhi, hence the +4                
    return (((int(cc / 2) * TOWER_IN_PHI) + 4) * CRYSTALS_IN_TOWER_PHI);
  }
}

// Given the RCT card number (0-35), get the crystal iEta of the "bottom left" corner (0-33*5)
int getCard_ref_iEta(int cc) {

  assert(isValidCard(cc));
  if ((cc % 2) == 1) {  // if cc is odd (positive eta)
    return ((16 * CRYSTALS_IN_TOWER_ETA) + 1);
  }
  else {   // if cc is even (negative eta)
    return (16 * CRYSTALS_IN_TOWER_ETA);
  }
}

// For a crystal with real (eta, phi) and falling in card cc, get its local iEta 
// relative to the bottom left corner of the card (possible local iEta ranges from 0 to 17 * 5,
// since in one card, there are 17 towers in eta, each with 5 crystals in eta.
int getCrystal_local_iEta(float hitEta, int cc) {
  assert(isValidCard(cc));

  // if ((cc % 2) == 1) {  // if cc is odd (positive eta)  
  //   return (getCrystal_iEta(hitEta) - getCard_ref_iEta(cc));
  // }
  // else { // if cc is even (negative eta)
  //   return (getCard_ref_iEta(cc) - getCrystal_iEta(hitEta));
  // }

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
  EBDetId id_;
  
public:
  // tool functions
  inline void setPt() { pt_ = (position_.mag2() > 0) ? energy_ * sin(position_.theta()) : 0; };
  inline void setEnergy(float et) { energy_ = et / sin(position_.theta()); };
  inline void setEt_uint(ap_uint<10> et_uint) { et_uint_ = et_uint; }
  inline void setPosition(const GlobalVector& pos) { position_ = pos; };
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
  ap_uint<10> totalE;
  ap_uint<10> crystalE[CRYSTALS_IN_TOWER_ETA][CRYSTALS_IN_TOWER_PHI];

public:
  // constructor                                                                                               
  linkECAL() { 
    totalE = 0;
    crystalE[CRYSTALS_IN_TOWER_ETA][CRYSTALS_IN_TOWER_PHI] = {};
  }

  // copy constructor
  linkECAL(const linkECAL &other) {
    totalE = other.totalE;
    std::cout << "Inside copy constructor: getting energy "<< totalE << std::endl;
    for (int i = 0; i < CRYSTALS_IN_TOWER_ETA; i++) {
      for (int j = 0; j < CRYSTALS_IN_TOWER_PHI; j++ ) {
        crystalE[i][j] = other.crystalE[i][j];
      }
    }
  }

  // overload operator= to use copy constructor
  linkECAL operator=(const linkECAL &other) {
    linkECAL newLink(other);
    return newLink;
  }
  
  // Set members
  inline void setCrystalE(int iEta, int iPhi, ap_uint<10> energy) { assert(iEta < 5); assert(iPhi < 5); crystalE[iEta][iPhi] = energy; }
  inline void addCrystalE(int iEta, int iPhi, ap_uint<10> energy) { 
    assert(iEta < 5); assert(iPhi < 5);
    crystalE[iEta][iPhi] += energy; 
    totalE += energy; }
  
  
  // Access members
  inline ap_uint<10> getCrystalE(int iEta, int iPhi) { assert(iEta < 5); assert(iPhi < 5); return crystalE[iEta][iPhi]; }
  inline ap_uint<10> getTotalE() { return totalE; }

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
  inline void setIdx(int idx) { idx_ = idx; };

  // get members
  inline float getIdx() const { return idx_; };
  inline linkECAL& getLinkECAL (int iEta, int iPhi) { return linksECAL[iEta][iPhi]; }
};

/*******************************************************************/

/* 
 * card class: represents one RCT card. Each card has five 3x4 regions.
 *             idx 0-35: odd values of cardIdx span eta = 0 to eta = 1.41 
 *                       even values of cardIdx span eta = -1.41 to eta = 0
 */

class card {
private:
  int idx_ = -1 ; 
  region3x4 card3x4Regions[5];
  
public:
  // constructor
  card() { 
    idx_ = -1; 
    for (int i = 0; i < 5; i++) {
      card3x4Regions[i].setIdx(i);
    }
  }
  
  // copy constructor
  card(const card& other) {
    idx_ = other.idx_;
    for (int i = 0; i < 5; i++) { card3x4Regions[i] = other.card3x4Regions[i]; }
  }

  // overload operator= to use copy constructor
  card operator=(const card& other) {
    card newCard(other);
    return newCard;
  }

  // set members
  inline void setIdx(int idx) { idx_ = idx; };

  // get members
  inline float getIdx() const { return idx_; };
  inline region3x4& getRegion3x4(int idx) { assert(idx < 5); return card3x4Regions[idx]; }

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


  // // Get all the HCAL hits
  // edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hbhecoll;
  // iEvent.getByToken(hcalTPToken_, hbhecoll);
  
  // for (const auto& hit : *hbhecoll.product()) {
  //   float et = decoder_->hcaletValue(hit.id(), hit.t0());
  //   if (et <= 0)
  //     continue;
    
  //   if (!(hcTopology_->validHT(hit.id()))) {
  //     LogError("EGammaCrystalsProducer")
  // 	<< " -- Hcal hit DetID not present in HCAL Geom: " << hit.id() << std::endl;
  //     throw cms::Exception("EGammaCrystalsProducer");
  //     continue;
  //   }
  //   const std::vector<HcalDetId>& hcId = theTrigTowerGeometry.detIds(hit.id());
  //   if (hcId.empty()) {
  //     LogError("EGammaCrystalsProducer")
  // 	<< "Cannot find any HCalDetId corresponding to " << hit.id() << std::endl;
  //     throw cms::Exception("EGammaCrystalsProducer");
  //     continue;
  //   }
  //   if (hcId[0].subdetId() > 1)
  //     continue;
  //   GlobalVector hcal_tp_position = GlobalVector(0., 0., 0.);
  //   for (const auto& hcId_i : hcId) {
  //     if (hcId_i.subdetId() > 1)
  //       continue;
  //     // get the first HCAL TP/ cell
  //     auto cell = hbGeometry->getGeometry(hcId_i);
  //     if (cell == nullptr)
  // 	continue;
  //     GlobalVector tmpVector = GlobalVector(cell->getPosition().x(), cell->getPosition().y(), cell->getPosition().z());
  //     hcal_tp_position = tmpVector;
      
  //     std::cout << "Found HCAL cell/TP with coordinates " << cell->getPosition().x() << ","
  //      		<< cell->getPosition().y() << ","
  //      		<< cell->getPosition().z() << " and ET (GeV) "
  // 		<< et << std::endl;
      
  //     break;
  //   }
  //   //    std::cout << "HCAL hit et: " << et << std::endl;
  // }

  //*******************************************************************
  //********************** Do RCT geometry  ***************************
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
	
	std::cout << "Card: " << cc << ", hit (Eta, phi): " 
		  << hit.position().eta() << ", " << hit.position().phi() << std::endl;

	

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
	
	// Get the crystal eta and phi index inside the 3x4 region (15x20)
	int inRegion_crystal_iEta = local_iEta % (TOWER_IN_ETA * CRYSTALS_IN_TOWER_ETA);
	int inRegion_crystal_iPhi = local_iPhi;

	// Get the crystal eta and phi index inside the link (5x5). E.g. crystal (7, 9) is crystal number (2, 4) in its link
	int inLink_crystal_iEta = (inRegion_crystal_iEta % CRYSTALS_IN_TOWER_ETA);
	int inLink_crystal_iPhi = (inRegion_crystal_iPhi % CRYSTALS_IN_TOWER_PHI);

	std::cout << "local_iEta, local_iPhi, regionNumber: "
		  << local_iEta << ", " 
		  << local_iPhi << ", "
		  << regionNumber << std::endl;
	std::cout << "inRegion_crystal_iEta, inRegion_crystal_iPhi (expecting 15x20), regionNumber (expecting 0-5): "
		  << inRegion_crystal_iEta << ", "
		  << inRegion_crystal_iPhi << ", "
		  << regionNumber << std::endl;
	
	// Within the region, figure out which crystal and link the hit falls into
	// Add the hit's energy to the right crystal/ right region
	
	if (regionNumber < 5) {
	  region3x4& myRegion = rctCard.getRegion3x4(regionNumber);
	  //	  myRegion.setIdx(regionNumber); // a bit redundant

	  std::cout << "inRegion_tower_iEta, inRegion_tower_iPhi (expecting 3x4): " << inRegion_tower_iEta << ", "
		    << inRegion_tower_iPhi << std::endl;
	  
	  // Get the link 
	  linkECAL& myLink = myRegion.getLinkECAL(inRegion_tower_iEta, inRegion_tower_iPhi);
	  
	  // Add the energy to the right 5x5 crystal 
	  std::cout << "inLink_crystal_iEta, inLink_crystal_iPhi (expecting 5x5): " 
		    << inLink_crystal_iEta << ", " << inLink_crystal_iPhi << std::endl;

	  float energyBefore = myLink.getCrystalE(inLink_crystal_iEta, inLink_crystal_iPhi);
	  //	  myLink.addCrystalE(inLink_crystal_iEta, inLink_crystal_iPhi, hit.energy());
	  myLink.addCrystalE(inLink_crystal_iEta, inLink_crystal_iPhi, hit.et_uint());
	  
	  float energy = myLink.getCrystalE(inLink_crystal_iEta, inLink_crystal_iPhi);                   
	  std::cout << "energy before/after " << energyBefore << " " << energy << std::endl; 
	  std::cout << "total energy in the link: " << myLink.getTotalE() << std::endl;


	  
	}
      } 
      
      
    } // end of loop over hits
    

    //*******************************************************************    
    //******* Within each region, read back the hits ********************
    //******************************************************************* 

    crystal temporary[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI];
    
    
    for (int idxRegion = 0; idxRegion < 5; idxRegion++) {

      if ((cc == 35) && (idxRegion == 0)) { // TEMP: only try one region: this should be equivalent to processInputLinks

	region3x4& myRegion = rctCard.getRegion3x4(idxRegion);
	
	std::cout << std::endl << "[----] DOING CARD  " << cc << " AND REGION IDX " << myRegion.getIdx() << std::endl;
	
	// In each 3x4 region, loop through the links (one link per tower)
	for (int iLinkEta = 0; iLinkEta < TOWER_IN_ETA; iLinkEta++) {
	  for (int iLinkPhi = 0; iLinkPhi < TOWER_IN_PHI; iLinkPhi++) {
	    
	    // Get the link (one link per tower) 
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

		std::cout << "Accessing temporary array: " << (ref_iEta + iEta) << ", " << (ref_iPhi + iPhi) 
			  << ", writing energy (uint:) " << uEnergy << std::endl;

		// Fill the 'temporary' array with a crystal object 
		temporary[ref_iEta + iEta][ref_iPhi + iPhi] = crystal(uEnergy);
	      }
	    }
	  }
	}

	// Use the 'temporary' array to fill an ecalRegion_t object 
	ecalRegion_t ecalRegion;
	ecalRegion = initStructure(temporary);

	

      } // end of "if"

      
      // 
      
      //      Cluster sort_clusterIn[NB_CLUSTERS];
      //      Cluster sort_clusterOut[NB_CLUSTERS];

      

    } // end of loop over regions  
  } // end of loop over cards
    
  std::cout << "I'm here!" << std::endl;
  
 
}

////////////////////////////////////////////////////////////////////////// 

//define this as a plug-in
DEFINE_FWK_MODULE(EGammaCrystalsProducer);

