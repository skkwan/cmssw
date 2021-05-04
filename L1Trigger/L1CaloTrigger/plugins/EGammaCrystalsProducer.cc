/* 
 * Description: Read crystal-level information
 */

// system include files
#include <array>
#include <cmath>
#include <cstdint>
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
static constexpr int n_towers_cardPhi = 4;    // new
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
  uint16_t et_uint_ = 0;
  GlobalVector position_;  // As opposed to GlobalPoint, so we can add them (for weighted average)
  EBDetId id_;
  
public:
  // tool functions
  inline void setPt() { pt_ = (position_.mag2() > 0) ? energy_ * sin(position_.theta()) : 0; };
  inline void setEnergy(float et) { energy_ = et / sin(position_.theta()); };
  inline void setEt_uint(int et) { et_uint_ = (uint16_t) et; }
  inline void setPosition(const GlobalVector& pos) { position_ = pos; };
  inline void setId(const EBDetId& id) { id_ = id; };

  inline float pt() const { return pt_; };
  inline float energy() const { return energy_; };
  inline uint16_t et_uint() const { return et_uint_; };
  inline const GlobalVector& position() const { return position_; };
  inline const EBDetId& id() const { return id_; };
};

/*******************************************************************/

/*
 * linkECAL class: represents one ECAL link (one tower: 5x5 crystals)
 */

class linkECAL {
private:
  uint16_t totalE = 0;
  uint16_t crystalE[CRYSTALS_IN_TOWER_ETA][CRYSTALS_IN_TOWER_PHI] = {};  // a 5x5 array

public:
  // constructor                                                                                               
  linkECAL() { }

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
  inline void setCrystalE(int iEta, int iPhi, uint16_t energy) { assert(iEta < 5); assert(iPhi < 5); crystalE[iEta][iPhi] = energy; }
  inline void addCrystalE(int iEta, int iPhi, uint16_t energy) { 
    assert(iEta < 5); assert(iPhi < 5);
    crystalE[iEta][iPhi] += energy; 
    totalE += energy; }
  
  
  // Access members
  inline uint16_t getCrystalE(int iEta, int iPhi) const { assert(iEta < 5); assert(iPhi < 5); return crystalE[iEta][iPhi]; }
  inline uint16_t getTotalE() const { return totalE; }

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
  uint16_t energy;   // formerly ap_uint<10>
  uint8_t  timing; // formerly ap_uint<4>

  crystal(){
    energy = 0;
    timing = 0;
  }

  crystal(uint16_t energy){  // To-do: add timing information
    this->energy = energy;
    this->timing = 0; 
  }

  crystal& operator=(const crystal& rhs){
    energy = rhs.energy;
    timing = rhs.timing;
    return *this;
  }
};

/*******************************************************************/

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
	ehit.setEt_uint(hit.encodedEt());  // also save the 10-bit Et
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

      if (idxRegion == 0) { // TEMP: only try one region

	region3x4& myRegion = rctCard.getRegion3x4(idxRegion);
	
	std::cout << std::endl << "[----] DOING REGION IDX " << myRegion.getIdx() << std::endl;
	
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
		uint16_t uEnergy = myLink.getCrystalE(iEta, iPhi);

		std::cout << "Accessing temporary array: " << (ref_iEta + iEta) << ", " 
			  << (ref_iPhi + iPhi) 
			  << ", writing energy (uint:) " << uEnergy << std::endl;

		// Fill the 'temporary' array with a crystal object 
		temporary[ref_iEta + iEta][ref_iPhi + iPhi] = crystal(uEnergy);

	      }
	    }
	  }
	}
      


      } // end of "if"
    } // end of loop over regions  
  } // end of loop over cards
    
  std::cout << "I'm here!" << std::endl;
  
 
}

////////////////////////////////////////////////////////////////////////// 

//define this as a plug-in
DEFINE_FWK_MODULE(EGammaCrystalsProducer);

