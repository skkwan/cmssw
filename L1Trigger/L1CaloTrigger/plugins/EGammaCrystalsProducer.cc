/* 
 * Description: Read crystal-level information
 */

// system include files
#include <memory>
#include <array>
#include <iostream>
#include <cmath>

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

static constexpr int n_crystals_towerEta = 5;
static constexpr int n_crystals_towerPhi = 5;
static constexpr int n_towers_Eta = 34;
static constexpr int n_towers_Phi = 72;
static constexpr int n_towers_halfPhi = 36;
static constexpr int n_towers_cardEta = 17;   // new: equivalent to n_towers_per_link
static constexpr int n_towers_cardPhi = 4;    // new
static constexpr int n_crystals_cardEta = (n_towers_Eta * n_towers_cardEta);
static constexpr int n_crystals_cardPhi = (n_towers_Phi * n_towers_cardPhi);
static constexpr int TOWER_IN_ETA = 3;
static constexpr int TOWER_IN_PHI = 4;
static constexpr float ECAL_eta_range = 1.4841;
static constexpr float cut_500_MeV = 0.5;

// Get crystal's iEta from real eta. (identical to getCrystal_etaID in L1EGammaCrystalsEmulatorProducer.cc)
int getCrystal_iEta(float eta) {
  float size_cell = 2 * ECAL_eta_range / (n_crystals_towerEta * n_towers_Eta);
  int iEta = int((eta + ECAL_eta_range) / size_cell);
  return iEta;
}

// Get crystal's iPhi from real phi. (identical to getCrystal_phiID in L1EGammaCrystalsEmulatorProducer.cc)
int getCrystal_iPhi(float phi) {
  float size_cell = 2 * M_PI / (n_crystals_towerPhi * n_towers_Phi);
  int iPhi = int((phi + M_PI) / size_cell);
  return iPhi;
}

// Card boundaries in eta (identical to getEtaMax_card)
int getCard_iEtaMax(int card) {
  int etamax = 0;
  if (card % 2 == 0)
    etamax = n_towers_cardEta* n_crystals_towerEta - 1;  // First eta half. 5 crystals in eta in 1 tower.
  else
    etamax = n_towers_Eta * n_crystals_towerEta - 1;
  return etamax;
}

int getCard_iEtaMin(int card) {
  int etamin = 0;
  if (card % 2 == 0)
    etamin = 0 * n_crystals_towerEta;  // First eta half. 5 crystals in eta in 1 tower.
  else
    etamin = n_towers_cardEta* n_crystals_towerEta;
  return etamin;
}

// Card boundaries in phi
int getCard_iPhiMax(int card) {
  int phimax = ((card / 2) + 1) * 4 * n_crystals_towerPhi - 1;
  return phimax;
}

int getCard_iPhiMin(int card) {
  int phimin = (card / 2) * 4 * n_crystals_towerPhi;
  return phimin;
}

// absolute Eta IDs range from 0-33 (Adapted from getTower_absoluteEtaID)
int getTower_absEtaID(float eta) {
  float size_cell = 2 * ECAL_eta_range / n_towers_Eta;
  int etaID = int((eta + ECAL_eta_range) / size_cell);
  return etaID;
}

// absolute Phi IDs range from 0-71 (Adapted from getTower_absolutePhiID)
int getTower_absPhiID(float phi) {
  float size_cell = 2 * M_PI / n_towers_Phi;
  int phiID = int((phi + M_PI) / size_cell);
  return phiID;
}

// given the RCT card number (0-35), get the TOWER iPhi of the "bottom left" corner (0-71)
int getCard_refTower_iPhi(int cc) {
  if ((cc < 0) || (cc > 35)) return -1;  // out of bounds!
      
  if ((cc % 2) == 1) { // if cc is odd
    return (int(cc / 2) * 4);
  }
  else {  // if cc is even
    return ((int(cc / 2) * 4) + 3);
  }
  
}

// given the RCT card number (0-35), get the TOWER iEta of the "bottom left" corner (0-33)
int getCard_refTower_iEta(int cc) {
  if ((cc < 0) || (cc > 35)) return -1;  // out of bounds!  

  if ((cc % 2) == 1) { // if cc is odd
    return 17;
  }
  else {  // if cc is even 
    return 16;
  }

}

// given the RCT card number (0-35), get the crystal iPhi of the "bottom left" corner (0- 71*5)
int getCard_ref_iPhi(int cc) {
  if ((cc < 0) || (cc > 35)) return -1;  // out of bounds!  
  if ((cc % 2) == 1) { // if cc is odd                                                                                     
    return (int(cc / 2) * TOWER_IN_PHI * n_crystals_towerPhi);
  }
  else {  // if cc is even, the bottom left corner is further in iPhi, hence the +4                
    return (((int(cc / 2) * TOWER_IN_PHI) + 4) * n_crystals_towerPhi);
  }
}

// given the RCT card number (0-35), get the crystal iEta of the "bottom left" corner (0-33*5)
int getCard_ref_iEta(int cc) {
  if ((cc < 0) || (cc > 35)) return -1;  // out of bounds! 
  if ((cc % 2) == 1) {  // if cc is odd
    return ((16 * n_crystals_towerEta) + 1);
  }
  else {   // if cc is even
    return (16 * n_crystals_towerEta);
  }
  
}


// get the RCT card region that a crystal is in, given the "local" iEta of the crystal 
// 0 is region closest to eta = 0. Regions 0, 1, 2, 3, 4 are in the barrel, Region 5 is in overlap
int getRegionNumber(int local_iEta) {
  int no = int(local_iEta / (TOWER_IN_ETA * n_crystals_towerEta));
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

// Declare the SimpleCaloHit class

class SimpleCaloHit {
private:
  float pt_ = 0;
  float energy_ = 0.;
  GlobalVector position_;  // As opposed to GlobalPoint, so we can add them (for weighted average)
  EBDetId id_;
  
public:
  // tool functions
  inline void setPt() { pt_ = (position_.mag2() > 0) ? energy_ * sin(position_.theta()) : 0; };
  inline void setEnergy(float et) { energy_ = et / sin(position_.theta()); };
  inline void setPosition(const GlobalVector& pos) { position_ = pos; };
  inline void setId(const EBDetId& id) { id_ = id; };

  inline float pt() const { return pt_; };
  inline float energy() const { return energy_; };
  inline const GlobalVector& position() const { return position_; };
  inline const EBDetId& id() const { return id_; };
};

/*******************************************************************/

/*
 * linkECAL class: represents one ECAL link (one tower: 5x5 crystals)
 */

class linkECAL {
private:
  float crystalE[n_crystals_towerEta][n_crystals_towerPhi] = {{0}};  // a 5x5 array

public:
  // Set members
  inline void addCrystalE(int iEta, int iPhi, float energy) { crystalE[iEta][iPhi] += energy; }

  // Access members
  inline float getCrystalE(int iEta, int iPhi) { return crystalE[iEta][iPhi]; }
  
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
  // set members
  inline void setIdx(int idx) { idx_ = idx; };

  // get members
  inline float getIdx() const { return idx_; };
  inline linkECAL getLinkECAL(int iEta, int iPhi) { return linksECAL[iEta][iPhi]; }
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
  // set members
  inline void setIdx(int idx) { idx_ = idx; };

  // get members
  inline float getIdx() const { return idx_; };
  inline region3x4 getRegion3x4(int idx) { assert(idx < 5); return card3x4Regions[idx]; }

};

/*******************************************************************/


//////////////////////////////////////////////////////////////////////////  

// Remaining EGammaCrystalsProducer initializer, destructor, and produce methods

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


  // 
  std::vector<card> rctCards; 
  for (int cc = 0; cc < n_towers_halfPhi; ++cc) {  // Loop over 36 L1 cards
  
    card rctCard;
    rctCard.setIdx(cc);
    for (const auto& hit : ecalhits) {
      // Check if the hit is in cards 0-35
      if (getCrystal_iPhi(hit.position().phi()) <= getCard_iPhiMax(cc) &&
  	  getCrystal_iPhi(hit.position().phi()) >= getCard_iPhiMin(cc) &&
  	  getCrystal_iEta(hit.position().eta()) <= getCard_iEtaMax(cc) &&
  	  getCrystal_iEta(hit.position().eta()) >= getCard_iEtaMin(cc)) {
	
	std::cout << "Card: " << cc << ", hit (Eta, phi): " 
		  << hit.position().eta() << ", " << hit.position().phi() << std::endl;

	

	// For now, only handle positive eta cards (i.e. cc is an odd integer)
	if ((cc % 2) == 1) {
	  // Get the crystal eta and phi, relative to the bottom left corner of the card 
	  // (again assuming positive eta) (0 up to 17*5, 0 up to 4*5) 
	  int inCard_crystal_iEta = getCrystal_iEta(hit.position().eta()) - getCard_ref_iEta(cc);
	  int inCard_crystal_iPhi = getCrystal_iPhi(hit.position().phi()) - getCard_ref_iPhi(cc);

	  // Once we have the iEta and iPhi of the crystal relative to the bottom left corner,
	  // everything else is determined.

	  // Figure out what region (0-5) the hit falls into 
	  // The region number (0-5) depends only on the local crystal iEta
	  int regionNumber = getRegionNumber(inCard_crystal_iEta);

	  // Get the tower eta and phi index inside the card (17x4)
	  int inCard_tower_iEta = int(inCard_crystal_iEta / n_crystals_towerEta); 
	  int inCard_tower_iPhi = int(inCard_crystal_iPhi / n_crystals_towerPhi);

	  // Get the tower eta and phi index inside the region (3x4)
	  int inRegion_tower_iEta = inCard_tower_iEta % TOWER_IN_ETA;
	  int inRegion_tower_iPhi = inCard_tower_iPhi % TOWER_IN_PHI;

	  // Get the crystal eta and phi index inside the region (15x20)
	  int inRegion_crystal_iEta = inCard_crystal_iEta % (TOWER_IN_ETA * n_crystals_towerEta);
	  int inRegion_crystal_iPhi = inCard_crystal_iPhi;

	  std::cout << "inCard_crystal_iEta, inCard_crystal_iPhi, regionNumber: "
		    << inCard_crystal_iEta << ", " 
		    << inCard_crystal_iPhi << ", "
		    << regionNumber << std::endl;
	  std::cout << "inRegion_crystal_iEta, inRegion_crystal_iPhi, regionNumber: "
                    << inRegion_crystal_iEta << ", "
                    << inRegion_crystal_iPhi << ", "
                    << regionNumber << std::endl;

	  // Within the region, figure out which crystal and link the hit falls into
	  // Add the hit's energy to the right crystal/ right region
	  
	  if (regionNumber < 5) {
	    region3x4 myRegion = rctCard.getRegion3x4(regionNumber);
	    linkECAL myLink = myRegion.getLinkECAL(inRegion_tower_iEta, inRegion_tower_iPhi);

	    // Add the energy
	    myLink.addCrystalE(inRegion_crystal_iEta, inRegion_crystal_iPhi, hit.energy());

	    // Read it back out
	    float energy = myLink.getCrystalE(inRegion_crystal_iEta, inRegion_crystal_iPhi);
	    std::cout << "energy " << energy << std::endl;
	  }
	}
	
      }
    }
    rctCards.push_back(rctCard);
  }

  // // Sanity check: print out
  // for (const auto& card : rctCards) {
  //   std::cout << card.getIdx() << std::endl;
  // }
    
  std::cout << "I'm here!" << std::endl;

 
}

////////////////////////////////////////////////////////////////////////// 

//define this as a plug-in
DEFINE_FWK_MODULE(EGammaCrystalsProducer);

