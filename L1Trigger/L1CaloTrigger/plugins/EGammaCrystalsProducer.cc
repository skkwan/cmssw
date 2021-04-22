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

////////////////////////////////////////////////////////////////////////// 

// Declare the class and its methods

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

// Class method definitions

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

  for (int cc = 0; cc < n_towers_halfPhi; ++cc) {  // Loop over 36 L1 cards
    
    for (const auto& hit : ecalhits) {
      if (getCrystal_iPhi(hit.position().phi()) <= getCard_iPhiMax(cc) &&
  	  getCrystal_iPhi(hit.position().phi()) >= getCard_iPhiMin(cc) &&
  	  getCrystal_iEta(hit.position().eta()) <= getCard_iEtaMax(cc) &&
  	  getCrystal_iEta(hit.position().eta()) >= getCard_iEtaMin(cc)) {
	
	std::cout << "Card: " << cc << ", hit (Eta, phi): " 
		  << hit.position().eta() << ", " << hit.position().phi() << std::endl;
      }
    }

  }
    
  std::cout << "I'm here!" << std::endl;

 
}

////////////////////////////////////////////////////////////////////////// 

//define this as a plug-in
DEFINE_FWK_MODULE(EGammaCrystalsProducer);

