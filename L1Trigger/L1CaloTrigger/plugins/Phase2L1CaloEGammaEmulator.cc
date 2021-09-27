/* 
 * Description: Phase 2 RCT Layer 1 emulator: create ECAL crystal collections
 */

// system include files
#include <ap_int.h>
#include <array>
#include <cmath>
// #include <cstdint>
#include <iostream>
#include <fstream>
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

#include "Phase2L1CaloEGammaEmulator.h"
#include "Phase2L1RCT.h"
#include "Phase2L1GCT.h"


// Declare the Phase2L1CaloEGammaEmulator class and its methods

class Phase2L1CaloEGammaEmulator : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1CaloEGammaEmulator(const edm::ParameterSet&);
  ~Phase2L1CaloEGammaEmulator() override;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  bool passes_ss(float pt, float ss);
  // bool passes_photon(float pt, float pss);
  // bool passes_iso(float pt, float iso);
  bool passes_looseTkss(float pt, float ss);
  // bool passes_looseTkiso(float pt, float iso);

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

// Phase2L1CaloEGammaEmulator initializer, destructor, and produce methods

Phase2L1CaloEGammaEmulator::Phase2L1CaloEGammaEmulator(const edm::ParameterSet & iConfig)
  : ecalTPEBToken_(consumes<EcalEBTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("ecalTPEB"))),
    hcalTPToken_(
		 consumes<edm::SortedCollection<HcalTriggerPrimitiveDigi> >(iConfig.getParameter<edm::InputTag>("hcalTP"))),
  calib_(iConfig.getParameter<edm::ParameterSet>("calib")) {
  produces<l1tp2::CaloCrystalClusterCollection>();
  produces<l1tp2::CaloTowerCollection>();
}

Phase2L1CaloEGammaEmulator::~Phase2L1CaloEGammaEmulator() {}



void Phase2L1CaloEGammaEmulator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
  // Declare output collections
  ///////////////////////////////////////////////////// 

  auto L1EGXtalClusters = std::make_unique<l1tp2::CaloCrystalClusterCollection>();
  auto L1CaloTowers = std::make_unique<l1tp2::CaloTowerCollection>();

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
	ehit.setEt_uint((ap_uint<10>) hit.encodedEt());  // also save the uint Et
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
      LogError("Phase2L1CaloEGammaEmulator")
  	<< " -- Hcal hit DetID not present in HCAL Geom: " << hit.id() << std::endl;
      throw cms::Exception("Phase2L1CaloEGammaEmulator");
      continue;
    }
    const std::vector<HcalDetId>& hcId = theTrigTowerGeometry.detIds(hit.id());
    if (hcId.empty()) {
      LogError("Phase2L1CaloEGammaEmulator")
  	<< "Cannot find any HCalDetId corresponding to " << hit.id() << std::endl;
      throw cms::Exception("Phase2L1CaloEGammaEmulator");
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
  //
  // L1 Outputs definition (using previous emulator by Cecile) - NOT firmware convention
  // Firmware convention -> CMSSW indexing conversion is done near the end.
  // 36 L1 cards, each with 4x17 towers. All using CMSSW indexing convention, NOT firmware convention
  //
  ap_uint<12> ECAL_tower_L1Card[n_links_card][n_towers_per_link][n_towers_halfPhi];
  ap_uint<12> HCAL_tower_L1Card[n_links_card][n_towers_per_link][n_towers_halfPhi];
  int iEta_tower_L1Card[n_links_card][n_towers_per_link][n_towers_halfPhi];  
  int iPhi_tower_L1Card[n_links_card][n_towers_per_link][n_towers_halfPhi];
  // 36 L1 cards with 4 links, each with 2 clusters (up to 8 per card)
  ap_uint<12> energy_cluster_L1Card[n_links_card][n_clusters_link][n_towers_halfPhi]; 
  int crystalID_cluster_L1Card[n_links_card][n_clusters_link][n_towers_halfPhi];  // range: [0, 5*5) 
  int towerID_cluster_L1Card[n_links_card][n_clusters_link][n_towers_halfPhi];    // range: [0, 17*4)

  //
  // L1 Outputs definition: Arrays that use firmware convention for indexing
  //
  tower_t towerHCALCard[n_towers_cardEta][n_towers_cardPhi][n_towers_halfPhi]; // 17x4x36 array (not to be confused with the 12x1 array of ap_uints, towerEtHCAL                                               
  tower_t towerECALCard[n_towers_cardEta][n_towers_cardPhi][n_towers_halfPhi];
  // There is one vector of clusters per card (up to 12 clusters before stitching across ECAL regions)
  std::vector<Cluster> cluster_list[n_towers_halfPhi];
  // After merging/stitching the clusters, we only take the 8 highest pt per card                      
  std::vector<Cluster> cluster_list_merged[n_towers_halfPhi];


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
	energy_cluster_L1Card[ii][jj][ll]    = 0;
	towerID_cluster_L1Card[ii][jj][ll]   = 0;
        crystalID_cluster_L1Card[ii][jj][ll] = 0;
	
      }
    }
  }
  //*******************************************************************
  //*************** Do RCT geometry (ECAL)  ***************************
  //*******************************************************************


  int theCard = 28; // Debugging only
  
  for (int cc = 0; cc < n_towers_halfPhi; ++cc) {  // Loop over 36 L1 cards
    // Debugging purposes only: TEMP: only do one card
    // if (cc != theCard) continue;
    
    // Initialize variables
    card rctCard;
    rctCard.setIdx(cc);
    
    // Debugging only: inject four hits at specific locations
    // int nHits = 4;
    // int arr_local_iEta[nHits] = {0, 1, 8, 9};
    // int arr_local_iPhi[nHits] = {0, 0, 10, 10};
    
    //    for (int j = 0; j < nHits; j++) {
    for (const auto& hit : ecalhits) {

      // Debugging only: do not use this block
      //  if (true) {

      // Check if the hit is in cards 0-35
      if ((getCrystal_iPhi(hit.position().phi()) <= getCard_iPhiMax(cc)) &&
      	  (getCrystal_iPhi(hit.position().phi()) >= getCard_iPhiMin(cc)) &&
      	  (getCrystal_iEta(hit.position().eta()) <= getCard_iEtaMax(cc)) &&
      	  (getCrystal_iEta(hit.position().eta()) >= getCard_iEtaMin(cc))) {
	
	// Get the crystal eta and phi, relative to the bottom left corner of the card 
	// (0 up to 17*5, 0 up to 4*5) 
	int local_iEta = getCrystal_local_iEta(hit.position().eta(), cc);
        int local_iPhi = getCrystal_local_iPhi(hit.position().phi(), cc);
      
	// int local_iEta = arr_local_iEta[j];
	// int local_iPhi = arr_local_iPhi[j];
	
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
	// 	    << regionNumber << std::endl;
	// std::cout << "inRegion_crystal_iEta, inRegion_crystal_iPhi (expecting 15x20), regionNumber (expecting 0-5): " << std::endl;
	// << inRegion_crystal_iEta << ", "
	// << inRegion_crystal_iPhi << ", "
	// << regionNumber << std::endl;
	  
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

	  // TEMP: add a dummy energy
	  // myLink.addCrystalE(inLink_crystal_iEta, inLink_crystal_iPhi, 10);
	  
	}

	if (cc == 28 && regionNumber == 3) {
	  std::cout << "Card: " << cc 
		    << ", hit (Eta, phi, et): "
		    << hit.position().eta() << ", " << hit.position().phi() << ", " << hit.et_uint() << ". "
		    << "(iEta, iPhi): " 
		    << getCrystal_iEta(hit.position().eta()) << ", "
		    << getCrystal_iPhi(hit.position().phi()) << ". "
		    << "(Crystal eta, phi): "
		    << getEta_fromCrystaliEta(getCrystal_iEta(hit.position().eta())) << ", "
		    << getPhi_fromCrystaliPhi(getCrystal_iPhi(hit.position().phi())) << ". "
		    // << "(Go back to iEta, iPhi in card: ) "
		    // << getCrystal_iEta_fromCardRegionInfo(cc, regionNumber, inRegion_tower_iEta, inLink_crystal_iEta) << ", "
		    // << getCrystal_iPhi_fromCardRegionInfo(cc, regionNumber, inRegion_tower_iPhi, inLink_crystal_iPhi) << "."
		  // << "local_iEta/iPhi: " << local_iEta << ", " << local_iPhi << ", " 
		  // << "inCard_tower_iEta/iPhi: " << inCard_tower_iEta     << ", " << inCard_tower_iPhi << ", "
		  // << "region/inRegion_tower_iEta/iPhi: " << regionNumber << ", " << inRegion_tower_iEta << ", "
		  // << inRegion_tower_iPhi << ", "
		  // << "inLink crystal iEta/iPhi: " << inLink_crystal_iEta << ", " << inLink_crystal_iPhi 
		    << std::endl;

	}
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
	// std::cout << "Real (eta, phi): " << eta << ", " << phi << "; " 
	// 	  << "iEta_tower_L1Card[" << jj << "][" << ii 
	// 	  << "][" << cc << "] = " << getTower_absEtaID(eta) << ", "
	// 	  << "iPhi_tower_L1Card[" << jj << "][" << ii
	// 	  << "][" << cc << "] = " << getTower_absPhiID(phi) << std::endl;
      }
    }
    

    //*******************************************************************
    //************* Do RCT geometry (HCAL) ******************************
    //*******************************************************************

    // Same idea as the ECAL RCT geometry, except we only care about the ET in towers 

    // Loop over hcal hits to get the HCAL towers.

    // TEMP: don't do HCAL hits
    //    for (int k = 0; k < 1; k++) {
    for (const auto& hit : hcalhits) {

      // TEMP: do not use this conditional
      //      if (true) {
      if (getCrystal_iPhi(hit.position().phi()) <= getCard_iPhiMax(cc) &&
          getCrystal_iPhi(hit.position().phi()) >= getCard_iPhiMin(cc) &&
          getCrystal_iEta(hit.position().eta()) <= getCard_iEtaMax(cc) &&
          getCrystal_iEta(hit.position().eta()) >= getCard_iEtaMin(cc) && hit.pt() > 0) {

      	// HCAL: Get the crystal eta and phi, relative to the bottom left corner of the card 
      	// (0 up to 17*5, 0 up to 4*5) 
      	int local_iEta = getCrystal_local_iEta(hit.position().eta(), cc);
      	int local_iPhi = getCrystal_local_iPhi(hit.position().phi(), cc);

	// int local_iEta = 0;
	// int local_iPhi = 0;
	
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

	// std::cout << "HCAL hit eta/phi : "  << hit.position().eta() << ", " << hit.position().phi() << ", region: " << regionNumber << ", "
	// 	  << "inRegion_tower_iEta and iPhi : " << inRegion_tower_iEta << ", " << inRegion_tower_iPhi 
	// 	  << ", energy " << hit.et_uint() 
	// 	  << std::endl;

	// Access the right HCAL region -> tower and increment the ET
	if (regionNumber < N_REGIONS_PER_CARD) {
	  towers3x4& myTowers3x4 = rctCard.getTowers3x4(regionNumber);
	  towerHCAL& myTower = myTowers3x4.getTowerHCAL(inRegion_tower_iEta, inRegion_tower_iPhi);

	  // TEMP: add a dummy energy
	  // myTower.addEt(1);
	  myTower.addEt(hit.et_uint());

	}
      }
    } // end of loop over hcal hits

    //*******************************************************************    
    //******* Within each ECAL region, read back the hits ***************
    //******************************************************************* 

    // Dummy for testing
    std::vector<Cluster> test_cluster;

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

		// if (uEnergy > 0) {
		//   std::cout << "energy>0: Accessing temporary array: " << (ref_iEta + iEta) << ", " << (ref_iPhi + iPhi) 
		// 	    << ", writing energy (uint:) " << uEnergy << std::endl;
		// }

		// Fill the 'temporary' array with a crystal object 
		temporary[ref_iEta + iEta][ref_iPhi + iPhi] = crystal(uEnergy);
	      }
	    } // end of loop over crystals

	    ///// Get HCAL tower ET
	    towerHCAL& myTower = myTowers.getTowerHCAL(iLinkEta, iLinkPhi);
	    std::cout << "towerEtHCAL element # " << (iLinkEta * TOWER_IN_PHI) + iLinkPhi << ": will fill with Et " << myTower.getEt() << std::endl;
	    towerEtHCAL[(iLinkEta * TOWER_IN_PHI) + iLinkPhi] = myTower.getEt();
	  }
	}
	
	// Iteratively find four clusters and remove them from 'temporary' as we go, and fill cluster_list
	for (int c = 0; c < N_CLUSTERS_PER_REGION; c++) {
	  Cluster newCluster = getClusterFromRegion3x4(temporary); // remove cluster from 'temporary' 
	  newCluster.setRegionIdx(idxRegion);                      // add the region number 
	  if (newCluster.clusterEnergy() > 0) {
	    cluster_list[cc].push_back(newCluster);                // do not push back 0-energy clusters
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
	  towerHCALCard[iEta][iPhi][cc] = tower_t(towerEtHCAL[i], 0, 0);
	  towerECALCard[iEta][iPhi][cc] = tower_t(towerEtECAL[i], 0, 0);
	}

      } // end of "if"

    } // end of loop over regions

    //-------------------------------------------// 
    // Stitching code (to be added)              //                                          
    //-------------------------------------------// 

    //-------------------------------------------//  
    // Sort the clusters:                        //
    // Take only the 8 highest ET clusters,      //
    // add leftover clusters back to ECAL towers //
    //-------------------------------------------//
    std::cout << "Card " << cc << ": unsorted ET: "; 
    for (int i = 0; i < 12; i++) {
      Cluster dummy = Cluster((ap_uint<12>) i+5, i%12, i%4, 0, 0, 0);
      test_cluster.push_back(dummy);
    }
    // START HERE: Replace cluster_list[cc] with test_cluster to do a test
    std::cout << "Up to twelve clusters going in: printing ET and tower (17x4):";
    if (cluster_list[cc].empty()) {
      std::cout << "The vector is empty; skip this step" << std::endl;
    }
    else {
      std::cout << "There will be " << cluster_list[cc].size() << " elements to loop over" << std::endl;
      for (unsigned int kk = 0; kk < cluster_list[cc].size(); ++kk) {
	std::cout << "Size is " << cluster_list[cc].size() << std::endl;
	Cluster c = cluster_list[cc][kk];
	std::cout << cluster_list[cc][kk].clusterEnergy() 
		  << " (" << cluster_list[cc][kk].towerEta() 
		  << ", " << cluster_list[cc][kk].towerPhi() << ") "; 
      }
      std::cout << std::endl;
      std::sort(cluster_list[cc].begin(), cluster_list[cc].end(), compareClusterET);
      std::cout << "NOMINALLY SORTED... ";
      for (unsigned int kk = 0; kk < cluster_list[cc].size(); ++kk) {
	std::cout << cluster_list[cc][kk].clusterEnergy() 
		  << " (" << cluster_list[cc][kk].towerEta()
		  << ", " << cluster_list[cc][kk].towerPhi() << ") ";
      }

      std::cout << "----" << std::endl;
      // If there are more than eight clusters, return the unused energy to the towers
      for (unsigned int kk = n_clusters_4link; kk < cluster_list[cc].size(); ++kk) {
	Cluster cExtra = cluster_list[cc][kk];
	if (cExtra.clusterEnergy() > 0) {
	  std::cout << "Extra cluster # " << kk << ": energy (GeV) " << cExtra.clusterEnergy()
		    << std::endl;
	  // Increment tower ET in towerECALCard[17][4] which uses RCT HLS geometry
	  // Get tower (eta, phi) (up to (17, 4)) in the RCT card
	  int whichTowerEtaInCard = ((cExtra.region() * TOWER_IN_ETA) + cExtra.towerEta());
	  int whichTowerPhiInCard = cExtra.towerPhi(); 
	  ap_uint<12> oldTowerEt = towerECALCard[whichTowerEtaInCard][whichTowerPhiInCard][cc].et();
	  ap_uint<12> newTowerEt = (oldTowerEt + cExtra.clusterEnergy());
	  ap_uint<3>  hoe        = towerECALCard[whichTowerEtaInCard][whichTowerPhiInCard][cc].hoe();
	  ap_uint<4>  satur      = towerECALCard[whichTowerEtaInCard][whichTowerPhiInCard][cc].satur();
	  towerECALCard[whichTowerEtaInCard][whichTowerPhiInCard][cc] = tower_t(newTowerEt, hoe, satur);
	  
	  std::cout << "... Adding to card eta (0-17) " << whichTowerEtaInCard
		    << ", card phi (0-4) " << whichTowerPhiInCard
		    << " old energy : " << oldTowerEt
		    << " new energy : " << newTowerEt
		    << std::endl;
	}
      }
      
      // Build the sorted cluster list IFF there were clusters 
      std::cout << "Building sorted cluster list: ";
      if (cluster_list[cc].empty()) {
	std::cout << "No clusters: do not build sorted cluster list" << std::endl;
      }
      else {
	// Save up to eight clusters: loop over cluster_list
	for (unsigned int kk = 0; kk < cluster_list[cc].size(); ++kk) {
	  if (kk >= n_clusters_4link) continue;
	  if (cluster_list[cc][kk].clusterEnergy() > 0) 
	    cluster_list_merged[cc].push_back(cluster_list[cc][kk]);
	}
	std::cout << "Sorted, up to 8 clusters: ";
	for (unsigned int kk = 0; kk < cluster_list_merged[cc].size(); ++kk) {
	  std::cout << cluster_list_merged[cc][kk].clusterEnergy() 
		    << " (" << cluster_list_merged[cc][kk].towerEta() 
		    << ", " << cluster_list_merged[cc][kk].towerPhi() << ") ";
	}
	std::cout << std::endl;
      }
    }
    // END HERE for str-replace testing

    //-------------------------------------------//
    // Write the L1 outputs                      //
    //-------------------------------------------//

    std::cout << "Sanity check: Card " << cc << ": SORTED ET: (if >0 clusters exist, print stuff)";
    
    // (May be deprecated:) If the cluster list isn't empty, write out the old CMSSW (by Cecile)
    // -style arrays energy_cluster_L1Card, crystalID_cluster_L1Card, and towerID_clusterL1Card
    // which treats (0, 0) as always the top left corner of the RCT card if you look at the
    // usual RCT diagram. (aka negative eta cards are NOT rotated)
    if (!cluster_list_merged[cc].empty()) {
      for (unsigned int jj = 0; jj < cluster_list_merged[cc].size(); ++jj) {
	Cluster c = cluster_list_merged[cc][jj];
	
	
	std::cout << c.clusterEnergy() << " (" << c.clusterEta() << ", " << c.clusterPhi() << ") "
		  << ", setting crystalID_cluster_L1Card[" << jj % n_links_card << "]["
		  << jj / n_links_card << "][" << cc << "] = " << getCrystalIDInTower_emulator(c.clusterEta(), c.clusterPhi(), cc) 
		  << ", setting towerID_cluster_L1Card[" << jj % n_links_card << "]["
		  << jj / n_links_card << "][" << cc << "] = " << getTowerID_emulator(c.towerEta(), c.towerPhi(), cc, c.region())
		  << std::endl;
	
	// Distribute (up to 12) clusters across 4 links
	energy_cluster_L1Card[jj % n_links_card][jj / n_links_card][cc] =
	  c.clusterEnergy();
	crystalID_cluster_L1Card[jj % n_links_card][jj / n_links_card][cc] = 
	  getCrystalIDInTower_emulator(c.clusterEta(), c.clusterPhi(), cc);
	towerID_cluster_L1Card[jj % n_links_card][jj / n_links_card][cc] =
	  getTowerID_emulator(c.towerEta(), c.towerPhi(), cc, c.region());
      }
    }
    
    std::cout << std::endl;

    // (May be deprecated:) Fill out ECAL_tower_L1Card and HCAL_tower_L1Card, which is in the style of the old CMSSW (by Cecile). The firmware 17x4 array treats the "bottom left" corner of the card as (0, 0)
    // (rotating the negative eta cards so that the endcap region is pointing up)
    // while the old CMSSW (by Cecile) treats the 4x17 array as always starting in the top left corner if we look
    // at the usual RCT diagram.
    std::cout << "Re-package towerECALCard into old CMSSW emulator geometry" << std::endl;
    for (int i = 0; i < n_towers_cardEta; i++) {
      for (int j = 0; j < n_towers_cardPhi; j++ ) {
        
    	// n.b. L1 output is 4*17*36, hence [j][i] instead of [i][j] on the L.H.S.    
    	if ((cc % 2) == 1) { // if cc is odd (positive eta)
    	  ECAL_tower_L1Card[j][i][cc] = towerECALCard[i][j][cc].et();
    	  HCAL_tower_L1Card[j][i][cc] = towerHCALCard[i][j][cc].et();
    	}
    	else {  // if cc is even (negative eta), we need to rotate the coordinates
    	  //	  std::cout << "writing to [" << (3-j) << "][" << (16-i) << "] ";
    	  ECAL_tower_L1Card[3-j][16-i][cc] = towerECALCard[i][j][cc].et();
    	  HCAL_tower_L1Card[3-j][16-i][cc] = towerHCALCard[i][j][cc].et();
    	}
      }
    }
    std::cout << std::endl;


    ///////////////////////////////////////////////////////////
    // Produce output collections for the event display
    //////////////////////////////////////////////////////////
    for (auto & c : cluster_list_merged[cc]) {                                                      
      std::cout << c.clusterEnergy() << " (" << c.clusterEta() << ", " << c.clusterPhi() << ") ";    

      float realEta = getEta_fromCrystaliEta( getCrystal_iEta_fromCardRegionInfo(cc, c.region(), c.towerEta(), c.clusterEta()) );
      float realPhi = getPhi_fromCrystaliPhi( getCrystal_iPhi_fromCardRegionInfo(cc, c.region(), c.towerPhi(), c.clusterPhi()) );
      
      std::cout << "Real eta, phi " << " (" << realEta << ", " << realPhi << ") ";
      reco::Candidate::PolarLorentzVector p4calibrated(c.clusterEnergy()/8.0,
						       realEta,
						       realPhi,
						       0.);
      
      // Constructor definition at: https://github.com/cms-l1t-offline/cmssw/blob/l1t-phase2-v3.3.11/DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h#L34
      // showerShape_cluster_L1Card: depends on passes_ss(cpt, 2x5/5x5)
      // equivalent to cshowershape_. equivalent to showerShape_cluster_L2Card
      // is_ss
      bool is_ss = passes_ss(c.clusterEnergy()/8.0,
			     (c.getEt2x5() / c.getEt5x5()) );
      // showerShapeLooseTk_cluster_L1Card: depends on passes_looseTkss(cpt, 2x5/5x5)
      // equivalent to cphotonshowershape_. equivalent to photonShowerShape_cluster_L2Card
      // is_looseTkss
      bool is_looseTkss = passes_looseTkss(c.clusterEnergy()/8.0,
					   (c.getEt2x5() / c.getEt5x5() ));
      std::cout << "is_ss: " << is_ss << ", is_looseTkss: " << is_looseTkss << std::endl;
      l1tp2::CaloCrystalCluster cluster(p4calibrated,
				        c.clusterEnergy()/8.0, // (convert to float)
					0,  // float h over e
					0,  // float iso
					0,  // DetId seedCrystal 
					0,  // puCorrPt
					c.getBrems(), // 0, 1, or 2 (as computed in firmware)
					0,            // et2x2 (not calculated)
					c.getEt2x5(), // et2x5 (as computed in firmware, save float)
					0,            // et3x5 (not calculated)
					c.getEt5x5(), // et5x5 (as computed in firmware, save float)
					is_ss, // standalone WP
					false, // electronWP98
					false, // photonWP80
					false, // electronWP90
					is_looseTkss, // looseL1TkMatchWP
					false  // stage2effMatch
					);
      
      L1EGXtalClusters->push_back(cluster);
    }  // end of loop over clusters                                                                                                           
    // Produce output tower collections
    for (int ii = 0; ii < n_towers_cardPhi; ++ii) { // 4 towers per card in phi
      for (int jj = 0; jj < n_towers_cardEta; ++jj) { // 17 towers per card in eta

	l1tp2::CaloTower l1CaloTower;
	// Divide ET by 8.0 to convert to GeV
	l1CaloTower.setEcalTowerEt(ECAL_tower_L1Card[ii][jj][cc]/8.0);
	// HCAL TPGs encoded ET: multiply by the LSB (0.5) to convert to GeV
	float hcalLSB = 0.5;
	l1CaloTower.setHcalTowerEt(HCAL_tower_L1Card[ii][jj][cc] * hcalLSB);
	l1CaloTower.setTowerIEta(getToweriEta_fromAbsID(iEta_tower_L1Card[ii][jj][cc]));
        l1CaloTower.setTowerIPhi(getToweriPhi_fromAbsID(iPhi_tower_L1Card[ii][jj][cc]));
        l1CaloTower.setTowerEta(getTowerEta_fromAbsID(iEta_tower_L1Card[ii][jj][cc]));
        l1CaloTower.setTowerPhi(getTowerPhi_fromAbsID(iPhi_tower_L1Card[ii][jj][cc]));
	
	L1CaloTowers->push_back(l1CaloTower);
      }
    }


  } // end of loop over cards

  iEvent.put(std::move(L1EGXtalClusters));
  iEvent.put(std::move(L1CaloTowers));
  std::cout << "Finished producing RCT emulator for this event" << std::endl;
  
}

bool Phase2L1CaloEGammaEmulator::passes_ss(float pt, float ss) {
  bool is_ss = ((c0_ss + c1_ss * std::exp(-c2_ss * pt)) <= ss);
  if (pt > 130)
    is_ss = true;
  return is_ss;
}

bool Phase2L1CaloEGammaEmulator::passes_looseTkss(float pt, float ss) {
  bool is_ss = ((e0_looseTkss - e1_looseTkss * std::exp(-e2_looseTkss * pt)) <= ss);
  if (pt > 130)
    is_ss = true;
  return is_ss;
}


////////////////////////////////////////////////////////////////////////// 

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2L1CaloEGammaEmulator);

