// -*- C++ -*-
//
// Package:    L1Trigger/L1CaloTrigger
// Class:      L1TCaloBarrelToCorrelator
//
/*
 Description: Creates digitized EGamma and ParticleFlow clusters to be sent to correlator. 

 Implementation: To be run together with Phase2L1CaloEGammaEmulator.
*/


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1TCalorimeterPhase2/interface/CaloPFCluster.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloPFDigiClusterToCorrLayer1.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/DigitizedClusterCorrelator.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/GCTBarrelDigiClusterToCorrLayer1.h"


#include <ap_int.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdio>
#include "L1Trigger/L1CaloTrigger/interface/Phase2L1CaloBarrelToCorrelator.h"
#include "L1Trigger/L1CaloTrigger/interface/Phase2L1CaloEGammaUtils.h"

//
// class declaration
//

class Phase2GCTBarrelToCorrelatorLayer1 : public edm::stream::EDProducer<> {
public:
  explicit Phase2GCTBarrelToCorrelatorLayer1(const edm::ParameterSet&);
  ~Phase2GCTBarrelToCorrelatorLayer1() override = default;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------
  const edm::EDGetTokenT<l1tp2::DigitizedClusterCorrelatorCollection> digiInputClusterToken_;
  const edm::EDGetTokenT<l1tp2::CaloPFClusterCollection> caloPFClustersSrc_;
};

//
// constructors and destructor
//
Phase2GCTBarrelToCorrelatorLayer1::Phase2GCTBarrelToCorrelatorLayer1(const edm::ParameterSet& iConfig) 
  // gctClustersSrc_(consumes<l1tp2::CaloCrystalClusterCollection >(cfg.getParameter<edm::InputTag>("gctClusters"))),
  : digiInputClusterToken_(consumes<l1tp2::DigitizedClusterCorrelatorCollection>(iConfig.getParameter<edm::InputTag>("gctDigiClustersInput"))),
    caloPFClustersSrc_(consumes<l1tp2::CaloPFClusterCollection>(iConfig.getParameter<edm::InputTag>("PFclusters")))
{
  
  produces<l1tp2::GCTBarrelDigiClusterToCorrLayer1CollectionFullDetector>("GCTBarrelDigiClustersToCorrLayer1");
  produces<l1tp2::CaloPFDigiClusterToCorrLayer1CollectionFullDetector>("CaloPFDigiClusterToCorrLayer1");
}

void Phase2GCTBarrelToCorrelatorLayer1::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace edm;

    //***************************************************//
    // Get the GCT digitized clusters and PF clusters
    //***************************************************//
    edm::Handle<l1tp2::DigitizedClusterCorrelatorCollection> inputGCTBarrelClusters;
    iEvent.getByToken(digiInputClusterToken_, inputGCTBarrelClusters);

    edm::Handle<l1tp2::CaloPFClusterCollection> inputPFClusters;
    iEvent.getByToken(caloPFClustersSrc_, inputPFClusters);
    
    //***************************************************//
    // Initialize outputs
    //***************************************************//

    // EGamma cluster output
    auto outputClustersFromBarrel = std::make_unique<l1tp2::GCTBarrelDigiClusterToCorrLayer1CollectionFullDetector>();
    // PF cluster output
    auto outputPFClusters = std::make_unique<l1tp2::CaloPFDigiClusterToCorrLayer1CollectionFullDetector>();

    // EG Clusters output by GCT SLR (duplicates included)
    l1tp2::GCTBarrelDigiClusterToCorrLayer1Collection out_eg_GCT1_SLR1_posEta;
    l1tp2::GCTBarrelDigiClusterToCorrLayer1Collection out_eg_GCT1_SLR1_negEta;
    l1tp2::GCTBarrelDigiClusterToCorrLayer1Collection out_eg_GCT1_SLR3_posEta;
    l1tp2::GCTBarrelDigiClusterToCorrLayer1Collection out_eg_GCT1_SLR3_negEta;
    l1tp2::GCTBarrelDigiClusterToCorrLayer1Collection out_eg_GCT2_SLR1_posEta;
    l1tp2::GCTBarrelDigiClusterToCorrLayer1Collection out_eg_GCT2_SLR1_negEta;
    l1tp2::GCTBarrelDigiClusterToCorrLayer1Collection out_eg_GCT2_SLR3_posEta;
    l1tp2::GCTBarrelDigiClusterToCorrLayer1Collection out_eg_GCT2_SLR3_negEta;
    l1tp2::GCTBarrelDigiClusterToCorrLayer1Collection out_eg_GCT3_SLR1_posEta;
    l1tp2::GCTBarrelDigiClusterToCorrLayer1Collection out_eg_GCT3_SLR1_negEta;
    l1tp2::GCTBarrelDigiClusterToCorrLayer1Collection out_eg_GCT3_SLR3_posEta;
    l1tp2::GCTBarrelDigiClusterToCorrLayer1Collection out_eg_GCT3_SLR3_negEta;

    // PF Clusters output by GCT SLR (duplicates included)
    l1tp2::CaloPFDigiClusterToCorrLayer1Collection out_pf_GCT1_SLR1_posEta;
    l1tp2::CaloPFDigiClusterToCorrLayer1Collection out_pf_GCT1_SLR1_negEta;

    l1tp2::CaloPFDigiClusterToCorrLayer1Collection out_pf_GCT1_SLR3_posEta;
    l1tp2::CaloPFDigiClusterToCorrLayer1Collection out_pf_GCT1_SLR3_negEta;

    l1tp2::CaloPFDigiClusterToCorrLayer1Collection out_pf_GCT2_SLR1_posEta;
    l1tp2::CaloPFDigiClusterToCorrLayer1Collection out_pf_GCT2_SLR1_negEta;

    l1tp2::CaloPFDigiClusterToCorrLayer1Collection out_pf_GCT2_SLR3_posEta;
    l1tp2::CaloPFDigiClusterToCorrLayer1Collection out_pf_GCT2_SLR3_negEta;

    l1tp2::CaloPFDigiClusterToCorrLayer1Collection out_pf_GCT3_SLR1_posEta;
    l1tp2::CaloPFDigiClusterToCorrLayer1Collection out_pf_GCT3_SLR1_negEta;

    l1tp2::CaloPFDigiClusterToCorrLayer1Collection out_pf_GCT3_SLR3_posEta;
    l1tp2::CaloPFDigiClusterToCorrLayer1Collection out_pf_GCT3_SLR3_negEta;

    //***************************************************//
    // Loop over the regions: in order: GCT1 SLR1, GCT1 SLR3, GCT2 SLR1, GCT2 SLR3, GCT3 SLR1, GCT3SLR3
    //***************************************************//

    int nRegions = 6;
    float regionCentersInDegrees[nRegions] = {10.0, 70.0, 130.0, -170.0, -110.0, -50.0};

    for (int i = 0; i < nRegions; i++) {

        // EG Clusters
        for (auto &clusterIn : *inputGCTBarrelClusters.product()) {
            
            l1tp2::GCTBarrelDigiClusterToCorrLayer1 clusterOut;

            // Check if this cluster falls into each card 
            float clusterRealPhiAsDegree = clusterIn.realPhi() * 180/M_PI; 
            float regionLowerPhiBound = regionCentersInDegrees[i] - 60;
            float regionUpperPhiBound = regionCentersInDegrees[i] + 60;
            if ((clusterRealPhiAsDegree > regionLowerPhiBound) && (clusterRealPhiAsDegree < regionUpperPhiBound)) {
                // Go from real phi to an index in the SLR
                // Calculate the distance in phi from the center of the region
                float phiDifference = clusterRealPhiAsDegree - regionCentersInDegrees[i];
                int iPhiCrystalDifference = (int) std::round(phiDifference);

                // For eta, the eta is already digitized, just needs to be converted from [0, +2*17*5) to [-17*5, +17*5)
                int temp_iEta_signed = clusterIn.eta() - (p2eg::CRYSTALS_IN_TOWER_ETA * p2eg::n_towers_per_link);
                // Default value: for clusters in positive eta, values go from 0, 1, 2, 3, 4
                int iEta = temp_iEta_signed; 
                // If cluster is in negative eta, instead of from -5, -4, -3, -2, -1, we want 4, 3, 2, 1, 0
                if (temp_iEta_signed < 0) {
                    // If in negative eta, convert to an absolute value, with 0 being the crystal nearest real eta = 0
                    iEta = std::abs(temp_iEta_signed + 1);
                }

                // Initialize the new cluster
                l1tp2::GCTBarrelDigiClusterToCorrLayer1 clusterOut = l1tp2::GCTBarrelDigiClusterToCorrLayer1(
                clusterIn.pt(),
                iEta,
                iPhiCrystalDifference,
                clusterIn.hoe(),
                clusterIn.hoeFlag(),
                clusterIn.iso(),
                clusterIn.isoFlags(),
                clusterIn.fb(),
                clusterIn.timing(),
                clusterIn.shapeFlags(),
                clusterIn.brems()
                );

                
                if (i == 0)      { if (temp_iEta_signed < 0) { out_eg_GCT1_SLR1_negEta.push_back(clusterOut); } else { out_eg_GCT1_SLR1_posEta.push_back(clusterOut);} }
                else if (i == 1) { if (temp_iEta_signed < 0) { out_eg_GCT1_SLR3_negEta.push_back(clusterOut); } else { out_eg_GCT1_SLR3_posEta.push_back(clusterOut);} }
                else if (i == 2) { if (temp_iEta_signed < 0) { out_eg_GCT2_SLR1_negEta.push_back(clusterOut); } else { out_eg_GCT2_SLR1_posEta.push_back(clusterOut);} }
                else if (i == 3) { if (temp_iEta_signed < 0) { out_eg_GCT2_SLR3_negEta.push_back(clusterOut); } else { out_eg_GCT2_SLR3_posEta.push_back(clusterOut);} }
                else if (i == 4) { if (temp_iEta_signed < 0) { out_eg_GCT3_SLR1_negEta.push_back(clusterOut); } else { out_eg_GCT3_SLR1_posEta.push_back(clusterOut);} }
                else if (i == 5) { if (temp_iEta_signed < 0) { out_eg_GCT3_SLR3_negEta.push_back(clusterOut); } else { out_eg_GCT3_SLR3_posEta.push_back(clusterOut);} }

            }
        }

        // Repeat for PF Clusters
        for (auto &pfIn : *inputPFClusters.product()) {
            l1tp2::CaloPFDigiClusterToCorrLayer1 pfOut;

            // Check if this cluster falls into each card 
            float clusterRealPhiAsDegree = pfIn.clusterPhi() * 180/M_PI; 
            float regionLowerPhiBound = regionCentersInDegrees[i] - 60;
            float regionUpperPhiBound = regionCentersInDegrees[i] + 60;
            if ((clusterRealPhiAsDegree > regionLowerPhiBound) && (clusterRealPhiAsDegree < regionUpperPhiBound)) {
                // Go from real phi to an index in the SLR
                // Calculate the distance in phi from the center of the region
                float phiDifference = clusterRealPhiAsDegree - regionCentersInDegrees[i];
                int iPhiCrystalDifference = (int) phiDifference;  // round down

                // For PFClusters, the method clusterEta returns a float, so we need to digitize this
                float eta_LSB = p2eg::ECAL_eta_range / (p2eg::N_GCTTOWERS_FIBER * p2eg::CRYSTALS_IN_TOWER_ETA);
                int temp_iEta_signed = pfIn.clusterEta() / eta_LSB;
                // Default value (for positive eta)
                temp iEta = temp_iEta_signed;
                // If cluster is in negative eta, instead of from -5, -4, -3, -2, -1, we want 4, 3, 2, 1, 0
                if (temp_iEta_signed < 0) {
                    // If in negative eta, convert to an absolute value, with 0 being the crystal nearest real eta = 0
                    iEta = std::abs(temp_iEta_signed + 1);
                }

                // Initialize the new cluster
                l1tp2::CaloPFDigiClusterToCorrLayer1 pfOut = l1tp2::CaloPFDigiClusterToCorrLayer1(
                    pfIn.clusterEt() / p2eg::ECAL_LSB,  // convert to integer
                    iEta,
                    iPhiCrystalDifference,
                    0  // no HoE value in PF Cluster
                );

                
                if (i == 0)      { if (temp_iEta_signed < 0) { out_pf_GCT1_SLR1_negEta.push_back(pfOut); } else { out_pf_GCT1_SLR1_posEta.push_back(pfOut); }}
                else if (i == 1) { if (temp_iEta_signed < 0) { out_pf_GCT1_SLR3_negEta.push_back(pfOut); } else { out_pf_GCT1_SLR3_posEta.push_back(pfOut); }}
                else if (i == 2) { if (temp_iEta_signed < 0) { out_pf_GCT2_SLR1_negEta.push_back(pfOut); } else { out_pf_GCT2_SLR1_posEta.push_back(pfOut); }}
                else if (i == 3) { if (temp_iEta_signed < 0) { out_pf_GCT2_SLR3_negEta.push_back(pfOut); } else { out_pf_GCT2_SLR3_posEta.push_back(pfOut); }}
                else if (i == 4) { if (temp_iEta_signed < 0) { out_pf_GCT3_SLR1_negEta.push_back(pfOut); } else { out_pf_GCT3_SLR1_posEta.push_back(pfOut); }}
                else if (i == 5) { if (temp_iEta_signed < 0) { out_pf_GCT3_SLR3_negEta.push_back(pfOut); } else { out_pf_GCT3_SLR3_posEta.push_back(pfOut); }}
            }
        }
    }

    // Need to push these back in a specific order
    outputClustersFromBarrel->push_back(out_eg_GCT1_SLR1_posEta);
    outputClustersFromBarrel->push_back(out_eg_GCT1_SLR1_negEta);
    outputClustersFromBarrel->push_back(out_eg_GCT1_SLR3_posEta);
    outputClustersFromBarrel->push_back(out_eg_GCT1_SLR3_negEta);
    outputClustersFromBarrel->push_back(out_eg_GCT2_SLR1_posEta);
    outputClustersFromBarrel->push_back(out_eg_GCT2_SLR1_negEta);
    outputClustersFromBarrel->push_back(out_eg_GCT2_SLR3_posEta);
    outputClustersFromBarrel->push_back(out_eg_GCT2_SLR3_negEta);
    outputClustersFromBarrel->push_back(out_eg_GCT3_SLR1_posEta);
    outputClustersFromBarrel->push_back(out_eg_GCT3_SLR1_negEta);
    outputClustersFromBarrel->push_back(out_eg_GCT3_SLR3_posEta);
    outputClustersFromBarrel->push_back(out_eg_GCT3_SLR3_negEta);
    
    outputPFClusters->push_back(out_pf_GCT1_SLR1_posEta);
    outputPFClusters->push_back(out_pf_GCT1_SLR1_negEta);
    outputPFClusters->push_back(out_pf_GCT1_SLR3_posEta);
    outputPFClusters->push_back(out_pf_GCT1_SLR3_negEta);
    outputPFClusters->push_back(out_pf_GCT2_SLR1_posEta);
    outputPFClusters->push_back(out_pf_GCT2_SLR1_negEta);
    outputPFClusters->push_back(out_pf_GCT2_SLR3_posEta);
    outputPFClusters->push_back(out_pf_GCT2_SLR3_negEta);
    outputPFClusters->push_back(out_pf_GCT3_SLR1_posEta);
    outputPFClusters->push_back(out_pf_GCT3_SLR1_negEta);
    outputPFClusters->push_back(out_pf_GCT3_SLR3_posEta);
    outputPFClusters->push_back(out_pf_GCT3_SLR3_negEta);

    iEvent.put(std::move(outputClustersFromBarrel), "GCTBarrelDigiClustersToCorrLayer1");
    iEvent.put(std::move(outputPFClusters), "CaloPFDigiClusterToCorrLayer1");
}

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2GCTBarrelToCorrelatorLayer1);
