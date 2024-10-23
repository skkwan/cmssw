#ifndef PHASE_2_L1_CALO_BARREL_TO_CORRELATOR
#define PHASE_2_L1_CALO_BARREL_TO_CORRELATOR

#include "DataFormats/L1TCalorimeterPhase2/interface/GCTBarrelDigiClusterToCorrLayer1.h"
#include "L1Trigger/L1CaloTrigger/interface/Phase2L1CaloEGammaUtils.h"

/* 
 * Comparators for sorting
 */
inline bool p2eg::compareBarrelDigiClusterCorrelatorET(const l1tp2::GCTBarrelDigiClusterToCorrLayer1& lhs, const l1tp2::GCTBarrelDigiClusterToCorrLayer1& rhs) {
    return (lhs.ptFloat() > rhs.ptFloat());
}

/*
 *
*/
void p2eg::sortAndPadSLR(l1tp2::GCTBarrelDigiClusterToCorrLayer1Collection &thisSLR) {
    // input is a vector and can be sorted
    std::sort(thisSLR.begin(), thisSLR.end(), p2eg::compareBarrelDigiClusterCorrelatorET);
    int nClusters = thisSLR.size();
    std::cout << ">>> p2eg::sortAndPadSLR: Total entries in SLR: " << nClusters << std::endl;
    // if fewer than designated number of clusters, pad with zeros
    if (nClusters < p2eg::N_EG_CLUSTERS_PER_SLR) {
        // do padding. if size == 2, push back four clusters
        for (int i = 0; i < (p2eg::N_EG_CLUSTERS_PER_SLR - nClusters); i++) {
            l1tp2::GCTBarrelDigiClusterToCorrLayer1 zeroCluster;
            thisSLR.push_back(zeroCluster);
        }
    }
    std::cout << ">>> p2eg::sortAndPadSLR: After push-back: Total entries in SLR: " << thisSLR.size() << std::endl;

}


#endif