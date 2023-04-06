#include <bitset>
#include <cassert>
#include "DataFormats/L1TCalorimeterPhase2/interface/DigitizedTowerCorrelator.h"

/*
 * Test basic constructor from digitized inputs
 */
int passItemizedDigitizedConstructorTest(void) {
    bool isFullyDigitizedInputs = true;
    ap_uint<10> pt = 0x1;
    ap_uint<4> hoe = 0x1;
    ap_uint<2> fb = 0x1;
    unsigned int nGCTCard = 1;
    unsigned int nFiber = 1;
    unsigned int nTower = 1;
    l1tp2::DigitizedTowerCorrelator tow = l1tp2::DigitizedTowerCorrelator(pt, hoe, fb, nGCTCard, nFiber, nTower, isFullyDigitizedInputs);
    return ((tow.et() == 0x1) &&
            (tow.hoe() == 0x1) && (tow.fb() == 0x1) && 
            (nGCTCard == 1) && (nFiber == 1) && (nTower == 1)
            );
}

/*
 * Positive saturated values constructor
 */
int passConstructorTest(void) {
    float pt_f = 512; // max value is (2^10 - 1) *  0.5 = 511.5
    ap_uint<4> hoe = 0xF;
    ap_uint<2> fb = 0x3;
    unsigned int nGCTCard = 2;
    unsigned int nFiber = 47;
    unsigned int nTower = 16;
    l1tp2::DigitizedTowerCorrelator tow = l1tp2::DigitizedTowerCorrelator(pt_f, hoe, fb, nGCTCard, nFiber, nTower);
    return ((tow.et() == 0x3FF) && (tow.hoe() = 0xF) && (tow.fb() == 0x3) && tow.hasValidIndices());
}

int main() {

    if (!passItemizedDigitizedConstructorTest()) {
        std::cout << "Fails itemized digitized constructor test!\n" << std::endl;
        return 1;
    }

    if (!passConstructorTest()) {
        std::cout << "Fails positive constructor test!\n" << std::endl;
        return 1;
    }

    return 0;
}