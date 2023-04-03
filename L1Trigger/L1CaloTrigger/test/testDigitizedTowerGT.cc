#include <bitset>
#include <cassert>
#include "DataFormats/L1TCalorimeterPhase2/interface/DigitizedTowerGT.h"

/*
 * Test basic constructor from digitized inputs
 */
int passItemizedDigitizedConstructorTest(void) {
    ap_uint<10> pt = 0x1;
    ap_uint<4> hoe = 0x1;
    ap_uint<2> fb = 0x1;
    unsigned int nGCTCard = 1;
    unsigned int nFiber = 1;
    unsigned int nTower = 1;
    l1tp2::DigitizedTowerGT towGT = l1tp2::DigitizedTowerGT(pt, hoe, fb, nGCTCard, nFiber, nTower);
    return ((towGT.et() == 0x1) &&
            (towGT.hoe() == 0x1) && (towGT.fb() == 0x1) && 
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
    l1tp2::DigitizedTowerGT towGT = l1tp2::DigitizedTowerGT(pt_f, hoe, fb, nGCTCard, nFiber, nTower);
    return ((towGT.et() == 0x3FF) && (towGT.hoe() = 0xF) && (towGT.fb() == 0x3) && towGT.hasValidIndices());
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