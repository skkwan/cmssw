#include <bitset>
#include <cassert>
#include "DataFormats/L1TCalorimeterPhase2/interface/DigitizedClusterCorrelator.h"


/*
 * Test constructor with digitized inputs
 */
bool passItemizedDigitizedConstructorTest(void) {
    l1tp2::DigitizedClusterCorrelator dc = l1tp2::DigitizedClusterCorrelator((ap_uint<12>) 0x1,
                                                                      (ap_uint<6>) 0x1, (ap_uint<5>) 0x1, 
                                                                      (ap_uint<3>) 0x1, (ap_uint<3>) 0x1, 
                                                                      (ap_uint<4>) 0x1, (ap_uint<3>) 0x1,
                                                                      (ap_uint<6>) 0x1, (ap_uint<5>) 0x1,
                                                                      (ap_uint<1>) 0x1, (ap_uint<1>) 0x1,
                                                                      (ap_uint<1>) 0x1);
    dc.print("in passItemizedDigitizedConstructorTest");
    return ((dc.pt() == 0x1) && (dc.eta() == 0x1) && (dc.phi() == 0x1) && 
            (dc.etaCr() == 0x1) && (dc.phiCr() == 0x1) && 
            (dc.hoe() == 0x1) && (dc.iso() == 0x1) && (dc.fb() == 0x1) && 
            (dc.timing() == 0x1) && (dc.shape() == 0x1) && (dc.brems() == 0x1) &&
            (dc.cardNumber() == 0x1) && 
            dc.passNullBitsCheck());

}

/*
 * Test constructor with positive float inputs
 */
bool passMaxValuesConstructorTest(void) {
    float maxPt = (4095 * 0.5) + 1; // test saturated 
    unsigned int maxiEta = 33;
    unsigned int maxiPhi = 23;  // max in GCT card unique area
    unsigned int maxiCrystal = 4; // same for eta and phi
    int nGCTCard = 2;
    l1tp2::DigitizedClusterCorrelator dc = l1tp2::DigitizedClusterCorrelator(maxPt, maxiEta, maxiPhi,
                                                                             maxiCrystal, maxiCrystal,
                                                                             0, true,
                                                                             0, 0,
                                                                             true, 0,
                                                                             nGCTCard);
    dc.print("in passFloatsConstructorTest");
    dc.printPt();
    dc.printEta();
    dc.printPhi();
    dc.printEtaCr();
    dc.printPhiCr();
    return ((dc.pt() == 0xFFF) && (dc.eta() == maxiEta) && (dc.phi() == maxiPhi) &&
            (dc.etaCr() == maxiCrystal) && (dc.phiCr() == maxiCrystal) &&
            (dc.iso() == true) &&
            (dc.shape() == true) &&
            dc.passNullBitsCheck() &&
            (dc.cardNumber() == nGCTCard)
            );
}



int main() {
    if (!passItemizedDigitizedConstructorTest()) {
        std::cout << "Fails itemized digitized constructor test!\n" << std::endl;
        return 1;
    }
    if (!passMaxValuesConstructorTest()) {
        std::cout << "Fails max values constructor test!\n" << std::endl;
        return 1;
    }
    return 0;
}