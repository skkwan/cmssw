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
                                                                      (ap_uint<1>) 0x1, (ap_uint<1>) 0x1);
    dc.print("in passItemizedDigitizedConstructorTest");
    return ((dc.pt() == 0x1) && (dc.eta() == 0x1) && (dc.phi() == 0x1) && 
            (dc.etaCr() == 0x1) && (dc.phiCr() == 0x1) && 
            (dc.hoe() == 0x1) && (dc.iso() == 0x1) && (dc.fb() == 0x1) && 
            (dc.timing() == 0x1) && (dc.shape() == 0x1) && (dc.brems() == 0x1) &&
             dc.passNullBitsCheck());

}

/*
 * Test constructor with float inputs
 */
bool passFloatsConstructorTest(void) {
    float dummyPt = 4.0;
    l1tp2::DigitizedClusterCorrelator dc = l1tp2::DigitizedClusterCorrelator(dummyPt);
    dc.print("in passFloatsConstructorTest");
    return (dc.pt() == 0x8);
    // && (dc.eta() == 0x1) && (dc.phi() == 0x1) && 
    //         (dc.etaCr() == 0x1) && (dc.phiCr() == 0x1) && 
    //         (dc.hoe() == 0x1) && (dc.iso() == 0x1) && (dc.fb() == 0x1) && 
    //         (dc.timing() == 0x1) && (dc.shape() == 0x1) && (dc.brems() == 0x1) &&
    //          dc.passNullBitsCheck());

}



int main() {
    if (!passItemizedDigitizedConstructorTest()) {
        std::cout << "Fails itemized digitized constructor test!\n" << std::endl;
        return 1;
    }
    if (!passFloatsConstructorTest()) {
        std::cout << "Fails floats constructor test!\n" << std::endl;
        return 1;
    }
    return 0;
}