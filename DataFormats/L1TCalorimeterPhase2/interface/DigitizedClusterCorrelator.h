#ifndef DataFormats_L1TCalorimeterPhase2_DigitizedClusterCorrelator_h
#define DataFormats_L1TCalorimeterPhase2_DigitizedClusterCorrelator_h

#include <ap_int.h>
#include <vector>


namespace l1tp2 {

    class DigitizedClusterCorrelator {

        private:
            ap_uint<64> clusterData;

            const float LSB_PT = 0.5; // 0.5 GeV
            const float LSB_ETA = 1.4/17; 
            const float LSB_PHI = M_PI/72;
            const float LSB_ETA_CR = 1.4/17/5;
            const float LSB_PHI_CR = M_PI/72/5; 

            // Private member functions to perform digitization 
            ap_uint<12> digitizePt(float pt_f) {
                return (ap_uint<12>) (pt_f / LSB_PT);
            }


        public:
            DigitizedClusterCorrelator() {
                clusterData = 0x0;
            }

            DigitizedClusterCorrelator(ap_uint<64> data) {
                clusterData = data;
            }

            // Constructor from digitized inputs
            DigitizedClusterCorrelator(ap_uint<12> pt, ap_uint<6> eta, ap_uint<5> phi,
                                       ap_uint<3> etaCr, ap_uint<3> phiCr,
                                       ap_uint<4> hoe, ap_uint<3> iso, ap_uint<6> fb, ap_uint<5> timing,
                                       ap_uint<1> shape, ap_uint<1> brems) {

                                        clusterData = ((ap_uint<64>) pt) | (((ap_uint<64>) eta) << 12) | (((ap_uint<64>) phi) << 18) |
                                                      (((ap_uint<64>) etaCr) << 23) | (((ap_uint<64>) phiCr) << 26) |
                                                      (((ap_uint<64>) hoe) << 29) | (((ap_uint<64>) iso << 33)) | 
                                                      (((ap_uint<64>) fb << 36)) | (((ap_uint<64>) timing << 42)) |
                                                      (((ap_uint<64>) shape << 47)) | (((ap_uint<64>) brems << 48));
                                       }

            // To-do: constructor from float inputs that will perform digitization
            DigitizedClusterCorrelator(float pt_f) {
                clusterData = ((ap_uint<64>) digitizePt(pt_f)); 
            }

            ap_uint<64> data() const { return clusterData; }

            // Other getters
            ap_uint<12> pt() { return (clusterData & 0xFFF); }
            ap_uint<6> eta() { return ((clusterData >> 12) & 0x3F); }  // (six 1's) 0b111111 = 0x3F
            ap_uint<5> phi() { return ((clusterData >> 18) & 0x1F); }  // (five 1's) 0b11111 = 0x1F
            ap_uint<3> etaCr() { return ((clusterData >> 23) & 0x7); }  // (three 1's) 0b111 = 0x7
            ap_uint<3> phiCr() { return ((clusterData >> 26) & 0x7); }
            ap_uint<4> hoe() { return ((clusterData >> 29) & 0xF); } // (four 1's) 0b1111 = 0xF 
            ap_uint<3> iso() { return ((clusterData >> 33) & 0x7); }
            ap_uint<6> fb() { return ((clusterData >> 36) & 0x3F); }
            ap_uint<5> timing() { return ((clusterData >> 42) & 0x1F); }
            ap_uint<1> shape() { return ((clusterData >> 47) & 0x1); }
            ap_uint<1> brems() { return ((clusterData >> 48) & 0x1); }

            const int unusedBitsStart() { return 49; } // unused bits start at bit 49

            // Prints
            void print(const std::string location);
            void printEta(void) {  std::cout << "eta: " << std::bitset<6>{eta()} << std::endl; };
            void printPhi(void) {  std::cout << "phi: " << std::bitset<5>{phi()} << std::endl; };
            void printEtaCr(void) { std::cout << "etaCr: " << std::bitset<3>{etaCr()} << std::endl; };
            void printPhiCr(void) { std::cout << "etaCr: " << std::bitset<3>{phiCr()} << std::endl; };
            void printHoE(void) { std::cout << "hoe: " << std::bitset<4>{hoe()} << std::endl; };
    

            // Other checks
            bool passNullBitsCheck(void)  {
                return ((data() >> unusedBitsStart()) == 0x0);
            }

    };

    // Collection typedef
    typedef std::vector<l1tp2::DigitizedClusterCorrelator> DigitizedClusterCorrelatorCollection;

} // namespace l1tp2

#endif