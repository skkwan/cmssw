#ifndef DataFormats_L1TCalorimeterPhase2_DigitizedClusterCorrelator_h
#define DataFormats_L1TCalorimeterPhase2_DigitizedClusterCorrelator_h

#include <ap_int.h>
#include <vector>


namespace l1tp2 {

    class DigitizedClusterCorrelator {

        private:
            // Data
            ap_uint<64> clusterData;
            int idxGCTCard;  // 0, 1, or 2

            // Constants
            const unsigned int n_towers_eta = 34; // in GCT card unique region
            const unsigned int n_towers_phi = 24;  // in GCT card unique region
            const unsigned int n_crystals_in_tower = 5;

            const float LSB_PT = 0.5; // 0.5 GeV

            const unsigned int n_bits_pt = 12; // 12 bits allocated for pt

            // Private member functions to perform digitization 
            ap_uint<12> digitizePt(float pt_f) {
                float maxPt_f = (std::pow(2, n_bits_pt) - 1) * LSB_PT;
                std::cout << "max pT: " << maxPt_f << std::endl;

                // If pT exceeds the maximum (extremely unlikely), saturate the value
                if (pt_f >= maxPt_f) {
                    return (ap_uint<12>) 0xFFF; 
                }

                return (ap_uint<12>) (pt_f / LSB_PT);
            }

            // 
            ap_uint<6> digitizeIEta(unsigned int iEta) {
                assert(iEta < n_towers_eta);
                return (ap_uint<6>) iEta;
            }

            // This is tower iPhi in the unique region of the GCT card, so this only goes from 0-23
            ap_uint<5> digitizeIPhi(unsigned int iPhi) {
                assert(iPhi < n_towers_phi);
                return (ap_uint<5>) iPhi;
            }

            ap_uint<3> digitizeIEtaCr(unsigned int iEtaCr) {
                assert(iEtaCr < n_crystals_in_tower);
                return (ap_uint<3>) iEtaCr;
            }

            ap_uint<3> digitizeIPhiCr(unsigned int iPhiCr) {
                assert(iPhiCr < n_crystals_in_tower);
                return (ap_uint<3>) iPhiCr;
            }

            // To-do: hoe information?
            ap_uint<4> digitizeHoE(unsigned int hoe) {
                return (ap_uint<4>) hoe;
            }

            // 
            ap_uint<3> digitizeIso(bool iso) {
                return (ap_uint<3>) iso;
            }

            // To-do: fb: no information yet
            ap_uint<6> digitizeFb(unsigned int fb) {
                return (ap_uint<6>) fb;
            }

            // To-do: timing: no information yet
            ap_uint<5> digitizeTiming(unsigned int timing) {
                return (ap_uint<5>) timing;
            }
            
            // Shape: shower shape working point
            ap_uint<1> digitizeShape(bool is_ss) {
                return (ap_uint<1>) is_ss;
            }

            // TO-DO: Brems: was brems applied (NOT STORED YET IN GCT)
            ap_uint<1> digitizeBrems(bool brems_applied) {
                return (ap_uint<1>) brems_applied;
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
                                       ap_uint<1> shape, ap_uint<1> brems,
                                       int iGCTCard) {

                                        clusterData = ((ap_uint<64>) pt) | (((ap_uint<64>) eta) << 12) | (((ap_uint<64>) phi) << 18) |
                                                      (((ap_uint<64>) etaCr) << 23) | (((ap_uint<64>) phiCr) << 26) |
                                                      (((ap_uint<64>) hoe) << 29) | (((ap_uint<64>) iso << 33)) | 
                                                      (((ap_uint<64>) fb << 36)) | (((ap_uint<64>) timing << 42)) |
                                                      (((ap_uint<64>) shape << 47)) | (((ap_uint<64>) brems << 48));
                                        idxGCTCard = iGCTCard;
                                       }

            // To-do: constructor from float inputs that will perform digitization
            DigitizedClusterCorrelator(float pt_f, unsigned int iEta, unsigned int iPhi,
                                       unsigned int iEtaCr, unsigned int iPhiCr,
                                       unsigned int hoe, bool iso,
                                       unsigned int fb, unsigned int timing,
                                       bool shape, unsigned int brems,
                                       int iGCTCard) {

                clusterData = (((ap_uint<64>) digitizePt(pt_f)) |
                               ((ap_uint<64>) digitizeIEta(iEta) << 12) | ((ap_uint<64>) digitizeIPhi(iPhi) << 18) | 
                               ((ap_uint<64>) digitizeIEtaCr(iEtaCr) << 23) | ((ap_uint<64>) digitizeIPhiCr(iPhiCr) << 26) |
                               ((ap_uint<64>) digitizeHoE(hoe) << 29) | ((ap_uint<64>) digitizeIso(iso) << 33) |
                               ((ap_uint<64>) digitizeFb(fb) << 36) | ((ap_uint<64>) digitizeTiming(timing) << 42) | 
                               ((ap_uint<64>) digitizeShape(shape) << 47) | ((ap_uint<64>) digitizeBrems(brems) << 48)
                               ); 
                idxGCTCard = iGCTCard;
            }

            ap_uint<64> data() const { return clusterData; }

            // Other getters
            float ptLSB() { return LSB_PT; }
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
            int cardNumber() { return idxGCTCard; }  // which GCT card (0, 1, or 2)

            const int unusedBitsStart() { return 49; } // unused bits start at bit 49

            // Prints
            void print(const std::string location);
            void printPt(void) {  std::cout << "pt: " << std::bitset<12>{pt()} << std::endl; };
            void printEta(void) {  std::cout << "eta: " << std::bitset<6>{eta()} << std::endl; };
            void printPhi(void) {  std::cout << "phi: " << std::bitset<5>{phi()} << std::endl; };
            void printEtaCr(void) { std::cout << "etaCr: " << std::bitset<3>{etaCr()} << std::endl; };
            void printPhiCr(void) { std::cout << "phiCr: " << std::bitset<3>{phiCr()} << std::endl; };
            void printHoE(void) { std::cout << "hoe: " << std::bitset<4>{hoe()} << std::endl; };
            void printCardNumber(void) { std::cout << "GCT card number: " << cardNumber() << std::endl; };

            // Other checks
            bool passNullBitsCheck(void)  {
                return ((data() >> unusedBitsStart()) == 0x0);
            }

    };

    // Collection typedef
    typedef std::vector<l1tp2::DigitizedClusterCorrelator> DigitizedClusterCorrelatorCollection;

} // namespace l1tp2

#endif