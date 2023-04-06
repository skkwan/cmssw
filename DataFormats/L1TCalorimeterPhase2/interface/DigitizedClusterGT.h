#ifndef DataFormats_L1TCalorimeterPhase2_DigitizedClusterGT_h
#define DataFormats_L1TCalorimeterPhase2_DigitizedClusterGT_h

#include <ap_int.h>
#include <bitset>
#include <vector>


namespace l1tp2 {

    class DigitizedClusterGT {

        private:
            // Data
            ap_uint<64> clusterData;

            // Constants
            static constexpr float LSB_PT = 0.03125; // 0.03125 GeV
            static constexpr unsigned int n_bits_eta_pi = 12; // 12 bits corresponds to pi in eta
            static constexpr unsigned int n_bits_phi_pi = 12; // 12 bits corresponds to pi in phi
            static constexpr unsigned int n_bits_pt = 16; // 12 bits allocated for pt
            static constexpr unsigned int n_bits_unused_start = 44; // unused bits start at bit number 44

            // Private member functions to perform digitization 
            ap_uint<1> digitizeIsValid(bool isValid) {
                return (ap_uint<1>) isValid;
            }

            ap_uint<16> digitizePt(float pt_f) {
                float maxPt_f = (std::pow(2, n_bits_pt) - 1) * LSB_PT;
                std::cout << "max pT: " << maxPt_f << std::endl;
                // If pT exceeds the maximum, saturate the value
                if (pt_f >= maxPt_f) {
                    return (ap_uint<16>) 0xFFFF; 
                }
                return (ap_uint<16>) (pt_f / LSB_PT);
            }

            // Use sign-magnitude convention
            ap_uint<13> digitizePhi(float phi_f) {
                ap_uint<1> sign = (phi_f >= 0) ? 0 : 1; 
                float phiMag_f = std::abs(phi_f);
                float lsb_phi = (M_PI / (std::pow(2, n_bits_phi_pi) - 1));
                ap_uint<12> phiMag_digitized = (phiMag_f / lsb_phi);
                ap_uint<13> phi_digitized = ((ap_uint<13>) sign) | ((ap_uint<13>) phiMag_digitized << 1);
                return phi_digitized;
            }


            // Use sign-magnitude convention
            ap_uint<14> digitizeEta(float eta_f) {
                ap_uint<1> sign = (eta_f >= 0) ? 0 : 1; 
                float etaMag_f = std::abs(eta_f);
                float lsb_eta = (M_PI / (std::pow(2, n_bits_eta_pi) - 1));
                ap_uint<13> etaMag_digitized = (etaMag_f / lsb_eta);
                ap_uint<14> eta_digitized = ((ap_uint<14>) sign) | ((ap_uint<14>) etaMag_digitized << 1);
                return eta_digitized;
            }


        public:
            DigitizedClusterGT() {
                clusterData = 0x0;
            }

            DigitizedClusterGT(ap_uint<64> data) {
                clusterData = data;
            }

            // Constructor from digitized inputs
            DigitizedClusterGT(ap_uint<1> isValid, ap_uint<16> pt, ap_uint<13> phi, ap_uint<14> eta) {
                clusterData = ((ap_uint<64>) isValid) | (((ap_uint<64>) pt) << 1) |
                              (((ap_uint<64>) phi) << 17) | (((ap_uint<64>) eta) << 30);
            }

            // To-do: constructor from float inputs that will perform digitization
            DigitizedClusterGT(bool isValid, float pt_f, float phi_f, float eta_f) {
                clusterData = (((ap_uint<64>) digitizeIsValid(isValid)) |
                               ((ap_uint<64>) digitizePt(pt_f) << 1) |
                               ((ap_uint<64>) digitizePhi(phi_f) << 17) | 
                               ((ap_uint<64>) digitizeEta(eta_f) << 30));
            }

            ap_uint<64> data() const { return clusterData; }

            // Other getters
            float ptLSB() const { return LSB_PT; }
            ap_uint<1> isValid() const { return (clusterData & 0x1); }
            ap_uint<16> pt() const { return ((clusterData >> 1) & 0xFFFF); } // 16 1's = 0xFFFF
            ap_uint<13> phi() const { return ((clusterData >> 17) & 0x1FFF); }  // (thirteen 1's)= 0x1FFF
            ap_uint<14> eta() const { return ((clusterData >> 30) & 0x3FFF); }  // (fourteen 1's) = 0x3FFF

            const int unusedBitsStart() const { return n_bits_unused_start; } // unused bits start at bit 44

            // Prints
            void print(const std::string location);
            void printIsValid(void) const { std::cout << "isValid: " << isValid() << std::endl; };
            void printPt(void) const {  std::cout << "pt: " << std::bitset<16>{pt()} << std::endl; };
            void printPhi(void) const {  std::cout << "phi: " << std::bitset<13>{phi()} << std::endl; };
            void printEta(void) const {  std::cout << "eta: " << std::bitset<14>{eta()} << std::endl; };

            // Other checks
            bool passNullBitsCheck(void) const {
                return ((data() >> unusedBitsStart()) == 0x0);
            }

    };

    // Collection typedef
    typedef std::vector<l1tp2::DigitizedClusterGT> DigitizedClusterGTCollection;

} // namespace l1tp2

#endif