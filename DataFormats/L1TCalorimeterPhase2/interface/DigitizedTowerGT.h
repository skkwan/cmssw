#ifndef DataFormats_L1TCalorimeterPhase2_DigitizedTowerGT_h
#define DataFormats_L1TCalorimeterPhase2_DigitizedTowerGT_h

#include <ap_int.h>
#include <vector>


namespace l1tp2 {

    class DigitizedTowerGT {

        private:
            // Data
            ap_uint<16> towerData;
            unsigned int idxCard; // 0, 1, or 2 (there are three GCT cards)
            unsigned int idxFiber; // 0 to 47 (there are 48 fibers in one GCT card)
            unsigned int idxTower; // 0 to 16 (there are 17 towers in one fiber)

            // Constants
            const float LSB_ET = 0.5; // 0.5 GeV, so max value is (2^10 - 1) * 0.5 = 511.5 GeV
            const unsigned int n_bits_pt = 10;
            const unsigned int n_towers_in_fiber = 17;
            const unsigned int n_fibers_in_card = 48;
            const unsigned int n_cards = 3;

            // Private member functions to perform digitization 
            ap_uint<10> digitizeEt(float et_f) {
                float maxEt_f = (std::pow(2, n_bits_pt) - 1) * LSB_ET;
                // If pT exceeds the maximum, saturate the value
                if (et_f >= maxEt_f) {
                    return (ap_uint<10>) 0x3FF;  
                }
                return (ap_uint<10>) (et_f / LSB_ET);
            }

            ap_uint<4> digitizeHoE(ap_uint<4> hoe) {
                return (ap_uint<4>) hoe;
            }

            // To-do: FB not implemented yet
            ap_uint<2> digitizeFB(ap_uint<2> fb) {
                return (ap_uint<2>) fb;
            }


        public:
            DigitizedTowerGT() {
                towerData = 0x0;
            }

            DigitizedTowerGT(ap_uint<16> data) {
                towerData = data;
            }

            // Constructor from digitized inputs
            DigitizedTowerGT(ap_uint<10> et, ap_uint<4> hoe, ap_uint<2> fb,
                            unsigned int indexCard, unsigned int indexFiber, unsigned int indexTower) {
                towerData = ((ap_uint<16>) et) | (((ap_uint<16>) hoe) << 10) |
                            (((ap_uint<16>) fb) << 14);
                idxCard = indexCard;
                idxFiber = indexFiber;
                idxTower = indexTower;
                assert(hasValidIndices());
            }

            // Constructor from float inputs
            DigitizedTowerGT(float et_f,  ap_uint<4> hoe, ap_uint<2> fb,
                             unsigned int indexCard, unsigned int indexFiber, unsigned int indexTower) {
                towerData = ((ap_uint<16>) digitizeEt(et_f)) | (((ap_uint<16>) hoe) << 10) |
                            (((ap_uint<16>) fb) << 14);
                idxCard = indexCard;
                idxFiber = indexFiber;
                idxTower = indexTower;
                assert(hasValidIndices());
            }

            ap_uint<16> data() const { return towerData; }

            // Other getters
            float etLSB() { return LSB_ET; }
            ap_uint<10> et() { return (towerData & 0x3FF); } // ten 1's = 0x3FF
            ap_uint<4> hoe() { return ((towerData >> 10) & 0xF); }  // four 1's= 0xF
            ap_uint<2> fb() { return ((towerData >> 14) & 0x3); }  // two 1's = 0x3
            unsigned int cardNumber() { return idxCard; } // GCT card number
            unsigned int fiberNumber() { return idxFiber; } // fiber number in card (hardware convention)
            unsigned int towerNumber() { return idxTower; } // tower number in fiber (hardware convention)

            // Prints
            void print(const std::string location);
            void printEt(void) {  std::cout << "et: " << std::bitset<10>{et()} << std::endl; };
            void printHoE(void) {  std::cout << "hoe: " << std::bitset<4>{hoe()} << std::endl; };
            void printFB(void) {  std::cout << "fb: " << std::bitset<2>{fb()} << std::endl; };
            void printInfo(void) { std::cout << "GCT card: " << cardNumber() << ", fiber number: " << fiberNumber() << ", tower number: " << towerNumber() << std::endl; };

            // Other checks
            bool hasValidIndices(void) { 
                return (idxTower < n_towers_in_fiber) && (idxFiber < n_fibers_in_card) && (idxCard < n_cards);
            }

    };

    // Collection typedef
    typedef std::vector<l1tp2::DigitizedTowerGT> DigitizedTowerGTCollection;

} // namespace l1tp2

#endif