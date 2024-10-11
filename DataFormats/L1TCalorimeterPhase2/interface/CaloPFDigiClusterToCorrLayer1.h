#ifndef DataFormats_L1TCalorimeterPhase2_CaloPFDigiClusterToCorrLayer1_h
#define DataFormats_L1TCalorimeterPhase2_CaloPFDigiClusterToCorrLayer1_h

#include <ap_int.h>
#include <vector>

namespace l1tp2 {

  class CaloPFDigiClusterToCorrLayer1 {
  private:
    // Data (note: positional information is entirely encoded in the location in the output array)
    unsigned long long int clusterData;

    // Constants
    static constexpr float LSB_PT = 0.5;                 // 0.5 GeV

    // start of the unused bits 
    static constexpr int n_bits_unused_start = 31; 

  public:
    CaloPFDigiClusterToCorrLayer1() { clusterData = 0x0; }

    CaloPFDigiClusterToCorrLayer1(ap_uint<64> data) { clusterData = data; }

    // Note types of the constructor
    CaloPFDigiClusterToCorrLayer1(
                               ap_uint<12> pt,
                               int etaCr,
                               int phiCr,
                               ap_uint<4> hoe) {
      // iEta is an unsigned quantity, in bits 12 through 18 (7 bits)
      ap_uint<64> temp_data_eta = 0x0; 
      temp_data_eta |= ((0x7F & abs(etaCr)) << 12);  

      // Repeat for phi, which is bits 20 through 26 (7 bits). The sign bit is 26, leaving 6 bits for the magnitude
      ap_uint<64> temp_data_phi = 0x0;
      if (phiCr > 0) temp_data_phi |= (0x1 << 25);   // set bit 26 to 1 if iPhi is positive 
      temp_data_phi |= ((0x3F & abs(phiCr)) << 19);   // 0x3F is 0b111111 (six 1's)

      clusterData = ((ap_uint<64>)pt) | temp_data_eta | temp_data_phi |
                    (((ap_uint<64>)hoe) << 26);
    }

    // Getters
    ap_uint<64> data() const { return clusterData; }

    // Other getters
    float ptLSB() const { return LSB_PT; }
    ap_uint<12> pt() const { return (clusterData & 0xFFF); }
    float ptFloat() const { return pt() * ptLSB(); }

    // crystal eta (unsigned quantity)
    int eta() const { return ((clusterData >> 12) & 0xFF);  }

    // crystal phi (signed quantity)
    int phi() const { 
      int signed_val; 
      // get the sign bit (the 26th bit). If it is 1, phi is positive. If it is 0, phi is negative
      // Magnitude is bits 20 through 25
      if (clusterData & (0b1 << 25)) signed_val = ((clusterData >> 19) & 0x3F);  // 0x3F is six 1's
      else signed_val = - ((clusterData >> 19) & 0x3F); 

      return signed_val;
    }  

    // HoE value and flag: not defined yet in the emulator
    ap_uint<4> hoe() const { return ((clusterData >> 26) & 0xF); }      // (four 1's) 0b1111 = 0xF
    const int unusedBitsStart() const { return n_bits_unused_start; }

    // Other checks
    bool passNullBitsCheck(void) const { return ((data() >> unusedBitsStart()) == 0x0); }

  };

  // Collection typedef
  typedef std::vector<l1tp2::CaloPFDigiClusterToCorrLayer1> CaloPFDigiClusterToCorrLayer1Collection;
  typedef std::vector<l1tp2::CaloPFDigiClusterToCorrLayer1Collection> CaloPFDigiClusterToCorrLayer1CollectionFullDetector;

}  // namespace l1tp2

#endif