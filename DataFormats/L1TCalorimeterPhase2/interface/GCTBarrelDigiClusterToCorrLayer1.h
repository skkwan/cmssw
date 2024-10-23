#ifndef DataFormats_L1TCalorimeterPhase2_GCTBarrelDigiClusterToCorrLayer1_h
#define DataFormats_L1TCalorimeterPhase2_GCTBarrelDigiClusterToCorrLayer1_h

#include <ap_int.h>
#include <vector>

namespace l1tp2 {

  class GCTBarrelDigiClusterToCorrLayer1 {
  private:
    // Data (note: positional information is entirely encoded in the location in the output array)
    unsigned long long int clusterData;

    // Constants
    static constexpr float LSB_PT = 0.5;                 // 0.5 GeV

    // start of the unused bits 
    static constexpr int n_bits_unused_start = 52; 

  public:
    GCTBarrelDigiClusterToCorrLayer1() { clusterData = 0x0; }

    GCTBarrelDigiClusterToCorrLayer1(ap_uint<64> data) { clusterData = data; }

    GCTBarrelDigiClusterToCorrLayer1(
                               ap_uint<12> pt,
                               int etaCr,
                               int phiCr,
                               ap_uint<4> hoe,
                               ap_uint<2> hoeFlag,
                               ap_uint<3> iso,
                               ap_uint<2> isoFlag,
                               ap_uint<6> fb,
                               ap_uint<5> timing,
                               ap_uint<2> shapeFlag,
                               ap_uint<2> brems) {
      // iEta is an unsigned quantity in bits 12 through 18 (7 bits)
      ap_uint<64> temp_data_eta = 0x0; 
      temp_data_eta |= ((ap_uint<64>) ((0x7F & abs(etaCr)) << 12));  

      // Repeat for phi, which is bits 19 through 25 (7 bits). The sign bit is 25, leaving 6 bits for the magnitude
      ap_uint<64> temp_data_phi = 0x0;
      if (phiCr > 0) temp_data_phi |= ((ap_uint<64>) (0x1 << 25));   // set bit 26 to 1 if iPhi is positive 
      temp_data_phi |= ((ap_uint<64>) ((0x3F & abs(phiCr)) << 19));   // 0x3F is 0b111111 (six 1's)

      clusterData = ((ap_uint<64>)pt) | temp_data_eta | temp_data_phi |
                    (((ap_uint<64>)hoe) << 26) | (((ap_uint<64>)hoeFlag) << 30) | (((ap_uint<64>)iso) << 32) |
                    (((ap_uint<64>)isoFlag) << 35) | (((ap_uint<64>)fb) << 37) | (((ap_uint<64>)timing) << 43) |
                    (((ap_uint<64>)shapeFlag << 48)) | (((ap_uint<64>)brems << 50));
    }

    // Getters
    ap_uint<64> data() const { return clusterData; }

    // Other getters
    float ptLSB() const { return LSB_PT; }
    ap_uint<12> pt() const { return (clusterData & 0xFFF); }
    float ptFloat() const { return pt() * ptLSB(); }

    // crystal eta (absolute value)
    int eta() const { return ((clusterData >> 12) & 0xFF); }  

    // crystal phi (signed quantity)
    int phi() const { 
      int signed_val; 
      // get the sign bit (the 25th bit). If it is 1, phi is positive. If it is 0, phi is negative
      // Magnitude is bits 19 through 24
      if (clusterData & (0b1 << 25)) signed_val = ((clusterData >> 19) & 0x3F);  // 0x3F is six 1's
      else signed_val = - ((clusterData >> 19) & 0x3F); 

      return signed_val;
    }  

    // HoE value and flag: not defined yet in the emulator
    ap_uint<4> hoe() const { return ((clusterData >> 26) & 0xF); }      // (four 1's) 0b1111 = 0xF
    ap_uint<2> hoeFlag() const { return ((clusterData >> 30) & 0x3); }  // (two 1's) 0b11 = 0x3

    // Raw isolation sum: not saved in the emulator
    ap_uint<3> iso() const { return ((clusterData >> 32) & 0x7); }

    // iso flag: two bits, least significant bit is the standalone WP (true or false), second bit is the looseTk WP (true or false)
    // e.g. 0b01 : standalone iso flag passed, loose Tk iso flag did not pass
    ap_uint<2> isoFlags() const { return ((clusterData >> 35) & 0x3); }  // (two 1's) 0b11 = 0x3
    bool passes_iso() const { return (isoFlags() & 0x1); }               // standalone iso WP
    bool passes_looseTkiso() const { return (isoFlags() & 0x2); }        // loose Tk iso WP

    // fb and timing: not saved in the current emulator
    ap_uint<6> fb() const { return ((clusterData >> 37) & 0x3F); }
    ap_uint<5> timing() const { return ((clusterData >> 43) & 0x1F); }

    // shower shape shape flag: two bits, least significant bit is the standalone WP, second bit is the looseTk WP
    // e.g. 0b01 : standalone shower shape flag passed, loose Tk shower shape flag did not pass
    ap_uint<2> shapeFlags() const { return ((clusterData >> 48) & 0x3); }

    bool passes_ss() const { return (shapeFlags() & 0x1); }         // standalone shower shape WP
    bool passes_looseTkss() const { return (shapeFlags() & 0x2); }  // loose Tk shower shape WP

    // brems: not saved in the current emulator
    ap_uint<2> brems() const { return ((clusterData >> 50) & 0x3); }

    const int unusedBitsStart() const { return n_bits_unused_start; }

    // Other checks
    bool passNullBitsCheck(void) const { return ((data() >> unusedBitsStart()) == 0x0); }

  };

  // Collection typedef
  typedef std::vector<l1tp2::GCTBarrelDigiClusterToCorrLayer1> GCTBarrelDigiClusterToCorrLayer1Collection;
  typedef std::vector<l1tp2::GCTBarrelDigiClusterToCorrLayer1Collection> GCTBarrelDigiClusterToCorrLayer1CollectionFullDetector;

}  // namespace l1tp2

#endif