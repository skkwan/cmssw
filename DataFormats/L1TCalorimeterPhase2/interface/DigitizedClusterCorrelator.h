#ifndef DataFormats_L1TCalorimeterPhase2_DigitizedClusterCorrelator_h
#define DataFormats_L1TCalorimeterPhase2_DigitizedClusterCorrelator_h

#include <ap_int.h>
#include <vector>

namespace l1tp2 {

  class DigitizedClusterCorrelator {
  private:
    // Data
    unsigned long long int clusterData;
    unsigned int idxGCTCard;  // 0, 1, or 2

    // Constants
    static constexpr unsigned int n_towers_eta = 34;  // in GCT card unique region
    static constexpr unsigned int n_towers_phi = 24;  // in GCT card unique region
    static constexpr unsigned int n_crystals_in_tower = 5;
    static constexpr float LSB_PT = 0.5;                     // 0.5 GeV
    static constexpr unsigned int n_bits_pt = 12;            // 12 bits allocated for pt
    static constexpr unsigned int n_bits_unused_start = 63;  // unused bits start at bit 63

    // Private member functions to perform digitization
    ap_uint<12> digitizePt(float pt_f) {
      float maxPt_f = (std::pow(2, n_bits_pt) - 1) * LSB_PT;
      // If pT exceeds the maximum (extremely unlikely), saturate the value
      if (pt_f >= maxPt_f) {
        return (ap_uint<12>)0xFFF;
      }

      return (ap_uint<12>)(pt_f / LSB_PT);
    }


    ap_uint<8> digitizeIEtaCr(unsigned int iEtaCr) {
      std::cout << "[inside digitizeIEtaCr: ] input is "<< iEtaCr << std::endl;
      return (ap_uint<8>) iEtaCr;
    }

    ap_uint<7> digitizeIPhiCr(unsigned int iPhiCr) {
      return (ap_uint<7>)iPhiCr;
    }

    // To-do: HoE is not defined for clusters
    ap_uint<4> digitizeHoE(unsigned int hoe) { return (ap_uint<4>) hoe; }
    ap_uint<2> digitizeHoeFlag(unsigned int hoeFlag) { return (ap_uint<2>) hoeFlag; }

    ap_uint<3> digitizeIso(unsigned int iso) { return (ap_uint<3>) iso; }
    ap_uint<2> digitizeIsoFlag(unsigned int isoFlag) { return (ap_uint<2>) isoFlag; }

    // To-do: fb: no information yet
    ap_uint<6> digitizeFb(unsigned int fb) { return (ap_uint<6>) fb; }

    // To-do: timing: no information yet
    ap_uint<5> digitizeTiming(unsigned int timing) { return (ap_uint<5>) timing; }

    // Shape: shower shape working point
    ap_uint<2> digitizeShapeFlag(unsigned int shapeFlag) { return (ap_uint<2>) shapeFlag; }

    // TO-DO: Brems: was brems applied (NOT STORED YET IN GCT)
    ap_uint<2> digitizeBrems(unsigned int brems) { return (ap_uint<2>) brems; }

  public:
    DigitizedClusterCorrelator() { clusterData = 0x0; }

    DigitizedClusterCorrelator(ap_uint<64> data) { clusterData = data; }

    // Constructor from digitized inputs
    DigitizedClusterCorrelator(ap_uint<12> pt,
                               ap_uint<8> etaCr,
                               ap_uint<7> phiCr,
                               ap_uint<4> hoe,
                               ap_uint<2> hoeFlag,
                               ap_uint<3> iso,
                               ap_uint<2> isoFlag,
                               ap_uint<6> fb,
                               ap_uint<5> timing,
                               ap_uint<2> shapeFlag,
                               ap_uint<2> brems,
                               int iGCTCard,
                               bool fullydigitizedInputs) {
      (void)fullydigitizedInputs;

      clusterData = ((ap_uint<64>)pt) |
                    (((ap_uint<64>)etaCr) << 20) | (((ap_uint<64>)phiCr) << 27) |
                    (((ap_uint<64>)hoe) << 31) | (((ap_uint<64>)hoeFlag) << 33) |
                    (((ap_uint<64>)iso) << 36) | (((ap_uint<64>)isoFlag) << 38) |
                    (((ap_uint<64>)fb) << 44) | (((ap_uint<64>)timing) << 49) |
                    (((ap_uint<64>)shapeFlag << 51)) | (((ap_uint<64>)brems << 53));
      idxGCTCard = iGCTCard;
    }

    // Constructor from float inputs
    DigitizedClusterCorrelator(float pt_f,
                               unsigned int iEtaCr,
                               unsigned int iPhiCr,
                               unsigned int hoe,
                               unsigned int hoeFlag,
                               unsigned int iso,
                               unsigned int isoFlag,
                               unsigned int fb,
                               unsigned int timing,
                               unsigned int shapeFlag,
                               unsigned int brems,
                               int iGCTCard) {
      std::cout << "input iEtaCr: " << iEtaCr << std::endl;
      clusterData = (((ap_uint<64>) digitizePt(pt_f)) | 
                     ((ap_uint<64>) digitizeIEtaCr(iEtaCr) << 20) |
                     ((ap_uint<64>) digitizeIPhiCr(iPhiCr) << 27) |
                     ((ap_uint<64>) digitizeHoE(hoe) << 31) |
                     ((ap_uint<64>) digitizeHoeFlag(hoeFlag) << 33) | 
                     ((ap_uint<64>) digitizeIso(iso) << 36) |
                     ((ap_uint<64>) digitizeIsoFlag(isoFlag) << 38) |
                     ((ap_uint<64>) digitizeFb(fb) << 44) |
                     ((ap_uint<64>) digitizeTiming(timing) << 49) |
                     ((ap_uint<64>) digitizeShapeFlag(shapeFlag) << 51) |
                     ((ap_uint<64>) digitizeBrems(brems) << 53));
      idxGCTCard = iGCTCard;
    }

    ap_uint<64> data() const { return clusterData; }

    // Other getters
    float ptLSB() const { return LSB_PT; }
    ap_uint<12> pt() const { return (clusterData & 0xFFF); } 

    // crystal eta in the correlator region (LSB: 2.8/170)
    ap_uint<8> eta() const { return ((clusterData >> 12) & 0xFF); }   // (eight 1's) 0b11111111 = 0xFF

    // crystal phi in the correlator region (LSB: 2pi/360)
    ap_uint<7> phi() const { return ((clusterData >> 20) & 0x7F); }   // (seven 1's) 0b1111111 = 0x7F

    // HoE value and flag: not defined yet in the emulator
    ap_uint<4> hoe() const { return ((clusterData >> 27) & 0xF); }  // (four 1's) 0b1111 = 0xF
    ap_uint<2> hoeFlag() const { return ((clusterData >> 31) & 0x3); }  // (two 1's) 0b11 = 0x3

    // Raw isolation sum: not saved in the emulator
    ap_uint<3> iso() const { return ((clusterData >> 33) & 0x7); }

    // iso flag: two bits, least significant bit is the standalone WP (true or false), second bit is the looseTk WP (true or false)
    // e.g. 0b01 : standalone iso flag passed, loose Tk iso flag did not pass
    ap_uint<2> isoFlags() const { return ((clusterData >> 36) & 0x3); } // (two 1's) 0b11 = 0x3

    // fb and timing: not saved in the current emulator
    ap_uint<6> fb() const { return ((clusterData >> 38) & 0x3F); }  
    ap_uint<5> timing() const { return ((clusterData >> 44) & 0x1F); }

    // shower shape shape flag: two bits, least significant bit is the standalone WP (true or false), second bit is the looseTk WP (true or false)
    // e.g. 0b01 : standalone shower shape flag passed, loose Tk shower shape flag did not pass
    ap_uint<2> shapeFlags() const { return ((clusterData >> 49) & 0x3); }

    // brems: not saved in the current emulator 
    ap_uint<2> brems() const { return ((clusterData >> 51) & 0x3); }

    // which GCT card (0, 1, or 2)
    unsigned int cardNumber() const { return idxGCTCard; } 

    const int unusedBitsStart() const { return n_bits_unused_start; }  

    // Other checks
    bool passNullBitsCheck(void) const { return ((data() >> unusedBitsStart()) == 0x0); }
  };

  // Collection typedef
  typedef std::vector<l1tp2::DigitizedClusterCorrelator> DigitizedClusterCorrelatorCollection;

}  // namespace l1tp2

#endif