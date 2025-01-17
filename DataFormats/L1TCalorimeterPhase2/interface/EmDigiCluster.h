#ifndef DataFormats_L1TCalorimeterPhase2_EmDigiCluster_h
#define DataFormats_L1TCalorimeterPhase2_EmDigiCluster_h

#include <ap_int.h>
#include <vector>

#include "DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/DigitizedClusterCorrelator.h"

namespace l1tp2 {

  class EmDigiCluster {
  private:
    // Data (note: positional information is entirely encoded in the location in the output array)
    unsigned long long int clusterData;

    // Constants
    static constexpr float LSB_PT = 0.5;                 // 0.5 GeV

    // start of the unused bits 
    static constexpr int n_bits_unused_start = 52; 

    // Reference to the original float cluster
    edm::Ref<l1tp2::CaloCrystalClusterCollection> clusterRef_;

    // reference to the original digitized cluster (before duplication in the output links)
    edm::Ref<l1tp2::DigitizedClusterCorrelatorCollection> digiClusterRef_;

  public:
    EmDigiCluster() { clusterData = 0; }

    EmDigiCluster(ap_uint<64> data) { clusterData = data; }

    EmDigiCluster(
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
      // To use .range() we need an ap class member
      ap_uint<64> temp_data; 
      ap_uint<7> etaCrDigitized = abs(etaCr);
      ap_int<7> phiCrDigitized = phiCr;

      temp_data.range(11, 0) = pt.range();
      temp_data.range(18, 12) = etaCrDigitized.range(); 
      temp_data.range(25, 19) = phiCrDigitized.range(); 
      temp_data.range(29, 26) = hoe.range();
      temp_data.range(31, 30) = hoeFlag.range();
      temp_data.range(34, 32) = iso.range();
      temp_data.range(36, 35) = isoFlag.range();
      temp_data.range(42, 37) = fb.range();
      temp_data.range(47, 43) = timing.range();
      temp_data.range(49, 48) = shapeFlag.range();
      temp_data.range(51, 50) = brems.range();

      clusterData = temp_data;
    }

    // Setters
    void setRef(const edm::Ref<l1tp2::CaloCrystalClusterCollection>& clusterRef) {
      clusterRef_ = clusterRef;
    }

    void setDigiRef(const edm::Ref<l1tp2::DigitizedClusterCorrelatorCollection>& digiClusterRef) {
      digiClusterRef_ = digiClusterRef;
    }

    // Getters
    ap_uint<64> data() const { return clusterData; }

    // Other getters
    float ptLSB() const { return LSB_PT; }
    ap_uint<12> pt() const { return data().range(11, 0); }
    float ptFloat() const { return pt() * ptLSB(); }

    // crystal eta (unsigned, 7 bits)
    int eta() const { return (ap_uint<7>) data().range(18, 12); }  

    // crystal phi (signed, 7 bits)
    int phi() const { return (ap_int<7>) data().range(25, 19); }

    // HoE value and flag: not defined yet in the emulator
    ap_uint<4> hoe() const { return data().range(29, 26); }     
    ap_uint<2> hoeFlag() const { return data().range(31, 30); }  

    // Raw isolation sum: not saved in the emulator
    ap_uint<3> iso() const { return data().range(34, 32); }

    // iso flag: two bits, least significant bit is the standalone WP (true or false), second bit is the looseTk WP (true or false)
    // e.g. 0b01 : standalone iso flag passed, loose Tk iso flag did not pass
    ap_uint<2> isoFlags() const { return data().range(36, 35); }  
    bool passes_iso() const { return (isoFlags() & 0x1); }               // standalone iso WP
    bool passes_looseTkiso() const { return (isoFlags() & 0x2); }        // loose Tk iso WP

    // fb and timing: not saved in the current emulator
    ap_uint<6> fb() const { return data().range(42, 37); }
    ap_uint<5> timing() const { return data().range(47, 43); }

    // shower shape shape flag: two bits, least significant bit is the standalone WP, second bit is the looseTk WP
    // e.g. 0b01 : standalone shower shape flag passed, loose Tk shower shape flag did not pass
    ap_uint<2> shapeFlags() const { return data().range(49, 48); }

    bool passes_ss() const { return (shapeFlags() & 0x1); }         // standalone shower shape WP
    bool passes_looseTkss() const { return (shapeFlags() & 0x2); }  // loose Tk shower shape WP

    // brems: not saved in the current emulator
    ap_uint<2> brems() const { return data().range(51, 50); }

    // Check that unused bits are zero
    const int unusedBitsStart() const { return n_bits_unused_start; }
    bool passNullBitsCheck(void) const { return ((data() >> unusedBitsStart()) == 0); }

    // Note: not possible to get real eta and phi without knowing the link

    // Get the underlying float cluster
    const edm::Ref<l1tp2::CaloCrystalClusterCollection>& clusterRef() const {
      return clusterRef_;
    }
    // Get the underlying digitized cluster (before duplication and zero-padding)
    const edm::Ref<l1tp2::DigitizedClusterCorrelatorCollection>& digiClusterRef() const {
      return digiClusterRef_;
    }
  };

  // Collection typedef
  typedef std::vector<l1tp2::EmDigiCluster> EmDigiClusterLink;
  typedef std::vector<l1tp2::EmDigiClusterLink> EmDigiClusterCollection;

}  // namespace l1tp2

#endif