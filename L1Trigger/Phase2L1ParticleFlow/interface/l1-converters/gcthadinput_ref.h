#ifndef L1Trigger_Phase2L1ParticleFlow_newfirmware_gcthadinput_ref_h
#define L1Trigger_Phase2L1ParticleFlow_newfirmware_gcthadinput_ref_h

#include "DataFormats/L1TParticleFlow/interface/layer1_emulator.h"

namespace edm {
  class ParameterSet;
  class ParameterSetDescription;
}  // namespace edm

namespace l1ct {
  class GctHadClusterDecoderEmulator {
    bool slim_;

  public:
    GctHadClusterDecoderEmulator(bool slim = false) : slim_{slim} {};
    GctHadClusterDecoderEmulator(const edm::ParameterSet &pset);

    ~GctHadClusterDecoderEmulator();

    double fracPart(const double total, const unsigned int hoe) const;

    static edm::ParameterSetDescription getParameterSetDescription();

    l1ct::HadCaloObjEmu decode(const ap_uint<64> &in) const;
  };
}  // namespace l1ct

#endif
