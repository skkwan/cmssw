#ifndef L1Trigger_Phase2L1ParticleFlow_newfirmware_gcteminput_ref_h
#define L1Trigger_Phase2L1ParticleFlow_newfirmware_gcteminput_ref_h

#include "DataFormats/L1TParticleFlow/interface/layer1_emulator.h"

namespace edm {
  class ParameterSet;
  class ParameterSetDescription;
}  // namespace edm

namespace l1ct {
  class GctEmClusterDecoderEmulator {

  public:
    GctEmClusterDecoderEmulator() {};
    GctEmClusterDecoderEmulator(const edm::ParameterSet &pset);

    ~GctEmClusterDecoderEmulator();

    static edm::ParameterSetDescription getParameterSetDescription();

    l1ct::EmCaloObjEmu decode(const ap_uint<64> &in) const;
  };
}  // namespace l1ct

#endif
