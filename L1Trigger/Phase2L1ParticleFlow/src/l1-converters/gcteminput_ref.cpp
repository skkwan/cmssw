#include "L1Trigger/Phase2L1ParticleFlow/interface/l1-converters/gcteminput_ref.h"

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

l1ct::GctEmClusterDecoderEmulator::GctEmClusterDecoderEmulator(const edm::ParameterSet &pset) {}

edm::ParameterSetDescription l1ct::GctEmClusterDecoderEmulator::getParameterSetDescription() {
  edm::ParameterSetDescription description;
  return description;
}
#endif

l1ct::GctEmClusterDecoderEmulator::~GctEmClusterDecoderEmulator() {}

l1ct::EmCaloObjEmu l1ct::GctEmClusterDecoderEmulator::decode(const ap_uint<64> &in) const {
  ap_uint<12> pt = in(11, 0);
  ap_int<7> eta = in(18, 12);
  ap_int<7> phi = in(25, 19);

  // need to add emid
  l1ct::EmCaloObjEmu out;
  out.clear();
  out.hwPt = pt * l1ct::pt_t(0.5);  // the LSB for GCT objects
  out.hwEta = eta * 4;
  out.hwPhi = phi * 4;

  // need to add emid

  return out;
}
