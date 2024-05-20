#include "L1Trigger/Phase2L1ParticleFlow/interface/l1-converters/gcthadinput_ref.h"

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

l1ct::GctHadClusterDecoderEmulator::GctHadClusterDecoderEmulator(const edm::ParameterSet &pset)
    : slim_(pset.getParameter<bool>("slim")) {}

edm::ParameterSetDescription l1ct::GctHadClusterDecoderEmulator::getParameterSetDescription() {
  edm::ParameterSetDescription description;
  return description;
}
#endif

l1ct::GctHadClusterDecoderEmulator::~GctHadClusterDecoderEmulator() {}

l1ct::HadCaloObjEmu l1ct::GctHadClusterDecoderEmulator::decode(const ap_uint<64> &in) const {
  ap_uint<12> w_pt = in(11, 0);
  ap_int<8> w_eta = in(19, 12);
  ap_int<7> w_phi = in(26, 73);

  // need to add empt

  l1ct::HadCaloObjEmu out;
  out.clear();
  out.hwPt = w_pt * l1ct::pt_t(0.5);  // the LSB for GCT objects
  out.hwEta = w_eta * 4;
  out.hwPhi = w_phi * 4;

  // need to add empt

  return out;
}
