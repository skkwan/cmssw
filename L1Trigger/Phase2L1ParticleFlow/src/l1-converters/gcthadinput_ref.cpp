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

double fracPart(double total, unsigned int hoe) {
  return total * std::pow(2.0, hoe) / (std::pow(2.0, hoe) + 1);
}

l1ct::HadCaloObjEmu l1ct::GctHadClusterDecoderEmulator::decode(const ap_uint<64> &in) const {
  ap_uint<12> pt = in(11, 0);
  ap_int<8> eta = in(19, 12);  // might need to shift over if this doesn't unify with em
  ap_int<7> phi = in(26, 20);

  l1ct::HadCaloObjEmu out;
  out.clear();
  out.hwPt = pt * l1ct::pt_t(0.5);  // the LSB for GCT objects
  out.hwEta = eta * 4;
  out.hwPhi = phi * 4;

  // need to add empt
  ap_uint<4> hoeVal = in(30, 27);
  // the lsb indicates what's bigger, EM or HAD
  auto isEMBigger = static_cast<bool>(hoeVal[0]);
  // This is not quite true. If HAD energy goes down to 0, then it flips and says that HAD is bigger
  ap_uint<3> hoe = hoeVal(3, 1);

  if (isEMBigger) {
    auto em = fracPart(out.hwPt.to_double(), hoe.to_uint());
    out.hwEmPt = em;
  } else {
    pt_t had = fracPart(out.hwPt.to_double(), hoe.to_uint());
    out.hwEmPt = out.hwPt - had;
  }


  return out;
}
