#include "L1Trigger/DemonstratorTools/interface/codecs/htsums.h"
#include "DataFormats/Math/interface/LorentzVector.h"

namespace l1t::demo::codecs {

  // Encode the missing HT from the missing HT module and the HT from the HT module
  // Currently the missing HT module input is not used, and only the HT scalar sum from the HT module is used
  ap_uint<64> encodeHtSum(const l1t::EtSum& htSum, const l1t::EtSum& htScalarSum) {
    l1tmhtemu::EtMiss htMiss;
    htMiss.Et = 0; // originally htSum.p4().energy();
    htMiss.Phi = 0; // originally htSum.hwPhi();
    ap_ufixed<l1thtemu::kScalarSumHTSize, l1thtemu::kScalarSumHTIntSize> HT; 
    HT.range() = htScalarSum.hwPt();  
    // originally: ap_uint<l1tmhtemu::kMHTSize> HT = htSum.hwPt();
    std::cout << "DemonstratorTools/src/codecs_htsums.cc: htScalarSum.hwPt(): " <<  htScalarSum.hwPt() << " after trying to cast to ap_ufixed: " << HT << std::endl;
    ap_uint<l1thtemu::kValidSize> valid = (htScalarSum.hwQual() > 0);
    ap_uint<l1thtemu::kUnassignedSize> unassigned = 0;
    ap_uint<64> htSumWord = (unassigned, HT.range(), htMiss.Phi, htMiss.Et.range(), valid);
    return htSumWord;
  }

  // Encodes htsum collection onto 1 output link
  std::array<std::vector<ap_uint<64>>, 1> encodeHtSums(const edm::View<l1t::EtSum>& htSums, const edm::View<l1t::EtSum>& htScalarSums) {
    std::vector<ap_uint<64>> htSumWords;

    if (htSums.size() != htScalarSums.size()) {
      throw cms::Exception("InvalidInput") << "htSums size: " << htSums.size() << " != htScalarSums size " << htScalarSums.size();
    }

    for (unsigned int i = 0; i < htSums.size(); i++) {
      htSumWords.push_back(encodeHtSum(htSums[i], htScalarSums[i]));
    }
   
    std::array<std::vector<ap_uint<64>>, 1> linkData;

    for (size_t i = 0; i < linkData.size(); i++) {
      // Pad etsum vectors -> full packet length (48 frames, but only 1 htsum max)
      htSumWords.resize(1, 0);
      linkData.at(i) = htSumWords;
    }

    return linkData;
  }

  // Decode HT sum
  std::vector<l1t::EtSum> decodeHtSums(const std::vector<ap_uint<64>>& frames) {
    std::vector<l1t::EtSum> htScalarSums;

    for (const auto& x : frames) {
      if (not x.test(0))
        break;

      math::XYZTLorentzVector v(0, 0, 0, l1thtemu::ht_t(x(l1thtemu::kVectorSumMSB, l1thtemu::kVectorSumLSB)).to_int());
      l1t::EtSum s(v,
                   l1t::EtSum::EtSumType::kTotalHt,
                   l1thtemu::ht_t(x(l1thtemu::kScalarSumHTMSB, l1thtemu::kScalarSumHTLSB)), // the only meaningful value
                   0,
                   0,
                   0);
      htScalarSums.push_back(s);
    }

    return htScalarSums;
  }


}  // namespace l1t::demo::codecs
