#include "L1Trigger/DemonstratorTools/interface/codecs/htsums.h"
#include "DataFormats/Math/interface/LorentzVector.h"

namespace l1t::demo::codecs {

  ap_uint<64> encodeHtSum(const l1t::EtSum& htSum) {
    l1tmhtemu::EtMiss htMiss;
    htMiss.Et = htSum.p4().energy();
    htMiss.Phi = htSum.hwPhi();
    ap_uint<l1tmhtemu::kMHTSize> HT = htSum.hwPt();
    ap_uint<l1tmhtemu::kValidSize> valid = (htSum.hwQual() > 0);
    ap_uint<l1tmhtemu::kUnassignedSize> unassigned = 0;
    ap_uint<64> htSumWord = (unassigned, HT, htMiss.Phi, htMiss.Et.range(), valid);
    return htSumWord;
  }

  ap_uint<64> encodeHtScalarSum(const l1t::EtSum& htSum) {
    l1thtemu::Ht ht; 
    ht.Et = htSum.p4().energy();
    ap_ufixed<l1thtemu::kPtSize, l1thtemu::kPtIntSize> HT = htSum.et(); 
    ap_uint<l1thtemu::kValidSize> valid = (htSum.hwQual() > 0); // TODO: check this is correct
    ap_uint<l1thtemu::kUnassignedSize> unassigned = 0;
    ap_uint<64> htScalarSumWord = (unassigned, HT.range(), valid);
    return htScalarSumWord;
  }


  // Encodes htsum collection onto 1 output link
  std::array<std::vector<ap_uint<64>>, 1> encodeHtSums(const edm::View<l1t::EtSum>& htSums) {
    std::vector<ap_uint<64>> htSumWords;

    for (const auto& htSum : htSums)
      htSumWords.push_back(encodeHtSum(htSum));

    std::array<std::vector<ap_uint<64>>, 1> linkData;

    for (size_t i = 0; i < linkData.size(); i++) {
      // Pad etsum vectors -> full packet length (48 frames, but only 1 htsum max)
      htSumWords.resize(1, 0);
      linkData.at(i) = htSumWords;
    }

    return linkData;
  }

  std::vector<l1t::EtSum> decodeHtSums(const std::vector<ap_uint<64>>& frames) {
    std::vector<l1t::EtSum> htSums;

    for (const auto& x : frames) {
      if (not x.test(0))
        break;

      math::XYZTLorentzVector v(0, 0, 0, l1tmhtemu::MHT_t(x(l1tmhtemu::kMHTMSB, l1tmhtemu::kMHTLSB)).to_int());
      l1t::EtSum s(v,
                   l1t::EtSum::EtSumType::kMissingHt,
                   l1tmhtemu::MHT_t(x(l1tmhtemu::kMHTMSB, l1tmhtemu::kMHTLSB)),
                   0,
                   l1tmhtemu::MHTphi_t(x(l1tmhtemu::kMHTPhiMSB, l1tmhtemu::kMHTPhiLSB)).to_int(),
                   0);
      htSums.push_back(s);
    }

    return htSums;
  }

  // Encodes htscalarsum collection onto 1 output link (htscalarsum: scalar sum only from adding 12 jets)
  std::array<std::vector<ap_uint<64>>, 1> encodeHtScalarSums(const edm::View<l1t::EtSum>& htScalarSums) {
    std::vector<ap_uint<64>> htScalarSumWords;

    for (const auto& htScalarSum : htScalarSums)
      htScalarSumWords.push_back(encodeHtScalarSum(htScalarSum));

    std::array<std::vector<ap_uint<64>>, 1> linkData;

    for (size_t i = 0; i < linkData.size(); i++) {
      // Pad etsum vectors -> full packet length (48 frames, but only 1 htsum max)
      htScalarSumWords.resize(1, 0);
      linkData.at(i) = htScalarSumWords;
    }

    return linkData;
  }

  std::vector<l1t::EtSum> decodeHtScalarSums(const std::vector<ap_uint<64>>& frames) {
    std::vector<l1t::EtSum> htScalarSums;

    for (const auto& x : frames) {
      if (not x.test(0))
        break;

      math::XYZTLorentzVector v(0, 0, 0, l1thtemu::ht_t(x(l1thtemu::kPtMSB, l1thtemu::kPtLSB)).to_int());
      l1t::EtSum s(v,
                   l1t::EtSum::EtSumType::kTotalHt,
                   l1thtemu::ht_t(x(l1thtemu::kPtMSB, l1thtemu::kPtLSB)),
                   0,
                   0,
                   0);
      htScalarSums.push_back(s);
    }

    return htScalarSums;
  }


}  // namespace l1t::demo::codecs
