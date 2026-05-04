
/**\class L1TrackerHTEmulatorProducer L1TrackerHTEmulatorProducer.cc
 L1Trigger/L1TTrackMatch/plugins/L1TrackerHTEmulatorProducer.cc
 Description: Takes L1TTkJets and performs a integer emulation of Track-based HT, outputting a collection of EtSum 
*/

// system include files
#include <memory>
#include <numeric>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/L1TCorrelator/interface/TkHT.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/L1Trigger/interface/TkJetWord.h"
#include "L1Trigger/L1TTrackMatch/interface/L1TkHTEmulatorProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace l1t;

class L1TkHTEmulatorProducer : public edm::stream::EDProducer<> {
public:
  explicit L1TkHTEmulatorProducer(const edm::ParameterSet&);
  ~L1TkHTEmulatorProducer() override = default;

private:
  virtual void beginJob();
  void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob();

  // ----------member data ---------------------------

  bool debug_ = false;
  bool displaced_;

  std::string L1HTCollectionName_;

  const edm::EDGetTokenT<TkJetWordCollection> jetToken_;
};

L1TkHTEmulatorProducer::L1TkHTEmulatorProducer(const edm::ParameterSet& iConfig)
    : jetToken_(consumes<TkJetWordCollection>(iConfig.getParameter<edm::InputTag>("L1TkJetEmulationInputTag"))) {
  debug_ = iConfig.getParameter<bool>("debug");
  displaced_ = iConfig.getParameter<bool>("displaced");

  // Name of output ED Product
  L1HTCollectionName_ = (std::string)iConfig.getParameter<std::string>("L1HTCollectionName");

  produces<std::vector<EtSum>>(L1HTCollectionName_);

  if (debug_) {
    edm::LogVerbatim("L1TrackerHTEmulatorProducer")
        << "-------------------------------------------------------------------------\n"
        << "====BITWIDTHS====\n"
        << "pt: " << l1t::TkJetWord::TkJetBitWidths::kPtSize << "\n"
        << "-------------------------------------------------------------------------\n";
  }
}

void L1TkHTEmulatorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  std::unique_ptr<std::vector<l1t::EtSum>> HTCollection(new std::vector<l1t::EtSum>(0));

  // L1 track-trigger jets
  edm::Handle<TkJetWordCollection> L1TkJetsHandle;
  iEvent.getByToken(jetToken_, L1TkJetsHandle);
  std::vector<TkJetWord>::const_iterator jetIter;

  if (!L1TkJetsHandle.isValid() && !displaced_) {
    LogError("TkHTEmulatorProducer") << "\nWarning: TkJetCollection not found in the event. Exit\n";
    return;
  }

  if (!L1TkJetsHandle.isValid() && displaced_) {
    LogError("TkHTEmulatorProducer") << "\nWarning: TkJetExtendedCollection not found in the event. Exit\n";
    return;
  }

  // floats used for debugging
  float HT_ = 0;

  l1thtemu::ht_t HT = 0;

  // loop over jets
  int jetn = 0;

  for (jetIter = L1TkJetsHandle->begin(); jetIter != L1TkJetsHandle->end(); ++jetIter) {

    //float tmp_jet_et_ = jetIter->pt();  // FIXME Get Et from the emulated jets
    float tmp_jet_pt_ = jetIter->pt();

    // bool tmp_jet_isDisplaced_ = jetIter->dispflag();

    l1thtemu::ht_t tmp_jet_pt = l1thtemu::digitizeSignedValue<l1thtemu::ht_t>(jetIter->pt(), l1thtemu::kPtSize, l1thtemu::kStepPt);
   
    jetn++;

    if (debug_) {
      edm::LogVerbatim("L1TrackerHTEmulatorProducer")
          << "****JET EMULATION" << jetn << "****\n"
          << "FLOATS ORIGINAL\n"
          << "PT: " << jetIter->pt() << "| ETA: " << jetIter->glbeta() << "| PHI: " << jetIter->glbphi()
          << "| NTRACKS: " << jetIter->nt() << "| COS(PHI): " << cos(jetIter->glbphi())
          << "| SIN(PHI): " << sin(jetIter->glbphi()) << "| Px: " << jetIter->pt() * cos(jetIter->glbphi())
          << "| Py: " << jetIter->pt() * sin(jetIter->glbphi()) << "\n"
          << "AP_INTS RAW\n"
          << "PT: " << jetIter->ptWord() << "| ETA: " << jetIter->glbEtaWord() << "| PHI: " << jetIter->glbPhiWord()
          << "| NTRACKS: " << jetIter->ntWord() << "\n"
          << "AP_INTS NEW\n"
          << "PT: " << tmp_jet_pt << "\n"
          << "AP_INTS NEW TO FLOATS\n"
          << "PT: " << (float)tmp_jet_pt * l1thtemu::kStepPt << "\n"
          << "-------------------------------------------------------------------------\n";
    }

    if (debug_) {
      HT_ += tmp_jet_pt_;
    }


    HT += tmp_jet_pt;

  }  // end jet loop

  // define missing HT

  if (debug_) {
    edm::LogVerbatim("L1TrackerHTEmulatorProducer")
        << "-------------------------------------------------------------------------\n"
        << "====HT FLOATS====\n"
        << "HT: " << HT_ 
        << "\n"
        << "====HT AP_INTS====\n"
        << "HT: " << HT 
        << "\n"
        << "====HT AP_INTS TO FLOATS====\n"
        << "HT: " << (float) HT * l1thtemu::kStepPt << "\n"
        << "-------------------------------------------------------------------------\n";
  }
  //rescale HT to correct output range
  HT = HT / (int)(1 / l1thtemu::kStepPt);

  math::XYZTLorentzVector vectorHt(0, 0, 0, HT);

  EtSum L1HTSum(vectorHt, EtSum::EtSumType::kTotalHt, (int) HT.range(), 0, 0, (int)jetn);

  HTCollection->push_back(L1HTSum);
  iEvent.put(std::move(HTCollection), L1HTCollectionName_);

}  //end producer

void L1TkHTEmulatorProducer::beginJob() {}

void L1TkHTEmulatorProducer::endJob() {}

DEFINE_FWK_MODULE(L1TkHTEmulatorProducer);
