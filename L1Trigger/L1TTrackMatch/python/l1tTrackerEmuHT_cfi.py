import FWCore.ParameterSet.Config as cms

l1tTrackerEmuHT = cms.EDProducer("L1TkHTEmulatorProducer",
    L1TkJetEmulationInputTag = cms.InputTag("l1tTrackJetsEmulation", "L1TrackJets"),
    L1HTCollectionName = cms.string("L1TrackerEmuHT"),
    jet_maxEta = cms.double(2.4),
    jet_minPt = cms.double(3.0),
    jet_minNtracksLowPt = cms.int32(0),
    jet_minNtracksHighPt = cms.int32(0),
    debug = cms.bool(True),
    displaced = cms.bool(False)
)

l1tTrackerEmuHTExtended = cms.EDProducer("L1TkHTEmulatorProducer",
    L1TkJetEmulationInputTag = cms.InputTag("l1tTrackJetsExtendedEmulation", "L1TrackJetsExtended"),
    L1HTCollectionName = cms.string("L1TrackerEmuHTExtended"),
    jet_maxEta = cms.double(2.4),
    jet_minPt = cms.double(3.0),
    jet_minNtracksLowPt = cms.int32(0),
    jet_minNtracksHighPt = cms.int32(0),
    debug = cms.bool(False),
    displaced = cms.bool(True)
)
