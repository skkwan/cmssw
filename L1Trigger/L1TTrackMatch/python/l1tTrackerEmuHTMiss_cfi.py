import FWCore.ParameterSet.Config as cms

l1tTrackerEmuHTMiss = cms.EDProducer("L1TkHTMissEmulatorProducer",
    L1TkJetEmulationInputTag = cms.InputTag("l1tTrackJetsEmulation", "L1TrackJets"),
    L1MHTCollectionName = cms.string("L1TrackerEmuHTMiss"),
    jet_maxEta = cms.double(2.4),
    jet_minPt = cms.double(3.0),
    jet_minNtracksLowPt = cms.int32(0),
    jet_minNtracksHighPt = cms.int32(0),
    debug = cms.bool(False),
    displaced = cms.bool(False),
    maxNJetsForHT = cms.int32(12),       # Maximum number of jets used for firmware-accurate TrackJets HT computation
)

l1tTrackerEmuHTMissExtended = cms.EDProducer("L1TkHTMissEmulatorProducer",
    L1TkJetEmulationInputTag = cms.InputTag("l1tTrackJetsExtendedEmulation", "L1TrackJetsExtended"),
    L1MHTCollectionName = cms.string("L1TrackerEmuHTMissExtended"),
    jet_maxEta = cms.double(2.4),
    jet_minPt = cms.double(3.0),
    jet_minNtracksLowPt = cms.int32(0),
    jet_minNtracksHighPt = cms.int32(0),
    debug = cms.bool(False),
    displaced = cms.bool(True),
    maxNJetsForHT = cms.int32(12),       # Maximum number of jets used for firmware-accurate TrackJets HT computation  
)
