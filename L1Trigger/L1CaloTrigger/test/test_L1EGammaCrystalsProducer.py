import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process("L1AlgoTest",eras.Phase2C4)

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )


#Search string: `dataset=/TT*/Phase2HLTTDRWinter20DIGI*/GEN-SIM-DIGI-RAW`
# Picking a PU 200 dataset: `dataset=/TT_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW`

# One file (10 GB)
#`xrdcp root://cmsxrootd.fnal.gov///store/mc/Phase2HLTTDRWinter20DIGI/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW/PU200_110X_mcRun4_realistic_v3-v2/110000/005E74D6-B50E-674E-89E6-EAA9A617B476.root .`

# file dataset=/SinglePion_PT0to200/Phase2HLTTDRWinter20DIGI-NoPU_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:/eos/user/s/skkwan/phase2RCTDevel/SinglePion_PT0to200_NoPU_04AC207E-AF58-C04A-9F90-746DDC628248.root'),
                            inputCommands = cms.untracked.vstring(
                                "keep *"
                            )
                        )


# --------------------------------------------------------------------------------------------                                                    
#                                                                                                                                                            
# ----   Run the relevant algorithms

# ---- Global Tag :
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')

# Choose a geometry
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

# Add HCAL Transcoder
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')


# --------------------------------------------------------------------------------------------
#
# ----    Produce the L1EGCrystal clusters using Emulator

process.load('L1Trigger.L1CaloTrigger.L1EGammaCrystalsEmulatorProducer_cfi')

process.pL1EG = cms.Path( process.L1EGammaClusterEmuProducer )

process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "l1egCrystalTest.root" ),
    outputCommands = cms.untracked.vstring(
                    "keep *_L1EGammaClusterEmuProducer_*_*",
                    "keep *_TriggerResults_*_*",
                    "keep *_simHcalTriggerPrimitiveDigis_*_*",
                    "keep *_EcalEBTrigPrimProducer_*_*"
                    )
)


process.end = cms.EndPath( process.Out )

process.schedule = cms.Schedule(process.pL1EG, process.end)

dump_file = open("dump_file.py", "w")
dump_file.write(process.dumpPython())
