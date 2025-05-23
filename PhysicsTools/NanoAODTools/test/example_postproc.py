#!/usr/bin/env python3
#
# Example of running the postprocessor to skim events with a cut, and 
# adding a new variable using a Module.
#
# Usage:
#   python3 example_postproc.py
#
from PhysicsTools.NanoAODTools.postprocessing.examples.exampleModule import *
# from PhysicsTools.NanoAODTools.postprocessing.examples.higgsinoSkimModule import *
from PhysicsTools.NanoAODTools.postprocessing.examples.oppositeFlavorModule import *

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

# DYJets
# /DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM
fnames = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/70000/B715A9DC-A458-3946-B3F6-34A0A8F44766.root"]

# TTTo2L2Nu:
# 2018: /TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM
# fnames = ["root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18NanoAODv9/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/130000/0804DEBA-97D5-BE46-BB9D-B1125570966E.root"]
# 2017: /TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM

# Signal:
# /SMS-TChiZH-mNLSP200To1500_mLSP0To600_TuneCP5_13TeV-madgraphMLM-pythia8RunIISummer20UL18NanoAODv9FSUL18_FSUL18_106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM
# fnames = []

# Data
fnames = ["root://cms-xrd-global.cern.ch//store/data/Run2018A/SingleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v2/2550000/36ED9511-D46A-0C4F-A485-C2DF1C874906.root"]

p = PostProcessor(outputDir=".",
                  inputFiles=fnames,
                  cut="(nMuon > 0) && (nElectron > 0) && (nJet >= 2) && (MET_pt > 50)",
                  modules=[oppositeFlavourModule()],
                  provenance=True,
                  maxEntries=1000, #just read the first maxEntries events
                  )
p.run()
