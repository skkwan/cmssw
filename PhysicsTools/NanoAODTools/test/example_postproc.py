#!/usr/bin/env python3
#
# Example of running the postprocessor to skim events with a cut, and 
# adding a new variable using a Module.
#
from PhysicsTools.NanoAODTools.postprocessing.examples.exampleModule import *
from PhysicsTools.NanoAODTools.postprocessing.examples.higgsinoSkimModule import *

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

# DY1Jets
# /DY1JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM
# fnames = ["root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18NanoAODv9/DY1JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/120000/13374A29-B61F-7443-AF78-C0C04D479595.root"]


# TTTo2L2Nu:
# 2018: /TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM
# fnames = []
# 2017: /TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM

# Signal:
# /SMS-TChiZH-mNLSP200To1500_mLSP0To600_TuneCP5_13TeV-madgraphMLM-pythia8RunIISummer20UL18NanoAODv9FSUL18_FSUL18_106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM
# fnames = []

# Data
fnames = ["root://cms-xrd-global.cern.ch///store/data/Run2018A/DoubleMuon/NANOAOD/UL2018_MiniAODv2_NanoAODv9-v1/270000/C489C20E-FD93-8B42-9F63-0AB2FB0F5C39.root"]

p = PostProcessor(outputDir=".",
                  inputFiles=fnames,
                  cut="((nMuon >= 2) || (nElectron >= 2)) && (nJet >= 2)",
                  modules=[higgsinoSkimModule()],
                  provenance=True,
                  maxEntries=1000, #just read the first maxEntries events
                  )
p.run()

