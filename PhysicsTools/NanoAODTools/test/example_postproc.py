#!/usr/bin/env python3
#
# Example of running the postprocessor to skim events with a cut, and 
# adding a new variable using a Module.
#
from PhysicsTools.NanoAODTools.postprocessing.examples.exampleModule import *

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

fnames = ["root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18NanoAODv9/DY1JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/120000/13374A29-B61F-7443-AF78-C0C04D479595.root"]

p = PostProcessor(outputDir=".",
                  inputFiles=fnames,
                  cut="((nMuon >= 2) || (nElectron >= 2))",
                  modules=[exampleModuleConstr()],
                  provenance=True,
                  maxEntries=10, #just read the first maxEntries events
                  )
p.run()

