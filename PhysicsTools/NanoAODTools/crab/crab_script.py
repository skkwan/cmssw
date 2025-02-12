#!/usr/bin/env python3
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import *

# this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.utils.crabhelper import inputFiles, runsAndLumis

from PhysicsTools.NanoAODTools.postprocessing.examples.exampleModule import *
from PhysicsTools.NanoAODTools.postprocessing.examples.higgsinoSkimModule import *


p = PostProcessor(".",
                  inputFiles(),
                  cut="((nMuon >= 2) || (nElectron >= 2)) && (nJet >= 2) && (MET_pt > 50)",
#                  modules=[higgsinoSkimModule()],
                  modules=[higgsinoSkimModule()],
                  provenance=True,
                  fwkJobReport=True,
                  jsonInput=runsAndLumis())
p.run()

print("DONE")
