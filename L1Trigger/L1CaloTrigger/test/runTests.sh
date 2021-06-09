#!/bin/bash

#-----------------------------------------------------#
# runTests.sh
#
# Script to run the CMSSW emulator and the Phase 2 firmware-based emulator
# and diff the .txt files that they write to.
#-----------------------------------------------------#

cmsenv

echo "runTests.sh: Running cmsRun on the two emulators..."
cmsRun test_L1EGammaCrystalsProducer.py &


cmsRun test_Phase2L1CaloEGammaEmulator.py &

wait

echo "runTests.sh: Running diff on the output files (<: in CMSSW, >: in firmware-based emulator)..."
diff cmsswProducerL1outputs.txt firmwareEmulatorL1outputs.txt

echo "runTests.sh: Done (<: in CMSSW, >: in firmware-based emulator)"
