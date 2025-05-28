# crab_template_cfg.py
# Variables like 'REQUEST_NAME' are replaced using a configuration .yml
# and the parseYaml.py wrapper

from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config

config = Configuration()

config.section_("General")
config.General.workArea     = 'crab_projects'
config.General.requestName  = 'REQUEST_NAME'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName   = 'PSet.py'
config.JobType.scriptExe  = 'crab_script.sh'
# config.JobType.scriptArgs = ['script=CRAB_SCRIPT_PYTHON_NAME']
config.JobType.inputFiles = ['crab_script.py'] # separated by list

config.section_("Data")
config.Data.inputDataset = 'DAS_NAME'
config.Data.inputDBS     = 'INPUT_DBS'
config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob  = 1
# config.Data.splitting    = 'Automatic'

config.Data.outLFNDirBase    = '/store/user/skkwan/HiggsinoSkimAllFlav'
config.Data.publication      = True
config.Data.outputDatasetTag = 'OUTPUT_TAG'

#config.JobType.maxJobRuntimeMin = 300
#config.JobType.numCores = 8
#config.JobType.maxMemoryMB = 9000

config.section_("Site")
config.Site.storageSite = "T2_US_Wisconsin"

