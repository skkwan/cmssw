import os

filesToSubmit = [
    "crabJobConfigs/2018/crab_DYJetsToLL_M-50_HT-70to100_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_DYJetsToLL_M-50_HT-100to200_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_DYJetsToLL_M-50_HT-200to400_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_DYJetsToLL_M-50_HT-400to600_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_DYJetsToLL_M-50_HT-600to800_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_DYJetsToLL_M-50_HT-800to1200_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_DYJetsToLL_M-50_HT-1200to2500_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_DYJetsToLL_M-50_HT-2500toInf_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_WJetsToLNu_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_TTTo2L2Nu_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_TTZToLLNuNu_M-10_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_TTWJetsToLNu_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_DoubleMuon_Run2018A_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_DoubleMuon_Run2018B_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_DoubleMuon_Run2018C_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_DoubleMuon_Run2018D_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_EGamma_Run2018A_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_EGamma_Run2018B_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_EGamma_Run2018C_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_EGamma_Run2018D_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_WminusH_HToBB_WToLNu_M-125_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_WplusH_HToBB_WToLNu_M-125_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_WWTo2L2Nu_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_WZTo2Q2L_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_WZTo3LNu_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_ZH_HToBB_ZToLL_M-125_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_ZZ_RunIISummer20UL18NanoAODv9_cfg.py",
    "crabJobConfigs/2018/crab_TChiZH_RunIISummer20UL18NanoAODv9_cfg.py",
]

for f in filesToSubmit:
    command = f"crab submit -c {f}"
    print(command)
    os.system(command)
