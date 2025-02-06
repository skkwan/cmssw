# README.md for CRAB NanoAODTools submission

## Basic steps

This code includes some `.yml` and wrapper `.py` so that we don't need to edit cmsRun and crab config files by hand for each dataset we want to run on.

1. Edit `datasetConfig2018.yml` which lists all the input DAS datasets in a Python dictionary.

   Edit `crab_template_cfg.py` which contains the parameters passed to CRAB, e.g.
   which T2 storage to use, multi-core options, and so on.

2. Now create one CRAB `_cfg.py` per dataset, by running this:
   ```bash
   # edit the top to point to the right dictionary entry in datasetConfig2018.yml
   python parseYaml.py
   ```
4. Set up
   ```bash
   cmsenv
   voms-proxy-init
   ```
5. Create the crab submit files and make sure they look OK:
   ```bash
   python3 parseYamlForBatchSubmit.py
   ```
   Check the resulting config file, e.g. `crabJobConfigs/2018/crab_TTTo2L2Nu_RunIISummer20UL18NanoAODv9_cfg.py`.
6. Submit manually:
   ```bash
   crab submit -c crabJobConfigs/2018/crab_TTTo2L2Nu_RunIISummer20UL18NanoAODv9_cfg.py
   # give it a few minutes and check status
   crab status -d crabJobConfigs/2018/crab_TTTo2L2Nu_RunIISummer20UL18NanoAODv9_cfg.py
   ```
7. (Optional, if running over many datasets, may take a long time) To submit multiple CRAB `_cfg.py` (hereafter referred to as "jobs": one job per `_cfg.py` file)
   edit `yamlForBatchSubmit.yml`. E.g. only run one dataset to see if it works first, by
   commenting out the other lines in this file.

   And make sure that the "string" at the top of `yamlForBatchSubmit.yml` is accurate.
   Submit one CRAB job per uncommented line in `yamlForBatchSubmit.yml`.
   ``` 
   python parseYamlForBatchSubmit.py --execute
   # May need to type GRID password again even if you did voms-proxy-init earlier
   ```

## Troubleshooting CRAB jobs

* Wall clock time exceeded
  * Could be due to broken paths in /hdfs (e.g. to outdated or duplicate datasets) which
    need to be marked as invalid.

* Multi-threading in CRAB config:
  ```
  config.JobType.maxJobRuntimeMin = 300
  config.JobType.numCores = 8
  config.JobType.maxMemoryMB = 9000
  ```