# parseYaml.py

import yaml
import os

config_file = open("datasetConfig2018.yml")

config = yaml.safe_load(config_file)


#-------------------------------------------------------#

for eraType in ["mc_2018"]:
   for d in config[eraType]["datasets"]:

      dir = "crabJobConfigs/" + str(config[eraType]["year"]) + "/" 
      os.system(f"mkdir -p {dir}")

      newCRABConfFile  = dir + "crab_" + d + "_" + config[eraType]["prodtag"] + "_cfg.py"  
      print(newCRABConfFile)

      requestName      = config["requestname"] + "_" + d + "_" + config[eraType]["prodtag"]  
      dasName          = config[eraType]["datasets"][d] # DAS name
      inputDBS         = config[eraType]["inputDBS"]    # inputDBS
      outputTag        = 'NanoSkimTest_' + d + "_" + config[eraType]["prodtag"]  

      # Make one new template file per dataset
      with open('crab_template_cfg.py', 'r') as templatefile:
         t1 = templatefile.read()
         with open(newCRABConfFile, 'w+') as writefile:

            t2 = t1.replace('REQUEST_NAME',            requestName)
            t3 = t2.replace('DAS_NAME',                dasName)
            t4 = t3.replace('INPUT_DBS',               inputDBS)
            t5 = t4.replace('OUTPUT_TAG',              outputTag)
            writefile.write(t5)

#-------------------------------------------------------#

