# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 09:41:26 2019
@author: GE
"""

import os
import subprocess
from joblib import Parallel, delayed
from tqdm import tqdm
from time import sleep
from datetime import datetime
import shutil
import pathlib

# Define functions =-=====================================
# Run LANDIS-II model function ----
def run_landisii(scenario_path):
    os.chdir(scenario_path)
    if os.name == "nt": # for windows
        identifiers = scenario_path.split('\\')
        command = 'landis-ii-7 scenario.txt | echo off'
    elif os.name == "posix": # for ubuntu, hopporin
        identifiers = scenario_path.split('/')
        command = 'dotnet /home/ge/landis/Core-Model-v7-LINUX-hopporin/build/Release/Landis.Console.dll scenario.txt | echo off'   # for Ubuntu, hopporin
    identifier = f'{identifiers[-3]}_{identifiers[-2]}_{identifiers[-1]}'
    print('Model run started: {}\n'.format(scenario_path))

    try:
        subprocess.call(command, shell = True)
        with open('Landis-log.txt', 'r') as f:
            lines = f.readlines()
        if "Model run is complete" in lines[-1]:
            print('Model run finished: {}\n'.format(scenario_path))
            shutil.copy('Landis-log.txt', os.path.join(db_path, 'Finished_' + identifier + '_landis-log.txt'))
        else:
            print('Error! See Landis-log.txt in {}\n'.format(scenario_path))
            shutil.copy('Landis-log.txt', os.path.join(db_path, 'ErrorLog_' + identifier + '_landis-log.txt'))
    except:
        print('Error in processing {} ...\n'.format(scenario_path))
        failpath = pathlib.Path(os.path.join(db_path, 'FaiedLog_' + identifier + '_landis-log.txt'))
        failpath.touch()
    try:
        # Remove non-necessary output files
        #  2  8 18 28 38 58 86 is needed
        biom_path = os.path.join(scenario_path, "OutputMaps", "biomass")
        biom_list = [p for p in os.listdir(biom_path) if ("-2.tif" not in p) and ("-8.tif" not in p) and ("-18.tif" not in p) and ("-28.tif" not in p) and ("-38.tif" not in p) and ("-58.tif" not in p) and ("-86.tif" not in p)]
        for f in biom_list:
            os.remove(os.path.join(biom_path, f))
        sppagebiom_path = os.path.join(scenario_path, "OutputMaps", "spp-biomass-by-age")
        sppagebiom_list = [p for p in os.listdir(sppagebiom_path) if ("-2.tif" not in p) and ("-8.tif" not in p) and ("-18.tif" not in p) and ("-28.tif" not in p) and ("-38.tif" not in p) and ("-58.tif" not in p) and ("-86.tif" not in p)]
        for f in sppagebiom_list:
            os.remove(os.path.join(sppagebiom_path, f))
    except:
        print('Error in deleting non-necessary data {} ...\n'.format(scenario_path))
        failpath = pathlib.Path(os.path.join(db_path, 'FaiedLog_' + identifier + '_delete_process.txt'))
        failpath.touch()

# Parallel control ------
def usejoblib(job, path_list):
    Parallel(n_jobs=job, verbose = 11)([delayed(run_landisii)(path) for path in path_list])


# Main ===================================================
if __name__ == "__main__":
    CORE_NUMBER = 10 # NUMBER OF CORES FOR MULTIPROCESSING
    scriptdir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(scriptdir); os.chdir('../')
    root_path = os.getcwd()
    # Dropbox path for log monitoring
    pcname = "landis-ws"
    dt_now = datetime.now()
    db_path = os.path.join(f'/home/ge/Dropbox/hg_tmp/Alpha_v3.1.0_{pcname}_{dt_now.strftime("%Y%m%d")}_{os.path.basename(root_path)}_log')
    if not os.path.exists(db_path):
        os.mkdir(db_path)
    # Create the list of scenarios
    os.chdir(root_path)
    mng_list = [path for path in os.listdir() if os.path.isdir(path) and path.startswith('BaU_')]
    clim_list = ['CSIRO-Mk3-6-0_RCP8.5', 'GFDL-CM3_RCP8.5', 'HadGEM2-ES_RCP8.5', 'MIROC5_RCP2.6', 'noCC']
    reg_list = [1]
    Niters = 5
    scenario_path = []
    for reg in reg_list:
        for clim_name in clim_list:
            for mng_name in mng_list:
                for i in range(Niters):
                    scenario_path.append(os.path.join(root_path, mng_name, clim_name, f'reg{reg}_iter{i}'))
                    print(os.path.join(root_path, mng_name, clim_name, 'reg{}'.format(reg)))
    usejoblib(CORE_NUMBER, scenario_path)
