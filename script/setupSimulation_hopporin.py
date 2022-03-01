# -*- coding: utf-8 -*-
"""
Created on 20/01/01
"""

import os
import pandas as pd


def MakeScenarioTxt():
    with open(os.path.join(rootdir, 'default_files', "scenario.txt"), "r", encoding = "utf-8") as fold:
        lines = fold.readlines()
        with open(os.path.join(scenariodir, 'scenario.txt'), 'w', encoding='utf-8') as fnew:
            for line in lines:
                if 'EcoregionsMap' in line:
                    fnew.write(f'EcoregionsMap\t../../../input/eco{regid}_2021-03-01.tif\n')
                elif 'Species' in line:
                    if mngname == 'nograss':
                        fnew.write('Species\t../../../input/species_donan_nograss_2021-08-10.txt\n')
                    else:
                        fnew.write(line)
                elif 'Output Biomass-by-Age' in line:
                    if mngname == 'nograss':
                        fnew.write('"Output Biomass-by-Age"\t../../../ini/output-biomass-by-age_nograss.txt\n')
                    else:
                        fnew.write(line)
                else:
                    fnew.write(line)


def MakeHarvestTxt():
    wt_y = 3 # year for windthrow event
    rm_ys = [4, 5, 6, 7, 17, 27, 37] # year for removin other cohorts
    if mngname == 'nograss':
        fname = os.path.join(rootdir, 'default_files', "BiomassHarvest_nograss_v3.1.0.txt")
    else:
        fname = os.path.join(rootdir, 'default_files', "BiomassHarvest_v3.1.0.txt")
    with open(fname, "r", encoding = "utf-8") as fold:
        lines = fold.readlines()
        with open(os.path.join(scenariodir, 'BiomassHarvest_v3.1.0.txt'), 'w', encoding='utf-8') as fnew:
            for line in lines:
                if 'ManagementAreas_replacehere' in line:
                    fnew.write('ManagementAreas\t../../../input/mng{}_2021-03-09.tif\n'.format(regid))
                elif 'Stands_replacehere' in line:
                    fnew.write('Stands\t../../../input/std{}_2021-03-09.tif\n'.format(regid))
                elif 'HarvestImplementations_replacehere' in line:
                    # for mngid, prescript in zip(mngids, prescripts):
                    for iter in range(len(mngids)):
                        mngid = mngids[iter]
                        fnew.write('>> setting for management area {}\n'.format(mngid))
                        if 'slpl' in mngname:
                            # Plantation setting ============
                            if mngid < 100:
                                # Windthrow + scarification + Plantation
                                fnew.write('{}\t{}\t100%\t{}\t{}\n'.format(mngid, plant_prescripts[iter], wt_y, wt_y))
                                # Remove other trees -----
                                for yr in rm_ys:
                                    fnew.write('{}\t{}\t100%\t{}\t{}\n'.format(mngid, remove_prescripts[iter], yr, yr))
                            # Natural forest setting ---------
                            elif mngid >= 100:
                                # Windthrow + Remove trees
                                fnew.write('{}\tslpl_naturalforest\t100%\t{}\t{}\n'.format(mngid, wt_y, wt_y))
                        else:
                            # Windthrow event setting ======
                            fnew.write('{}\tCLwindthrow\t100%\t{}\t{}\n'.format(mngid, wt_y, wt_y))
                        # Sasa reduction by climate change
                elif 'SasaReduction_replacehere' in line:
                    if climcase != 'noCC':
                        # import pdb; pdb.set_trace()
                        for itr in range(len(mngids)):
                            mngid = mngids[itr]
                            if mngid in [15, 16, 17, 101]:
                                eco = 1
                            elif mngid in [21, 22, 25, 26, 111]:
                                eco = 2
                            elif mngid in [31, 35, 36, 121]:
                                eco = 3
                            fnew.write('>> setting for sasareduction in area {}\n'.format(mngid))
                            for yr in range(2021, 2101):
                                if climcase == 'HadGEM2-ES RCP8.5' and yr == 2100:
                                    clim_gcm_df = clim_df.query('identifier == @climcase and year == 2099 and clusterID == @eco')
                                else:
                                    clim_gcm_df = clim_df.query('identifier == @climcase and year == @yr and clusterID == @eco')
                                if clim_gcm_df['sasaReduct_cut'].values[0] > 0:
                                    sasared = min(int(clim_gcm_df['sasaReduct_cut'].values[0] * 100), 99)
                                    fnew.write('{}\tsasareduct{}\t100%\t{}\t{}\n'.format(mngid, sasared, yr - 2014 + 1, yr - 2014 + 1))
                else:
                    fnew.write(line)


def MakeNECNsuccessionTxt():
    if mngname == 'nograss':
        fname = os.path.join(rootdir, 'default_files', "NECN-succession_donan_nograss_v3.1.0.txt")
    else:
        fname = os.path.join(rootdir, 'default_files', "NECN-succession_donan_v3.1.0.txt")
    with open(fname, "r", encoding = "utf-8") as fold:
        lines = fold.readlines()
        with open(os.path.join(scenariodir, 'NECN-succession_donan_v3.1.0.txt'), 'w', encoding='utf-8') as fnew:
            for line in lines:
                if 'InitialCommunitiesMap_replacehere' in line:
                    fnew.write('InitialCommunitiesMap\t../../../input/ic_reg{}.tif\n'.format(regid))
                elif 'ClimateConfigFile_replacehere' in line:
                    fnew.write('ClimateConfigFile\t../../../ini/clim-gen_{}.txt\n'.format(climname))
                elif 'REGIDHERE' in line:
                    fnew.write(line.replace('REGIDHERE', str(regid)))
                elif 'HarvestReductionParameters_replacehere' in line:
                    # Plant or removing prescripts ------
                    if len(plant_prescripts) > 0:
                        for prescript in plant_prescripts:
                            fnew.write('{}\t{}\t0.0\t0.0\t0.0\t0.0\n'.format(prescript, rmfrac))
                    if len(remove_prescripts) > 0:
                        for prescript in remove_prescripts:
                            fnew.write('{}\t{}\t0.0\t0.0\t0.0\t0.0\n'.format(prescript, rmfrac))
                    if 'cl' in mngname:
                        fnew.write('CLwindthrow\t{}\t0.0\t0.0\t0.0\t0.0\n'.format(rmfrac))
                else:
                    fnew.write(line)


def ReturnHarvestParams():
    if mngname == 'slpl1':
        rmfrac = 1.0
        plant_prescripts = ['slpl1_larikaem', 'slpl1_abiesach', 'slpl1_crypjapo', 'slpl1_betuerma', 'slpl1_fagucren', 'slpl1_larikaem', 'slpl1_abiesach', 'slpl1_betuerma', 'slpl1_larikaem', 'slpl1_abiesach']
        remove_prescripts = ['cutothers_larikaem', 'cutothers_abiesach', 'cutothers_crypjapo', 'cutothers_betuerma', 'cutothers_fagucren', 'cutothers_larikaem', 'cutothers_abiesach', 'cutothers_betuerma', 'cutothers_larikaem', 'cutothers_abiesach']
    elif mngname == 'slpl2':
        rmfrac = 1.0
        plant_prescripts = ['slpl2' for _ in range(10)]
        remove_prescripts = ['cutothers_crypjapo' for _ in range(10)]
    elif mngname == 'cl':
        rmfrac = 0
        plant_prescripts = []
        remove_prescripts = []
    elif mngname == 'nograss':
        rmfrac = 0
        plant_prescripts = []
        remove_prescripts = []
    return [rmfrac, plant_prescripts, remove_prescripts]


# Main ==================================================
if __name__ == "__main__":
    # set root directory
    scriptdir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(scriptdir); os.chdir('../')
    rootdir = os.getcwd()
    # Set params
    mngnames = ['slpl1', 'slpl2', 'cl']
    climnames = ['CSIRO-Mk3-6-0_RCP8.5', 'GFDL-CM3_RCP8.5', 'HadGEM2-ES_RCP8.5', 'MIROC5_RCP2.6', 'noCC'] # 2021.03.09 updated
    regids = [1]
    Niters = 5
    mngids = [15, 16, 17, 21, 22, 25, 26, 31, 35, 36, 101, 111, 121] # 2021.03.09 updated
    # sasa reduction
    clim_df = pd.read_csv(os.path.join(rootdir, 'data', 'climate', 'climate_yearly_selected.csv'))
    print('Root directory: %s'%rootdir)
    for mngname in mngnames:
        mngdir = os.path.join(rootdir, 'BaU_{}'.format(mngname))
        if not os.path.isdir(mngdir):
            os.mkdir(mngdir)
        rmfrac, plant_prescripts, remove_prescripts = ReturnHarvestParams()
        for climname in climnames:
            climdir = os.path.join(mngdir, climname)
            climcase = climname.replace('_', ' ')
            if not os.path.isdir(climdir):
                os.mkdir(climdir)
            for regid in regids:
                for i in range(Niters):
                    scenariodir = os.path.join(climdir, f'reg{regid}_iter{i}')
                    print(scenariodir)
                    if not os.path.isdir(scenariodir):
                        os.mkdir(scenariodir)
                    # Process 0. Make scenario.txt ==============================
                    try:
                        MakeScenarioTxt()
                    except:
                        print('failed to make scenario.txt')

                    # Replace 1: BiomassHarvest.txt ==================================
                    try:
                        MakeHarvestTxt()
                    except:
                        print('failed to make BiomassHarvest.txt')

                    # Replace 2: NECN-succession.txt ==================================
                    try:
                        MakeNECNsuccessionTxt()
                    except:
                        print('failed to make NECN-succession.txt')

