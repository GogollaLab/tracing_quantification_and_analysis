# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 15:22:01 2019

@author: gehrlach
"""

import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import filedialog
import os



def read_and_clean_data(csv):
    """ Read in the CSV file that was created by python script 'summary_RV'. It adds all countings from 
    all sections together.
    This function reads the csv, drops unnecassary columns, remove entries where the count was zero 
    """
    
    df = pd.DataFrame()
    df = pd.read_csv(csv)
    df = df[df["Count"] != 0]
    df = df.drop(df.columns[0],1)
    
    #add metadata here
    df['genotype'] = genotype
    df['tracing_from'] = tracing_from
    df['animal'] = animal_number
    
    return df


def remove_starter_volume(df, start_Bregma, end_Bregma):
    """ 
    cleaning of dataframe. Makes values positive, removes starter volume, adds injection center and 
    percent total output to df 
    """
    
    #bring everything into the positive numbers
    df['Bregma'] += 2
    start_Bregma += 2
    end_Bregma +=2
    injection_center = start_Bregma - (start_Bregma - end_Bregma)/2 -2
    
    #remove starter cells from insula and claustrum at from start to end bregma here
    starter_volume = df[(df['Label'].isin(['AIP','AID', 'GI', 'DI', 'AIV', 'AI', 'claustrum'])) & 
                        (df['Bregma'] < start_Bregma) & (df['Bregma'] > end_Bregma)]
    cleaned_df = pd.concat([df,starter_volume,starter_volume]).drop_duplicates(keep=False)
    cleaned_df['Bregma'] -= 2
    cleaned_df['injection center'] = injection_center
    
    cleaned_df['percent_total_input'] = cleaned_df['Count'].apply(lambda x: 100*(x/cleaned_df['Count'].sum()))
    
    cleaned_df = cleaned_df.assign(cell_density = lambda x: x['Count'].div(x['Total Area']))
    return cleaned_df


def exclude_regions(df):
    """ 
    excludes regions with which there is no real connectivity, but due to off target signal shows
    artifical signals. manual curation of Striatum, CeA and Cerebellum. 
    """
    
    regions = ['CPu','NAcc', 'Nacc', 'NAc', 'Nac','NAcsh', 'Nacsh', 'NaSh','IPAC','Tu', 
               'CB','CPu', 'CPU', 'Cpu','2Cb', 'Int', 'IntA', 'PM','CeA', 'CeM', 'CeC', 'CeL',]
    
    excluded_idx = df[df['Label'].isin(regions)].index
    df.drop(excluded_idx, inplace=True)
            
    return df

def cell_count_per_ROI(df):
    """
    total cell count across all bregmas per ROI
    """
    
    df = df.groupby(['Label'])
    cell_count = df['Count'].sum()
    cell_count_percent = (cell_count / cell_count.sum() *100.00)
    return cell_count, cell_count_percent

def cell_count_per_Bregma(df):
    """
    total cell count across all ROIs but this time per Bregma
    """
    
    df = df.groupby(['Bregma'])
    cell_count = df['Count'].sum()
    return cell_count

def pivot_counts(df):
    """
    First, this function delets every column from the dataframe, except Bregma and Cell count, then it creates 
    a pivot table summing the counts and filling missing values with 0
    
    shows number of cells per ROI sorted by bregma 
    """
    df = df.drop('Total Area', 1)
    df = df.drop('Average Size',1)
    df = df.drop('Section',1)
    df = df.drop('animal',1)
    df = df.drop('genotype',1)
    df = df.drop('tracing_from',1)
    df = df.drop('cell_density',1)
    df = df.drop('percent_total_input',1)
    df = df.drop('injection center',1)
    df = df.pivot_table(df, index=["Bregma"], columns=["Label"], aggfunc=[np.sum], fill_value=0 )
    
    return df

def pivot_percent(df):
    """
    First, this function delets every column from the dataframe, except Bregma and Percent Input, then it creates 
    a pivot table summing the counts and filling missing values with 0
    
    shows number of cells per ROI sorted by bregma 
    """
    
    df = df.drop('Total Area', 1)
    df = df.drop('Average Size',1)
    df = df.drop('Section',1)
    df = df.drop('animal',1)
    df = df.drop('genotype',1)
    df = df.drop('tracing_from',1)
    df = df.drop('Count',1)
    df = df.drop('cell_density',1)
    df = df.drop('injection center',1)
    df = df.pivot_table(df, index=["Bregma"], columns=["Label"], aggfunc=[np.sum], fill_value=0 )
    
    return df


def pivot_density(df):
    """
    First, this function delets every column from the dataframe, except Bregma and density, then it creates 
    a pivot table summing the counts and filling missing values with 0
    
    shows number of cells per ROI sorted by bregma 
    """
    df = df.drop('Total Area', 1)
    df = df.drop('Average Size',1)
    df = df.drop('Section',1)
    df = df.drop('animal',1)
    df = df.drop('genotype',1)
    df = df.drop('tracing_from',1)
    df = df.drop('Count',1)
    df = df.drop('percent_total_input',1)
    df = df.drop('injection center',1)
    df = df.pivot_table(df, index=["Bregma"], columns=["Label"], aggfunc=[np.sum], fill_value=0 )
    
    return df


def cell_density_ROI(df):
    """
    groups dataframe by ROI name and calculates the mean average cell density
    """
    
    grouped_ROI_density = df.groupby(['Label']) 
    average_cell_density = grouped_ROI_density['cell_density'].mean()
    return average_cell_density


def cell_density_Bregma(df):
    """
    groups dataframe by Bregma level name and calculates the mean average cell density
    """
    
    grouped_Bregma_density = df.groupby(['Bregma']) 
    average_cell_density_bregma = grouped_Bregma_density['cell_density'].mean()
    return average_cell_density_bregma


def sort_into_hierarchies(density,strength):
    """ 
    creation of dictionaries for clustering of ROI names into higher and lower hierarchies. 
    """
    
    higher_hierarchy_nomenclature = {'Insular Cortex': ['AI', 'AID', 'AIV', 'GI', 'DI', 'AIP'],
                                     'Sensory Cortex': ['S1', 'S2', 'Au1', 'AuD', 'Auditory cortex', 'auditory cortex', 'AuV', 'V1', 'V2L', 'Visual cortex', 'Visual cortex-1', 'Visual Cortex', 'visual cortex', 'V1B', 'V1M', 'VM-1', 'V2MM'],
                                     'Motor Cortex': ['M1', 'M2'],
                                     'Prefrontal Cortex':['FrA', 'MO', 'VO' , 'DLO', 'LO', 'A32', 'A24a', 'cg', 'Cg1', 'Cg2', 'A25', 'A24b', 'Fr3', 'DP', 'IL', 'PrL'],
                                     'Association Cortex': ['parietal cortex', 'Ect','PRh', 'Prh', 'DLEnt', 'DIEnt', 'VIEnt', 'VlEnt', 'MEnt', 'TeA','A29c', 'A29', 'A30', 'a30', 'A29a', 'A29b', 'CEnt', 'DlEnt', 'VLEnt', 'VIEnt', 'Vlent', 'LEnt'],
                                     'Piriform Cortex':['Pir'],
                                     'Olfactory areas': ['CxA','olfactory bulb', 'olfacotry bulb', 'AON', 'tenia tecta', 'taenia tecta', 'teania tectum', 'taenia tectum', 'endopiriform nucleus', 'endopirform nucleus', 'endopiriform nuncleus', 'EndoPir','Endo Pir', 'EnoPir', 'VTT', 'LOT', 'DTT'],
                                     'Claustrum': ['claustrum'],
                                     'Striatum':['NAcc', 'NAc', 'Nac', 'Nacc', 'Nacsh','NAcsh', 'NaSh', 'IPAC','Tu','CPu', 'CPU', 'Cpu', 'LSS', 'CB', 'LSD'],
                                     'Pallidum':['VP','GP','BNST', 'BST', 'MS'],
                                     'Basal Forebrain':['SIB', 'diagonal band', 'HDB', 'LDB', 'B'],
                                     'Hippocampus': ['vHip', 'vHIp', 'AHiP', 'DG', 'PaS', 'DS', 'POST', 'Post', 'PROST', 'PrS', 'Py', 'SubC', 'VS', 'STr', 'Str'],
                                     'Amygdala': ['CeA', 'CeM', 'CeC', 'CeL', 'PMCo', 'PLCo','AA', 'BLP', 'LA', 'BMA', 'BMP', 'EA','BLA', 'BLS', 'BLv', 'BLV', 'CeA', 'MeA', 'MeAD', 'ACo','AHiAL', 'AHiPL', 'I', 'IM', 'MeP', 'STIA', 'APir','RAPir', 'RaPir'],    
                                     'Hypothalamus':['Arc', 'AHA', 'AH','LH', 'MHb', 'STh', 'PSTh', 'VMH', 'MPOM', 'MPA', 'hypothalamus', 'F', 'dorsal hypothalamus','AHC', 'LPO' ,'MnPO', 'MCLH', 'ZI', 'ZIR', 'LHb', 'LHB', 'RML', 'PH', 'SubI', 'Subl', 'Sub1', 'PHD', 'RI', 'RMM', 'PHA'],
                                     'Thalamus':['DM','AD','MD','MGN', 'MGm', 'MG', 'MGV', 'MGv', 'medial geniculate', 'VPL', 'VPM', 'VL', 'VM', 'VPPC', 'VPLpc', 'VPLPC', 'VPMpc', 'VPLPc', 'VPMpv', 'CM', 'AM', 'AV', 'Po', 'PO', 'PV', 'PVA', 'PaF', 'PIL', 'PoT', 'LD', 'IAM', 'Re', 'Rh', 'SPFPC', 'SPFp', 'SPF', 'medial geniculate', 'CM', 'REth', 'Rt', 'Sub', 'SubD', 'PC', 'SG', 'sm', 'Gus', 'DLG', 'Eth', 'PrG', 'SO', 'SubV', 'CL', 'DT', 'DTg', 'LP', 'LVe', 'RRe', 'VA'],
                                     'Midbrain':['substsantia nigra', 'Substantia nigra', 'SN', 'cp', 'SNc', 'SNC', 'SNr', 'SNR', 'dorsal raphe', 'median raphe', 'medial raphe', 'isRt', 'CIC', 'PAG', 'DLPAG', 'VTA', 'VTAR', 'RLi', 'Rli', 'mRt','PVG', 'EW', 'LDTgV', 'PrEW', 'VLPAG', 'p1Rt', 'PBP', 'DA8', 'DAB', 'RPC', 'RP', 'PR', 'APT', 'APTD', 'BIC', 'ECIC','DpG', 'DR', 'IP', 'IPD', 'IRt', 'DCIC', 'MnR', 'RMC', 'RMg', 'RRF', 'CI', 'CLi', 'PN', 'Pa4', 'PaR', 'RPa', 'SubB', 'InC', 'PMnR', 'IF', 'InG', 'PIF', 'RM', 'IPR', 'DpWh', 'PMnT', 'SuG', 'In', 'DpMe'],
                                     'Hindbrain':['PBN', 'parabrachial nucleus', 'NTS', 'MPB', 'LPB', 'PCRtA', 'PTg', 'PrCnF', 'SPTg', 'PnO', 'PnV', 'MPL', 'PL', 'KF', '2Cb', 'SubCV', '5', 'VeCb', 'CAT', 'CnF', 'DMTg', 'MiTg', 'PBG', 'LSO', 'PLV', 'PPTg', 'Pn', 'PnC', 'PNV', 'Pr5', 'RtTg', 'SPO', 'Tz', 'VLL', '7N', 'B9', 'Bar', 'CG', 'DPGi', 'Gi', 'GiA', 'IRtA', 'Int', 'IntA', 'LC', 'LPGi', 'LPGiA', 'MSO', 'MVe', 'MVeMC', 'MVePC', 'P7', 'PM', 'Ve', 'Pa5', 'Pr', 'Su5', 'VTg', 'InCo', 'MTg', 'Sp5', 'SuVe'] }
    
    
    lower_order_hierarchy = {
    
           #Sensory Cortex        
            'Somatosensory cortex, primary': ['S1'],
            'Somatosensory cortex, secondary': ['S2'],
            'Auditory cortex, primary': ['Au1', 'Auditory cortex', 'auditory cortex'],
            'Auditory cortex, secondary or higher': ['AuD', 'AuV'],
            'Visual cortex, primary': ['V1', 'V1B', 'V1M', 'VM-1', 'Visual cortex', 'Visual cortex-1', 'Visual Cortex', 'visual cortex'],
            'Visual cortex, secondary or higher': ['V2L', 'V2ML','V2MM'],
            'Piriform Area':['Pir'],
            
            #Motor Cortex
            'Motor cortex, primary': ['M1', 'Fr3'],
            'Motor cortex, secondary': ['M2'],
            
            #Prefrontal Areas
            'Frontal Association Cortex':['FrA'],
            'Medial Prefrontal Cortex': ['A24a', 'A25', 'A32', 'IL', 'PrL', 'DP', 'Cg2'],
            'Orbitofrontal Cortex': ['VO', 'LO', 'MO', 'DLO'],
            'Cingulate Cortex': ['A24b', 'cg', 'Cg1'],
            
            #Association Cortex
            'Parietal Cortex':['LPtA', 'MPtA', 'PtA', 'parietal cortex'],
            'Ectorhinal Cortex':['Ect'],
            'Perirhinal Cortex': ['PRh', 'Prh'],
            'Temporal Association Cortex':['TeA'],
            'Rertosplenial Cortex': ['A30', 'a30', 'A29c', 'A29b', 'A29a', 'A29', 'RSA', 'RSG', 'RSGa', 'RSGb' ],
            'Enthorhinal Cortex': ['DLEnt', 'Dlent', 'DIEnt', 'DlEnt', 'VIEnt', 'Vlent', 'VlEnt', 'MEnt', 'CEnt', 'LEnt'],
            
            #Insular Ctx
            'Agranular Insular Cortex': ['AI'],
            'Agranular Insular Cortex, dorsal': ['AID'],
            'Agranular Insular Cortex, ventral': ['AIV'],
            'Agranular Insular Cortex, posterior': ['AIP'],
            'Dysgranular Insular Cortex': ['DI'],
            'Granular Insular Cortex': ['GI'],
            
            #Olfactory Areas
            'Anterior olfactory nucleus': ['AON', 'AOP', 'AOM', 'AOL', 'AOD'],
            'Piriform-amygdalar area':['CxA'],
            'Endopiriform Nucleus': ['DEn', 'VEn', 'endopiriform nucleus', 'endopirform nucleus', 'endopiriform nuncleus', 'EndoPir', 'Endo Pir', 'EnoPir'],
            'tenia tecta': ['DTT', 'VTT', 'tenia tecta', 'taenia tecta', 'teania tectum', 'taenia tecta', 'taenia tectum'],
            'Nucleus of the lateral olfactory tract':['LOT'],
            'Main olfactory bulb':['olfactory bulb', 'olfacotry bulb'],
                        
            #Claustrum
            'Claustrum': ['claustrum'],
            
            #Striatum
            'Caudate Putamen':['CPu', 'CPU', 'Cpu', 'LSS', 'LSD'],
            'Nucleus Accumbens, core':['NAcc', 'Nacc', 'NAc', 'Nac'],
            'Nucleus Accumbens, shell': ['NAcsh', 'Nacsh', 'NaSh'],
            'Fundus of Striatum (IPAC)':['IPAC'],
            'olfactory tubercle':['Tu', 'CB'],
            
            #Pallidum
            'Ventral Pallidum':['VP'],
            'Bed Nuclei of the Stria Terminalis (BNST)': ['BNST', 'BST'],
            'Globus Pallidus': ['GP', 'LGP', 'MGP'],
            'Medial Pallidum': ['MS'],
            
            #Hippocampus
            'Hippocampus (CA1,CA2,CA3)': ['CA1', 'CA2', 'CA3', 'fi', 'LMol', 'Or', 'Py', 'Rad', 'STr', 'Str', 'SLu', 'alv', 'vHip', 'vHIp', 'AHiP', 'PaS', 'DS', 'POST', 'PROST', 'Post', 'PrS', 'SubC', 'VS'],
            'Hippocampus(dentate gyrus)': ['DG', 'Mol'],
            
            #Amygdala
            'Lateral Amygdala' :['LA'],
            'Basolateral nucleus': ['BLA','BLP','BL', 'BLv', 'BLV', 'BLS'],
            'Basomedial nucleus': ['BMA', 'BMP', 'IM'],
            'Cortical Amygdala': ['PMCo', 'PLCo', 'ACo'],
            'Anterior Amygdaloid Area':['AA', 'AAD', 'AAV', 'I'],
            'Amygdalopiriform Transition Area': ['RAPir', 'RaPir', 'APir', 'AHiAL', 'AHiPL'],
            'Central Nucleus of Amygdala':['CeC', 'CeL', 'CeA', 'CeM', 'CeMAD', 'CeMAV', 'Ce', 'STIA'],            
            'Medial Amygdala': ['MeA', 'MeAD', 'MeAV', 'MePD', 'MeP'],
            'Extended Amygdala': ['EA'],
                        
            #Basal Forebrain
            'Basal Forebrain':['MDB', 'LDB', 'diagonal band', 'VDB', 'HDB', 'SIB', 'SI/B', 'SI', 'B'],
            
            #Hypothalamus
            'Hypothalamic lateral zone': ['Arc', 'LH', 'LHA', 'STh', 'PSTh', 'MCLH', 'ZI', 'ZIR', 'F', 'hypothalamus', 'dorsal hypothalamus', 'LPO', 'RML', 'SO', 'SubI', 'Subl','Sub1', 'RMM'],
            'Hypothalamic medial zone' :['MnPO', 'MPO', 'MPA', 'MPOM', 'AHC', 'MHb', 'PH', 'AHA', 'VMH', 'PHD', 'RI', 'PHA', 'AH'],
                        
            #Thalamus
            'Habenular nuclei': ['LHb', 'LHB'],
            'Ventral posterolateral complex ':['VPL','VPM','VL','VM', 'VPPC', 'VPLpc', 'VPMpc' ,'VPLPC', 'VPLPc', 'Gus', 'RRe'],
            'Medial geniculate nuclei': ['MGN','MGM','MGv', 'MGV', 'MGm','medial geniculate', 'MG'],            
            'Posterior complex of the thalamus':['Po', 'PO', 'PoT', 'PIL', 'REth', 'SPFPC', 'SPFp', 'SG', 'Eth', 'SPF', 'DT', 'DTg'],
            'Mediodorsal thalamic nuclei': ['DM', 'MD', 'MDC','MDL','MDM', 'Sub', 'SubD', 'LP'],
            'Anterodorsal thalamic nucleus': ['AD', 'AM', 'AV', 'LD', 'DLG', 'PrG'],
            'Intralaminar nuclei of the dorsal thalamus': ['Re', 'Rh', 'CM', 'PF', 'IAM', 'PaF', 'Rt', 'PC', 'SubV', 'CL', 'VA'],
            'Paraventricular thalamic nucleus': ['PVT', 'PVA', 'PV', 'PVG', 'sm'],
                                             
            #Midbrain
            'Periaqueductal gray':['PAG', 'EW', 'PrEW', 'VLPAG', 'DLPAG', 'InC', 'In'],     
            'Substantia Nigra': ['SNc', 'SNC', 'SNr', 'SNR','SN', 'substantia nigra', 'Substantia nigra', 'cp', 'PR', 'IP', 'IPD', 'RMC', 'CI', 'PN'],
            'Ventral tegmental area': ['VTA', 'PBP', 'Pa4', 'PIF', 'VTAR'],
            'Raphe Nuclei': ['dorsal raphe', 'DR', 'median raphe', 'medial raphe', 'RLi', 'Rli', 'MnR', 'PMnR', 'PMnT', 'RMg', 'CLi', 'RPa', 'IF', 'RM', 'IPR'],
            'Midbrain reticular nucleus':['isRt', 'mRt', 'p1Rt', 'RPC', 'RP', 'DA8', 'DAB', 'APT', 'APTD', 'BIC', 'IRt', 'RRF', 'PaR', 'SubB', 'DpMe'],
            'Superior Colliculus':['InG', 'DpWh', 'SuG', 'DpG'],
            'Inferior colliculus':['CIC', 'ECIC', 'DCIC'],
            
            #Hindbrain
            'Parabrachial nucleus': ['PBN', 'LPB', 'MPB', 'parabrachial nucleus'],
            'Pedunculopontine nucleus': ['PTg', 'LDTgV', 'SPTg', 'MPL', 'PL', 'KF', 'PrCnF', 'PnO', 'PnC', 'InCo', 'Pn', 'PnV', 'VeCb', 'CnF', 'CAT', 'DMTg', 'MiTg', 'MTg', 'PBG', 'LSO', 'PLV', 'PPTg', 'VTg', 'PNV', 'Pr5', 'RtTg', 'SPO', 'Tz', 'VLL', 'B9', 'Bar', 'CG', 'Gi', 'GiA', 'LC', 'MSO', 'Su5'],
            'Medulla': ['NTS', 'PCRtA', '7N', 'DPGi', 'IRtA', 'LPGi', 'LPGiA', 'LVe', 'MVe', 'MVeMC', 'MVePC', 'P7', 'Pa5', 'Pr', 'Ve', 'Sp5', 'SuVe'],
            'Cerebellum':['2Cb', 'Int', 'IntA', 'PM'],
            'Motor trigeminal nucleus':['5', 'SubCV']
            }
    
    #Check if there are ROIs in the data, that don´t fit the terms above
    missing_rois_lower = []
    for rois in density.index:
        if rois in [j for i in list(lower_order_hierarchy.values()) for j in i]:
            continue
        else:
            missing_rois_lower.append(rois)
            
    missing_rois_higher = []
    for rois in density.index:
        if rois in [j for i in list(higher_hierarchy_nomenclature.values()) for j in i]:
            continue
        else:
            missing_rois_higher.append(rois)
    
    
    counter = pd.Series().astype('float64') 

    # average of ROI densities, if they are grouped in lower order hierarchy
    for lower_region, ROI_name in lower_order_hierarchy.items():
        counter[lower_region] = 0.0
        
    index_counter = 0
    for lower_region, ROI_name in lower_order_hierarchy.items():
    
        for index in density.index:
            if index in ROI_name:
                print(index)
                counter[lower_region] = counter[lower_region] + density[index]
                print(density[index])
                index_counter += 1
        
        if index_counter == 0:
            counter[lower_region] = 0.0
        else:
            counter[lower_region] = (counter[lower_region] / index_counter) 
        index_counter=0
    counter.rename('average cell density')
    counter_density = counter.copy()
   
       
    #Sort projection strength (% of total output pixel) into lower order hierarchy
    for lower_region, ROI_name in lower_order_hierarchy.items():
        counter[lower_region] = 0

    index_counter = 0
    for lower_region, ROI_name in lower_order_hierarchy.items():
    
        for index in strength.index:
            if index in ROI_name:
                counter[lower_region] = counter[lower_region] + strength[index]
                index_counter += 1
        if index_counter == 0:
            counter[lower_region] = 0.0
       
        index_counter=0 
    
    counter.rename('percent of total input')
    counter_strength = counter.copy()
    
    
    dataframe = pd.concat([counter_density, counter_strength], axis=1)
    column_names = ['average cell density [cells/um²]','% of total inputs']
    dataframe.columns= column_names
    
    return dataframe, missing_rois_lower, missing_rois_higher



def cell_count_ROI_higher_hierarchy_in_percent(cell_count_ROI):
    """
    clusters ROIs into highest abstraction (e.g. amygdala, sensory cortex,etc.)
    """
    
    counter = pd.Series().astype('float64')
    total_cells = cell_count_ROI.sum()      
    higher_hierarchy_nomenclature = {'Insular Cortex': ['AI', 'AID', 'AIV', 'GI', 'DI', 'AIP'],
                                     'Sensory Cortex': ['S1', 'S2', 'Au1', 'AuD', 'Auditory cortex', 'auditory cortex', 'AuV', 'V1', 'V2L', 'Visual cortex', 'Visual cortex-1', 'Visual Cortex', 'visual cortex', 'V1B', 'V1M', 'VM-1', 'V2MM'],
                                     'Motor Cortex': ['M1', 'M2'],
                                     'Prefrontal Cortex':['FrA', 'MO', 'VO' , 'DLO', 'LO', 'A32', 'A24a', 'cg', 'Cg1', 'Cg2', 'A25', 'A24b', 'Fr3', 'DP', 'IL', 'PrL'],
                                     'Association Cortex': ['parietal cortex', 'Ect','PRh', 'Prh', 'DLEnt', 'DIEnt', 'VIEnt', 'VlEnt', 'MEnt', 'TeA','A29c', 'A29', 'A30', 'a30', 'A29a', 'A29b', 'CEnt', 'DlEnt', 'VLEnt', 'VIEnt', 'Vlent', 'LEnt'],
                                     'Piriform Cortex':['Pir'],
                                     'Olfactory areas': ['CxA','olfactory bulb', 'olfacotry bulb', 'AON', 'tenia tecta', 'taenia tecta', 'teania tectum', 'taenia tectum', 'endopiriform nucleus', 'endopirform nucleus', 'endopiriform nuncleus', 'EndoPir','Endo Pir', 'EnoPir', 'VTT', 'LOT', 'DTT'],
                                     'Claustrum': ['claustrum'],
                                     'Striatum':['NAcc', 'NAc', 'Nac', 'Nacc', 'Nacsh','NAcsh', 'NaSh', 'IPAC','Tu','CPu', 'CPU', 'Cpu', 'LSS', 'CB', 'LSD'],
                                     'Pallidum':['VP','GP','BNST', 'BST', 'MS'],
                                     'Basal Forebrain':['SIB', 'diagonal band', 'HDB', 'LDB', 'B'],
                                     'Hippocampus': ['vHip', 'vHIp', 'AHiP', 'DG', 'PaS', 'DS', 'POST', 'Post', 'PROST', 'PrS', 'Py', 'SubC', 'VS', 'STr', 'Str'],
                                     'Amygdala': ['CeA', 'CeM', 'CeC', 'CeL', 'PMCo', 'PLCo','AA', 'BLP', 'LA', 'BMA', 'BMP', 'EA','BLA', 'BLS', 'BLv', 'BLV', 'CeA', 'MeA', 'MeAD', 'ACo','AHiAL', 'AHiPL', 'I', 'IM', 'MeP', 'STIA', 'APir','RAPir', 'RaPir'],    
                                     'Hypothalamus':['Arc', 'AHA', 'AH','LH', 'MHb', 'STh', 'PSTh', 'VMH', 'MPOM', 'MPA', 'hypothalamus', 'F', 'dorsal hypothalamus','AHC', 'LPO' ,'MnPO', 'MCLH', 'ZI', 'ZIR', 'LHb', 'LHB', 'RML', 'PH', 'SubI', 'Subl', 'Sub1', 'PHD', 'RI', 'RMM', 'PHA'],
                                     'Thalamus':['DM','AD','MD','MGN', 'MGm', 'MG', 'MGV', 'MGv', 'medial geniculate', 'VPL', 'VPM', 'VL', 'VM', 'VPPC', 'VPLpc', 'VPLPC', 'VPMpc', 'VPLPc', 'VPMpv', 'CM', 'AM', 'AV', 'Po', 'PO', 'PV', 'PVA', 'PaF', 'PIL', 'PoT', 'LD', 'IAM', 'Re', 'Rh', 'SPFPC', 'SPFp', 'SPF', 'medial geniculate', 'CM', 'REth', 'Rt', 'Sub', 'SubD', 'PC', 'SG', 'sm', 'Gus', 'DLG', 'Eth', 'PrG', 'SO', 'SubV', 'CL', 'DT', 'DTg', 'LP', 'LVe', 'RRe', 'VA'],
                                     'Midbrain':['substsantia nigra', 'Substantia nigra', 'SN', 'cp', 'SNc', 'SNC', 'SNr', 'SNR', 'dorsal raphe', 'median raphe', 'medial raphe', 'isRt', 'CIC', 'PAG', 'DLPAG', 'VTA', 'VTAR', 'RLi', 'Rli', 'mRt','PVG', 'EW', 'LDTgV', 'PrEW', 'VLPAG', 'p1Rt', 'PBP', 'DA8', 'DAB', 'RPC', 'RP', 'PR', 'APT', 'APTD', 'BIC', 'ECIC','DpG', 'DR', 'IP', 'IPD', 'IRt', 'DCIC', 'MnR', 'RMC', 'RMg', 'RRF', 'CI', 'CLi', 'PN', 'Pa4', 'PaR', 'RPa', 'SubB', 'InC', 'PMnR', 'IF', 'InG', 'PIF', 'RM', 'IPR', 'DpWh', 'PMnT', 'SuG', 'In', 'DpMe'],
                                     'Hindbrain':['PBN', 'parabrachial nucleus', 'NTS', 'MPB', 'LPB', 'PCRtA', 'PTg', 'PrCnF', 'SPTg', 'PnO', 'PnV', 'MPL', 'PL', 'KF', '2Cb', 'SubCV', '5', 'VeCb', 'CAT', 'CnF', 'DMTg', 'MiTg', 'PBG', 'LSO', 'PLV', 'PPTg', 'Pn', 'PnC', 'PNV', 'Pr5', 'RtTg', 'SPO', 'Tz', 'VLL', '7N', 'B9', 'Bar', 'CG', 'DPGi', 'Gi', 'GiA', 'IRtA', 'Int', 'IntA', 'LC', 'LPGi', 'LPGiA', 'MSO', 'MVe', 'MVeMC', 'MVePC', 'P7', 'PM', 'Ve', 'Pa5', 'Pr', 'Su5', 'VTg', 'InCo', 'MTg', 'Sp5', 'SuVe'] }
   
    for higher_region, sub_region in higher_hierarchy_nomenclature.items():
        counter[higher_region] = 0.0

    for higher_region, sub_region in higher_hierarchy_nomenclature.items():
    
        for index in cell_count_ROI.index:
            if index in sub_region:
                counter[higher_region] = counter[higher_region] + cell_count_ROI[index]
        counter[higher_region] = (counter[higher_region] / total_cells) *100.00
    return counter    


def amygdala_only(df):
    """
    creates df including only regions of the amygdala 
    """
    
    amygdala = ['LA','BLA','BLP','BLv', 'BLV','BMA','BMP', 'IM','BMP','PMCo', 'PLCo', 
                'ACo','AA', 'AAD', 'AAV', 'I','RAPir', 'RaPir', 'APir','CeC','CeL','CeM',
                'MeA', 'MeAD', 'MeAV', 'MePD', 'MeP','EA']
    df=df[df['Label'].isin(amygdala)]
    df['percent_total_input_amygdala'] = df['Count'].apply(lambda x: 100*(x/df['Count'].sum()))
    return df


def sort_amygdala_pivot(df):
    """
    adds a column to the df according to which lower hierarchy the labels belong to 
    (CeA CeM CeC into "central amygdala",etc.)
    """
    lower_order_hierarchy = {
            'Lateral Amygdala' :['LA'],
            'Basolateral amygdalar nucleus, anterior (BLA)': ['BLA'],
            'Basolateral amygdalar nucleus, posterior (BLP)': ['BLP'],
            'Basolateral amygdalar nucleus, ventral (BLv)': ['BLv', 'BLV'],
            'Basomedial amygdalar nucleus, anterior': ['BMA'],
            'Basomedial amygdalar nucleus, posterior': ['BMP'],
            'Cortical amygdalar area': ['PMCo', 'PLCo', 'ACo'],
            'Anterior amygdaloid area':['AA', 'AAD', 'AAV', 'I'],
            'Amygdalopirifrom transition area':['RAPir', 'RaPir', 'APir'],
            'Central Nucleus of Amygdala, capsular part':['CeC'],
            'Central Nucleus of Amygdala, lateral part':['CeL'],
            'Central Nucleus of Amygdala, medial part':['CeM'],
            'Medial Amygdala': ['MeA', 'MeAD', 'MeAV', 'MePD', 'MeP'],
            'Extended Amygdala': ['EA']}
      
    temp1 = []        
    for i in df.Label:
        for key, value in lower_order_hierarchy.items():
            if i in value :
                temp1.append(key)
    return temp1
        

def sort_amygdala(df,df2):
    """
    adds a column to the df according to which lower hierarchy the labels belong to 
    (CeA CeM CeC into "central amygdala",etc.)
    sorts percent of total input into hierarchy
    """
    
    lower_order_hierarchy = {
            'Lateral Amygdala' :['LA'],
            'Basolateral amygdalar nucleus, anterior (BLA)': ['BLA'],
            'Basolateral amygdalar nucleus, posterior (BLP)': ['BLP'],
            'Basolateral amygdalar nucleus, ventral (BLv)': ['BLv', 'BLV'],
            'Basomedial amygdalar nucleus, anterior': ['BMA'],
            'Basomedial amygdalar nucleus, posterior': ['BMP'],
            'Cortical amygdalar area': ['PMCo', 'PLCo', 'ACo'],
            'Anterior amygdaloid area':['AA', 'AAD', 'AAV', 'I'],
            'Amygdalopirifrom transition area':['RAPir', 'RaPir', 'APir'],
            'Central Nucleus of Amygdala, capsular part':['CeC'],
            'Central Nucleus of Amygdala, lateral part':['CeL'],
            'Central Nucleus of Amygdala, medial part':['CeM'],
            'Medial Amygdala': ['MeA', 'MeAD', 'MeAV', 'MePD', 'MeP'],
            'Extended Amygdala': ['EA']}
    
    counter = pd.Series().astype('float64') 
                   
    #Sort percent of total input strength (% of total output pixel) into lower order hierarchy
    for lower_region, ROI_name in lower_order_hierarchy.items():
        counter[lower_region] = 0

    index_counter = 0
    for lower_region, ROI_name in lower_order_hierarchy.items():
    
        for index in df.index:
            if index in ROI_name:
                counter[lower_region] = counter[lower_region] + df[index]
                index_counter += 1
        if index_counter == 0:
            counter[lower_region] = 0.0
       
        index_counter=0 
    
    counter.rename('total cell count')
    counter_cells = counter.copy()
    counter_cells_percent = counter_cells.apply(lambda x: 100*(x/df2.Count.sum()))
    
    return counter_cells, counter_cells_percent



def th_only(df):
    """
    creates df including only regions of the thalamus 
    """
    thalamus = ['VPL','VPM','VPPC', 'VPLpc', 'VPMpc' ,'VPLPC', 'VPLPc', 'Gus', 'RRe','VM','VL','MGN',
                'MGM','MGv', 'MGV', 'MGm','medial geniculate', 'MG','DM', 'MD', 'MDC','MDL','MDM', 
                'Sub', 'SubD', 'LP','Po', 'PO', 'PoT', 'PIL', 'REth', 'SPFPC', 'SPFp', 'SG', 'Eth', 
                'SPF', 'DT', 'DTg','AD', 'AM', 'AV', 'LD', 'DLG', 'PrG','CM','CL','PF', 'IAM', 'PaF',
                'Rh','PVT', 'PVA', 'PV', 'PVG', 'sm']
    df=df[df['Label'].isin(thalamus)]
    df['percent_total_input_thalamus'] = df['Count'].apply(lambda x: 100*(x/df['Count'].sum()))
    return df


def sort_th_pivot(df):
    """
    adds a column to the df according to which lower hierarchy the  labels belong to 
    (CeA CeM CeC into "central amygdala",etc.)
    """
    
    lower_order_hierarchy = {
            'Ventral Posterior Complex' :['VPL','VPM','VPPC', 'VPLpc', 'VPMpc' ,'VPLPC', 'VPLPc', 'Gus', 'RRe'],
            'Ventral medial nucleus': ['VM'],
            'Ventral anterior-lateral complex': ['VL'],
            'Medial geniculate nuclei': ['MGN','MGM','MGv', 'MGV', 'MGm','medial geniculate', 'MG'],
            'Medio-dorsal nucleus': ['DM', 'MD', 'MDC','MDL','MDM', 'Sub', 'SubD', 'LP'],
            'Posterior complex': ['Po', 'PO', 'PoT', 'PIL', 'REth', 'SPFPC', 'SPFp', 'SG', 'Eth', 'SPF', 
                                  'DT', 'DTg'],
            'anterio-dorsal nuclei':['AD', 'AM', 'AV', 'LD', 'DLG', 'PrG'],
            'centro-medial nucleus':['CM','CL'],
            'Parafascicular nucleus':['PF', 'IAM', 'PaF'],
            'Rhomboid nucleus':['Rh'],
            'Paraventricular Thalamus':['PVT', 'PVA', 'PV', 'PVG', 'sm']
            }
      
    temp1 = []        
    for i in df.Label:
        for key, value in lower_order_hierarchy.items():
            if i in value :
                temp1.append(key)
    return temp1
        

def sort_thalamus(df, df2):
    """
    adds a column to the df according to which lower hierarchy the  labels belong to 
    (CeA CeM CeC into "central amygdala",etc.)
    """
    
    lower_order_hierarchy = {
            'Ventral Posterior Complex' :['VPL','VPM', 'VPPC', 'VPLpc', 'VPMpc' ,'VPLPC', 'VPLPc', 'Gus', 'RRe'],
            'Ventral medial nucleus': ['VM'],
            'Ventral anterior-lateral complex': ['VL'],
            'Medial geniculate nuclei': ['MGN','MGM','MGv', 'MGV', 'MGm','medial geniculate', 'MG'],
            'Medio-dorsal nucleus': ['DM', 'MD', 'MDC','MDL','MDM', 'Sub', 'SubD', 'LP'],
            'Posterior complex': ['Po', 'PO', 'PoT', 'PIL', 'REth', 'SPFPC', 'SPFp', 'SG', 'Eth', 'SPF', 
                                  'DT', 'DTg'],
            'anterio-dorsal nuclei':['AD', 'AM', 'AV', 'LD', 'DLG', 'PrG'],
            'centro-medial nucleus':['CM','CL'],
            'Parafascicular nucleus':['PF', 'IAM', 'PaF'],
            'Rhomboid nucleus':['Rh'],
            'Paraventricular Thalamus':['PVT', 'PVA', 'PV', 'PVG', 'sm']
            }

    counter = pd.Series().astype('float64') 
   
    #Sort percent of total input strength (% of total output pixel) into lower order hierarchy
    for lower_region, ROI_name in lower_order_hierarchy.items():
        counter[lower_region] = 0

    index_counter = 0
    for lower_region, ROI_name in lower_order_hierarchy.items():
    
        for index in df.index:
            if index in ROI_name:
                counter[lower_region] = counter[lower_region] + df[index]
                index_counter += 1
        if index_counter == 0:
            counter[lower_region] = 0.0
       
        index_counter=0 
    
    counter.rename('total cell count')
    counter_cells = counter.copy()
    counter_cells_percent = counter_cells.apply(lambda x: 100*(x/df2.Count.sum()))
    
    return counter_cells, counter_cells_percent
    

def striatum_only(df):
    """
    creates df including only regions of the thalamus 
    """
    
    striatum = ['NAcc', 'NAc', 'Nac', 'Nacc', 'Nacsh','NAcsh', 'NaSh', 'IPAC','Tu','CPu', 
                'CPU', 'Cpu', 'LSS', 'CB', 'LSD']
    df=df[df['Label'].isin(striatum)]
    return df


def sort_striatum_pivot(df):
    """
    adds a column to the df according to which lower hierarchy the  labels belong to
    """
    
    lower_order_hierarchy = {
            'Nucleus Accumbens, core':['NAcc', 'Nacc', 'NAc', 'Nac'],
            'Nucleus Accumbens, shell': ['NAcsh', 'Nacsh', 'NaSh'],
            'Fundus of Striatum (IPAC)':['IPAC'],
            'olfactory tubercle':['Tu', 'CB'],
            'Caudate Putamen':['CPu', 'CPU', 'Cpu', 'LSS', 'LSD']}
      
    temp2 = []        
    for i in df.Label:
        for key, value in lower_order_hierarchy.items():
            if i in value :
                temp2.append(key)
    

    return temp2


def sort_striatum(df,df2):
    """
    adds a column to the df according to which lower hierarchy the labels belong to 
    (CeA CeM CeC into "central amygdala",etc.)
    sorts percent of total input into hierarchy
    """
    
    lower_order_hierarchy = {
            'Nucleus Accumbens, core':['NAcc', 'Nacc', 'NAc', 'Nac'],
            'Nucleus Accumbens, shell': ['NAcsh', 'Nacsh', 'NaSh'],
            'Fundus of Striatum (IPAC)':['IPAC'],
            'olfactory tubercle':['Tu', 'CB'],
            'Caudate Putamen':['CPu', 'CPU', 'Cpu', 'LSS', 'LSD']}
    
    counter = pd.Series().astype('float64') 
                   
    #Sort percent of total input strength (% of total output pixel) into lower order hierarchy
    for lower_region, ROI_name in lower_order_hierarchy.items():
        counter[lower_region] = 0

    index_counter = 0
    for lower_region, ROI_name in lower_order_hierarchy.items():
    
        for index in df.index:
            if index in ROI_name:
                counter[lower_region] = counter[lower_region] + df[index]
                index_counter += 1
        if index_counter == 0:
            counter[lower_region] = 0.0
       
        index_counter=0 
    
    counter.rename('total cell count')
    counter_cells = counter.copy()
    counter_cells_percent = counter_cells.apply(lambda x: 100*(x/df2.Count.sum()))
    
    return counter_cells, counter_cells_percent



###############################################################################################

"""
1) loop through all files in a folder
2) process data as before and create the csv files in a folder with the name of the animal
3) at the end, append the raw summery data to a big dataframe, where all the data is contained
"""

###############################################################################################

"""
make dictionary including sample number and Injection site Start and End

injection_sites = {
        'sample number' :[start,end],    
        
        }
"""

df_big = pd.DataFrame()

root = tk.Tk()
root.withdraw()
directory = filedialog.askdirectory()

for filename in os.listdir(directory):
    
    #parse filename out, also without .csv ending
    file_name_with_csv = os.path.basename(filename)
    file_name = file_name_with_csv[:-4]
    file_name_split = file_name.split('_')

    #parse animal #, and condition out here
    animal_number = file_name_split[1]
    condition = file_name_split[0]
    condition_split = condition.split('.')
    genotype = condition_split[1]
    tracing_from = condition_split[2]
    

    #get Bregma levels of Starter Cell Start and End
    for keys,values in injection_sites.items():
            if keys == animal_number:
                start_Bregma = values[0]
                end_Bregma = values[1]
   
    
    #load in the csv, remove starter volume and exclude regions that cannot provide input to IC
    df = read_and_clean_data(filename)
    df = remove_starter_volume(df, start_Bregma, end_Bregma)
    df = exclude_regions(df)

    #get data grouped by ROIs
    cell_count_ROI, cell_count_ROI_percent = cell_count_per_ROI(df)
    cell_density = cell_density_ROI(df)
    
    #sort into lower and higher order hierarchy
    lower_order_sorted_df, missing_ROIs_lower, missing_ROIs_higher = sort_into_hierarchies(cell_density,cell_count_ROI_percent)
    
    #cutoff of 0.03 % here
    lower_order_sorted_df[lower_order_sorted_df['% of total inputs'] < 0.03] = 0
    
    
    #additional calculations:
    cell_count_Bregma = cell_count_per_Bregma(df)
    
    #create pivot tables to plot along rostro-caudal axis:
    pivot_table_counts = pivot_counts(df)
    pivot_table_percent = pivot_percent(df)
    pivot_table_density = pivot_density(df)
    
    cell_density_bregma = cell_density_Bregma(df)
    
    cell_count_ROI_percent_grouped_in_hierarchy = cell_count_ROI_higher_hierarchy_in_percent(cell_count_ROI)
    
    
    
######################################################################################################
    
    #AMYGDALA ANALYSIS
    
    amygdala_only_cellnumber, amygdala_only_percent = sort_amygdala(cell_count_ROI,df)
    amygdala_only_df = pd.concat([amygdala_only_cellnumber, amygdala_only_percent],  axis=1)
    amygdala_only_df.columns =['cell_count', 'percent']
    amygdala_only_df[amygdala_only_df['percent'] < 0.03] = 0
    
    #remove all ROIs that are not amygdala from the df
    df_amygdala = df.copy()
    df_amygdala = amygdala_only(df_amygdala)
    
    #relabel amygdala lowest lever ROIS, like AA, AAD, AAV into "anterior amygdala"
    df_amygdala = df_amygdala.assign(Label=sort_amygdala_pivot(df_amygdala))
    
    #create pivot tables to plot Anterior-Posterior graphs
    pivot_table_counts_amy = pivot_counts(df_amygdala)
    pivot_table_percent_amygdala = pivot_percent(df_amygdala)
    pivot_table_density_amy = pivot_density(df_amygdala)
    
    
    
########################################################################################################
    
    
    #THALAMUS ANALYSIS
    
    # sort and process thalamus ROIs
    thalamus_only_cellnumber, thalamus_only_percent = sort_thalamus(cell_count_ROI,df)

    thalamus_only_df = pd.concat([thalamus_only_cellnumber, thalamus_only_percent],  axis=1)
    thalamus_only_df.columns =['cell_count', 'percent']
    thalamus_only_df[thalamus_only_df['percent'] < 0.03] = 0
    
    
    #remove all ROIs that are not thalamus from the df
    df_thalamus = df.copy() 
    df_thalamus = th_only(df)

    #relabel thalamus lowest lever ROIS
    df_thalamus = df_thalamus.assign(Label=sort_th_pivot(df_thalamus))

    #create pivot tables to plot Anterior-Posterior graphs
    pivot_table_counts_th = pivot_counts(df)
    pivot_table_percent_th = pivot_percent(df)
    pivot_table_density_th = pivot_density(df)
    

##########################################################################################################    
   
    
    #STRIATUM ANALYSIS
    
    # sort and process striatum ROIs
    striatum_only_cellnumber, striatum_only_percent = sort_striatum(cell_count_ROI,df)

    striatum_only_df = pd.concat([striatum_only_cellnumber, striatum_only_percent],  axis=1)
    striatum_only_df.columns =['cell_count', 'percent']
    striatum_only_df[striatum_only_df['percent'] < 0.03] = 0
    
    
    #remove all ROIs that are not striatum from the df
    df_striatum = df.copy() 
    df_striatum = striatum_only(df)

    #relabel striatum lowest lever ROIS
    df_striatum = df_striatum.assign(Label=sort_striatum_pivot(df_striatum))

    #create pivot tables to plot Anterior-Posterior graphs
    pivot_table_counts_th = pivot_counts(df)
    pivot_table_percent_th = pivot_percent(df)
    pivot_table_density_th = pivot_density(df)
    

###########################################################################################################  
    
#save everything as .csv files
     
    save_path = file_name 
    
    #create target directory 
    if not os.path.exists(save_path):
        os.mkdir(save_path)
        print("Directory " , save_path ,  " Created ")
    else:    
        print("Directory " , save_path ,  " already exists")
    
    #global data
    df.to_csv("{0}\{1}_summary_raw.csv".format(save_path, file_name))    
    cell_count_Bregma.to_csv("{0}\{1}_cell_count_Bregma.csv".format(save_path, file_name))
    cell_count_ROI.to_csv("{0}\{1}_cell_count_every_ROIs.csv".format(save_path, file_name))
    lower_order_sorted_df.to_csv("{0}\{1}_input_and_density_hierarchy.csv".format(save_path, file_name))
    cell_count_ROI_percent_grouped_in_hierarchy.to_csv("{0}\{1}_cell_count_ROIs_percent_grouped_hierarchy.csv".format(save_path, file_name))
    cell_density_bregma.to_csv("{0}\{1}_average_cell_density_Bregma.csv".format(save_path, file_name))
    pivot_table_percent.to_csv("{0}\{1}_pivot_percent_total_input.csv".format(save_path, file_name))
    pivot_table_counts.to_csv("{0}\{1}_pivot_raw_counts.csv".format(save_path, file_name))
    pivot_table_density.to_csv("{0}\{1}_pivot_density.csv".format(save_path, file_name))
    
    #Amygdala
    amygdala_only_df.to_csv("{0}\{1}_amygdala_count_and_percent.csv".format(save_path, file_name))
    pivot_table_percent_amygdala.to_csv("{0}\{1}_pivot_percent_of_amygdala_input.csv".format(save_path, file_name))
    pivot_table_counts_amy.to_csv("{0}\{1}_pivot_raw_counts_amygdala.csv".format(save_path, file_name))
    pivot_table_density_amy.to_csv("{0}\{1}_pivot_density_amygdala.csv".format(save_path, file_name))
    
    #Thalamus
    thalamus_only_df.to_csv("{0}\{1}_thalamus_count_and_percent.csv".format(save_path, file_name))
    pivot_table_percent_th.to_csv("{0}\{1}_pivot_percent_of_thalamus_input.csv".format(save_path, file_name))
    pivot_table_counts_th.to_csv("{0}\{1}_pivot_raw_counts_thalamus.csv".format(save_path, file_name))
    pivot_table_density_th.to_csv("{0}\{1}_pivot_density_thalamus.csv".format(save_path, file_name))
    
    #Striatum
    striatum_only_df.to_csv("{0}\{1}_striatum_count_and_percent.csv".format(save_path, file_name))
    pivot_table_percent_striatum.to_csv("{0}\{1}_pivot_percent_of_striatums_input.csv".format(save_path, file_name))
    pivot_table_counts_striatum.to_csv("{0}\{1}_pivot_raw_counts_striatum.csv".format(save_path, file_name))
    pivot_table_density_striatum.to_csv("{0}\{1}_pivot_density_striatum.csv".format(save_path, file_name))
    
###########################################################################################################################################
    
    #concatenating everything
    
    df_big = df_big.append(df, ignore_index=True, sort=False)
    df_big.to_csv("all rabies tracings combined.csv")