# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 09:14:53 2019

@author: gehrlach
"""

import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import filedialog
import os


def read_and_clean_data(csv):
    """ Read in the CSV file that was created by python script 'summary_AAV'. It adds all countings from all sections together.
    
    This function reads the csv, drops unnecassary columns, removes entries where the count was zero. 
    """
    
    df = pd.DataFrame()
    df = pd.read_csv(csv, encoding='latin1')
    df = df.drop(df.columns[0],1)
    
    #add metadata here
    df['genotype'] = genotype
    df['tracing_from'] = tracing_from
    df['animal'] = animal_number
    
    df['Bregma'] = df['Bregma'].round(1)
              
        
    return df


def remove_starter_volume(df, start_Bregma, end_Bregma):
    """ 
    cleaning of dataframe. Makes values positive, removes starter volume, adds injection center and percent total output to df 
    """
    
    #bring everything into the positive numbers
    df['Bregma'] += 2
    start_Bregma += 2
    end_Bregma +=2
    injection_center = start_Bregma - (start_Bregma - end_Bregma)/2 -2
    
    #remove starter cell volume from insula and claustrum at from start to end bregma here
    #add center of injection site to dataframe
    starter_volume = df[(df['Label'].isin(['AIP','AID', 'GI', 'DI', 'AIV', 'AI', 'claustrum'])) & (df['Bregma'] <= start_Bregma) & (df['Bregma'] >= end_Bregma)]
    cleaned_df = pd.concat([df,starter_volume,starter_volume]).drop_duplicates(keep=False)
    cleaned_df['Bregma'] -= 2
    cleaned_df['injection center'] = injection_center
    
    #add percent_total_output to dataframe
    cleaned_df['percent_total_output'] = cleaned_df['Count'].apply(lambda x: 100*(x/cleaned_df['Count'].sum()))
    cleaned_df = percent_density(cleaned_df)
    
    return cleaned_df


def exclude_regions(df):
    """ 
    excludes regions with which there is no real connectivity, but due to off target signal shows artifical signals. manual curation of striatum, cea and cerebellum. 
    """
   
    regions = ['0001-3014-3317', 'Classified image', 'Lo', '0001-2452-1470', '0001-2927-2944', '0001-2761-2481', 'V2ML', 'fmi' , 'ic', 'ml', 'mp', 'ns', 'rmx', 'rs', 'scp', 'corpus callosum', 'cc', 'LL', 'Op', 'PVP', 'EP', 'LDVL', 'LPM', 'acp', 'Com', 'II','LDDM', 'LBh', 'MBh', 'IPL', 'MM', 'Or','SNL','py']
    set_to_zero = ['septum']
    
    excluded_idx = df[df['Label'].isin(regions)].index
    df.drop(excluded_idx, inplace=True)
    
    zero_idx = df[df['Label'].isin(set_to_zero)].index
    df['Count'].loc[zero_idx]=0         
    return df


def assign_level_1(df):
    """ 
    creation of dictionaries for clustering of ROI names into lower and higher hierarchies. 
    """
    
    lower_order_hierarchy = {
    
            #Sensory Cortex        
            'Somatosensory Cortex, primary': ['S1', 'S1J'],
            'Somatosensory Cortex, secondary': ['S2'],
            'Auditory Cortex, primary': ['Au1', 'Auditory cortex', 'auditory cortex'],
            'Auditory Cortex, secondary or higher': ['AuD', 'AuV'],
            'Visual Cortex, primary': ['V1', 'V1B', 'V1M', 'VM-1', 'Visual cortex', 'Visual cortex-1', 'Visual Cortex', 'visual cortex'],
            'Visual Cortex, secondary or higher': ['V2L', 'V2ML','V2MM'],
            'Piriform Area':['Pir'],
            
            #Motor Cortex
            'Motor Cortex, primary': ['M1', 'Fr3'],
            'Motor Cortex, secondary': ['M2'],
            
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
            'Retrosplenial Cortex': ['A30', 'a30', 'A29c', 'A29b', 'A29a', 'A29', 'RSA', 'RSG', 'RSGa', 'RSGb' ],
            'Entorhinal Cortex': ['DLEnt', 'Dlent', 'DIEnt', 'DlEnt', 'VIEnt', 'Vlent', 'VlEnt', 'MEnt', 'CEnt', 'LEnt'],
            
            #Insular Ctx
            'Agranular Insular Cortex': ['AI'],
            'Agranular Insular Cortex, dorsal': ['AID'],
            'Agranular Insular Cortex, ventral': ['AIV'],
            'Agranular Insular Cortex, posterior': ['AIP'],
            'Dysgranular Insular Cortex': ['DI'],
            'Granular Insular Cortex': ['GI'],
            
            #Olfactory Areas
            'Anterior Olfactory Nucleus': ['AON', 'AOP', 'AOM', 'AOL', 'AOD'],
            'Piriform-Amygdalar Area':['CxA'],
            'Endopiriform Nucleus': ['DEn', 'VEn', 'endopiriform nucleus', 'endopirform nucleus', 'endopiriform nuncleus', 'EndoPir', 'Endo Pir', 'EnoPir'],
            'Tenia Tecta': ['DTT', 'VTT', 'tenia tecta', 'taenia tecta', 'teania tectum', 'taenia tecta', 'taenia tectum'],
            'Nucleus of the Lateral Olfactory Tract':['LOT'],
            'Main Olfactory Bulb':['olfactory bulb', 'olfacotry bulb'],
                        
            #Claustrum
            'Claustrum': ['claustrum'],
            
            #Striatum
            'Caudate Putamen':['CPu', 'CPU', 'Cpu', 'LSS', 'LSD'],
            'Nucleus Accumbens, core':['NAcc', 'Nacc', 'NAc', 'Nac', 'Nacc-1'],
            'Nucleus Accumbens, shell': ['NAcsh', 'Nacsh', 'NaSh'],
            'Fundus of Striatum (IPAC)':['IPAC'],
            'Olfactory Tubercle':['Tu', 'CB', 'TU'],
                        
            #Pallidum
            'Ventral Pallidum':['VP'],
            'Bed Nuclei of the Stria Terminalis (BNST)': ['BNST', 'BST'],
            'Globus Pallidus': ['GP', 'LGP', 'MGP'],
            'Medial Septum': ['MS'],
            'Lateral Septum': ['septum'],
            
            #Hippocampus
            'Hippocampus (CA1,CA2,CA3)': ['CA1', 'CA2', 'CA3', 'fi', 'LMol', 'Or', 'Py', 'Rad', 'STr', 'Str', 'SLu', 'alv', 'vHip', 'vHIp', 'AHiP', 'PaS', 'DS', 'POST', 'PROST', 'Post', 'PrS', 'SubC', 'VS'],
            'Hippocampus (Dentate Gyrus)': ['DG', 'Mol'],
            
            #Amygdala
            'Lateral Amygdala' :['LA'],
            'Basolateral Nucleus': ['BLA','BLP','BL', 'BLv', 'BLV', 'BLS'],
            'Basomedial Nucleus': ['BMA', 'BMP', 'IM'],
            'Cortical Amygdala': ['PMCo', 'PLCo', 'ACo'],
            'Anterior Amygdaloid Area':['AA', 'AAD', 'AAV', 'I'],
            'Amygdalopiriform Transition Area': ['RAPir', 'RaPir', 'APir', 'AHiAL', 'AHiPL'],
            'Central Nucleus of Amygdala':['CeC', 'CeL', 'CeA', 'CeM', 'CeMAD', 'CeMAV', 'Ce', 'STIA'],            
            'Medial Amygdala': ['MeA', 'MeAD', 'MeAV', 'MePD', 'MeP'],
            'Extended Amygdala': ['EA'],
                        
            #Basal Forebrain
            'Basal Forebrain':['MDB', 'LDB', 'diagonal band', 'VDB', 'HDB', 'SIB', 'SI/B', 'SI', 'B'],
            
            #Hypothalamus
            'Hypothalamic Lateral Zone': ['Arc', 'LH', 'LHA', 'STh', 'PSTh', 'MCLH', 'ZI', 'ZIR', 'F', 'hypothalamus', 'dorsal hypothalamus', 'LPO', 'RML', 'SO', 'SubI', 'Subl','Sub1', 'RMM'],
            'Hypothalamic Medial Zone' :['MnPO', 'MPO', 'MPA', 'MPOM', 'AHC', 'MHb', 'PH', 'AHA', 'VMH', 'PHD', 'RI', 'PHA', 'AH'],
                        
            #Thalamus
            'Habenular Nuclei': ['LHb', 'LHB'],
            'Ventral Group of the Dorsal Thalamus':['VPL','VPM','VL','VM', 'VPPC', 'VPLpc', 'VPMpc' ,'VPLPC', 'VPLPc', 'Gus', 'RRe'],
            'Medial Geniculate Nuclei': ['MGN','MGM','MGv', 'MGV', 'MGm','medial geniculate', 'MG'],            
            'Posterior Complex of the Thalamus':['Po', 'PO', 'PoT', 'PIL', 'REth', 'SPFPC', 'SPFp', 'SG', 'Eth', 'SPF', 'DT', 'DTg'],
            'Mediodorsal Thalamic Nuclei': ['DM', 'MD', 'MDC','MDL','MDM', 'Sub', 'SubD', 'LP'],
            'Anterodorsal Thalamic Nucleus': ['AD', 'AM', 'AV', 'LD', 'DLG', 'PrG'],
            'Intralaminar Nuclei of the Dorsal Thalamus': ['Re', 'Rh', 'CM', 'PF', 'IAM', 'PaF', 'Rt', 'PC', 'SubV', 'CL', 'VA'],
            'Paraventricular Thalamic Nucleus': ['PVT', 'PVA', 'PV', 'PVG', 'sm'],
                        
            #Midbrain
            'Periaqueductal Gray':['PAG', 'EW', 'PrEW', 'VLPAG', 'DLPAG', 'InC', 'In'],     
            'Substantia Nigra': ['SNc', 'SNC', 'SNr', 'SNR','SN', 'substantia nigra', 'Substantia nigra', 'cp', 'PR', 'IP', 'IPD', 'RMC', 'CI', 'PN'],
            'Ventral Tegmental Area': ['VTA', 'PBP', 'Pa4', 'PIF', 'VTAR'],
            'Raphe Nuclei': ['dorsal raphe', 'DR', 'median raphe', 'medial raphe', 'RLi', 'Rli', 'MnR', 'PMnR', 'PMnT', 'RMg', 'CLi', 'RPa', 'IF', 'RM', 'IPR'],
            'Midbrain Reticular Nucleus':['isRt', 'mRt', 'p1Rt', 'RPC', 'RP', 'DA8', 'DAB', 'APT', 'APTD', 'BIC', 'IRt', 'RRF', 'PaR', 'SubB', 'DpMe'],
            'Superior Colliculus':['InG', 'DpWh', 'SuG', 'DpG'],
            'Inferior Colliculus':['CIC', 'ECIC', 'DCIC'],
            
            #Hindbrain
            'Parabrachial Nucleus': ['PBN', 'LPB', 'MPB', 'parabrachial nucleus'],
            'Pedunculopontine Nucleus': ['PTg', 'LDTgV', 'SPTg', 'MPL', 'PL', 'KF', 'PrCnF', 'PnO', 'PnC', 'InCo', 'Pn', 'PnV', 'VeCb', 'CnF', 'CAT', 'DMTg', 'MiTg', 'MTg', 'PBG', 'LSO', 'PLV', 'PPTg', 'VTg', 'PNV', 'Pr5', 'RtTg', 'SPO', 'Tz', 'VLL', 'B9', 'Bar', 'CG', 'Gi', 'GiA', 'LC', 'MSO', 'Su5'],
            'Medulla': ['NTS', 'PCRtA', '7N', 'DPGi', 'IRtA', 'LPGi', 'LPGiA', 'LVe', 'MVe', 'MVeMC', 'MVePC', 'P7', 'Pa5', 'Pr', 'Ve', 'Sp5', 'SuVe'],
            'Cerebellum':['2Cb', 'Int', 'IntA', 'PM'],
            'Motor Trigeminal Nucleus':['5', 'SubCV']
            }
    
    level1_label=[]
    for i in df.Label:
        for key, value in lower_order_hierarchy.items():
            if i in value:
                level1_label.append(key)
    
    print(len(level1_label))
    print(len(df))
    df['Level_1'] = level1_label
    return df


def assign_level_2(df):
    higher_hierarchy_nomenclature = {'Insular Cortex': ['AI', 'AID', 'AIV', 'GI', 'DI', 'AIP'],
                                     'Sensory Cortex': ['S1','S1J', 'S2', 'Au1', 'AuD', 'Auditory cortex', 'auditory cortex', 'AuV', 'V1', 'V2L', 'Visual cortex', 'Visual cortex-1', 'Visual Cortex', 'visual cortex', 'V1B', 'V1M', 'VM-1', 'V2MM'],
                                     'Motor Cortex': ['M1', 'M2'],
                                     'Prefrontal Cortex':['FrA', 'MO', 'VO' , 'DLO', 'LO', 'A32', 'A24a', 'cg', 'Cg1', 'Cg2', 'A25', 'A24b', 'Fr3', 'DP', 'IL', 'PrL'],
                                     'Association Cortex': ['parietal cortex', 'Ect','PRh', 'Prh', 'DLEnt', 'DIEnt', 'VIEnt', 'VlEnt', 'MEnt', 'TeA','A29c', 'A29', 'A30', 'a30', 'A29a', 'A29b', 'CEnt', 'DlEnt', 'VLEnt', 'VIEnt', 'Vlent', 'LEnt'],
                                     'Piriform Cortex':['Pir'],
                                     'Olfactory areas': ['CxA','olfactory bulb', 'olfacotry bulb', 'AON', 'tenia tecta', 'taenia tecta', 'teania tectum', 'taenia tectum', 'endopiriform nucleus', 'endopirform nucleus', 'endopiriform nuncleus', 'EndoPir','Endo Pir', 'EnoPir', 'VTT', 'LOT', 'DTT'],
                                     'Claustrum': ['claustrum'],
                                     'Striatum':['NAcc', 'NAc', 'Nac', 'Nacc', 'Nacsh','NAcsh', 'NaSh', 'IPAC','Tu','TU','CPu', 'CPU', 'Cpu', 'LSS', 'CB', 'LSD', 'Nacc-1'],
                                     'Pallidum':['VP','GP','BNST', 'BST', 'MS', 'septum'],
                                     'Basal Forebrain':['SIB', 'diagonal band', 'HDB', 'LDB', 'B'],
                                     'Hippocampus': ['vHip', 'vHIp', 'AHiP', 'DG', 'PaS', 'DS', 'POST', 'Post', 'PROST', 'PrS', 'Py', 'SubC', 'VS', 'STr', 'Str'],
                                     'Amygdala': ['CeA', 'CeM', 'CeC', 'CeL', 'PMCo', 'PLCo','AA', 'BLP', 'LA', 'BMA', 'BMP', 'EA','BLA', 'BLS', 'BLv', 'BLV', 'CeA', 'MeA', 'MeAD', 'ACo','AHiAL', 'AHiPL', 'I', 'IM', 'MeP', 'STIA', 'APir','RAPir', 'RaPir'],    
                                     'Hypothalamus':['Arc', 'AHA', 'AH','LH', 'MHb', 'STh', 'PSTh', 'VMH', 'MPOM', 'MPA', 'hypothalamus', 'F', 'dorsal hypothalamus','AHC', 'LPO' ,'MnPO', 'MCLH', 'ZI', 'ZIR', 'LHb', 'LHB', 'RML', 'PH', 'SubI', 'Subl', 'Sub1', 'PHD', 'RI', 'RMM', 'PHA'],
                                     'Thalamus':['DM','AD','MD','MGN', 'MGm', 'MG', 'MGV', 'MGv', 'medial geniculate', 'VPL', 'VPM', 'VL', 'VM', 'VPPC', 'VPLpc', 'VPLPC', 'VPMpc', 'VPLPc', 'VPMpv', 'CM', 'AM', 'AV', 'Po', 'PO', 'PV', 'PVA', 'PaF', 'PIL', 'PoT', 'LD', 'IAM', 'Re', 'Rh', 'SPFPC', 'SPFp', 'SPF', 'medial geniculate', 'CM', 'REth', 'Rt', 'Sub', 'SubD', 'PC', 'SG', 'sm', 'Gus', 'DLG', 'Eth', 'PrG', 'SO', 'SubV', 'CL', 'DT', 'DTg', 'LP', 'LVe', 'RRe', 'VA'],
                                     'Midbrain':['substsantia nigra', 'Substantia nigra', 'SN', 'cp', 'SNc', 'SNC', 'SNr', 'SNR', 'dorsal raphe', 'median raphe', 'medial raphe', 'isRt', 'CIC', 'PAG', 'DLPAG', 'VTA', 'VTAR', 'RLi', 'Rli', 'mRt','PVG', 'EW', 'LDTgV', 'PrEW', 'VLPAG', 'p1Rt', 'PBP', 'DA8', 'DAB', 'RPC', 'RP', 'PR', 'APT', 'APTD', 'BIC', 'ECIC','DpG', 'DR', 'IP', 'IPD', 'IRt', 'DCIC', 'MnR', 'RMC', 'RMg', 'RRF', 'CI', 'CLi', 'PN', 'Pa4', 'PaR', 'RPa', 'SubB', 'InC', 'PMnR', 'IF', 'InG', 'PIF', 'RM', 'IPR', 'DpWh', 'PMnT', 'SuG', 'In', 'DpMe'],
                                     'Hindbrain':['PBN', 'parabrachial nucleus', 'NTS', 'MPB', 'LPB', 'PCRtA', 'PTg', 'PrCnF', 'SPTg', 'PnO', 'PnV', 'MPL', 'PL', 'KF', '2Cb', 'SubCV', '5', 'VeCb', 'CAT', 'CnF', 'DMTg', 'MiTg', 'PBG', 'LSO', 'PLV', 'PPTg', 'Pn', 'PnC', 'PNV', 'Pr5', 'RtTg', 'SPO', 'Tz', 'VLL', '7N', 'B9', 'Bar', 'CG', 'DPGi', 'Gi', 'GiA', 'IRtA', 'Int', 'IntA', 'LC', 'LPGi', 'LPGiA', 'MSO', 'MVe', 'MVeMC', 'MVePC', 'P7', 'PM', 'Ve', 'Pa5', 'Pr', 'Su5', 'VTg', 'InCo', 'MTg', 'Sp5', 'SuVe'] }
    
    level2_label=[]
    for i in df.Label:
        for key, value in higher_hierarchy_nomenclature.items():
            if i in value:
                level2_label.append(key)
    
    df['Level_2'] = level2_label
    return df


def sort_striatum(df_striatum, df):
    """
    adds a column to the df according to which lower hierarchy the labels belong to (CeA CeM CeC into "central amygdala",etc.)
    """
    
    lower_order_hierarchy = {
            'Caudate Putamen':['CPu', 'CPU', 'Cpu', 'LSS', 'LSD'],
            'Nucleus Accumbens, core':['NAcc', 'Nacc', 'NAc', 'Nac', 'Nacc-1'],
            'Nucleus Accumbens, shell': ['NAcsh', 'Nacsh', 'NaSh'],
            'Fundus of Striatum (IPAC)':['IPAC'],
            'Olfactory Tubercle':['Tu', 'CB','TU'],
            }

    counter = pd.Series().astype('float64') 
              
    #Sort percent of total input strength (% of total output pixel) into lower order hierarchy
    for lower_region, ROI_name in lower_order_hierarchy.items():
        counter[lower_region] = 0

    index_counter = 0
    for lower_region, ROI_name in lower_order_hierarchy.items():
    
        for index in df_striatum.index:
            if index in ROI_name:
                counter[lower_region] = counter[lower_region] + df_striatum[index]
                index_counter += 1
        if index_counter == 0:
            counter[lower_region] = 0.0
       
        index_counter=0 
    
    counter.rename('total cell count')
    counter_cells = counter.copy()
    counter_cells_percent = counter_cells.apply(lambda x: 100*(x/df.Count.sum()))
        
    return counter_cells, counter_cells_percent


def striatum_only(df):
    """
    creates df including only regions of the striatum 
    """
    
    striatum = ['NAcc', 'NAc', 'Nac', 'Nacc', 'Nacsh','NAcsh', 'NaSh', 'IPAC','Tu','TU','CPu', 'CPU', 'Cpu', 'LSS', 'CB', 'LSD', 'Nacc-1']
    df=df[df['Label'].isin(striatum)]
    return df


def sort_striatum_pivot(df):
    """
    adds a column to the df according to which lower hierarchy the labels belong to
    """
    
    lower_order_hierarchy = {
            'Caudate Putamen':['CPu', 'CPU', 'Cpu', 'LSS', 'LSD'],
            'Nucleus Accumbens, core':['NAcc', 'Nacc', 'NAc', 'Nac', 'Nacc-1'],
            'Nucleus Accumbens, shell': ['NAcsh', 'Nacsh', 'NaSh'],
            'Fundus of Striatum (IPAC)':['IPAC'],
            'Olfactory Tubercle':['Tu', 'CB','TU'],
            }
      
    temp2 = []        
    for i in df.Label:
        for key, value in lower_order_hierarchy.items():
            if i in value :
                temp2.append(key)
    
    return temp2


def th_only(df):
    """
    creates df including only regions of the thalamus 
    """
    
    thalamus = ['VPL','VPM','VPPC', 'VPLpc', 'VPMpc' ,'VPLPC', 'VPLPc', 'Gus', 'RRe','VM','VL','MGN','MGM','MGv', 'MGV', 'MGm','medial geniculate', 'MG','DM', 'MD', 'MDC','MDL','MDM', 'Sub', 'SubD', 'LP','Po', 'PO', 'PoT', 'PIL', 'REth', 'SPFPC', 'SPFp', 'SG', 'Eth', 'SPF', 'DT', 'DTg','AD', 'AM', 'AV', 'LD', 'DLG', 'PrG','CM','CL','PF', 'IAM', 'PaF','Rh','PVT', 'PVA', 'PV', 'PVG', 'sm']
    df=df[df['Label'].isin(amygdala)]
    df['percent_total_output_thalamus'] = df['Count'].apply(lambda x: 100*(x/df['Count'].sum()))
    return df

def sort_thalamus(df_thalamus, df):
    """
    adds a column to the df according to which lower hierarchy the labels belong to (CeA CeM CeC into "central amygdala",etc.)
    """
    
    lower_order_hierarchy = {
            'Ventral Posterior Complex' :['VPL','VPM', 'VPPC', 'VPLpc', 'VPMpc' ,'VPLPC', 'VPLPc', 'Gus', 'RRe'],
            'Ventral medial nucleus': ['VM'],
            'Ventral anterior-lateral complex': ['VL'],
            'Medial geniculate nuclei': ['MGN','MGM','MGv', 'MGV', 'MGm','medial geniculate', 'MG'],
            'Medio-dorsal nucleus': ['DM', 'MD', 'MDC','MDL','MDM', 'Sub', 'SubD', 'LP'],
            'Posterior complex': ['Po', 'PO', 'PoT', 'PIL', 'REth', 'SPFPC', 'SPFp', 'SG', 'Eth', 'SPF', 'DT', 'DTg'],
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
    
        for index in df_thalamus.index:
            if index in ROI_name:
                counter[lower_region] = counter[lower_region] + df_thalamus[index]
                index_counter += 1
        if index_counter == 0:
            counter[lower_region] = 0.0
       
        index_counter=0 
    
    counter.rename('total cell count')
    counter_cells = counter.copy()
    counter_cells_percent = counter_cells.apply(lambda x: 100*(x/df.Count.sum()))
   
    return counter_cells, counter_cells_percent


def sort_th_pivot(df):
    """
    adds a column to the df according to which lower hierarchy the  labels belong to (CeA CeM CeC into "central amygdala",etc.)
    """
    
    lower_order_hierarchy = {
            'Ventral Posterior Complex' :['VPL','VPM','VPPC', 'VPLpc', 'VPMpc' ,'VPLPC', 'VPLPc', 'Gus', 'RRe'],
            'Ventral medial nucleus': ['VM'],
            'Ventral anterior-lateral complex': ['VL'],
            'Medial geniculate nuclei': ['MGN','MGM','MGv', 'MGV', 'MGm','medial geniculate', 'MG'],
            'Medio-dorsal nucleus': ['DM', 'MD', 'MDC','MDL','MDM', 'Sub', 'SubD', 'LP'],
            'Posterior complex': ['Po', 'PO', 'PoT', 'PIL', 'REth', 'SPFPC', 'SPFp', 'SG', 'Eth', 'SPF', 'DT', 'DTg'],
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

def amygdala_only(df):
    """
    creates df including only regions of the amygdala 
    """
    
    amygdala = ['LA','BLA','BLP','BLv', 'BLV','BMA','BMP', 'BMP','PMCo', 'PLCo', 'ACo','AA', 'AAD', 'AAV', 'RAPir', 'RaPir', 'APir','CeC','CeL','CeM','MeA', 'MeAD', 'MeAV', 'MePD', 'MeP','EA']
    df=df[df['Label'].isin(amygdala)]
    df['percent_total_output_amygdala'] = df['Count'].apply(lambda x: 100*(x/df['Count'].sum()))
    return df

def sort_amygdala(df_amygdala, df):
    """
    adds a column to the df according to which lower hierarchy the  labels belong to (CeA CeM CeC into "central amygdala",etc.)
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
    
        for index in df_amygdala.index:
            if index in ROI_name:
                counter[lower_region] = counter[lower_region] + df_amygdala[index]
                index_counter += 1
        if index_counter == 0:
            counter[lower_region] = 0.0
       
        index_counter=0 
    
    counter.rename('total cell count')
    counter_cells = counter.copy()
    counter_cells_percent = counter_cells.apply(lambda x: 100*(x/df.Count.sum()))

    return counter_cells, counter_cells_percent


def sort_amygdala_pivot(df):
    """
    adds a column to the df according to which lower hierarchy the  labels belong to (CeA CeM CeC into "central amygdala",etc.)
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


def percent_density(df):
    """
    Calculates how many percent of possible pixels are axons. Pixel Density normalized to the maximum possible value. 
    for 1024x1024 confocal : 0.509 ; for 512x512: 0.109; for 839 its 0.857
    """
    
    # normalizes the percent density to the strongest projection (usually thalamus, AID, etc.)  
    df['percent_density'] = df['Pixel Density [Pixel/um²]'].apply(lambda x: 100*(x/df['Pixel Density [Pixel/um²]'].max()))
    
    return df


def average_density_per_ROI(df):
    """
    groups dataframe by ROI name and calculates the meaan average pixel density
    """
    
    grouped_ROI_density = df.groupby(['Label']) 
    average_pixel_density = grouped_ROI_density['percent_density'].mean()
    
    return average_pixel_density


def pixel_count_per_ROI(df):
    """
    total pixel count across all bregmas per ROI
    """
    
    df = df.groupby(['Label'])
    pixel_count = df['Count'].sum()
    pixel_count_percent = (pixel_count / pixel_count.sum() *100.00)
    return pixel_count, pixel_count_percent


def pixel_count_per_Bregma(df):
    """
    total pixel count across all ROIs but this time per Bregma
    """
    
    df = df.groupby(['Bregma'])
    pixel_count = df['Count'].sum()
    return pixel_count


def pivot_density(df):
    """
    First, this function deletes every column from the dataframe, except Bregma and axon densities, 
    then it creates a pivot table
    summing the counts and filling missing values with 0
    
    shows number of pixels per ROI sorted by bregma 
    """
    
    df = df.drop('Count',1)
    df = df.drop('Total Area', 1)
    df = df.drop('Pixel Density [Pixel/um²]',1)
    df = df.drop('Section',1)
    df = df.drop('genotype',1)
    df = df.drop('tracing_from',1)
    df = df.drop('animal',1)
    df = df.drop('injection center',1)
    df = df.drop('percent_total_output',1) 
    df = df.drop('Level_1',1)
    df = df.drop('Level_2',1)
    df = df.pivot_table(df, index=["Bregma"], columns=["Label"], aggfunc=[np.sum], fill_value=0 )
    
    
    return df


def pivot_percent(df):
    """
    First, this function delets every column from the dataframe, except Bregma and Cell count, then it creates a pivot table
    summing the counts and filling missing values with 0
    
    shows number of pixels per ROI sorted by bregma 
    """
    
    df = df.drop('Count',1)
    df = df.drop('Total Area', 1)
    df = df.drop('Pixel Density [Pixel/um²]',1)
    df = df.drop('Section',1)
    df = df.drop('animal',1)
    df = df.drop('genotype',1)
    df = df.drop('tracing_from',1)
    df = df.drop('injection center',1)
    df = df.drop('percent_total_output',1)
    df = df.drop('percent_density',1)
    df = df.drop('Level_1',1)
    df = df.drop('Level_2',1)
    df = df.pivot_table(df, index=["Bregma"], columns=["Label"], aggfunc=[np.sum], fill_value=0 )
    
    return df


def pixel_density_ROI(df):
    """
    groups dataframe by ROI name and calculates the meaan average cell density
    """
    
    grouped_ROI_density = df.groupby(['Label']) 
    average_pixel_density = grouped_ROI_density['Pixel Density [Pixel/um²]'].mean()
    return average_pixel_density


def pixel_density_Bregma(df):
    """
    groups dataframe by Bregma level name and calculates the meaan average cell density
    """
    
    grouped_Bregma_density = df.groupby(['Bregma']) 
    average_cell_density_bregma = grouped_Bregma_density['Pixel Density [Pixel/um²]'].mean()
    return average_cell_density_bregma


def cell_count_ROI_higher_hierarchy_in_percent(cell_count_ROI):
    """
    clusters ROIs into highest abstraction (e.g. amygdala, sensory cortex,etc.)
    """
    
    counter = pd.Series() 
    total_cells = cell_count_ROI.sum()      
    higher_hierarchy_nomenclature = {'Insular Cortex': ['AI', 'AID', 'AIV', 'GI', 'DI', 'AIP'],
                                     'Sensory Cortex': ['S1','S1J', 'S2', 'Au1', 'AuD', 'Auditory cortex', 'auditory cortex', 'AuV', 'V1', 'V2L', 'Visual cortex', 'Visual cortex-1', 'Visual Cortex', 'visual cortex', 'V1B', 'V1M', 'VM-1', 'V2MM'],
                                     'Motor Cortex': ['M1', 'M2'],
                                     'Prefrontal Cortex':['FrA', 'MO', 'VO' , 'DLO', 'LO', 'A32', 'A24a', 'cg', 'Cg1', 'Cg2', 'A25', 'A24b', 'Fr3', 'DP', 'IL', 'PrL'],
                                     'Association Cortex': ['parietal cortex', 'Ect','PRh', 'Prh', 'DLEnt', 'DIEnt', 'VIEnt', 'VlEnt', 'MEnt', 'TeA','A29c', 'A29', 'A30', 'a30', 'A29a', 'A29b', 'CEnt', 'DlEnt', 'VLEnt', 'VIEnt', 'Vlent', 'LEnt'],
                                     'Piriform Cortex':['Pir'],
                                     'Olfactory areas': ['CxA','olfactory bulb', 'olfacotry bulb', 'AON', 'tenia tecta', 'taenia tecta', 'teania tectum', 'taenia tectum', 'endopiriform nucleus', 'endopirform nucleus', 'endopiriform nuncleus', 'EndoPir','Endo Pir', 'EnoPir', 'VTT', 'LOT', 'DTT'],
                                     'Claustrum': ['claustrum'],
                                     'Striatum':['NAcc', 'NAc', 'Nac', 'Nacc', 'Nacsh','NAcsh', 'NaSh','TU', 'IPAC','Tu','CPu', 'CPU', 'Cpu', 'LSS', 'CB', 'LSD', 'Nacc-1'],
                                     'Pallidum':['VP','GP','BNST', 'BST', 'MS', 'septum'],
                                     'Basal Forebrain':['SIB', 'diagonal band', 'HDB', 'LDB', 'B'],
                                     'Hippocampus': ['vHip', 'vHIp', 'AHiP', 'DG', 'PaS', 'DS', 'POST', 'Post', 'PROST', 'PrS', 'Py', 'SubC', 'VS', 'STr', 'Str'],
                                     'Amygdala': ['CeA', 'CeM', 'CeC', 'CeL', 'PMCo', 'PLCo','AA', 'BLP', 'LA', 'BMA', 'BMP', 'EA','BLA', 'BLS', 'BLv', 'BLV', 'CeA', 'MeA', 'MeAD', 'ACo','AHiAL', 'AHiPL', 'I', 'IM', 'MeP', 'STIA', 'APir','RAPir', 'RaPir'],    
                                     'Hypothalamus':['Arc', 'AHA', 'AH','LH', 'MHb', 'STh', 'PSTh', 'VMH', 'MPOM', 'MPA', 'hypothalamus', 'F', 'dorsal hypothalamus','AHC', 'LPO' ,'MnPO', 'MCLH', 'ZI', 'ZIR', 'LHb', 'LHB', 'RML', 'PH', 'SubI', 'Subl', 'Sub1', 'PHD', 'RI', 'RMM', 'PHA'],
                                     'Thalamus':['DM','AD','MD','MGN', 'MGm', 'MG', 'MGV', 'MGv', 'medial geniculate', 'VPL', 'VPM', 'VL', 'VM', 'VPPC', 'VPLpc', 'VPLPC', 'VPMpc', 'VPLPc', 'VPMpv', 'CM', 'AM', 'AV', 'Po', 'PO', 'PV', 'PVA', 'PaF', 'PIL', 'PoT', 'LD', 'IAM', 'Re', 'Rh', 'SPFPC', 'SPFp', 'SPF', 'medial geniculate', 'CM', 'REth', 'Rt', 'Sub', 'SubD', 'PC', 'SG', 'sm', 'Gus', 'DLG', 'Eth', 'PrG', 'SO', 'SubV', 'CL', 'DT', 'DTg', 'LP', 'LVe', 'RRe', 'VA'],
                                     'Midbrain':['substsantia nigra', 'Substantia nigra', 'SN', 'cp', 'SNc', 'SNC', 'SNr', 'SNR', 'dorsal raphe', 'median raphe', 'medial raphe', 'isRt', 'CIC', 'PAG', 'DLPAG', 'VTA', 'VTAR', 'RLi', 'Rli', 'mRt','PVG', 'EW', 'LDTgV', 'PrEW', 'VLPAG', 'p1Rt', 'PBP', 'DA8', 'DAB', 'RPC', 'RP', 'PR', 'APT', 'APTD', 'BIC', 'ECIC','DpG', 'DR', 'IP', 'IPD', 'IRt', 'DCIC', 'MnR', 'RMC', 'RMg', 'RRF', 'CI', 'CLi', 'PN', 'Pa4', 'PaR', 'RPa', 'SubB', 'InC', 'PMnR', 'IF', 'InG', 'PIF', 'RM', 'IPR', 'DpWh', 'PMnT', 'SuG', 'In', 'DpMe'],
                                     'Hindbrain':['PBN', 'parabrachial nucleus', 'NTS', 'MPB', 'LPB', 'PCRtA', 'PTg', 'PrCnF', 'SPTg', 'PnO', 'PnV', 'MPL', 'PL', 'KF', '2Cb', 'SubCV', '5', 'VeCb', 'CAT', 'CnF', 'DMTg', 'MiTg', 'PBG', 'LSO', 'PLV', 'PPTg', 'Pn', 'PnC', 'PNV', 'Pr5', 'RtTg', 'SPO', 'Tz', 'VLL', '7N', 'B9', 'Bar', 'CG', 'DPGi', 'Gi', 'GiA', 'IRtA', 'Int', 'IntA', 'LC', 'LPGi', 'LPGiA', 'MSO', 'MVe', 'MVeMC', 'MVePC', 'P7', 'PM', 'Ve', 'Pa5', 'Pr', 'Su5', 'VTg', 'InCo', 'MTg', 'Sp5', 'SuVe'] }
    
    for higher_region, sub_region in higher_hierarchy_nomenclature.items():
        counter[higher_region] = 0

    for higher_region, sub_region in higher_hierarchy_nomenclature.items():
    
        for index in cell_count_ROI.index:
            if index in sub_region:
                counter[higher_region] = counter[higher_region] + cell_count_ROI[index]
        counter[higher_region] = (counter[higher_region] / total_cells) *100.00
    return counter    


def sort_into_higher_hierarchy_in_percent(cell_count_ROI):
    """
    clusters ROIs into highest abstraction (e.g. amygdala, sensory cortex,etc.)
    """
    
    counter = pd.Series().astype('float64')
    total_cells = cell_count_ROI.sum()      
    higher_hierarchy_nomenclature = {'Insular Cortex': ['AI', 'AID', 'AIV', 'GI', 'DI', 'AIP'],
                                     'Sensory Cortex': ['S1','S1J', 'S2', 'Au1', 'AuD', 'Auditory cortex', 'auditory cortex', 'AuV', 'V1', 'V2L', 'Visual cortex', 'Visual cortex-1', 'Visual Cortex', 'visual cortex', 'V1B', 'V1M', 'VM-1', 'V2MM'],
                                     'Motor Cortex': ['M1', 'M2'],
                                     'Prefrontal Cortex':['FrA', 'MO', 'VO' , 'DLO', 'LO', 'A32', 'A24a', 'cg', 'Cg1', 'Cg2', 'A25', 'A24b', 'Fr3', 'DP', 'IL', 'PrL'],
                                     'Association Cortex': ['parietal cortex', 'Ect','PRh', 'Prh', 'DLEnt', 'DIEnt', 'VIEnt', 'VlEnt', 'MEnt', 'TeA','A29c', 'A29', 'A30', 'a30', 'A29a', 'A29b', 'CEnt', 'DlEnt', 'VLEnt', 'VIEnt', 'Vlent', 'LEnt'],
                                     'Piriform Cortex':['Pir'],
                                     'Olfactory areas': ['CxA','olfactory bulb', 'olfacotry bulb', 'AON', 'tenia tecta', 'taenia tecta', 'teania tectum', 'taenia tectum', 'endopiriform nucleus', 'endopirform nucleus', 'endopiriform nuncleus', 'EndoPir','Endo Pir', 'EnoPir', 'VTT', 'LOT', 'DTT'],
                                     'Claustrum': ['claustrum'],
                                     'Striatum':['NAcc', 'NAc', 'Nac', 'Nacc', 'Nacsh','NAcsh', 'NaSh', 'IPAC','TU','Tu','CPu', 'CPU', 'Cpu', 'LSS', 'CB', 'LSD', 'Nacc-1'],
                                     'Pallidum':['VP','GP','BNST', 'BST', 'MS', 'septum'],
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


def sort_into_hierarchies(density,strength):
    
    higher_hierarchy_nomenclature = {'Insular Cortex': ['AI', 'AID', 'AIV', 'GI', 'DI', 'AIP'],
                                     'Sensory Cortex': ['S1','S1J', 'S2', 'Au1', 'AuD', 'Auditory cortex', 'auditory cortex', 'AuV', 'V1', 'V2L', 'Visual cortex', 'Visual cortex-1', 'Visual Cortex', 'visual cortex', 'V1B', 'V1M', 'VM-1', 'V2MM'],
                                     'Motor Cortex': ['M1', 'M2'],
                                     'Prefrontal Cortex':['FrA', 'MO', 'VO' , 'DLO', 'LO', 'A32', 'A24a', 'cg', 'Cg1', 'Cg2', 'A25', 'A24b', 'Fr3', 'DP', 'IL', 'PrL'],
                                     'Association Cortex': ['parietal cortex', 'Ect','PRh', 'Prh', 'DLEnt', 'DIEnt', 'VIEnt', 'VlEnt', 'MEnt', 'TeA','A29c', 'A29', 'A30', 'a30', 'A29a', 'A29b', 'CEnt', 'DlEnt', 'VLEnt', 'VIEnt', 'Vlent', 'LEnt'],
                                     'Piriform Cortex':['Pir'],
                                     'Olfactory areas': ['CxA','olfactory bulb', 'olfacotry bulb', 'AON', 'tenia tecta', 'taenia tecta', 'teania tectum', 'taenia tectum', 'endopiriform nucleus', 'endopirform nucleus', 'endopiriform nuncleus', 'EndoPir','Endo Pir', 'EnoPir', 'VTT', 'LOT', 'DTT'],
                                     'Claustrum': ['claustrum'],
                                     'Striatum':['NAcc', 'NAc', 'Nac', 'Nacc', 'Nacsh','NAcsh', 'NaSh', 'IPAC','TU','Tu','CPu', 'CPU', 'Cpu', 'LSS', 'CB', 'LSD', 'Nacc-1'],
                                     'Pallidum':['VP','GP','BNST', 'BST', 'MS','septum'],
                                     'Basal Forebrain':['SIB', 'diagonal band', 'HDB', 'LDB', 'B'],
                                     'Hippocampus': ['vHip', 'vHIp', 'AHiP', 'DG', 'PaS', 'DS', 'POST', 'Post', 'PROST', 'PrS', 'Py', 'SubC', 'VS', 'STr', 'Str'],
                                     'Amygdala': ['CeA', 'CeM', 'CeC', 'CeL', 'PMCo', 'PLCo','AA', 'BLP', 'LA', 'BMA', 'BMP', 'EA','BLA', 'BLS', 'BLv', 'BLV', 'CeA', 'MeA', 'MeAD', 'ACo','AHiAL', 'AHiPL', 'I', 'IM', 'MeP', 'STIA', 'APir','RAPir', 'RaPir'],    
                                     'Hypothalamus':['Arc', 'AHA', 'AH','LH', 'MHb', 'STh', 'PSTh', 'VMH', 'MPOM', 'MPA', 'hypothalamus', 'F', 'dorsal hypothalamus','AHC', 'LPO' ,'MnPO', 'MCLH', 'ZI', 'ZIR', 'LHb', 'LHB', 'RML', 'PH', 'SubI', 'Subl', 'Sub1', 'PHD', 'RI', 'RMM', 'PHA'],
                                     'Thalamus':['DM','AD','MD','MGN', 'MGm', 'MG', 'MGV', 'MGv', 'medial geniculate', 'VPL', 'VPM', 'VL', 'VM', 'VPPC', 'VPLpc', 'VPLPC', 'VPMpc', 'VPLPc', 'VPMpv', 'CM', 'AM', 'AV', 'Po', 'PO', 'PV', 'PVA', 'PaF', 'PIL', 'PoT', 'LD', 'IAM', 'Re', 'Rh', 'SPFPC', 'SPFp', 'SPF', 'medial geniculate', 'CM', 'REth', 'Rt', 'Sub', 'SubD', 'PC', 'SG', 'sm', 'Gus', 'DLG', 'Eth', 'PrG', 'SO', 'SubV', 'CL', 'DT', 'DTg', 'LP', 'LVe', 'RRe', 'VA'],
                                     'Midbrain':['substsantia nigra', 'Substantia nigra', 'SN', 'cp', 'SNc', 'SNC', 'SNr', 'SNR', 'dorsal raphe', 'median raphe', 'medial raphe', 'isRt', 'CIC', 'PAG', 'DLPAG', 'VTA', 'VTAR', 'RLi', 'Rli', 'mRt','PVG', 'EW', 'LDTgV', 'PrEW', 'VLPAG', 'p1Rt', 'PBP', 'DA8', 'DAB', 'RPC', 'RP', 'PR', 'APT', 'APTD', 'BIC', 'ECIC','DpG', 'DR', 'IP', 'IPD', 'IRt', 'DCIC', 'MnR', 'RMC', 'RMg', 'RRF', 'CI', 'CLi', 'PN', 'Pa4', 'PaR', 'RPa', 'SubB', 'InC', 'PMnR', 'IF', 'InG', 'PIF', 'RM', 'IPR', 'DpWh', 'PMnT', 'SuG', 'In', 'DpMe'],
                                     'Hindbrain':['PBN', 'parabrachial nucleus', 'NTS', 'MPB', 'LPB', 'PCRtA', 'PTg', 'PrCnF', 'SPTg', 'PnO', 'PnV', 'MPL', 'PL', 'KF', '2Cb', 'SubCV', '5', 'VeCb', 'CAT', 'CnF', 'DMTg', 'MiTg', 'PBG', 'LSO', 'PLV', 'PPTg', 'Pn', 'PnC', 'PNV', 'Pr5', 'RtTg', 'SPO', 'Tz', 'VLL', '7N', 'B9', 'Bar', 'CG', 'DPGi', 'Gi', 'GiA', 'IRtA', 'Int', 'IntA', 'LC', 'LPGi', 'LPGiA', 'MSO', 'MVe', 'MVeMC', 'MVePC', 'P7', 'PM', 'Ve', 'Pa5', 'Pr', 'Su5', 'VTg', 'InCo', 'MTg', 'Sp5', 'SuVe'] }
    
    lower_order_hierarchy = {
    
            #Sensory Cortex        
            'Somatosensory Cortex, primary': ['S1','S1J'],
            'Somatosensory Cortex, secondary': ['S2'],
            'Auditory Cortex, primary': ['Au1', 'Auditory cortex', 'auditory cortex'],
            'Auditory Cortex, secondary or higher': ['AuD', 'AuV'],
            'Visual Cortex, primary': ['V1', 'V1B', 'V1M', 'VM-1', 'Visual cortex', 'Visual cortex-1', 'Visual Cortex', 'visual cortex'],
            'Visual Cortex, secondary or higher': ['V2L', 'V2ML','V2MM'],
            'Piriform Area':['Pir'],
            
            #Motor Cortex
            'Motor Cortex, primary': ['M1', 'Fr3'],
            'Motor Cortex, secondary': ['M2'],
            
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
            'Retrosplenial Cortex': ['A30', 'a30', 'A29c', 'A29b', 'A29a', 'A29', 'RSA', 'RSG', 'RSGa', 'RSGb' ],
            'Entorhinal Cortex': ['DLEnt', 'Dlent', 'DIEnt', 'DlEnt', 'VIEnt', 'Vlent', 'VlEnt', 'MEnt', 'CEnt', 'LEnt'],
            
            #Insular Ctx
            'Agranular Insular Cortex': ['AI'],
            'Agranular Insular Cortex, dorsal': ['AID'],
            'Agranular Insular Cortex, ventral': ['AIV'],
            'Agranular Insular Cortex, posterior': ['AIP'],
            'Dysgranular Insular Cortex': ['DI'],
            'Granular Insular Cortex': ['GI'],
            
            #Olfactory Areas
            'Anterior Olfactory Nucleus': ['AON', 'AOP', 'AOM', 'AOL', 'AOD'],
            'Piriform-Amygdalar Area':['CxA'],
            'Endopiriform Nucleus': ['DEn', 'VEn', 'endopiriform nucleus', 'endopirform nucleus', 'endopiriform nuncleus', 'EndoPir', 'Endo Pir', 'EnoPir'],
            'Tenia Tecta': ['DTT', 'VTT', 'tenia tecta', 'taenia tecta', 'teania tectum', 'taenia tecta', 'taenia tectum'],
            'Nucleus of the Lateral Olfactory Tract':['LOT'],
            'Main Olfactory Bulb':['olfactory bulb', 'olfacotry bulb'],
                        
            #Claustrum
            'Claustrum': ['claustrum'],
            
            #Striatum
            'Caudate Putamen':['CPu', 'CPU', 'Cpu', 'LSS', 'LSD'],
            'Nucleus Accumbens, core':['NAcc', 'Nacc', 'NAc', 'Nac', 'Nacc-1'],
            'Nucleus Accumbens, shell': ['NAcsh', 'Nacsh', 'NaSh'],
            'Fundus of Striatum (IPAC)':['IPAC'],
            'Olfactory Tubercle':['Tu', 'CB','TU'],
            
            #Pallidum
            'Ventral Pallidum':['VP'],
            'Bed Nuclei of the Stria Terminalis (BNST)': ['BNST', 'BST'],
            'Globus Pallidus': ['GP', 'LGP', 'MGP'],
            'Medial Pallidum': ['MS'],
            
            #Hippocampus
            'Hippocampus (CA1,CA2,CA3)': ['CA1', 'CA2', 'CA3', 'fi', 'LMol', 'Or', 'Py', 'Rad', 'STr', 'Str', 'SLu', 'alv', 'vHip', 'vHIp', 'AHiP', 'PaS', 'DS', 'POST', 'PROST', 'Post', 'PrS', 'SubC', 'VS'],
            'Hippocampus (Dentate Gyrus)': ['DG', 'Mol'],
            
            #Amygdala
            'Lateral Amygdala' :['LA'],
            'Basolateral Nucleus': ['BLA','BLP','BL', 'BLv', 'BLV', 'BLS'],
            'Basomedial Nucleus': ['BMA', 'BMP', 'IM'],
            'Cortical Amygdala': ['PMCo', 'PLCo', 'ACo'],
            'Anterior Amygdaloid Area':['AA', 'AAD', 'AAV', 'I'],
            'Amygdalopiriform Transition Area': ['RAPir', 'RaPir', 'APir', 'AHiAL', 'AHiPL'],
            'Central Nucleus of Amygdala':['CeC', 'CeL', 'CeA', 'CeM', 'CeMAD', 'CeMAV', 'Ce', 'STIA'],            
            'Medial Amygdala': ['MeA', 'MeAD', 'MeAV', 'MePD', 'MeP'],
            'Extended Amygdala': ['EA'],
            
            #Basal Forebrain
            'Basal Forebrain':['MDB', 'LDB', 'diagonal band', 'VDB', 'HDB', 'SIB', 'SI/B', 'SI', 'B'],
            
            #Hypothalamus
            'Hypothalamic Lateral Zone': ['Arc', 'LH', 'LHA', 'STh', 'PSTh', 'MCLH', 'ZI', 'ZIR', 'F', 'hypothalamus', 'dorsal hypothalamus', 'LPO', 'RML', 'SO', 'SubI', 'Subl','Sub1', 'RMM'],
            'Hypothalamic Medial Zone' :['MnPO', 'MPO', 'MPA', 'MPOM', 'AHC', 'MHb', 'PH', 'AHA', 'VMH', 'PHD', 'RI', 'PHA', 'AH'],
                        
            #Thalamus
            'Habenular Nuclei': ['LHb', 'LHB'],
            'Ventral Group of the Dorsal Thalamus':['VPL','VPM','VL','VM', 'VPPC', 'VPLpc', 'VPMpc' ,'VPLPC', 'VPLPc', 'Gus', 'RRe'],
            'Medial Geniculate Nuclei': ['MGN','MGM','MGv', 'MGV', 'MGm','medial geniculate', 'MG'],            
            'Posterior Complex of the Thalamus':['Po', 'PO', 'PoT', 'PIL', 'REth', 'SPFPC', 'SPFp', 'SG', 'Eth', 'SPF', 'DT', 'DTg'],
            'Mediodorsal Thalamic Nuclei': ['DM', 'MD', 'MDC','MDL','MDM', 'Sub', 'SubD', 'LP'],
            'Anterodorsal Thalamic Nucleus': ['AD', 'AM', 'AV', 'LD', 'DLG', 'PrG'],
            'Intralaminar Nuclei of the Dorsal Thalamus': ['Re', 'Rh', 'CM', 'PF', 'IAM', 'PaF', 'Rt', 'PC', 'SubV', 'CL', 'VA'],
            'Paraventricular Thalamic Nucleus': ['PVT', 'PVA', 'PV', 'PVG', 'sm'],
            
            #Midbrain
            'Periaqueductal Gray':['PAG', 'EW', 'PrEW', 'VLPAG', 'DLPAG', 'InC', 'In'],     
            'Substantia Nigra': ['SNc', 'SNC', 'SNr', 'SNR','SN', 'substantia nigra', 'Substantia nigra', 'cp', 'PR', 'IP', 'IPD', 'RMC', 'CI', 'PN'],
            'Ventral Tegmental Area': ['VTA', 'PBP', 'Pa4', 'PIF', 'VTAR'],
            'Raphe Nuclei': ['dorsal raphe', 'DR', 'median raphe', 'medial raphe', 'RLi', 'Rli', 'MnR', 'PMnR', 'PMnT', 'RMg', 'CLi', 'RPa', 'IF', 'RM', 'IPR'],
            'Midbrain Reticular Nucleus':['isRt', 'mRt', 'p1Rt', 'RPC', 'RP', 'DA8', 'DAB', 'APT', 'APTD', 'BIC', 'IRt', 'RRF', 'PaR', 'SubB', 'DpMe'],
            'Superior Colliculus':['InG', 'DpWh', 'SuG', 'DpG'],
            'Inferior Colliculus':['CIC', 'ECIC', 'DCIC'],
            
            #Hindbrain
            'Parabrachial Nucleus': ['PBN', 'LPB', 'MPB', 'parabrachial nucleus'],
            'Pedunculopontine Nucleus': ['PTg', 'LDTgV', 'SPTg', 'MPL', 'PL', 'KF', 'PrCnF', 'PnO', 'PnC', 'InCo', 'Pn', 'PnV', 'VeCb', 'CnF', 'CAT', 'DMTg', 'MiTg', 'MTg', 'PBG', 'LSO', 'PLV', 'PPTg', 'VTg', 'PNV', 'Pr5', 'RtTg', 'SPO', 'Tz', 'VLL', 'B9', 'Bar', 'CG', 'Gi', 'GiA', 'LC', 'MSO', 'Su5'],
            'Medulla': ['NTS', 'PCRtA', '7N', 'DPGi', 'IRtA', 'LPGi', 'LPGiA', 'LVe', 'MVe', 'MVeMC', 'MVePC', 'P7', 'Pa5', 'Pr', 'Ve', 'Sp5', 'SuVe'],
            'Cerebellum':['2Cb', 'Int', 'IntA', 'PM'],
            'Motor Trigeminal Nucleus':['5', 'SubCV']
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

    for lower_region, ROI_name in lower_order_hierarchy.items():
        counter[lower_region] = 0.0
        
    # average of ROI densities, if they are grouped in lower order hierarchy
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
    counter.rename('average pixel density')
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
         # average of ROI densities, if they are grouped in lower order hierarchy
        index_counter=0
    
    counter.rename('average projection strength')
    counter_strength = counter.copy()
    
    
    dataframe = pd.concat([counter_density, counter_strength], axis=1)
    dataframe.fillna(value=0, inplace=True)
    column_names = ['average pixel density [%]','average projection strength [%]']
    dataframe.columns= column_names
    
    return dataframe, missing_rois_lower, missing_rois_higher



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


for filename in os.listdir('.'):
    if "AAV" in filename and ".csv" in filename:
    
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
        df = exclude_regions(df)
        df = remove_starter_volume(df, start_Bregma, end_Bregma)
        
        # get projection strength in absolute and relative values for all ROI labels
        pixel_count_ROI, pixel_count_ROI_percent = pixel_count_per_ROI(df)
        
        #get averaged pixel densities in percent for all ROI labels
        average_pixel_density = average_density_per_ROI(df)
        
        #pixel count per bregma
        pixel_count_bregma = pixel_count_per_Bregma(df)
        
        #sort values into lower and higher order categories
        lower_hierarchy_df, missing_ROIs_lower, missing_ROIs_higher = sort_into_hierarchies(average_pixel_density, pixel_count_ROI_percent)
        lower_hierarchy_df[lower_hierarchy_df['average projection strength [%]'] < 0.03] = 0

        higher_hierarchy= sort_into_higher_hierarchy_in_percent(pixel_count_ROI)
        
        #assign lower order hierachry and higher order hierarchy labels to df and add as separate columns
        df = assign_level_1(df)
        df = assign_level_2(df)
        
        
        
############################################################################################################
        
        #AMYGDALA ANALYSIS
        
        amygdala_only_pixelnumber, amygdala_only_percent = sort_amygdala(pixel_count_ROI,df)
        amygdala_only_df = pd.concat([amygdala_only_pixelnumber, amygdala_only_percent],  axis=1)
        amygdala_only_df.columns =['pixel_count', 'percent']
        amygdala_only_df[amygdala_only_df['percent'] < 0.03] = 0
        
        #remove all ROIs that are not amygdala from the df
        df_amygdala = df.copy()
        df_amygdala = amygdala_only(df_amygdala)
        
        #relabel amygdala lowest lever ROIS, like AA, AAD, AAV into "anterior amygdala"
        df_amygdala = df_amygdala.assign(Label=sort_amygdala_pivot(df_amygdala))
        
        #create pivot tables to plot Anterior-Posterior graphs
        #pivot_table_counts_amy = pivot_counts(df_amygdala)
        pivot_table_percent_amy = pivot_percent(df_amygdala)
        pivot_table_density_amy = pivot_density(df_amygdala)
        
############################################################################################################
        
        #THALAMUS ANALYSIS
        
        # sort and process thalamus ROIs
        thalamus_only_pixelnumber, thalamus_only_percent = sort_thalamus(pixel_count_ROI,df)
    
        thalamus_only_df = pd.concat([thalamus_only_pixelnumber, thalamus_only_percent],  axis=1)
        thalamus_only_df.columns =['cell_count', 'percent']
        thalamus_only_df[thalamus_only_df['percent'] < 0.03] = 0
        
        
        #remove all ROIs that are not amygdala from the df
        df_thalamus = df.copy() 
        df_thalamus = th_only(df)
    
        #relabel amygdala lowest lever ROIS, like AA, AAD, AAV into "anterior amygdala"
        df_thalamus = df_thalamus.assign(Label=sort_th_pivot(df_thalamus))
    
        #create pivot tables to plot Anterior-Posterior graphs
        #pivot_table_counts_th = pivot_counts(df)
        pivot_table_percent_th = pivot_percent(df_thalamus)
        pivot_table_density_th = pivot_density(df_thalamus)
        
############################################################################################################
                
        #STRIATUM ANALYSIS
        # sort and process thalamus ROIs
        striatum_only_pixelnumber, striatum_only_percent = sort_striatum(pixel_count_ROI,df)
    
        striatum_only_df = pd.concat([striatum_only_pixelnumber, striatum_only_percent],  axis=1)
        striatum_only_df.columns =['cell_count', 'percent']
        striatum_only_df[striatum_only_df['percent'] < 0.03] = 0
        
        
        #remove all ROIs that are not amygdala from the df
        df_striatum = df.copy() 
        df_striatum = striatum_only(df)
    
        #relabel amygdala lowest lever ROIS, like AA, AAD, AAV into "anterior amygdala"
        df_striatum = df_striatum.assign(Label=sort_striatum_pivot(df_striatum))
    
        #create pivot tables to plot Anterior-Posterior graphs
        #pivot_table_counts_str = pivot_counts(df)
        pivot_table_percent_str = pivot_percent(df_striatum)
        pivot_table_density_str = pivot_density(df_striatum)
        
#############################################################################################################
        
        save_path = file_name 
        
        # Create target Directory if don't exist
        if not os.path.exists(save_path):
            os.mkdir(save_path)
            print("Directory " , save_path ,  " Created ")
        else:    
            print("Directory " , save_path ,  " already exists")

        # save everything as .csv files
        
        #global data
        df.to_csv("{0}\{1}_summary_raw.csv".format(save_path, file_name))
        lower_hierarchy_df.to_csv("{0}\{1}_sorted_lower_hierarchy.csv".format(save_path, file_name))
        higher_hierarchy.to_csv("{0}\{1}_pixel_count_highest_hierarchy.csv".format(save_path, file_name))
        pixel_count_ROI.to_csv("{0}\{1}_pixel_count_ROIs.csv".format(save_path, file_name))
        pixel_count_ROI_percent.to_csv("{0}\{1}_pixel_count_ROIs_percent.csv".format(save_path, file_name))
        pixel_count_bregma.to_csv("{0}\{1}_pixel_count_bregma.csv".format(save_path, file_name))
        
        #amygdala analysis
        amygdala_only_df.to_csv("{0}\{1}_amygdala_count_and_percent.csv".format(save_path, file_name))
        pivot_table_percent_amy.to_csv("{0}\{1}_pivot_percent_of_amygdala_input.csv".format(save_path, file_name))
        #pivot_table_counts_amy.to_csv("{0}\{1}_pivot_raw_counts_amygdala.csv".format(save_path, file_name))
        pivot_table_density_amy.to_csv("{0}\{1}_pivot_density_amygdala.csv".format(save_path, file_name))
        
        
        #thalamus analysis
        thalamus_only_df.to_csv("{0}\{1}_thalamus_count_and_percent.csv".format(save_path, file_name))
        pivot_table_percent_th.to_csv("{0}\{1}_pivot_percent_of_thalamus_input.csv".format(save_path, file_name))
        #pivot_table_counts_th.to_csv("{0}\{1}_pivot_raw_counts_thalamus.csv".format(save_path, file_name))
        pivot_table_density_th.to_csv("{0}\{1}_pivot_density_thalamus.csv".format(save_path, file_name))
        
        #striatum analysis
        
        striatum_only_df.to_csv("{0}\{1}_striatum_count_and_percent.csv".format(save_path, file_name))
        pivot_table_percent_str.to_csv("{0}\{1}_pivot_percent_of_striatum_input.csv".format(save_path, file_name))
        #pivot_table_counts_str.to_csv("{0}\{1}_pivot_raw_counts_striatums.csv".format(save_path, file_name))
        pivot_table_density_str.to_csv("{0}\{1}_pivot_density_striatum.csv".format(save_path, file_name))
        
        df_big = df_big.append(df, ignore_index=True)
    else:
        continue


df_big.to_csv("all axonal tracings combined.csv")
        
        