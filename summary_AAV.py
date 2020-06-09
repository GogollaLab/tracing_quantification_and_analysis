# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 22:04:26 2017

@author: gehrlach

script finds .xls files that contain summary in their names.

then reads out total cell count and area, supdates it into a dataframe
using distance from bregma and the ROI name

"""

import pandas as pd
import os

df2 = pd.DataFrame()
df = pd.DataFrame()


for file in os.listdir('.'):
    filename = file[:-4]
    
    if 'axon analysis' in file:
        
        # parse bregma and ROI out of filename
        bregma = file.split('_')[0]
        section = file.split('_')[1]
        
        #read table
        df = pd.read_csv(file, sep ='\t' , encoding='latin1')
        
        #add bregma and roi column
        df['Section'] = section
        df['Bregma'] = bregma
        
        #append to summary dataframe
        df2= df2.append(df)
        
df2 = df2.drop(' ',1)    
df2.to_csv('filename.csv')

##############


