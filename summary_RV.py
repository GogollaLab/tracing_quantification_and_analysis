# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 22:04:26 2017

@author: gehrlach
script should find .xls files that contain summary in their names.

then read out total cell count and area and supdate it into a dataframe
using distance from bregma and the ROI name

"""

import pandas as pd
import os

df = pd.DataFrame()


for file in os.listdir('.'):
    filename = file[:-4]
   
    if 'summary' in file:
        
        # parse bregma and ROI out of filename
        bregma = file.split('_')[0]
        section = file.split('_')[1]
        roi = file.split('_')[3][:-4]
                
        #read table
        df = pd.read_table(file)
        
        #add bregma and roi column
        df['Section'] = section
        df['Bregma'] = bregma
        df['ROI'] = roi
        
        
        #remove slice column
        #df = df.drop('Slice',1)
        df = df.drop('%Area', 1)
       
        #calculate cell density and append to dataframe
        
        
        #append to summary dataframe
        d2= df2.append(df)
        
df.to_csv('all cells of all ROI analysis.csv')

##############


