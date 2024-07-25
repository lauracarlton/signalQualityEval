#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run this from a folder that contains the following files:
    - SQE_metrics
    - generateFullReport
    - 10-5-System_Mastoids_EGI129.xlsx

rootDir = path to the BIDS compliant data
saveFig = 1 if you want to save the generated figures
excluded = specify any subjects you want to exclude
remove_short = 1 if you want to remove the short channels from the plot
    
@author: lauracarlton
"""

from generateFullReport import generateFullReport
import pandas as pd 

## replace with your root directory to BIDS data
rootDir = '/Users/lauracarlton/Downloads/Co-Location Study 2/'

resultDict = generateFullReport(rootDir, saveFig=1, excluded=None, remove_short=0)

#%% extract SNR

task = 'task-Simon'

snr = resultDict[task + '_snr']

# to save this as a csv
df = pd.DataFrame(resultDict)
df.to_csv(rootDir + '_DQR_output.csv', index=None)


