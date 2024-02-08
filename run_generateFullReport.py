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

## replace with your root directory to BIDS data
rootDir = '/Users/lauracarlton/Downloads/Co-Location Study 2/'

generateFullReport(rootDir, saveFig=1, excluded=None, remove_short=0)





