#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate a DQR for each run for a subject in the motion artefact study 
@author: lauracarlton
"""

from pysnirf2 import Snirf
from SQE_metrics import generate_report
import os

#%%

subj = '4'
rootDir = '/Users/lauracarlton/Library/CloudStorage/GoogleDrive-lcarlton@bu.edu/My Drive/fNIRS/Data/motionArtefactStudy/subj-' + subj + '/'
saveDir = rootDir  + 'DQR_export/'

if not os.path.exists(saveDir):
    os.makedirs(saveDir)
    
#%%

tasks = [#['RS',['01','02']],
         # ['WM', ['01', '02', '03', '04']],
         ['MA', ['01']],
         # ['MAaudio', ['01']],
          ['squats', ['01']],
         ['Wcont', ['01']],
         ['WSalt', ['01']]  ]


for task in tasks:
    
    for run in task[1]:
        
        runPath = 'subj-' + subj + '_task-' + task[0] + '_run-' + run 
        fileName = rootDir + runPath + '.snirf'
        saveName = saveDir + runPath + '_DQR.jpeg'
        
        snirf_obj = Snirf(fileName, 'r+')
        
        generate_report(snirf_obj, title = runPath, remove_short=1, savePath = saveName)
                
        





