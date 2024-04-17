#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make full report for each run and group

@author: lauracarlton
"""

from SQE_metrics import generate_report, snr, GVTD, sci_psp
import os
from mne_bids import get_entity_vals, BIDSPath, read_raw_bids
import numpy as np 
from pysnirf2 import Snirf
import pathlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def generateFullReport(rootDir, excluded = None, saveFig=0, remove_short=0):
    
    # get list of all the subjects
    all_subjs = get_entity_vals(rootDir,'subject')
    
    # if there are subjects to exclude, remove them here 
    if excluded != None:
        for e in excluded: 
            all_subjs.remove(str(e))
    
    # get list of all the sessions
    sessions = get_entity_vals(rootDir, 'session')
    
    # get list of all the tasks
    tasks = get_entity_vals(rootDir,'task')
    
    # get list of all the runs 
    runs = get_entity_vals(rootDir, 'run')
    
    # if any of the entities returns an empty list, set it equal to [None]
    if len(all_subjs) == 0 :
        all_subjs = [None]
    
    if len(sessions) == 0: 
        sessions  = [None]

    if len(tasks) == 0:
        tasks = [None]
    
    if len(runs) == 0:
        runs = [None]
    
    # if derivatives/DQE/ directory does not exist, then create it 
    if not os.path.exists(rootDir+'derivatives/DQE'):
        os.makedirs(rootDir+'derivatives/DQE')

    ## loop through all sessions, runs, tasks and subjects ##
    for ses in sessions:
        
        for run in runs:
            
            for task in tasks:
                
                # initialize empty arrays 
                paths = []
                all_snirfs = []
                mean_snr_runs = []
                median_GVTD_runs = []
                mean_sci_runs = []
                mean_psp_runs = []
                mean_percent_runs = []
                run_names = []
                
                for sub in all_subjs:
                    
                    # get the bids path from the subject, task, session and run 
                    bids_path = BIDSPath(root=rootDir, subject=sub, task=task, session=ses, run=run, datatype='nirs', suffix='nirs', extension='.snirf')
                    savePath = bids_path.root.as_posix() + '/derivatives/DQE/'
                    
                    # try reading the path to ensure it is valid 
                    try:
                        snirfPath = read_raw_bids(bids_path)
                    except:
                        print('could not load ' + sub)
                        continue
                    
                    # save path in list 
                    paths.append(bids_path.fpath.as_posix())
                    
                    # load the snirf file 
                    snirf = Snirf(bids_path.fpath.as_posix())
                    all_snirfs.append(snirf)
                    
                    # get the SNR for both wavelengths 
                    snr_valsa = snr(snirf, lam=0, ax=None, remove_short=remove_short)
                    snr_valsb= snr(snirf, lam=1, ax=None, remove_short=remove_short)
                    snr_vals = np.hstack([snr_valsa, snr_valsb])
                    mean_snr = np.mean(snr_vals) # calculate the mean across both wavelengths
                    mean_snr_runs.append(mean_snr)
                    
                    # calculate GVTD 
                    GVTD_vals = GVTD(snirf, ax=None)
                    median_GVTD = np.median(GVTD_vals) # calculate the median 
                    median_GVTD_runs.append(median_GVTD)
                    
                    # calculate sci and psp and get the channel percentage 
                    sci, psp, chan_percent =  sci_psp(snirf, mode=None, ax=None, remove_short=remove_short)
                    mean_sci = np.mean(sci, axis=(0,1)) # calculate mean sci
                    mean_psp = np.mean(psp, axis=(0,1)) # calculate mean psp
                    mean_sci_runs.append(mean_sci)
                    mean_psp_runs.append(mean_psp)
                
                    mean_percent = np.mean(chan_percent)*100 # calculate mean percent of good channel
                    mean_percent_runs.append(mean_percent)
                  
                    
                    name = 'sub-' + sub
                    run_names.append(sub)

                    savePath_scan = savePath + name
                    title = 'Data Quality Evaluation: sub-' + sub + '; task-' + task
                    fileName = '/DQE_sub-' + sub + '_task-' + task 
                    
                    if ses != None:
                        title += '; session-' +  ses  
                        fileName += '_ses-' + ses
                        savePath_scan += 'ses-'+ses
                        
                    if run != None:
                        title += '; run-' + run
                        fileName += '_run-' + run            
                    
                    if not os.path.exists(savePath_scan):
                            os.makedirs(savePath_scan)
                            
                    generate_report(snirf, savePath = savePath_scan+fileName, title = title, remove_short=remove_short)
                   
                
                ## compile mean sci, psp, perc, snr, GVTD for table ##
                sci_allruns = np.array([round(a, 2) for a in mean_sci_runs])
                psp_allruns = np.array([round(a, 2) for a in mean_psp_runs])
                perc_allruns = np.array([round(a, 2) for a in mean_percent_runs])
                snr_allruns = np.array([round(a, 2) for a in mean_snr_runs])
                GVTD_allruns = np.array([round(a, 2) for a in median_GVTD_runs])
                
                ## set text for the body of the table
                cellText = np.vstack([sci_allruns, psp_allruns, snr_allruns, GVTD_allruns,perc_allruns]).T
                
                # initialize the figure
                plt.rcParams.update({'font.size': 10})
                fig = plt.figure(figsize=(20,10))
                
                # set the title and savePath depending on task, session and run 
                title = 'Group level report: task-' + task 
                savePath_group = savePath + 'DQEgroup_task-' + task
                
                if ses != None:
                    title += '; session-' +  ses  
                    savePath_group += '_session-' + ses
                if run != None:
                    title += '; run-' + run
                    savePath_group += '_run-' + run
                
                
                savePath_group += '.pdf'
                
                ## generate the figure 
                fig.suptitle(title)
                gs = gridspec.GridSpec(3, 2, figure=fig)
                
                ax1 = fig.add_subplot(gs[0,0])
                ax1.bar(run_names, mean_percent_runs)
                ax1.set_title('% of run above SCIxPSP threshold')
                plt.xticks(rotation=30)

                ax2 = fig.add_subplot(gs[1,0])
                ax2.bar(run_names, snr_allruns)
                ax2.set_title('Mean SNR')
                plt.xticks(rotation=30)

                ax3 = fig.add_subplot(gs[2,0])
                ax3.bar(run_names, GVTD_allruns)
                ax3.set_title('Median GVTD')
                plt.xticks(rotation=30)
                
                ax4 = fig.add_subplot(gs[:,1])
                rowColours = ['y']*len(all_subjs)
                colColours = ['y']*4 + ['r']
                the_table = ax4.table(cellText, rowLabels= run_names, colLabels = ["mean SCI", "mean PSP", "mean SNR", "median GVTD", "% above threshold"], loc='center', rowColours = rowColours, colColours = colColours)
                the_table.scale(1.1,3)
                the_table.auto_set_font_size(False)
                ax4.axis('off')
                plt.tight_layout()
                
                plt.savefig(savePath_group, dpi=1200)
      

                resultDict = {}
                fieldName = 'task-' + task 
                        
                if ses != None:
                    fieldName += '_session-' + ses
                if run != None:
                    fieldName += '_run-' + run
                
                resultDict[fieldName + '_snr'] = snr_allruns
                resultDict[fieldName + '_GVTD'] = GVTD_allruns
                resultDict[fieldName + '_mean_percent_runs'] = mean_percent_runs
                resultDict[fieldName + '_psp_allruns'] = psp_allruns
                resultDict[fieldName + '_sci_allruns'] = sci_allruns
                
    
    return resultDict
