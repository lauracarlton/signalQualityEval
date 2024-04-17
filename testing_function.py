#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:44:26 2023

@author: lauracarlton
"""

'''n
things to do:
    - figure out best way to display channel names in SCIxPSP
    - put legend on the top of the figure 
    
'''
import numpy as np
from pysnirf2 import Snirf
from SQE_metrics import generate_report  , GVTD, sci_psp, snr, sort_Data, filter_data
import matplotlib.pyplot as plt
# from group_level_report import groupLevelReport
import mne 
from mne.preprocessing.nirs import optical_density, _validate_nirs_info
from mne_nirs.preprocessing import scalp_coupling_index_windowed, peak_power
from mne_nirs.visualisation import plot_timechannel_quality_metric
from itertools import compress
from scipy import signal

#%% test scan level report



# fileName_walkingRun = '/Users/lauracarlton/Library/CloudStorage/GoogleDrive-lcarlton@bu.edu/.shortcut-targets-by-id/1Qja7XXGiJbseXCER2cQ-c7FviFH3swGF/NN22 Data Transfer/231027/ninjaNIRS2022_2023-10-27-14-56-01.snirf'
# savePath_walkingRun = '/Users/lauracarlton/Library/CloudStorage/GoogleDrive-lcarlton@bu.edu/.shortcut-targets-by-id/1Qja7XXGiJbseXCER2cQ-c7FviFH3swGF/NN22 Data Transfer/231027/ninjaNIRS2022_2023-10-27-14-56-01_DQR.jpeg'

# fileName_motionRun = '/Users/lauracarlton/Library/CloudStorage/GoogleDrive-lcarlton@bu.edu/.shortcut-targets-by-id/1Qja7XXGiJbseXCER2cQ-c7FviFH3swGF/NN22 Data Transfer/231027/ninjaNIRS2022_2023-10-27-14-40-33.snirf'
# savePath_motionRun = '/Users/lauracarlton/Library/CloudStorage/GoogleDrive-lcarlton@bu.edu/.shortcut-targets-by-id/1Qja7XXGiJbseXCER2cQ-c7FviFH3swGF/NN22 Data Transfer/231027/ninjaNIRS2022_2023-10-27-14-30-33_DQR.jpeg'


# fileName1 = '/Users/lauracarlton/Library/CloudStorage/GoogleDrive-lcarlton@bu.edu/.shortcut-targets-by-id/1Qja7XXGiJbseXCER2cQ-c7FviFH3swGF/NN22 Data Transfer/230911/ninjaNIRS2022_2023-09-11-13-50-06.snirf'
# savePath1 = '/Users/lauracarlton/Library/CloudStorage/GoogleDrive-lcarlton@bu.edu/My Drive/fNIRS/Data/NN22/230911/ninjaNIRS2022_2023-09-11-13-50-06_SQE.jpeg'

# fileName1 = '/Users/lauracarlton/Google Drive/My Drive/fNIRS/Data/NN22 Data Transfer/230919/ninjaNIRS2022_2023-09-19-18-26-16.snirf'
# savePath1 = '/Users/lauracarlton/Google Drive/My Drive/fNIRS/Data/NN22 Data Transfer/230919/ninjaNIRS2022_2023-09-19-18-26-16'
# fileName1 = '/Users/lauracarlton/fNIRS/Data/NN22/2023-08-02/ninjaNIRS2022_2023-08-03-17-38-32_fingerTap1.snirf'
# savePath1 = '/Users/lauracarlton/fNIRS/Data/NN22/2023-08-02/export/ninjaNIRS2022_2023-08-03-17-38-32_fingerTap1_SQE.jpeg'

# fileName2 = '/Users/lauracarlton/fNIRS/Data/NN22/2023-08-02/ninjaNIRS2022_2023-08-03-17-45-38_fingerTap2.snirf'
# fileName = '/Users/lauracarlton/fNIRS/Data/NN22/2023-08-02/ninjaNIRS2022_2023-08-03-17-45-38_fingerTap2.snirf'
# savePath2 = '/Users/lauracarlton/fNIRS/Data/NN22/2023-08-02/export/ninjaNIRS2022_2023-08-03-17-45-38_fingerTap2_SQE.jpeg'

# title = ''
# title = 'walk/talk/stand paradigm test 2'
# load the snirf object 
# fig, ax = plt.subplots(1,1)
fileName = '/Users/lauracarlton/Library/CloudStorage/GoogleDrive-lcarlton@bu.edu/My Drive/fNIRS/Data/motionArtefactStudy/subj-3/subj-3_task-MA_run-01.snirf'

# snr_lam2 = snr(snirf_obj, lam=1,ax=ax, savePath=savePath1)

snirf_obj = Snirf(fileName, 'r+')


 #%%
fig2,ax = plt.subplots(1,1)
gvtd = GVTD(snirf_obj, ax=ax)


#%%
title = ''
generate_report(snirf_obj, title =title, remove_short=1)

# sortedData, channels = sort_Data(snirf_obj)
# time = nirs_obj_data[0].time
# nChannels = 567 
# cardiac_data = filter_data(sortedData, time, nChannels)

#%%
nirs_obj_data = snirf_obj.nirs[0].data

data = nirs_obj_data[0].dataTimeSeries
OD = data/np.mean(data,axis=0)
#%%
# fig, ax = plt.subplots(1,1)
# ax.plot(cardiac_data[82,:,0])
# ax.plot(cardiac_data[82,:,1])
#%%
# fig = plt.figure()
# ax = fig.add_subplot()
# nirs_obj_data = snirf_obj.nirs[0].data
# time = nirs_obj_data[0].time

# sortedData, channels = sort_Data(snirf_obj)
# snr_lam1 = snr(snirf_obj, lam=1,ax=None)

# snr_bad = np.where(snr_lam1 < 5)[0]

# fig, ax = plt.subplots(1,1)
# ax.plot(time, np.transpose(sortedData[snr_bad[10:20],:,0]))


#%%

# GVTD(snirf_obj, ax)
# sci_psp(snirf_obj, mode='thresholded', ax=ax,savePath = None)
# savePath = '/Users/lauracarlton/fNIRS/Data/jessie_lowDensityExample_report.png'
# savePath = '/Users/lauracarlton/fNIRS/Data/example_report_sub-1_task-MW_run-1_nirs.png'

#%% find number of channels with NaNs
chansNaN_lam1 = 0
chansNaN_lam2 = 0

for c in range(nChannels):
    
    if np.isnan(sortedData[c,:,0]).any():
        chansNaN_lam1 +=1

    if np.isnan(sortedData[c,:,1]).any():
        chansNaN_lam2 +=1


#%%
nirs_obj_data = snirf_obj.nirs[0].data
data = nirs_obj_data[0].dataTimeSeries
nMeasurements = np.shape(data)[1]
measurements = nirs_obj_data[0].measurementList
time = nirs_obj_data[0].time

# for iChan in range(nMeasurements):
#     chanData = data[:,iChan]
#     x = np.arange(len(chanData))
#     xp = np.arange(len(chanData))[np.isnan(chanData) == False]
#     fp = chanData[np.isnan(chanData) == False]
#     data[:,iChan] = np.interp(x, xp, fp)

fig, ax = plt.subplots(1,1)
ax.plot(time,data[:,545], label='S46-D76')
ax.plot(time,data[:,547], label='S46-D88')
ax.legend()
ax.set_xlabel('Time')

    #%%
measurements = nirs_obj_data[0].measurementList
aux = snirf_obj.nirs[0].aux[0].dataTimeSeries
aux_peaks = signal.find_peaks(aux)

lam=0
sortedData, channelNames = sort_Data(snirf_obj)
data_lam = sortedData[:,:,lam]
#%%

# A = np.isnan(data_lam)
mean_intensity = []

for c in range(len(measurements)//2):
    chan_data = data_lam[c,:][np.logical_not(np.isnan(data_lam[c,:]))]
    mean_intensity.append(np.mean(chan_data))
    

#%%
# nirs_obj_data = snirf_obj.nirs[0].data
# 
# data = nirs_obj_data[0].dataTimeSeries
# time = nirs_obj_data[0].time
# measurements = nirs_obj_data[0].measurementList
# nMeasures=len(measurements)
# mean_intensity = np.mean(data_lam, axis=1)

# A[np.isnan(A)] = 0

signal_threshold = 6*(10**6)
saturated = [1 if chan>signal_threshold else 0 for chan in mean_intensity]

   



# generate_report(snirf_obj, savePath = None)
#%%
fig, ax = plt.subplots(1,1)
aux_dict = {}
height_list = []
ax.plot(time, aux)
for i in range(len(aux_peaks[0])):
    # ax.axvline(time[aux_peaks[0][i]], color'r')
    # 
    
    height = aux[aux_peaks[0][i]]
    
    if str(height) not in aux_dict:
        aux_dict[str(height)] = [aux_peaks[0][i]]
    else:
        aux_dict[str(height)].append(aux_peaks[0][i])
    
colours = ['b', 'r', 'c', 'g', 'y']
keys = aux_dict.keys()

for i,key in enumerate(keys):
    print(key)
    triggers = aux_dict[key]
    for j in range(len(triggers)):
        onset = time[triggers[j]]
        plt.axvline(x=onset, color=colours[i], lw=1)

#%%

import matplotlib.gridspec as gridspec

plt.rcParams.update({'font.size': 10})
fig = plt.figure(figsize=(20,10))
# if title == None:
#     fig.suptitle('Signal Quality Evaluation')
# else:
fig.suptitle('Signal Quality Evaluation')

gs = gridspec.GridSpec(2, 2, height_ratios=[8,3],figure=fig)

ax1 = fig.add_subplot(gs[0,0])
snr(snirf_obj,lam=0, ax=ax1)

ax2 = fig.add_subplot(gs[0,1])
snr(snirf_obj,lam=1, ax=ax2)

# ax3 = fig.add_subplot(gs[0,2])
# sci_psp(snirf_obj,mode='montage', ax=ax3)

# ax4 = fig.add_subplot(gs[1,:])
# sci_psp(snirf_obj,mode='thresholded', ax=ax4)

ax3 = fig.add_subplot(gs[1,:])
GVTD(snirf_obj, ax=ax3)

plt.tight_layout()
plt.savefig(savePath)

# if savePath != None:
#      plt.savefig(savePath)




#%% test group level report

root_path = '/Users/lauracarlton/fNIRS/Data/MovieWatching/'
sess=1
task = 'MW'
excluded = [1, 14]
groupLevelReport(root_path, task, sess, excluded=excluded, saveFig=1)


#%% test with MNE


raw_intensity = mne.io.read_raw_snirf(fileName)

raw_intensity.load_data()
raw_od = optical_density(raw_intensity)

picks = _validate_nirs_info(raw_od.info, fnirs='od', which='Scalp coupling index')
raw = raw_od.copy().pick(picks).load_data()
raw_copy = raw.copy()

raw_haemo = mne.preprocessing.nirs.beer_lambert_law(raw)


sci = mne.preprocessing.nirs.scalp_coupling_index(raw)
fig, ax = plt.subplots()
ax.hist(sci)
ax.set(xlabel='Scalp Coupling Index', ylabel='Count', xlim=[0, 1])
plt.show()

# nBadChannels_raw = len(list(compress(raw.ch_names, sci < 0.7)))
# raw.info['bads'] = list(compress(raw.ch_names, sci < 0.7))

_, scores_sci, times = scalp_coupling_index_windowed(raw_copy, time_window=5)
plot_timechannel_quality_metric(raw_haemo, scores_sci, times, threshold=0.7, title="Scalp Coupling Index Quality Evaluation")

