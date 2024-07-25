#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 10:34:10 2024

@author: lauracarlton
"""

from SQE_metrics import generate_report, sci_psp
from snirf import Snirf
import pdb 
import matplotlib.pyplot as plt 

#%%
snirfPath = '/Users/lauracarlton/Documents/DATA/MAFC_raw/sub-01/sub-01_task-MA_run-01_nirs_wStim.snirf'
#%%
snirf_obj = Snirf(snirfPath)
#%%

# pdb.set_trace()
generate_report(snirf_obj)
#%%


pdb.set_trace()
sci, psp, chan_perc = sci_psp(snirf_obj, mode='montage', ax=None)

#%%

from cedalion.sigproc.quality import sci, psp 
import cedalion.io as io 
from cedalion import Quantity, units
from cedalion.sigproc.data_quality_report import plot_metrics_on_probe
import matplotlib.pyplot as plt

elements = io.read_snirf(snirfPath)
#%%
import xarray as xr
SCI, SCI_mask = sci(elements[0].data[0], 5*units.s, 0.7)
PSP, PSP_mask = psp(elements[0].data[0], 5*units.s, 0.7)

scixpsp_mask = xr.where(SCI_mask & PSP_mask == False, SCI_mask, PSP_mask)

perc_channel = 1 - scixpsp_mask.sum('time')/scixpsp_mask.shape[1]
#%%
fig, ax = plt.subplots(1,1)
plot_metrics_on_probe(elements[0], SCI.mean('time'), ax=ax, colormap=plt.cm.jet)

fig, ax = plt.subplots(1,1)
plot_metrics_on_probe(elements[0], PSP.mean('time'), ax=ax, colormap=plt.cm.jet, vmax=0.2)

#%%
fig, ax = plt.subplots(1,1)
plot_metrics_on_probe(elements[0], perc_channel, ax=ax, colormap=plt.cm.jet)
