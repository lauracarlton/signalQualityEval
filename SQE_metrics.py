#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GENERATE SQE PLOTS
@author: lauracarlton
"""
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats
import matplotlib.colors as clrs
from scipy import signal 
import seaborn as sns

#%% Data preprocessing

def sort_Data(snirf_obj, interpolate=True):
    
    '''
    Sort the data from the snirf object
    turn nMeasurements x nTimePoints to nChannels x nTimePoints x nWavelengths
    returns channels in order of measurementList
    '''
    
    nirs_obj_data = snirf_obj.nirs[0].data
    data = nirs_obj_data[0].dataTimeSeries

    time = nirs_obj_data[0].time
    measurements = nirs_obj_data[0].measurementList
    channels = []
    nMeasures = len(measurements)
    nChannels = nMeasures//2
    nSamples = len(time)
    nWLs = 2
    
    ### interpolate NaNs
    if interpolate == True:
        for iChan in range(nMeasures):
            chanData = data[:,iChan]
            x = np.arange(len(chanData))
            xp = np.arange(len(chanData))[np.isnan(chanData) == False]
            fp = chanData[np.isnan(chanData) == False]
            data[:,iChan] = np.interp(x, xp, fp)
            
    ### extract the channel names based on the measurement list ###
    for u in range(nMeasures):
        sourceIndex =  measurements[u].sourceIndex
        detectorIndex = measurements[u].detectorIndex
        name_temp = 'S'+str(sourceIndex)+'_D'+str(detectorIndex)
        if name_temp not in channels:
            channels.append(name_temp)
        
    #### sort channels and data to make sure they are in the correct order ###
    sortedData = np.zeros([nChannels, nSamples, nWLs])
    i=0
    for u in range(nChannels):
        channel1 = [measurements[u].sourceIndex, measurements[u].detectorIndex]
        wl1 = measurements[u].wavelengthIndex
        for j in range(nChannels,nMeasures):
            channel2 = [measurements[j].sourceIndex, measurements[j].detectorIndex]
            wl2 = measurements[j].wavelengthIndex
            if channel1 == channel2 and wl1 != wl2:
                sortedData[i,:,0] = data[:,u] 
                sortedData[i,:,1] = data[:,j]
                i+=1
    return sortedData, channels

def filter_data(sortedData, time, nChannels, fcut_min = 0.5, fcut_max = 2.5):
    
    '''
    use a butterworth bandpass filter between 0.5 and 2.4
    normalized frequency is fcut/(2/fs)
    normalize by the cardiac data 
    -> divide by standard deviation of the filtered signal per channel
    if there is an NaN at the end of the data, replaces it with average of the last 50 points
    '''

    T = time[1]-time[0]
    fs = 1/T
    b,a = signal.butter(4,[fcut_min*(2/fs), fcut_max*(2/fs)], btype='bandpass')
    
    cardiac_data = np.zeros(np.shape(sortedData))
    
    # for c in range(nChannels):
    #     chan_data = sortedData[c,:,0][np.logical_not(np.isnan(sortedData[c,:,0]))]
    #     mean_intensity = np.mean(chan_data[-50:-1])
    #     sortedData[c,:,0][np.isnan(sortedData[c,:,0])] = mean_intensity
        
    #     chan_data = sortedData[c,:,1][np.logical_not(np.isnan(sortedData[c,:,1]))]
    #     mean_intensity = np.mean(chan_data[-50:-1])
    #     sortedData[c,:,1][np.isnan(sortedData[c,:,1])] = mean_intensity
    
    
    cardiac_data[:,:,0] = signal.filtfilt(b,a,sortedData[:,:,0])
    cardiac_data[:,:,0] = np.array([cardiac_data[i,:,0]/np.std(cardiac_data[i,:,0]) for i in range(nChannels)])
    
    cardiac_data[:,:,1] = signal.filtfilt(b,a,sortedData[:,:,1])
    cardiac_data[:,:,1] = np.array([cardiac_data[i,:,1]/np.std(cardiac_data[i,:,1]) for i in range(nChannels)])
    
    return cardiac_data

#%% Plot 2D montage 

def SQE_2Dplot_func(snirfObj, metric, ax, colormap=plt.cm.jet, title='DQR', threshold_ind = None, threshold_col = None, saturation=None, vmin=0, vmax=1, savePath = None, remove_short=0):
    '''
    CREATE A 2D MONTAGE OF OPTODES WITH CHANNELS COLOURED ACCORDING TO A GIVEN METRIC
    
    Parameters:
        snirfObj -> valid snirfObj to extract the measurement list
        ax -> axis object to plot the montage
        metric -> metric to plot onto the channels
        colormap -> colormap to use to color the channels (default = jet)
        title -> title for the plot (default = SQE)
        threshold -> metrics values below threshold will be plotted as dotted lines (default = 0)
        vmin -> minimum value for the colorbar (default = 0)
        vmax -> maximum value for the colorbar (default = 1)
    '''

    def cart2sph(x, y, z):
        hxy = np.hypot(x, y)
        r = np.hypot(hxy, z)
        el = np.arctan2(z, hxy)
        az = np.arctan2(y, x)
        return az, el, r
    
    def pol2cart(theta, rho):
        x = rho * np.cos(theta)
        y = rho * np.sin(theta)
        return x, y
    
    def convert_optodePos3D_to_circular2D(pos, tranformation_matrix, norm_factor):
        pos = np.append(pos, np.ones((pos.shape[0],1)), axis=1)
        pos_sphere = np.matmul(pos,tranformation_matrix)
        pos_sphere_norm = np.sqrt(np.sum(np.square(pos_sphere), axis=1))
        pos_sphere_norm= pos_sphere_norm.reshape(-1,1)
        pos_sphere = np.divide(pos_sphere,pos_sphere_norm)
        azimuth, elevation, r = cart2sph(pos_sphere[:,0], pos_sphere[:,1], pos_sphere[:,2])
        elevation = math.pi/2-elevation;
        x, y = pol2cart(azimuth,elevation)
        x = x/norm_factor
        y = y/norm_factor
        return x, y
    
    channels_df = pd.read_excel('10-5-System_Mastoids_EGI129.xlsx') 
    probe_landmark_pos3D = []
    circular_landmark_pos3D = []
    
    #### find the landmarks in the probe ####
    for u in range(len(snirfObj.nirs[0].probe.landmarkLabels)):
        idx_list = channels_df.index[channels_df['Label']==snirfObj.nirs[0].probe.landmarkLabels[u]].tolist()
        if idx_list:
            circular_landmark_pos3D.append([channels_df['X'][idx_list[0]],channels_df['Y'][idx_list[0]], channels_df['Z'][idx_list[0]]])
            landmark_pos3D = snirfObj.nirs[0].probe.landmarkPos3D[u,0:3].tolist()
            landmark_pos3D.extend([1])
            probe_landmark_pos3D.append(landmark_pos3D)
            
    landmarks2D_theta = (channels_df['Theta']*2*math.pi/360).to_numpy()
    landmarks2D_phi = ((90-channels_df['Phi'])*2*math.pi/360).to_numpy()
    x,y = pol2cart(landmarks2D_theta, landmarks2D_phi)
    
    norm_factor = max(np.sqrt(np.add(np.square(x),np.square(y))))
    temp = np.linalg.inv(np.matmul(np.transpose(probe_landmark_pos3D),probe_landmark_pos3D))
    tranformation_matrix = np.matmul(temp,np.matmul(np.transpose(probe_landmark_pos3D),circular_landmark_pos3D))        
    tranformation_matrix = tranformation_matrix
    # tranformation_matrix = np.linalg.lstsq(probe_landmark_pos3D, circular_landmark_pos3D, rcond=None)
        
    measurements = snirfObj.nirs[0].data[0].measurementList
    skipped_channels = []
    skipped_detectors = []
    skipped_metrics = []
    
    if remove_short == 1: # then remove any channels that are less than 10mm 
        
        for u in range(len(measurements)//2):
            sourceIndex =  measurements[u].sourceIndex
            detectorIndex =  measurements[u].detectorIndex
            x = snirfObj.nirs[0].probe.sourcePos3D[sourceIndex-1]
            y= snirfObj.nirs[0].probe.detectorPos3D[detectorIndex-1]
            dist = math.dist(x,y)
                
            if dist < 10:
                    skipped_channels.append([sourceIndex, detectorIndex])
                    skipped_detectors.append(detectorIndex)
                    skipped_metrics.append(u)
    
    # if the metrics/threshold_col given include those for short channels, remove them from the array 
    if len(metric) == len(measurements)//2:
        metric = np.delete(metric,skipped_metrics)
    
    if type(threshold_col) == list:
        if len(threshold_col) == len(measurements)//2:
            threshold_col = np.delete(threshold_col,skipped_metrics)

    #### scale indices #####
    sourcePos2DX , sourcePos2DY = convert_optodePos3D_to_circular2D(snirfObj.nirs[0].probe.sourcePos3D, tranformation_matrix, norm_factor)
    detectorPos2DX , detectorPos2DY = convert_optodePos3D_to_circular2D(snirfObj.nirs[0].probe.detectorPos3D, tranformation_matrix, norm_factor)
    
    # sourcePos2D = snirfObj.nirs[0].probe.sourcePos2D
    # sourcePos2DX = sourcePos2D[:,0]
    # sourcePos2DY = sourcePos2D[:,1]
    # detectorPos2D = snirfObj.nirs[0].probe.detectorPos2D
    # detectorPos2DX = detectorPos2D[:,0]
    # detectorPos2DY = detectorPos2D[:,1]
    
    
    scale = 1.3
    sourcePos2DX = sourcePos2DX*scale
    detectorPos2DX = detectorPos2DX*scale
    sourcePos2DY = sourcePos2DY*scale
    detectorPos2DY = detectorPos2DY*scale
        
    #### plot the positions on the unit circle ####
    t = np.linspace(0, 2 * np.pi, 100)
    head_x = [math.sin(i) for i in t]
    head_y = [math.cos(i) for i in t]
        
    
    #### plot the channels according to the metrics ####
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = matplotlib.cm.ScalarMappable(cmap=colormap,norm=norm)
    fontDict_src = dict(color='r', fontweight = 'bold', fontstretch= 'expanded',fontsize = 7)
    fontDict_det = dict(color='b', fontweight = 'bold', fontstretch= 'expanded',fontsize = 7)
    
    i =0
    for u in range(len(measurements)//2):
        sourceIndex =  measurements[u].sourceIndex
        detectorIndex =  measurements[u].detectorIndex
        
        # skip the short_channels 
        if [sourceIndex, detectorIndex] in skipped_channels:
            continue
        
        
        x = [sourcePos2DX[sourceIndex-1], detectorPos2DX[detectorIndex-1]]
        y = [sourcePos2DY[sourceIndex-1], detectorPos2DY[detectorIndex-1]]
        

        try:
            assert(threshold_col == None)
        except:
            if threshold_col[i] == 1: #metric[u] < threshold: 
                linestyle = '-'
                alpha = 0.4
            else:
                linestyle = '-'
                alpha = 1
        else:
            linestyle = '-'
            alpha=1
        
        try:
            assert(saturation == None)
        except:
            if saturation[i] == 1:
                color = '0.7'
                alpha = 1
            else:
                color = colormap(norm(metric[i]))
        else:
            color = colormap(norm(metric[i]))
            
        ax.plot(x,y, color=color,linestyle=linestyle, linewidth = 2, alpha=alpha)
        # ax.text(sourcePos2DX[sourceIndex-1], sourcePos2DY[sourceIndex-1],str(sourceIndex),fontdict=fontDict_src) # bbox=dict(color = 'r',boxstyle = "round, pad=0.3", alpha=0.05))
        # ax.text(detectorPos2DX[detectorIndex-1], detectorPos2DY[detectorIndex-1], str(detectorIndex),fontdict=fontDict_det) # bbox=dict(color='b',boxstyle = "round, pad=0.3", alpha=0.05))
        i+=1
    
    ax.plot(head_x,head_y,'k')
    for u in range(len(sourcePos2DX)):
        ax.plot(sourcePos2DX[u] , sourcePos2DY[u], 'r.', markersize=8)
        
    for u in range(len(detectorPos2DX)):
        if u+1 in skipped_detectors:
            continue
        ax.plot(detectorPos2DX[u] , detectorPos2DY[u], 'b.',markersize=8)
     
    if threshold_ind != None:
        ticks = [vmin, (vmin+vmax)//2, threshold_ind, vmax]
        ticks = np.squeeze(ticks)

    else:   
        ticks = [vmin, (vmin+vmax)//2, vmax]
        ticks = np.squeeze(ticks)
        
    ax.plot(0, 1 , marker="^",markersize=16)
    # plt.colorbar(sm, shrink =0.6, ticks=ticks)
    ax.set_title(title)
    plt.tight_layout()
    plt.axis('equal')
    plt.axis('off')
    
    if savePath is not None: 
        plt.savefig(savePath, dpi=1200)
    


#%% SNR

def snr(snirf_obj,lam, ax, threshold=5, savePath=None, remove_short=0):
    '''
    Signal-to-noise ratio 
    lam -> determines which wavelength to evaluate - 0 or 1
    ax -> axes object to plot the figure
    threhold -> metrics below threshold plotted as dotted line 
    '''
    
    wavelength = int(snirf_obj.nirs[0].probe.wavelengths[lam])

    sortedData, channelNames = sort_Data(snirf_obj)
    data_lam = sortedData[:,:,lam]
    nChans = np.shape(data_lam)[0]
    mean = np.zeros(nChans)
    std = np.zeros(nChans)
    
    ### calculate the mean and standard deviation - exclude NaNs ###
    for c in range(nChans):
        chan_data = data_lam[c,:][np.logical_not(np.isnan(data_lam[c,:]))]
        mean[c] = np.mean(chan_data)
        std[c] = np.std(chan_data)

    
    ### take the ratio between the mean and standard deviation for each channel ###
    SNR = mean/std
    SNR[np.isnan(SNR)] = 0 
    SNR[np.isneginf(SNR)] = 0 
    SNR[np.isinf(SNR)] = 0 

    if ax == None:
        return SNR
    
    saturation = saturationCheck(snirf_obj, lam=lam)
    
    threshold_col = [1 if SNR[u] < threshold else 0 for u in range(len(SNR))]

    colors=["red","yellow","green"]
    nodes = [0.0, 0.5, 1.0]
    cmap = clrs.LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))        

    if savePath != None: 
        savePath_snr = savePath+'SNR_' + str(wavelength) + '_montage.jpeg'
        SQE_2Dplot_func(snirf_obj, ax=ax, metric=SNR, colormap=cmap, vmin=0, vmax=25, threshold_col = threshold_col, threshold_ind = threshold, saturation=saturation, savePath=savePath_snr, title='SNR: $\lambda$ = ' + str(wavelength), remove_short=remove_short)
    else: 
        SQE_2Dplot_func(snirf_obj, ax=ax, metric=SNR, colormap=cmap, vmin=0, vmax=25, threshold_col = threshold_col, threshold_ind = threshold, saturation=saturation, savePath=savePath, title='SNR: $\lambda$ = ' + str(wavelength), remove_short=remove_short)


#%% SCI x PSP

def sci_channel(snirf_obj, ax, savePath=None, remove_short=0):
    '''
    generate histogram of sci per channel
    '''
    nirs_obj_data = snirf_obj.nirs[0].data
    time = nirs_obj_data[0].time
    measurements = nirs_obj_data[0].measurementList
    nMeasures = len(measurements)
    nChannels = nMeasures//2
    
    sortedData, channels = sort_Data(snirf_obj, measurements, nChannels)
    cardiac_data = filter_data(sortedData, time, nChannels)
    
    sci_channel = np.zeros(nChannels)
    
    for ii in range(nChannels):
        c = np.corrcoef(cardiac_data[:,ii,0], cardiac_data[:,ii,1])[0][1]
        if not np.isfinite(c):  
            c = 0
        sci_channel[ii] = c
    
    ax.hist(sci_channel)
    if savePath != None: 
        SQE_2Dplot_func(snirf_obj,sci_channel,ax, savePath=savePath, title='WHole Channel SCI', remove_short=remove_short)
        plt.savefig(savePath+'sci_hist.png')
        

def sci_psp(snirf_obj, mode, ax, window = 5, threshold_sci=0.7,threshold_psp=0.1, savePath=None, remove_short=0):
    '''
    calculate windowed scalp coupling index and peak spectral power 
    Parameters:
        mode -> possible figures to generate, options: heatmap, threholded, montage
        ax -> axes object to plot figure 
        window -> window length in samples 
        threshold_sci -> threshold for sci (default = 0.7)
        threhold_psp -> threshold for psp (default = 0.1)
    '''
    nirs_obj_data = snirf_obj.nirs[0].data
    time = nirs_obj_data[0].time
    stim = snirf_obj.nirs[0].stim    
    T = time[1]-time[0]
    fs = 1/T
    nSamples = len(time)
    measurements = nirs_obj_data[0].measurementList
    nMeasures = len(measurements)
    nChannels = nMeasures//2
    
    window_samples = int(np.ceil(window *fs))
    n_windows = int(np.floor(nSamples/window_samples))
    
    sortedData, channels = sort_Data(snirf_obj)
    # cardiac_data = sortedData
    try:
        cardiac_data = filter_data(sortedData, time, nChannels)
    except:
        cardiac_data = np.zeros(sortedData.shape)
        
    cardiac_data[np.isnan(cardiac_data)] = 0
    
    sci = np.zeros([nChannels, n_windows])
    psp = np.zeros([nChannels, n_windows])
    times = []
    
    for w in range(n_windows):    
        
        start_sample = int(w * window_samples)
        end_sample = start_sample + window_samples
        end_sample = np.min([end_sample, nSamples - 1])

        t_start = time[start_sample]
        t_stop = time[end_sample]
        times.append((t_start, t_stop))
        
        for ii in range(0, nChannels):
            c1 = cardiac_data[ii, start_sample:end_sample, 0]
            # c1[np.isnan(c1)] = 0
            c2 = cardiac_data[ii, start_sample:end_sample, 1]
            # c2[np.isnan(c2)] = 0

            similarity = signal.correlate(c1, c2, 'full')
            # similarity = window_samples*similarity/np.sqrt(np.sum(abs(c1)**2 * np.sum(abs(c2)**2)))


            sci[ii, w] = stats.pearsonr(c1,c2)[0] #similarity[lags==0] #s

            [f, pxx] = signal.periodogram(similarity, fs=fs, window=np.hamming(len(similarity)), nfft=len(similarity), scaling='spectrum')
            
            psp[ii, w] = max(pxx)
            
    # sci  = (sci - np.nanmin(sci))/(np.nanmax(sci) - np.nanmin(sci))
    # psp  = (psp - np.nanmin(psp))/(np.nanmax(psp) - np.nanmin(psp))
    # sci[np.isnan(sci)] = 0
    # sci[np.isinf(sci)] = 0
    # sci[np.isneginf(sci)] = 0
    # psp[np.isnan(psp)] = 0
    # psp[np.isinf(psp)] = 0
    # psp[np.isneginf(psp)] = 0
    
    cols = [np.round(t[0]) for t in times]
    sci_to_plot = pd.DataFrame(data=sci,columns=pd.Index(cols, name='Time (s)'), index=pd.Index(channels, name='Channel'))
    psp_to_plot = pd.DataFrame(data=psp,columns=pd.Index(cols, name='Time (s)'), index=pd.Index(channels, name='Channel'))
    
    sci_binary = np.ones(np.shape(sci))
    psp_binary = np.ones(np.shape(psp))
    scixpsp_binary = np.ones(np.shape(psp))
    chan_percent = np.zeros(nChannels)

    for chan in range(nChannels):
        
        for w in range(n_windows):
            if sci[chan,w] < threshold_sci:
                sci_binary[chan, w] = 0
            
            if psp[chan,w] < threshold_psp:
                psp_binary[chan,w] = 0
            
            if sci[chan,w] < threshold_sci or psp[chan,w] < threshold_psp:
                scixpsp_binary[chan,w] = 0
        
        
        chan_percent[chan] = sum(scixpsp_binary[chan,:])/n_windows
        
    if ax == None:
        return sci, psp, chan_percent
        
    if mode=='heatmaps':                              

        sns.heatmap(data=sci_to_plot,vmin=0,vmax=1, cmap='Reds_r',
                     cbar_kws=dict(label='Score'), ax=ax[0])
        ax[0].set_title('Scalp Coupling Index')
        ax[0].set_yticks(ticks = np.arange(nChannels),labels=channels, fontsize=8)
        sns.heatmap(data=psp_to_plot,vmin=0,vmax=1, cmap='Reds_r',
                     cbar_kws=dict(label='Score'), ax=ax[1])
        ax[1].set_title('Peak Spectral Power')
        ax[1].set_yticks(ticks = np.arange(nChannels),labels=channels, fontsize=8)
        plt.tight_layout()
        
        if savePath != None: 
            plt.savefig(savePath+'scixpsp_heatmap.png')
        return 
    
            
    elif mode=='thresholded':
 
        # sci_b_to_plot = pd.DataFrame(data=sci_binary,columns=pd.Index(cols, name='Time (s)'), index=pd.Index(channels, name='Channel'))
        # psp_b_to_plot = pd.DataFrame(data=psp_binary,columns=pd.Index(cols, name='Time (s)'), index=pd.Index(channels, name='Channel'))
                                        
        # fig,ax = plt.subplots(3,1,figsize=(20,8))
        # sns.heatmap(data=sci_b_to_plot, cmap='Reds_r',
        #              cbar_kws=dict(label='Score',ticks=[0,1]), ax=ax[0])
        # ax[0].set_title('Scalp Coupling Index')
        # ax[0].set_yticks(ticks = np.arange(nChannels),labels=channels, fontsize=8)
        # sns.heatmap(data=psp_b_to_plot, cmap='Reds_r',
        #              cbar_kws=dict(label='Score',ticks=[0,1]), ax=ax[1],vmin=0,vmax=1)
        # ax[1].set_title('Peak Spectral Power')
        # ax[1].set_yticks(ticks = np.arange(nChannels),labels=channels, fontsize=8)
        
        scixpsp_b_to_plot = pd.DataFrame(data=scixpsp_binary,columns=pd.Index(cols, name='Time (s)'), index=pd.Index(channels, name='Channel'))

        colors=["black","white"]
        # nodes = [0,1]
        cmap = clrs.ListedColormap(colors)
        times = np.linspace(0,time[-1], num=n_windows)
        c = ax.pcolor(times,np.arange(nChannels), scixpsp_b_to_plot, cmap=cmap)
        ax.set_xlim([0,time[-1]])
        ax.set_yticks([])
        ax.set_ylabel('Channels')
        plt.colorbar(c, location='bottom', ticks = [0,1], shrink = 0.1, pad = 0.2)
        ax.set_title('SCI X PSP')
        # ax.set_yticks(ticks = np.arange(nChannels),labels=channels, fontsize=8)
        
        ### add markers for events
        colours = ['b', 'r', 'g', 'y', 'c']

        for i,stimChan in enumerate(stim):
            stimData = stimChan.data
            nStims = np.shape(stimData)[0]
            for n in range(nStims):
                onset = stimData[n,0]
                plt.axvline(x=onset, color=colours[i], lw=1)

        plt.tight_layout()
        if savePath != None: 
            plt.savefig(savePath+'scixpsp_thresholded.png')
        return
    
    elif mode == 'montage':
        threshold = 0.6
        threshold_col = [1 if chan_percent[u] < threshold else 0 for u in range(len(chan_percent))]
        
        saturation1 = saturationCheck(snirf_obj, lam=0, signal_threshold = 6*(10**6))
        saturation2= saturationCheck(snirf_obj, lam=1, signal_threshold = 6*(10**6))
        saturation = []
        
        for i in range(len(saturation1)):
            if saturation1[i] == 1 or saturation2[i] == 1:
                saturation.append(1)
            else:
                saturation.append(0)

        colors=["red","yellow","green"]
        nodes = [0.0, 0.5, 1.0]
        cmap = clrs.LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))        
    
        if savePath != None: 
            plt.savefig(savePath+'sci_psp_thresh.png')
            SQE_2Dplot_func(snirf_obj, metric=chan_percent, ax=ax, colormap = cmap, threshold_col = threshold_col, threshold_ind = threshold, saturation=saturation, savePath=savePath+'sci_psp_percent.png', title='% of time below thresholds', remove_short=remove_short)
        else :   
            SQE_2Dplot_func(snirf_obj, metric=chan_percent, ax=ax, colormap = cmap, threshold_col = threshold_col, threshold_ind = threshold, saturation=saturation, savePath=savePath, title='% of time above thresholds', remove_short=remove_short)
        return
    
    
    
#%% GVTD 

def GVTD(snirf_obj, ax, savePath=None):
    
    nirs_obj_data = snirf_obj.nirs[0].data
    

    data = nirs_obj_data[0].dataTimeSeries

    # convert to OD
    OD = - np.log( data / np.mean(data, axis=0) )
    OD[np.isnan(OD)] = 0
    OD[np.isinf(OD)] = 0
    
    # filter OD
    fny = (1/(snirf_obj.nirs[0].data[0].time[1] - snirf_obj.nirs[0].data[0].time[0])) / 2
    fmin = 0.01/fny
    fmax = 0.5/fny
    b, a = signal.butter(4, (fmin, fmax), "bandpass")
    OD_filt = signal.filtfilt(b,a,OD.T)
    
    time = nirs_obj_data[0].time
    stim = snirf_obj.nirs[0].stim

    # Step 1: Find the matrix of the temporal derivatives
    dataDiff = OD_filt - np.roll(OD_filt, shift=(0, -1), axis=(0, 1))
    
    # Step 2: Find the RMS across the channels for each time-point of dataDiff
    gvtdTimeTrace = np.sqrt(np.mean(dataDiff[:, :dataDiff.shape[1]-1]**2, axis=0))
    
    # Step 3: Add a zero in the beginning for GVTD to have the same number of time-points as your original dataMatrix
    gvtdTimeTrace = np.concatenate(([0], gvtdTimeTrace))
    
    if ax == None:
        return gvtdTimeTrace
    
    ax.plot(time,gvtdTimeTrace,'k')
    # ax.set_ylim([0,1])
    ax.set(xlabel='Time (s)',ylabel='GVTD',title='GVTD')
    # ax.set_ylim([min(GVTD_z), 1])
    
    ### add markers for events
    colours = ['b', 'r', 'g', 'y', 'c']
    
    for i,stimChan in enumerate(stim):
        stimData = stimChan.data
        nStims = np.shape(stimData)[0]
        for n in range(nStims):
            onset = stimData[n,0]
            plt.axvline(x=onset, color=colours[i], lw=1)
            
    if savePath != None:
        plt.savefig(savePath+'GVTD_timecourse.png')
        
        
#%% check saturation 

def saturationCheck(snirf_obj, lam, signal_threshold=6*(10**6)):
    
    sortedData, channelNames = sort_Data(snirf_obj)
    data_lam = sortedData[:,:,lam]
    mean_intensity = []
    nChans = np.shape(data_lam)[0]
    
    for c in range(nChans):
        chan_data = data_lam[c,:][np.logical_not(np.isnan(data_lam[c,:]))]
        mean_intensity.append(np.mean(chan_data))    
    
    saturated = [1 if chan>signal_threshold else 0 for chan in mean_intensity]
   
    return saturated

#%%  generate large plots 
def generate_figures(snirf_obj):
    '''
    generate individual figures
    '''
    snr(snirf_obj)
    sci_channel(snirf_obj)
    sci_psp(snirf_obj)
    GVTD(snirf_obj)

def generate_report(snirf_obj, savePath = None, title=None, remove_short=0):
    import matplotlib.gridspec as gridspec
    import os
    
    plt.rcParams.update({'font.size': 10})
    fig = plt.figure(figsize=(20,10))
    if title == None:
        fig.suptitle('Signal Quality Evaluation')
    else:
        fig.suptitle(title)
    
    gs = gridspec.GridSpec(3, 3, height_ratios=[8,2,1],figure=fig)
    
    ax1 = fig.add_subplot(gs[0,0])
    snr(snirf_obj,lam=0, ax=ax1,remove_short=remove_short)
    
    ax2 = fig.add_subplot(gs[0,1])
    snr(snirf_obj,lam=1, ax=ax2, remove_short=remove_short)

    ax3 = fig.add_subplot(gs[0,2])
    sci_psp(snirf_obj,mode='montage', ax=ax3, remove_short=remove_short)

    ax4 = fig.add_subplot(gs[1,:])
    sci_psp(snirf_obj,mode='thresholded', ax=ax4, remove_short=remove_short)
    
    ax5 = fig.add_subplot(gs[2,:], sharex=ax4)
    GVTD(snirf_obj, ax=ax5)
    
    plt.tight_layout()
    if savePath != None:
         dirname = os.path.dirname(savePath)
         if not os.path.exists(dirname):
             os.makedirs(dirname)
         plt.savefig(savePath)
    
