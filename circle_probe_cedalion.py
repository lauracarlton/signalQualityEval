#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
montage plotting for cedalion 

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
import cedalion 
import cedalion.nirs
import cedalion.xrutils as xrutils

#%%
def plot_circle_probe(snirfObj, metric, ax, colormap=plt.cm.bwr, title=None, threshold_ind = None, threshold_col = None, saturation=None, vmin=0, vmax=1, savePath = None, remove_short=0):
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
    geo3d = snirfObj.geo3d
    
    landmarks = geo3d.loc[geo3d.type == cedalion.io.snirf.PointType.LANDMARK]
    sources = geo3d.loc[geo3d.type == cedalion.io.snirf.PointType.SOURCE]
    detectors = geo3d.loc[geo3d.type == cedalion.io.snirf.PointType.DETECTOR]
    
    
    #### find the landmarks in the probe ####
    for u in range(len(landmarks)):
        idx_list = channels_df.index[channels_df['Label']==landmarks.label[u]].tolist()
        if idx_list:
            circular_landmark_pos3D.append([channels_df['X'][idx_list[0]],channels_df['Y'][idx_list[0]], channels_df['Z'][idx_list[0]]])
            landmark_pos3D = landmarks[u,0:3].to_numpy().tolist()
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
        
    # measurements = snirfObj.nirs[0].data[0].measurementList
    skipped_channels = []
    skipped_detectors = []
    skipped_metrics = []
    data = snirfObj.data[0]
    nMeas = len(data.channel)
    
    if remove_short == 1: # then remove any channels that are less than 10mm 
        
    
        for u in range(nMeas):
            
            sourceIndex =  data.source[u]
            detectorIndex =  data.detector[u]
            
            dist = xrutils.norm(geo3d.loc[data.source[u]] - geo3d.loc[data.detector[u]], dim="pos")


            # x = snirfObj.nirs[0].probe.sourcePos3D[sourceIndex-1]
            # y= snirfObj.nirs[0].probe.detectorPos3D[detectorIndex-1]
            # dist = math.dist(x,y)
                
            if dist < 10:
                    skipped_channels.append([sourceIndex, detectorIndex])
                    skipped_detectors.append(detectorIndex)
                    skipped_metrics.append(u)
    
    # if the metrics/threshold_col given include those for short channels, remove them from the array 
    if len(metric) == nMeas//2:
        metric = np.delete(metric,skipped_metrics)
    
    if type(threshold_col) == list:
        if len(threshold_col) == nMeas//2:
            threshold_col = np.delete(threshold_col,skipped_metrics)

    #### scale indices #####
    sourcePos2DX , sourcePos2DY = convert_optodePos3D_to_circular2D(sources, tranformation_matrix, norm_factor)
    detectorPos2DX , detectorPos2DY = convert_optodePos3D_to_circular2D(detectors, tranformation_matrix, norm_factor)
    
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
    for u in range(nMeas):
        sourceIndex =  data.source[u]
        detectorIndex =  data.detector[u]
        
        # skip the short_channels 
        if [sourceIndex, detectorIndex] in skipped_channels:
            continue
        
        
        iS = int(sourceIndex.to_numpy().tolist()[1:])
        iD = int(detectorIndex.to_numpy().tolist()[1:])
        x = [sourcePos2DX[iS-1], detectorPos2DX[iD-1]]
        y = [sourcePos2DY[iS-1], detectorPos2DY[iD-1]]
        

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
    else:   
        ticks = [vmin, (vmin+vmax)//2, vmax]
        
    ax.plot(0, 1 , marker="^",markersize=16)
    plt.colorbar(sm,shrink =0.6, ticks=ticks)
    ax.set_title(title)
    plt.tight_layout()
    plt.axis('equal')
    plt.axis('off')
    
    if savePath is not None: 
        plt.savefig(savePath, dpi=1200)
    






