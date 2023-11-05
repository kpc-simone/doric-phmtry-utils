from lmfit import minimize, Parameters, minimize, Parameters, Parameter, report_fit, printfuncs
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

import sys, os
import h5py

### file I/O ###

def traverse_datasets(hdf_file):
    def h5py_dataset_iterator(g, prefix=''):
        for key in g.keys():
            item = g[key]
            path = f'{prefix}/{key}'
            if isinstance(item, h5py.Dataset): # test for dataset
                yield (path, item)
            elif isinstance(item, h5py.Group): # test for group (go down)
                yield from h5py_dataset_iterator(item, path)

    for path, _ in h5py_dataset_iterator(hdf_file):
        yield path

def load_phmtry_raw_doric(filename):
    data = {}
    with h5py.File(filename, 'r') as f:
        for dset in traverse_datasets(f):
            path_info = dset.split('/')
            print(path_info)
            if path_info[1] == 'DataAcquisition':
                signal = '{}/{}'.format(path_info[-2],path_info[-1])
                data[signal] = f[dset][:].ravel()
            elif dset.split('/')[1] == 'Traces':
                signal = path_info[-1]
                data[signal] = f[dset][:].ravel()
            print('loading {} ... '.format(signal) )
            print('Path:', dset)
            print('Shape:', f[dset].shape)
            print('Data type:', f[dset].dtype)
    return data

def load_phmtry_raw_csv(ph_path):
    phtdf_raw = pd.read_csv(ph_path,skiprows=1)
    
    phdf = pd.DataFrame( )
    phdf['time'] = phtdf_raw['Time(s)']
    phdf['F-465'] = phtdf_raw['AIn-1 - Dem (AOut-1)']
    phdf['AF-405'] = phtdf_raw['AIn-2 - Dem (AOut-2)']

    return phdf

### photobleaching correction ###
def exp_func(t,tau,alpha,beta):
    return alpha * np.exp(-t/tau) + beta

def correct_photobleaching(ts,ys,pre_interval,post_interval):
    
    idx_pre_start = np.argmin( [np.abs(pre_interval[0]-t) for t in ts] )
    idx_pre_end = np.argmin( [np.abs(pre_interval[1]-t) for t in ts] )
    idx_post_start = np.argmin( [np.abs(post_interval[0]-t) for t in ts] )
    idx_post_end = np.argmin( [np.abs(post_interval[1]-t) for t in ts] )
    
    ts_sub = np.hstack([ ts[idx_pre_start:idx_pre_end],ts[idx_post_start:idx_post_end] ])
    ys_sub = np.hstack([ ys[idx_pre_start:idx_pre_end],ys[idx_post_start:idx_post_end] ])

    popt,pcov = curve_fit(exp_func,ts_sub,ys_sub,
    
                        # initial guess of fitting parameters
                        p0 = ( 300., 0.1, np.mean(ys_sub) ) ,
                        maxfev = 10000)
    fs = exp_func(ts, *popt)
    ys_pb_corrected = ys - fs
    
    return ys_pb_corrected

### motion artifact correction ###
def compute_residual(ps, data, template_to_align):
    
    a = ps['a']
    b = ps['b']
    
    template_aligned = a * template_to_align + b
    
    return data.flatten() - template_aligned.flatten()

def correct_motion(data, control):
    params = Parameters()
    params.add('a',value = 5., min = 0.1, max = 10.)
    params.add('b',value = 0., min = -1, max = 1)

    result = minimize(compute_residual, params, args=(data, control), method='leastsq')
    print(result.params)
    
    # to recover the aligned control, subtract residuals from channel 1
    control_aligned = data - result.residual.reshape(control.shape)
    signal_corrected = data - control_aligned
    
    return control_aligned, signal_corrected
    
### z-score transformation ###
def transform_to_zscore(ts, ys, baseline_interval):
    idx_bl_start = np.argmin( [np.abs(baseline_interval[0]-t) for t in ts] )
    idx_bl_end = np.argmin( [np.abs(baseline_interval[1]-t) for t in ts] )
    
    bl_mean = np.mean( ys[idx_bl_start:idx_bl_end] )
    bl_std = np.std( ys[idx_bl_start:idx_bl_end] )
    
    ys_zscore = ( ys - bl_mean ) / bl_std
    
    return ys_zscore
    
### resampling ###
def resample_signal(ts,ys,SR_desired):
    ys_func = interp1d(ts,ys,bounds_error=False)
    ts_new = np.arange( 0., round(ts.max(),2), 1/SR_desired )
    ys_new = ys_func(ts_new)
    return ts_new,ys_new

### extracting segments ###
def extract_segments(ts,ys,seg_times,labels,bl_times=(0.,5.)):
    # TODO: check that all seg_times are the same duration

    segs_df = pd.DataFrame( )
    
    #fig,(ax1,ax2) = plt.subplots(2,1,sharex=True)
    #ax1.plot(ts,ys,color='dimgray')
    
    # if bl_times is a nested list, then distinct baseline_intervals are being passed for each segment
    if any(isinstance(i, list) for i in bl_times):
        for seg_time,bl_time,label in zip(seg_times,bl_times,labels):
            # I would like to recycle the transform_to_zscore function
            # I could do this by first passing the all timeseries data to the z-score function
            # and then simply cutting out everything but that associated with the segment interval
            
            t_start = min(bl_time[0],seg_time[0])
            t_end = max(bl_time[1],seg_time[1])
            idx_start = np.argmin( [ np.abs( t_start - t ) for t in ts ] )
            idx_end = np.argmin( [ np.abs( t_end - t ) for t in ts ] )
            
            # step 1: don't align in time
            ts_seg = ts[idx_start:idx_end]
            ys_seg = ys[idx_start:idx_end]
            ys_zsc = transform_to_zscore(ts_seg,ys_seg,baseline_interval=bl_time)
            
            # step 2: cut out intervals outside the desired segment
            idx_start = np.argmin( [ np.abs( seg_time[0] - t ) for t in ts_seg ] )
            if segs_df.shape[1] == 0:
                idx_end = np.argmin( [ np.abs( seg_time[1] - t ) for t in ts_seg ] )
            else:
                idx_end = idx_start + segs_df.shape[0]
            
            ys_zsc = ys_zsc[idx_start:idx_end]
            segs_df[label] = ys_zsc.tolist()
            
            # for debugging
            ts_zsc = ts_seg[idx_start:idx_end]
            # ax1.plot(ts_seg,ys_seg)
            # ax2.plot(ts_zsc,ys_zsc)
        
            # for ax in (ax1,ax2):
                # ax.axvspan(bl_time[0],bl_time[1],color='k',alpha=0.2)
                # ax.axvline(seg_time[0],color='k',linestyle='--')
                # ax.axvline(seg_time[1],color='k',linestyle='--')
            
    else:
        for seg_time,label in zip(seg_times,labels):
            idx_start = np.argmin( [ np.abs( seg_time[0] - t ) for t in ts ] )
            
            if segs_df.shape[1] == 0:
                idx_end = np.argmin( [ np.abs( seg_time[1] - t ) for t in ts ] )
            else:
                idx_end = idx_start + segs_df.shape[0]
                
            ts_seg = ts[idx_start:idx_end] - ts[idx_start]
            ys_seg = ys[idx_start:idx_end]
            
            ys_zsc = transform_to_zscore(ts_seg,ys_seg,baseline_interval=bl_times)
            #ax1.plot(ts[idx_start:idx_end],ys_seg)
            #ax2.plot(ts[idx_start:idx_end],ys_zsc)
            
            segs_df[label] = ys_zsc.tolist()
    #ax2.axhline(0,color='k',linestyle='--')
    #plt.show(block=False)
    
    return segs_df