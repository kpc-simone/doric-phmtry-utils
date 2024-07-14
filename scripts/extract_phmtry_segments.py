import matplotlib.pyplot as plt
from tkinter.filedialog import askdirectory,askopenfilename
import pandas as pd
import numpy as np
import argparse
import sys,os

sys.path.insert(0,'../utils')
from phmtry_tools import resample_signal, extract_segments

parser = argparse.ArgumentParser()

parser.add_argument('--segstart_colname', required = True,
                        help = 'The name of the column defining the start of the segment of interest.')
parser.add_argument('--baseline_end', required = True,
                        help = 'The name of the column to take as the end of the baseline period, for z-score normalization.')
parser.add_argument('--outcome_colname', required = True,
                        help = 'The name of the column to associate with some experimental outcome.')
                        
parser.add_argument('--pre_interval', type = float, nargs = '?', const = -5, default = -5,
                        help = 'Beginning of extraction, relative to segment start. Default is -5 seconds.')
parser.add_argument('--event_dur', type = float, nargs = '?', const = 8., default = 8.,
                        help = 'Duration of each event. Default is 8 seconds.')
parser.add_argument('--post_interval', type = float, nargs = '?', const = 12., default = 12.,
                        help = 'Time to keep, in seconds, after event has concluded. Default is 12 seconds.')
parser.add_argument('--baseline_pre_interval', type = float, nargs = '?', const = -5, default = -5.,
                        help = 'Beginning of baseline period, relative to its end. Default is -5 seconds.')

parser.add_argument('--SR_desired', type = int, nargs = '?', const = 100, default = 100,
                        help = 'The desired final sampling rate of the timeseries data. Default is 20 samples/second.')


### Script parameters ###
#SR_desired = 20

#pre_interval = -5
#segstart_colname = 'TPoI-2'
#event_dur = 8.
#post_interval = 12.

#baseline_pre_interval = -5
#baseline_end = 'TPoI-1'

#outcome_colname = 'observation'

if __name__ == '__main__':
    args = parser.parse_args()
    SR_desired = args.SR_desired
    pre_interval = args.pre_interval
    segstart_colname = args.segstart_colname
    event_dur = args.event_dur
    post_interval = args.post_interval
    baseline_pre_interval = args.baseline_pre_interval
    baseline_end = args.baseline_end
    outcome_colname = args.outcome_colname

    data_dir = askdirectory(title = 'Select folder containing processed photometry recordings, from which to extract all events. ')
    ec_file = askopenfilename(title = 'Select experimental conditions (ec) file, that defines the start and end times of the events. ')
    out_data_dir = askdirectory(title = 'Select directory to save the ouput data. ')

    phmtry_files = [f for f in os.listdir(data_dir) if '.csv' in f]

    outfilename = 'events-alignedto-{}-zscore.csv'.format(segstart_colname)
    ecdf = pd.read_csv(ec_file)

    ecdf = ecdf[ecdf[segstart_colname].notna()]

    ts_out = np.arange(start = pre_interval,stop = event_dur+post_interval,step = 1/SR_desired )
    out_df = pd.DataFrame( data = { 'time' : ts_out } )

    for phmtry_file in phmtry_files:

        animal = phmtry_file.split('-')[0]
        print(animal)
        an_ecdf = ecdf[ecdf['animal'] == animal]
        print(an_ecdf.head(10))
        
        if len(an_ecdf) > 0:
            bl_times = [ [row[baseline_end]+baseline_pre_interval,row[baseline_end]] for r,row in an_ecdf.iterrows() ]
            seg_times = [ [ row[segstart_colname]+pre_interval,row[segstart_colname]+event_dur+post_interval ] for r,row in an_ecdf.iterrows() ]
            
            labels = [ '{}-t{}-{}'.format(row['animal'],row['trial'],row[outcome_colname]) for r,row in an_ecdf.iterrows() ]
            
            phdf = pd.read_csv(os.path.join(data_dir,phmtry_file))
            ts_res,ys_res = resample_signal(ts = phdf['time'], ys = phdf['fluo465-zsc'], SR_desired = SR_desired)
                    
            '''
            Probem: I need to generalize the z-score normalization procedure to have arbitrary baseline intervals
            Solution: Pass a list of baseline interval times as a separate argument to extract_segments
            '''

            segs_df = extract_segments(ts = ts_res, 
                                            ys = ys_res,
                                            seg_times = seg_times,
                                            labels = labels,
                                            bl_times = bl_times,
                                            )
            
            out_df = pd.concat( [out_df,segs_df], axis = 1 )
            
            # plot 1: verify segment extraction
            fig,(ax1,ax2) = plt.subplots(1,2,figsize=(9,3),gridspec_kw={'width_ratios':[2,1]})
            ax1.plot(ts_res,ys_res,color='dimgray')
            for trial,(seg_start,seg_end) in zip(segs_df.columns,seg_times):
                ax1.plot(out_df['time']+seg_start-pre_interval,out_df[ trial ])
                ax1.axvline(seg_start,color='k',linestyle='--')
                ax1.axvline(seg_end,color='k',linestyle='--')
            
            # plot 2: verify z-score normalization
            ax2.plot(out_df['time'],out_df[segs_df.columns])
            ax2.axvline(0.,color='k',linestyle='--')
            ax2.axvline(event_dur,color='k',linestyle='--')
            
            for ax in (ax1,ax2):
                ax.set_xlabel('Time (s)')
                ax.set_ylabel('z-score')
            
            fig.suptitle(animal)
            fig.tight_layout()
            plt.show()

    out_df.to_csv(os.path.join(out_data_dir,outfilename))

    fig,ax = plt.subplots(1,1)

    zsc_mean = out_df.iloc[:,1:].mean(axis=1)
    zsc_std = out_df.iloc[:,1:].std(axis=1)
    ax.plot(out_df['time'],zsc_mean, color='tab:green')
    ax.fill_between(out_df['time'],zsc_mean-zsc_std,zsc_mean+zsc_std,color='tab:green',alpha=0.3)
    ax.set_xlabel('Time (s) from event onset')
    ax.set_ylabel('Z-Score')
    ax.axvline(0.,color='k',linestyle='--')
    ax.axvline(event_dur,color='k',linestyle='--')
    ax.set_xlim(pre_interval,event_dur+post_interval)
    plt.show()