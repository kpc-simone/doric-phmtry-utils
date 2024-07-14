import matplotlib.pyplot as plt
from tkinter.filedialog import askdirectory
import pandas as pd
import argparse
import sys,os

'''
Processes raw photometry files given a directory. 
Performs photobleaching correction, motion-artifact correction, and transformations to z-score. 

'''

parser = argparse.ArgumentParser()
parser.add_argument('--baseline_pre_start', type = float, nargs = '?', const = 100., default = 100.,
                    help = 'Start of the pre-baseline interval for photobleaching correction, relative to the start of the recording. Default is t=100 seconds. ')
parser.add_argument('--baseline_pre_end', type = float, nargs = '?', const = 600., default = 600.,
                    help = 'End of the pre-baseline interval for photobleaching correction, relative to the start of the recording. Default is t=600 seconds. ')
parser.add_argument('--baseline_post_start', type = float, nargs = '?', const = 500., default = 500.,
                    help = 'Start of the pre-baseline interval for photobleaching correction, relative to the end of the recording. Default is t=500 seconds. ')
parser.add_argument('--baseline_post_end', type = float, nargs = '?', const = 0., default = 0.,
                    help = 'End of the pre-baseline interval for photobleaching correction, relative to the end of the recording. Default is t=0 seconds. ')

sys.path.insert(0,'../utils')
from phmtry_tools import load_phmtry_raw_doric, correct_photobleaching, correct_motion, transform_to_zscore

# the organization of doric files can take multiple forms , depending on the particular version of doric studio installed
# we need to homogenize the names of signals so that we can process them
signals_map = {
        'Console_time(s)'       : 'time',
        'AIn-1 - Dem (AOut-1)'  : 'F-465',
        'AIn-2 - Dem (AOut-2)'  : 'AF-405',
        'DI--O-3'               : 'DI/03',
        'DI--O-4'               : 'DI/04',
        
        'AIN01xAOUT01-LockIn/Time'        : 'time',
        'AIN01xAOUT01-LockIn/Values'      : 'F-465',
        'AIN02xAOUT02-LockIn/Values'      : 'AF-405',
        'DigitalIO/DIO01'                 : 'DI/01',
        'DigitalIO/DIO02'                 : 'DI/02',
        'DigitalIO/DIO03'                 : 'DI/03',
        }

if __name__ == '__main__':
    args = parser.parse_args()
    baseline_pre_start = args.baseline_pre_start
    baseline_pre_end = args.baseline_pre_end
    baseline_post_start = args.baseline_post_start
    baseline_post_end = args.baseline_post_end
    plot = True

    data_dir = askdirectory( title = "Select folder containing raw photometry data" )
    phmtry_files = [f for f in os.listdir(data_dir) if '.doric' in f]
    
    out_data_dir = os.path.join(data_dir,'../phmtry-processed')
    if not os.path.exists(out_data_dir):
        os.mkdir(out_data_dir)
    for phmtry_file in phmtry_files:

        animal = phmtry_file.split('_')[0]
        print(animal)

        filepath = os.path.join(os.path.dirname(__file__),data_dir,phmtry_file)

        # load the data from doric file
        data = load_phmtry_raw_doric(filepath)
        phdf = pd.concat([pd.DataFrame(v, columns=[k]) for k, v in data.items()], axis=1).rename( columns = signals_map )
        phdf = phdf[['time','F-465','AF-405']].dropna()
        # phdf = pd.DataFrame( data = data ).rename( columns = signals_map ).dropna()
        print(phdf.head())

        pre_interval_ = (baseline_pre_start,baseline_pre_end)
        post_interval_ = ( phdf['time'].max()-baseline_post_start,phdf['time'].max()-baseline_post_end )
        print(animal,pre_interval_,post_interval_)
        
        # correct for photobleaching
        phdf['fluo465-pbc'] = correct_photobleaching( ts = phdf['time'], ys = phdf['F-465'],
                                                            pre_interval = pre_interval_, 
                                                            post_interval = post_interval_ 
                                                            )
        phdf['fluo405-pbc'] = correct_photobleaching( ts = phdf['time'], ys = phdf['AF-405'],
                                                            pre_interval = pre_interval_,
                                                            post_interval = post_interval_
                                                            )
        
        # correct for motion-artifacts
        phdf['fluo405-maf'], phdf['fluo465-mac'] = correct_motion(phdf['fluo465-pbc'].to_numpy(),phdf['fluo405-pbc'].to_numpy())
        
        # transform to z-score
        phdf['fluo465-zsc'] = transform_to_zscore( ts = phdf['time'], 
                                            ys = phdf['fluo465-mac'], 
                                            baseline_interval = (120,480)
                                            )

        fig,(ax1,ax2,ax3,ax4) = plt.subplots(4,1,sharex=True,figsize=(6,10))
        ax1.plot(phdf['time'],phdf['F-465'],label='465 channel')
        ax1.plot(phdf['time'],phdf['AF-405'],label='405 channel')
        ax1.axvspan(xmin=pre_interval_[0],xmax=pre_interval_[1],color='k',alpha=0.2)
        ax1.axvspan(xmin=post_interval_[0],xmax=post_interval_[1],color='k',alpha=0.2)
        ax1.set_title('Raw')
        ax1.legend(loc='upper right')
        
        ax2.plot(phdf['time'],phdf['fluo465-pbc'],label='465 channel')
        ax2.plot(phdf['time'],phdf['fluo405-pbc'],label='405 channel')
        ax2.set_title('Photobleaching corrected')
        ax2.legend(loc='lower right')

        ax3.plot(phdf['time'],phdf['fluo465-pbc'],label='465 channel, before')
        ax3.plot(phdf['time'],phdf['fluo405-maf'],label='405 channel, fit')
        ax3.plot(phdf['time'],phdf['fluo465-mac'],label='465 channel, after')
        ax3.set_title('Motion artifact correction')
        ax3.legend(loc='lower right')
        
        ax4.plot(phdf['time'],phdf['fluo465-zsc'],label='465 channel')
        ax4.set_title('Z-score normalization')
        ax4.legend(loc='lower right')
        
        fig.suptitle(phmtry_file)
        fig.tight_layout()
        plt.savefig(os.path.join(out_data_dir,'{}-phmtry.pdf'.format(animal)),format='pdf')
        plt.show()
        
        phdf.to_csv(os.path.join(out_data_dir,'{}-phmtry.csv'.format(animal)))