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

import os
from pathlib import Path
from typing import Dict, Optional, Tuple, Union, List

# Default mapping to homogenize Doric signal names across versions
DEFAULT_SIGNALS_MAP: Dict[str, str] = {
    "Console_time(s)": "time",
    "AIn-1 - Dem (AOut-1)": "F-465",
    "AIn-2 - Dem (AOut-2)": "AF-405",
    "DI--O-3": "DI/03",
    "DI--O-4": "DI/04",
    "AIN01xAOUT01-LockIn/Time": "time",
    "AIN01xAOUT01-LockIn/Values": "F-465",
    "AIN02xAOUT02-LockIn/Values": "AF-405",
    "DigitalIO/DIO01": "DI/01",
    "DigitalIO/DIO02": "DI/02",
    "DigitalIO/DIO03": "DI/03",
}


def process_photometry_directory(
    data_dir: Union[str, Path, None] = None,
    *,
    baseline_pre_start: float = 100.0,
    baseline_pre_end: float = 400.0,
    baseline_post_start: float = 300.0,
    baseline_post_end: float = 0.0,
    zscore_baseline_interval: Tuple[float, float] = (120.0, 480.0),
    signals_map: Optional[Dict[str, str]] = None,
    out_dir: Union[str, Path, None] = None,
    make_plots: bool = True,
    show_plots: bool = True,
    save_plots: bool = True,
    save_csv: bool = True,
    overwrite: bool = True,
    file_suffix: str = ".doric",
) -> List[Path]:
    """
    Processes raw Doric photometry files in a directory:
      1) loads Doric files
      2) photobleaching correction (pre + post baseline windows)
      3) motion artifact correction
      4) z-score normalization

    Parameters
    ----------
    data_dir:
        Folder containing raw photometry .doric files. If None, opens a folder picker dialog.
    baseline_pre_start, baseline_pre_end:
        Pre-baseline interval (seconds) relative to recording start for photobleaching correction.
    baseline_post_start, baseline_post_end:
        Post-baseline interval (seconds) relative to recording end for photobleaching correction.
        Example: post_interval = (tmax - baseline_post_start, tmax - baseline_post_end)
    zscore_baseline_interval:
        Baseline interval (seconds) used for z-score normalization.
    signals_map:
        Optional mapping from Doric signal names to canonical names. Defaults to DEFAULT_SIGNALS_MAP.
    out_dir:
        Output directory for processed CSVs and PDFs. Defaults to ../phmtry-processed relative to data_dir.
    make_plots / show_plots / save_plots:
        Controls for generating, displaying, and saving the QC plots.
    save_csv:
        Whether to save processed dataframes as CSV.
    overwrite:
        Whether to overwrite existing outputs.
    file_suffix:
        File suffix to include (default ".doric").

    Returns
    -------
    outputs:
        List of paths to generated outputs (CSVs and/or PDFs).
    """
    # Lazy import so the function works in non-GUI contexts if data_dir is provided
    if data_dir is None:
        from tkinter.filedialog import askdirectory  # local import to avoid Tk dependency if unused

        data_dir = askdirectory(title="Select folder containing raw photometry data")

    data_dir = Path(data_dir).expanduser().resolve()
    if not data_dir.exists():
        raise FileNotFoundError(f"data_dir does not exist: {data_dir}")

    if signals_map is None:
        signals_map = DEFAULT_SIGNALS_MAP

    # Default output directory mirrors your script: ../phmtry-processed relative to data_dir
    if out_dir is None:
        out_dir = (data_dir / ".." / "phmtry-processed").resolve()
    else:
        out_dir = Path(out_dir).expanduser().resolve()

    out_dir.mkdir(parents=True, exist_ok=True)

    # Find Doric files
    phmtry_files = sorted([p for p in data_dir.iterdir() if p.is_file() and p.name.endswith(file_suffix)])
    if len(phmtry_files) == 0:
        raise FileNotFoundError(f"No files ending with '{file_suffix}' found in: {data_dir}")

    outputs: List[Path] = []

    for filepath in phmtry_files:
        phmtry_file = filepath.name
        animal = phmtry_file.split("_")[0]
        print(f"\nProcessing: {phmtry_file} (animal={animal})")

        # load the data from doric file
        data = load_phmtry_raw_doric(str(filepath))

        # build dataframe and normalize column names
        phdf = (
            pd.concat([pd.DataFrame(v, columns=[k]) for k, v in data.items()], axis=1)
            .rename(columns=signals_map)
        )

        # keep canonical columns and drop missing
        required = ["time", "F-465", "AF-405"]
        missing = [c for c in required if c not in phdf.columns]
        if missing:
            raise KeyError(
                f"Missing required columns after rename for file '{phmtry_file}': {missing}\n"
                f"Available columns: {list(phdf.columns)}"
            )

        phdf = phdf[required].dropna()
        if phdf.empty:
            raise ValueError(f"After selecting {required} and dropping NaNs, dataframe is empty for: {phmtry_file}")

        tmax = float(phdf["time"].max())
        pre_interval = (baseline_pre_start, baseline_pre_end)
        post_interval = (tmax - baseline_post_start, tmax - baseline_post_end)
        print(f"  pre_interval  = {pre_interval}")
        print(f"  post_interval = {post_interval} (tmax={tmax:.3f})")

        # photobleaching correction
        phdf["fluo465-pbc"] = correct_photobleaching(
            ts=phdf["time"],
            ys=phdf["F-465"],
            pre_interval=pre_interval,
            post_interval=post_interval,
        )
        phdf["fluo405-pbc"] = correct_photobleaching(
            ts=phdf["time"],
            ys=phdf["AF-405"],
            pre_interval=pre_interval,
            post_interval=post_interval,
        )

        # motion artifact correction
        phdf["fluo405-maf"], phdf["fluo465-mac"] = correct_motion(
            phdf["fluo465-pbc"].to_numpy(),
            phdf["fluo405-pbc"].to_numpy(),
        )

        # z-score normalization
        phdf["fluo465-zsc"] = transform_to_zscore(
            ts=phdf["time"],
            ys=phdf["fluo465-mac"],
            baseline_interval=zscore_baseline_interval,
        )

        # Outputs
        pdf_path = out_dir / f"{animal}-phmtry.pdf"
        csv_path = out_dir / f"{animal}-phmtry.csv"

        # Save CSV
        if save_csv:
            if csv_path.exists() and not overwrite:
                print(f"  Skipping CSV (exists, overwrite=False): {csv_path}")
            else:
                phdf.to_csv(csv_path, index=False)
                outputs.append(csv_path)
                print(f"  Wrote CSV: {csv_path}")

        # Make/save plots
        if make_plots:
            fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, figsize=(6, 10))

            ax1.plot(phdf["time"], phdf["F-465"], label="465 channel")
            ax1.plot(phdf["time"], phdf["AF-405"], label="405 channel")
            ax1.axvspan(xmin=pre_interval[0], xmax=pre_interval[1], color="k", alpha=0.2)
            ax1.axvspan(xmin=post_interval[0], xmax=post_interval[1], color="k", alpha=0.2)
            ax1.set_title("Raw")
            ax1.legend(loc="upper right")

            ax2.plot(phdf["time"], phdf["fluo465-pbc"], label="465 channel")
            ax2.plot(phdf["time"], phdf["fluo405-pbc"], label="405 channel")
            ax2.set_title("Photobleaching corrected")
            ax2.legend(loc="lower right")

            ax3.plot(phdf["time"], phdf["fluo465-pbc"], label="465 channel, before")
            ax3.plot(phdf["time"], phdf["fluo405-maf"], label="405 channel, fit")
            ax3.plot(phdf["time"], phdf["fluo465-mac"], label="465 channel, after")
            ax3.set_title("Motion artifact correction")
            ax3.legend(loc="lower right")

            ax4.plot(phdf["time"], phdf["fluo465-zsc"], label="465 channel")
            ax4.set_title("Z-score normalization")
            ax4.legend(loc="lower right")

            fig.suptitle(phmtry_file)
            fig.tight_layout()

            if save_plots:
                if pdf_path.exists() and not overwrite:
                    print(f"  Skipping PDF (exists, overwrite=False): {pdf_path}")
                else:
                    fig.savefig(pdf_path, format="pdf")
                    outputs.append(pdf_path)
                    print(f"  Wrote PDF: {pdf_path}")

            if show_plots:
                plt.show()
            else:
                plt.close(fig)

    return outputs


# Example usage (optional):
# outputs = process_photometry_directory(r"C:\path\to\raw_photometry")
# print(outputs)

def extract_event_aligned_photometry(
    data_dir: Union[str, Path, None] = None,
    ec_file: Union[str, Path, None] = None,
    out_data_dir: Union[str, Path, None] = None,
    *,
    segstart_colname: str = "TPoI-2",
    baseline_end: str = "TPoI-1",
    outcome_colname: str = "observation",
    pre_interval: float = -5.0,
    event_dur: float = 8.0,
    post_interval: float = 12.0,
    baseline_pre_interval: float = -5.0,
    SR_desired: int = 100,
    signal_col: str = "fluo465-zsc",
    time_col: str = "time",
    animal_col_in_ec: str = "animal",
    trial_col_in_ec: str = "trial",
    outfilename: Optional[str] = None,
    make_trial_qc_plots: bool = True,
    show_plots: bool = True,
    save_out_csv: bool = True,
    overwrite: bool = True,
) -> Tuple[pd.DataFrame, Optional[Path]]:
    """
    Extract trial-aligned segments from processed photometry CSVs, using an experimental-conditions (EC) CSV
    that provides event times per trial.

    Workflow:
      1) Load EC file, drop rows missing segstart_colname
      2) For each photometry CSV in data_dir:
         - match EC rows by animal id
         - resample z-scored signal to SR_desired
         - extract segments for each trial, with per-trial baseline windows
         - concatenate into a single wide dataframe: [time, trial1, trial2, ...]
         - optional QC plots per animal
      3) Save combined dataframe to out_data_dir
      4) Plot overall mean ± std across trials

    Returns
    -------
    out_df:
        DataFrame with 'time' column and one column per extracted trial.
    out_path:
        Path to saved CSV (or None if not saved).
    """
    # Lazy GUI imports so function works headless if args are provided
    if data_dir is None or ec_file is None or out_data_dir is None:
        from tkinter.filedialog import askdirectory, askopenfilename  # local import

        if data_dir is None:
            data_dir = askdirectory(
                title="Select folder containing processed photometry recordings, from which to extract all events."
            )
        if ec_file is None:
            ec_file = askopenfilename(
                title="Select experimental conditions (ec) file, that defines the start and end times of the events."
            )
        if out_data_dir is None:
            out_data_dir = askdirectory(title="Select directory to save the output data.")

    data_dir = Path(data_dir).expanduser().resolve()
    ec_file = Path(ec_file).expanduser().resolve()
    out_data_dir = Path(out_data_dir).expanduser().resolve()
    out_data_dir.mkdir(parents=True, exist_ok=True)

    if not data_dir.exists():
        raise FileNotFoundError(f"data_dir does not exist: {data_dir}")
    if not ec_file.exists():
        raise FileNotFoundError(f"ec_file does not exist: {ec_file}")

    # Load EC and filter rows that define the event start
    ecdf = pd.read_csv(ec_file)
    if segstart_colname not in ecdf.columns:
        raise KeyError(f"EC file missing segstart_colname column '{segstart_colname}'. Columns: {list(ecdf.columns)}")
    if baseline_end not in ecdf.columns:
        raise KeyError(f"EC file missing baseline_end column '{baseline_end}'. Columns: {list(ecdf.columns)}")
    if outcome_colname not in ecdf.columns:
        raise KeyError(f"EC file missing outcome_colname column '{outcome_colname}'. Columns: {list(ecdf.columns)}")
    if animal_col_in_ec not in ecdf.columns:
        raise KeyError(f"EC file missing animal column '{animal_col_in_ec}'. Columns: {list(ecdf.columns)}")
    if trial_col_in_ec not in ecdf.columns:
        raise KeyError(f"EC file missing trial column '{trial_col_in_ec}'. Columns: {list(ecdf.columns)}")

    ecdf = ecdf[ecdf[segstart_colname].notna()].copy()

    # Build common output timebase
    ts_out = np.arange(start=pre_interval, stop=event_dur + post_interval, step=1 / SR_desired)
    out_df = pd.DataFrame({"time": ts_out})

    # Default output filename
    if outfilename is None:
        outfilename = f"events-alignedto-{segstart_colname}-zscore.csv"
    out_path = out_data_dir / outfilename

    # Gather photometry CSVs
    phmtry_files = sorted([p for p in data_dir.iterdir() if p.is_file() and p.suffix.lower() == ".csv"])
    if len(phmtry_files) == 0:
        raise FileNotFoundError(f"No .csv files found in: {data_dir}")

    for ph_path in phmtry_files:
        phmtry_file = ph_path.name
        animal = phmtry_file.split("-")[0]  # preserves your original naming convention
        print(f"\nAnimal: {animal} | file: {phmtry_file}")

        an_ecdf = ecdf[ecdf[animal_col_in_ec] == animal]
        if len(an_ecdf) == 0:
            print(f"  No EC rows found for animal '{animal}'. Skipping.")
            continue

        # Per-trial baseline and segment windows
        bl_times: List[List[float]] = [
            [float(row[baseline_end]) + baseline_pre_interval, float(row[baseline_end])]
            for _, row in an_ecdf.iterrows()
        ]
        seg_times: List[List[float]] = [
            [float(row[segstart_colname]) + pre_interval, float(row[segstart_colname]) + event_dur + post_interval]
            for _, row in an_ecdf.iterrows()
        ]

        labels: List[str] = [
            f"{row[animal_col_in_ec]}-t{row[trial_col_in_ec]}-{row[outcome_colname]}"
            for _, row in an_ecdf.iterrows()
        ]

        # Load photometry
        phdf = pd.read_csv(ph_path)
        if time_col not in phdf.columns:
            raise KeyError(f"{phmtry_file} missing time column '{time_col}'. Columns: {list(phdf.columns)}")
        if signal_col not in phdf.columns:
            raise KeyError(f"{phmtry_file} missing signal column '{signal_col}'. Columns: {list(phdf.columns)}")

        # Resample
        ts_res, ys_res = resample_signal(ts=phdf[time_col], ys=phdf[signal_col], SR_desired=SR_desired)

        # Extract segments (assumes extract_segments handles baseline normalization when bl_times provided)
        segs_df = extract_segments(
            ts=ts_res,
            ys=ys_res,
            seg_times=seg_times,
            labels=labels,
            bl_times=bl_times,
        )

        # Merge into global output df (aligned on implicit row order of ts_out)
        out_df = pd.concat([out_df, segs_df], axis=1)

        # QC plots per animal
        if make_trial_qc_plots:
            fig, (ax1, ax2) = plt.subplots(
                1, 2, figsize=(9, 3), gridspec_kw={"width_ratios": [2, 1]}
            )

            # Plot 1: verify extraction in original time coordinates
            ax1.plot(ts_res, ys_res, color="dimgray")
            for trial_col, (seg_start, seg_end) in zip(segs_df.columns, seg_times):
                ax1.plot(out_df["time"] + seg_start - pre_interval, out_df[trial_col])
                ax1.axvline(seg_start, color="k", linestyle="--")
                ax1.axvline(seg_end, color="k", linestyle="--")

            # Plot 2: verify normalization / alignment in event time coords
            ax2.plot(out_df["time"], out_df[segs_df.columns])
            ax2.axvline(0.0, color="k", linestyle="--")
            ax2.axvline(event_dur, color="k", linestyle="--")

            for ax in (ax1, ax2):
                ax.set_xlabel("Time (s)")
                ax.set_ylabel("z-score")

            fig.suptitle(animal)
            fig.tight_layout()

            if show_plots:
                plt.show()
            else:
                plt.close(fig)

    # Save combined output
    saved_path: Optional[Path] = None
    if save_out_csv:
        if out_path.exists() and not overwrite:
            print(f"Output exists and overwrite=False, not saving: {out_path}")
        else:
            out_df.to_csv(out_path, index=False)
            saved_path = out_path
            print(f"\nWrote: {out_path}")

    # Summary plot: mean ± std across trials (if any trials exist)
    trial_cols = list(out_df.columns[1:])
    if len(trial_cols) > 0:
        fig, ax = plt.subplots(1, 1)

        zsc_mean = out_df[trial_cols].mean(axis=1)
        zsc_std = out_df[trial_cols].std(axis=1)

        ax.plot(out_df["time"], zsc_mean)  # (no explicit color)
        ax.fill_between(out_df["time"], zsc_mean - zsc_std, zsc_mean + zsc_std, alpha=0.3)

        ax.set_xlabel("Time (s) from event onset")
        ax.set_ylabel("Z-Score")
        ax.axvline(0.0, color="k", linestyle="--")
        ax.axvline(event_dur, color="k", linestyle="--")
        ax.set_xlim(pre_interval, event_dur + post_interval)
        fig.tight_layout()

        if show_plots:
            plt.show()
        else:
            plt.close(fig)
    else:
        print("No trial columns were extracted; skipping summary plot.")

    return out_df, saved_path


# Example usage:
# out_df, out_path = extract_event_aligned_photometry(
#     data_dir=r"C:\data\phmtry-processed",
#     ec_file=r"C:\data\conditions.csv",
#     out_data_dir=r"C:\data\outputs",
#     segstart_colname="TPoI-2",
#     baseline_end="TPoI-1",
#     outcome_colname="observation",
# )
# print(out_path)
