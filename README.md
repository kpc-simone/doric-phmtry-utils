# doric-phmtry-utils: Tools for processing photometry timeseries data
This repository contains scripts for common processing/analysis steps of photometry data steps in behavioral/systems neuroscience.

Scripts Included:

- [x] Convert from .doric format to .csv and apply basic processing steps.
	- [x] Apply photobleaching correction
	- [x] Apply motion artifact correction
	- [x] Convert to z-score
- [x] Extract segments of timeseries data given a spreadsheet

`
git clone https://github.com/kpc-simone/doric-phmtry-utils.git
cd doric-phmtry-utils
conda env create -f environment.yml --name doric-phmtry-utils-env 
conda activate doric-phmtry-utils
`

# Usage
The scripts folder contains tools for processing photometry data in the .doric file format.

### 1. Convert from .doric format to .csv and apply basic processing steps.

From the /script folder, run this command with the following optional arguments:

`
python process_doric_phmtry.py --baseline_pre_start [baseline_pre_start] --baseline_pre_end [baseline_pre_end] --baseline_post_start [baseline_post_start] --baseline_post_end [baseline_post_end]
`

Meaning of parameters:

`baseline_pre_start`: Start of the pre-baseline interval for photobleaching correction, relative to the start of the recording. Default is t=100 seconds.

`baseline_pre_end`: End of the pre-baseline interval for photobleaching correction, relative to the start of the recording. Default is t=600 seconds.

`baseline_post_start`: Start of the pre-baseline interval for photobleaching correction, relative to the end of the recording. Default is t=500 seconds.

`baseline_post_end`: End of the pre-baseline interval for photobleaching correction, relative to the end of the recording. Default is t=0 seconds.
