# Purpose

This repo gathers the required data and code to reproduce the neural stitching results in the published LFADS paper, specifically those results presented in Figure 4.

Pandarinath C, O’Shea DJ, Collins J, Jozefowicz R, Stavisky SD, Kao JC, et al. Inferring single-trial neural population dynamics using sequential auto-encoders. Nat Methods. 2018;15: 805–815. doi:10.1038/s41592-018-0109-9
https://www.nature.com/articles/s41592-018-0109-9

![Figure 4](https://raw.githubusercontent.com/nplcode/lfads-neural-stitching-reproduce/master/fig4_raster.png "Figure 4")

This repo was prepared by Daniel O'Shea, for questions or issues, it would be preferable to file an issue on Github directly.

# What you'll need

* *Git LFS* (required to download the larger portions of this repo): https://git-lfs.github.com
* *MATLAB*, particularly a relatively modern version. I'm using 2019a but 2017 onwards should work.
* *LFADS Run Manager*: https://github.com/lfads/lfads-run-manager (see docs at https://lfads.github.io/lfads-run-manager/), if you want to prepare an LFADS run and load the results into MATLAB easily
* *LFADS Python / TensorFlow*, if you want to retrain the model: see https://lfads.github.io/lfads-run-manager/install/ for installation instructions

# Details

If you just want the data, the data fed into LFADS are located in `export_v05_broadbandRethreshNonSorted_filtered` as mat files 
and the outputs (LFADS posterior means) have been exported as HDF5 files into posterior_means_export for all of the single
dataset runs and the stitched model.

If you'd like to rerun LFADS yourself, follow along in the notebook provided `neural_stitching_walkthrough.ipynb`. 
Then you can train LFADS with instructions provided at https://lfads.github.io/lfads-run-manager/.

The decoding and analysis scripts are too interwoven into my personal analysis toolkit to extract easily, 
but I've included the code regardless for your reference. The saved outputs from these scripts are included in the 
`results` folder. Some of the scripts inside `+PierreEricLFADS` reference these saved outputs to generate individual figure
panels directly, and thus do not depend on the kinematic decoding. The primary result decoding kinematics from LFADS factors 
is done in `kinematicDecodeFromFactors.m`.

# Questions?

Please file an issue on Github first unless you need to contact me directly at djoshea at stanford. Thanks!

