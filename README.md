# biPersonal_pipeline

biPersonal_pipeline consists in a set of Matlab functions for the processing of dual-EEG data (hyperscanning).
Each of the functions performs a step in the processing of EEG hyperscanning by using EEGLAB functions arranged in neccesary sequential order.
These functions were tailor-maded to work on data from the experiment reported at:
Pérez et al. (2019) Differential brain-to-brain entrainment while speaking and listening in native and foreign languages. Cortex, 111, 303-315.
doi: 10.1016/j.cortex.2018.11.026;
but it could be easily adapted to other EEG formats and designs.

The functions comprised in the biPersonal_Pipeline are:

biPer_preproc_v2.m
   The following steps are performed over the raw EEG hyperscanning data (same order)
   1.- Import
   2.- Downsampling
   3.- High-pass filtering
   4.- Line noise removal
   5.- Selecting data recorded during the experiment (i.e. removing EEG recorded pre and post experiment)
   6.- Splitings data into the two individual recordings
   ---- and then, for each individual recording (i.e. twice)----
   7.- Find noisy channels
   8.- Removing bad channles
   9.- Interpolation
  10.- Re-referencing to the average
  11.- Referencing with REST
  12.- ASR cleaning
  13.- ICA decomposition (AMICA)
  14.- Dipole fitting
  
biPer_cleaning_v2.m
   To be used after function biPer_preproc(). Assumes the existence inside folder (xxx) of two EEG recordings (xxx_A.set and xxx_B.set) .
   These datasets should be decomposed by ICA and with fitted dipoles.
   The following steps are performed over the individual data (same order)
   1.- estimating symetrically constrained bilateral dipoles
   2.- keeping ICs most probably containing neural activity
