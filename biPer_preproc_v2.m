function biPer_preproc_v2(filename, montage_name)
% biPer_prepoc() % Performs the pre-processing of dual EEG recordings using
%                  EEGLAB functions.
%%INPUT
% filename        [string] File name (including extension) of the dual EEG
%                 recording. e.g. 'eegxxx.vhdr'
%                 If the  EEG file is not in your current folder, the
%                 filepath should be also included in the string.
%                 e.g. 'C:\Alejandro\biper_second\eegxxx.vhdr'
% montage_file    [string] File name (including the extension) of the
%                 channels locations file. The file has to be in a
%                 folder included in the Matlab's path.
%                 e.g. 'montage_biPer_27.ced'.
%%OUTPUT
%                 A folder named as the raw EEG file (eegxx) containing several steps of
%                 pre-processed EEG data in EEGLAB format (.set). Last step is named
%                 eegxxx_A_ICA_dipfit.set and eegxxx_B_ICA_dipfit.set
% 
% The following EEGLAB plugins are used:
%   1.- Fileio20200820
%   2.- REST_reference_v1.1_20190818
%   3.- AMICA(1.5.1)
%   4.- PrepPipeline (v0.55.4)
%   5.- clean_rawdata (v2.2)
%   6.- Fieldtrip-lite20200820
%   Also calls the function hyper4leadfield()
%
% This function is part of BiPersonal_pipeline set of functions and it have been tested 
% using EEGLAB versions v2019.1 and v2020.0

% Author: Alejandro Perez
% version 1: University of Toronto, UTSC, June 7th 2018.
% version 2: MRC-CBU, University of Cambridge, March 21, 2020.

% adding functions from biPersonal_pipeline to the path
addpath(genpath('U:\Alejandro\Matlab_functions\biPersonal_pipeline'));

% Checking out if the path was provided
if contains(filename,filesep)
    foldersep = strfind(filename,filesep);
    filepath = filename(1:foldersep(end));
    cd(filepath);
    filename = filename(foldersep(end)+1:end); clear foldersep;
end

% Obtaining the file name (without the extension)
point = strfind(filename,'.');
name = filename(1:point-1); clear point;

% creating a new directory to save data
mkdir(name);
% path to save the data
saving_path  = [pwd filesep name];
% getting the path to the montage
montage_biPer = [fileparts(which(montage_name)) filesep montage_name];

% Variables to be used in the dipole fitting (adding folder to the path)
path_dipfit = what('dipfit'); path_dipfit = path_dipfit.path;
addpath(genpath(path_dipfit));
templateChannelFilePath       = [fileparts(which('standard_1005.elc')) filesep 'standard_1005.elc'];
hdmFilePath                   = [fileparts(which('standard_vol.mat')) filesep 'standard_vol.mat'];
MRIfile                       = [fileparts(which('standard_mri.mat')) filesep 'standard_mri.mat'];

%%% Importing the raw EEG data using Fileio plugin
EEG = pop_fileio(filename);

%%% Resampling data
if EEG.srate>256
    EEG = pop_resample(EEG, 250); % reduces storage and unnecessary high-freq info
end

% Estimating filter order
filter_order = pop_firwsord('blackman', EEG.srate, 1.5);

%%% High-pass filtering
lw = 1; % low cutoff or highpass frequency (use 2 Hz to obtain very stable ICAs)
EEG = pop_firws(EEG, 'fcutoff', lw, 'ftype', 'highpass', 'wtype', 'blackman', 'forder', filter_order);
clear lw filter_order;

%%% Removing line noise 
% Option 1: using cleanLineNoise() included in Prep Pipeline
EEG = cleanLineNoise(EEG, struct('lineNoiseChannels', [1:EEG.nbchan], ...
       'lineFrequencies', [50  100], 'Fs', EEG.srate, 'p', 0.01, ...
       'fScanBandWidth', 2, 'taperBandWidth', 2, 'taperWindowSize', 4, ...
       'pad', 0, 'taperWindowStep', 1, 'fPassBand', [0  EEG.srate/2], ...
       'tau', 100, 'maximumIterations', 10));

% % Option 2: Spectrum Interpolation using FieldTrip (untested).
% path_Fieldtrip = what('Fieldtrip-lite20200820'); % add Fieldtrip to path
% path_Fieldtrip = path_Fieldtrip.path;
% addpath(genpath(path_Fieldtrip));
% Fline = [50 100]; % line noise frequency (and harmonics) in Hz
% [temp_data] = ft_preproc_dftfilter(EEG.data, EEG.srate, Fline, ...
%     'Flreplace','neighbour', 'Flwidth',[1 2], 'Neighwidth',[2 2]);
% EEG.data = temp_data; clea temp_data;
% % removing from path to avoid potential conflicts
% rmpath(genpath(path_Fieldtrip));

% % Cutting signal close to the onset/offset of first/last stimuli
% In the BVA format, event(1) corresponds to the marker 'New Segment'.
lat1 = EEG.event(2).latency;  
lat2 = EEG.event(end).latency;
% % Option 1: if there is a mark at the end of the trial
% howmuch = 3; % time in seconds before the first and after the last marker.
% EEG = pop_select( EEG,'point',[(lat1 - howmuch*EEG.srate) (lat2 + howmuch*EEG.srate)]);
% clear howmuch lat1 lat2;
% % Option 2: If you don't mark the end of the trial
trial_length = 120;             % Trial length (sec) e.g. trial_lenght = 2;
before = 3;                     % Time (sec) to keep before the first mark appears
after  = before + trial_length; % Time (sec) to keep after the last mark appears
EEG = pop_select( EEG,'point',[(lat1 - before*EEG.srate) (lat2 + after*EEG.srate)]);
clear before after lat1 lat2;

n_ch = size(EEG.data,1)/2;  % Number of channels per participant

% Loop for the two EEG recordings
for subj=[1 2]
    switch(subj)
        case {1}
            AorB = 'A';
            % data subject A (e.g. channel 1:32)
            EEG_temp = pop_select(EEG,'channel',[1:n_ch]);
            % adding channel locations
            EEG_temp = pop_editset(EEG_temp, 'chanlocs', montage_name);
            % removing the EOG and EMG channels but keeping it in the EEG struct
            EEG_temp.emg = EEG_temp.data(end-4:end,:);
            EEG_temp = pop_select( EEG_temp, 'channel',[1:27]); % {'cheek' 'HEOL' 'boca_up' 'boca_down' 'HEOR'}
            EEG_temp = eeg_checkset( EEG_temp );
            
            % Calculating the leadfield using REST (once)
            [LeadField, coordinateTransformParameters] = hyper4leadfield(EEG_temp);
            save([saving_path filesep 'leadfield_' name '.mat'], 'LeadField');
            
        case {2}
            AorB = 'B';
            % data subject B (e.g. channel 33:64)
            EEG_temp = pop_select(EEG,'channel',[n_ch+1:2*n_ch]);
            % adding channel locations
            EEG_temp = pop_editset(EEG_temp, 'chanlocs', montage_name);
            % removing the EOG and EMG channels but keeping it in the EEG struct
            EEG_temp.emg = EEG_temp.data(end-4:end,:);
            EEG_temp = pop_select( EEG_temp, 'channel',[1:27]); %{'M_B' 'HEOL_B' 'boca_upB' 'boca_aB' 'HEOR_B'}
            EEG_temp = eeg_checkset( EEG_temp );
            load([saving_path filesep 'leadfield_' name '.mat'], 'LeadField');
    end
    
    %%% Apply average reference after adding initial reference
    EEG_temp.nbchan = EEG_temp.nbchan+1;
    EEG_temp.data(end+1,:) = zeros(1, EEG_temp.pnts);
    EEG_temp.chanlocs(1,EEG_temp.nbchan).labels = 'initialReference';
    EEG_temp = pop_reref(EEG_temp, []);
    EEG_temp = pop_select( EEG_temp,'nochannel',{'initialReference'});
    
    %%% Referencing with the Reference Electrode Standardization Technique (REST)
    [data_z] = rest_refer(EEG_temp.data,LeadField);
    EEG_temp.data = data_z; clear data_z;
    
    %%% Bad channel detection using PREP_pipeline
    reference = findNoisyChannels(EEG_temp); % using the PREP_pipeline function
    EEG_temp.bad_chans = reference.noisyChannels; clear reference;
    
    % saving data first time
    pop_saveset( EEG_temp, 'filename',[name '_' AorB '.set'],'filepath',saving_path);
    
    %%% rejecting bad channels
    ch4interp = EEG_temp.chanlocs; % var for later interpolation
    
    % Channels Fp1,Fp2,F7,F8 not consider as bad if detected by correlation method
    b_ch = EEG_temp.bad_chans.badChannelsFromCorrelation;
    ch = find(b_ch == 1 | b_ch == 2 | b_ch== 11 | b_ch==12); % index depends on your montage
    b_ch(ch) = [];
    b_ch = [b_ch, EEG_temp.bad_chans.badChannelsFromLowSNR, ...
        EEG_temp.bad_chans.badChannelsFromHFNoise, ...
        EEG_temp.bad_chans.badChannelsFromDeviation];
    b_ch = unique(b_ch);
    if ~isempty(ch)
    EEG_temp = pop_select( EEG_temp, 'nochannel', b_ch);
    end
    
    %%% interpolating removed channels
    if ~isempty(ch)
        EEG_temp = pop_interp( EEG_temp, ch4interp, 'spherical');
    end
    clear ch b_ch;
    
    % estimating data rank
    tmprank = rank(EEG_temp.data);
    covarianceMatrix = cov(EEG_temp.data', 1);
    [~, D] = eig (covarianceMatrix);
    rankTolerance = 1e-7;
    tmprank2=sum (diag (D) > rankTolerance);
    if tmprank ~= tmprank2
        tmprank2 = min(tmprank, tmprank2);
    end; clear D tmprank covarianceMatrix rankTolerance;
    
    % Discard channels to make the data full ranked.
    channelSubset = loc_subsets(EEG_temp.chanlocs, tmprank2); clear tmprank2;
    EEG_temp = pop_select( EEG_temp,'channel', channelSubset{1});
    chans_temp = pop_chancenter(EEG_temp.chanlocs, [],[]);
    EEG_temp = pop_editset(EEG_temp, 'chanlocs', chans_temp);
    
    %%% ASR algorithm
    cleaned_eeg = clean_artifacts(EEG_temp, 'FlatlineCriterion','off','ChannelCriterion','off', ...
        'LineNoiseCriterion','off','Highpass','off','BurstCriterion',10,'WindowCriterion',0.4, ...
        'BurstRejection','off','Distance','Euclidian','WindowCriterionTolerances',[-inf 7]);

    % getting irrecoverable time windows (i.e. missing data)
    mask = cleaned_eeg.etc.clean_sample_mask;
    
    % visualization of the cleaned data (comment if runing several dyads serially)
    vis_artifacts(cleaned_eeg,EEG_temp);
    pause; close all; clear cleaned_eeg;
    
    EEG_temp = clean_artifacts(EEG_temp, 'FlatlineCriterion','off','ChannelCriterion','off', ...
        'LineNoiseCriterion','off','Highpass','off','BurstCriterion',10,'WindowCriterion','off', ...
        'BurstRejection','off','Distance','Euclidian','WindowCriterionTolerances',[-inf 7]);
    EEG_temp.mask = mask;
    
    % saving data ASR cleaned
    pop_saveset( EEG_temp, 'filename',[name '_ASR_' AorB '.set'],'filepath',saving_path);
    
    % Adding the mask as a channel
    EEG_temp.nbchan = EEG_temp.nbchan+1;
    EEG_temp.data(end+1,:) = mask;
    EEG_temp.chanlocs(1,EEG_temp.nbchan).labels = 'mask';
    
    %%% Epoching data, including the mask from ASR
    % Markers and epoching parameters depends on your experiment
    EEG_temp = pop_epoch( EEG_temp, { 'S113' 'S213' 'S114' 'S214' 'S115' ...
        'S215' 'S13' 'S14' 'S15' 'S21' 'S22' 'S26' 'S121' 'S221' 'S122' ...
        'S222' 'S126' 'S226'},[-1  121], 'newname', ' resampled epochs', ...
        'epochinfo', 'yes');
    
    % Including the epoched mask
    mask_epoched = EEG_temp.data(end,:,:);
    EEG_temp.mask_epoched = mask_epoched;
    
    % Removing the mask channel
    EEG_temp = pop_select( EEG_temp,'nochannel',{'mask'});
    
    % Saving epoched data
    pop_saveset( EEG_temp, 'filename',[name '_epoched_' AorB '.set'],'filepath',saving_path);
    
    %%% Performing ICA with AMICA. (time consuming)
    % Calculating data rank for epoched data
    Temp = EEG_temp.data; Temp = reshape(Temp,size(Temp,1),[],1);
    dataRank = rank(double(Temp')); clear Temp;
    
    runamica15(EEG_temp.data, 'num_chans', EEG_temp.nbchan,...
        'outdir', [saving_path filesep 'AMICA_' AorB],...
        'pcakeep', dataRank, 'num_models', 1, 'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1);
    EEG_temp.etc.amica   = loadmodout15([saving_path filesep 'AMICA_' AorB]);
    EEG_temp.icaweights  = EEG_temp.etc.amica.W;
    EEG_temp.icasphere   = EEG_temp.etc.amica.S;
    EEG_temp = eeg_checkset(EEG_temp, 'ica');
    pop_saveset( EEG_temp, 'filename', [name '_' AorB '_ICA.set'],'filepath', saving_path);
    
    % Fitting dipoles to ICs
    [~,coordinateTransformParameters] = coregister(EEG_temp.chanlocs, templateChannelFilePath, 'warp', 'auto', 'manual', 'off');
    EEG_temp = pop_dipfit_settings( EEG_temp, 'hdmfile', hdmFilePath, 'coordformat', 'MNI', 'mrifile', MRIfile, ...
        'chanfile', templateChannelFilePath, 'coord_transform', coordinateTransformParameters, 'chansel', 1:EEG_temp.nbchan);
    EEG_temp = pop_multifit(EEG_temp, 1:EEG_temp.nbchan,'threshold', 100, 'dipplot','off','plotopt',{'normlen' 'on'}); % 'dipplot','on',
    pop_saveset( EEG_temp, 'filename', [name '_' AorB '_ICA_dipfit.set'],'filepath', saving_path);
    clear EEG_temp;
end
end

