% EEG Data Processing Script for Secondment at Tecnalia
%
% - Experimental Design: Closed-loop digit span test.
% - Closed-Loop Control: Volume (amplitude) modulation based on EEG alpha-band
%   spectral power during interdigit intervals.
%
% - Dataset: 16 participants, each with two EEG recordings:
%   1) Brain Recorder Software: Single continuous recording.
%   2) BCI2000: Four separate recordings, each corresponding to a block.
%
% Software to be used:
%                    Matlab R2022.b
%                    eeglab2023.1
%                    EEGLAB plugins:
%                                   "AMICA" v1.7
%                                   "BCI2000import" v0.36
%                                   "Biosig" v3.8.1
%                                   "Cleanline" v2.00
%                                   "Fieldtrip-lite" v20230926
%                                   "Fileio" v20230926
%                                   "ICLabel" v1.4
%                                   "Viewprops" v1.5.4
%                                   "bva-io" v1.72
%                                   "clean_rawdata" v2.91
%                                   "dipfit" v5.3
%                                   "firfilt" v2.7.1
%                                   "fitTwoDipoles" v1.00

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Alejandro Perez. University of Surrey.
% v1.0: 15/05/2023
% v2.0: 20/09/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SECTION DEPRECATED
% % This section (until line 156) have been deprecated.
%
% % ensuring EEGLAB functions are added to the path
% eeglab;
% % clearing unnecessary vars and closing pop up GUI
% clear; close all;
%
% % path to EEG data at Surrey's OneDrive (pre-defined)
% data_path = 'C:\Users\ap0082\OneDrive - University of Surrey\Alejandro\Experiment3_closed-loop\';
% cd(data_path);
%
% % getting folders corresponding to each subject and organising data paths
% A = dir ('subj*');
% A = struct2table(A);
% A.BrainRecorder_folder = append(A.folder, filesep, A.name, '\Experiment\BrainRecorder\');
% A.BCI2000_folder = append(A.folder, filesep, A.name, '\Experiment\BCI2000\');
%
% for i=1:height(A) % loop across the participants (16)
%
%     % Pre-processing will be performed on the Brain Recorder data which is
%     % continuous and better suited for ICA. However, markers are only in
%     % the BCI2000 recordings. Thus we need to synchronise recordings (find
%     % the equivalence).
%
%     %%%% Brain Recorder data %%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Only one recording is expected. More, will produce an error.
%
%     % participant's path to data recorded w/ Brain Recorder
%     BV_path = char(A.BrainRecorder_folder(i));
%     cd(BV_path);
%     BV_file = dir('*.vhdr');
%     EEG = pop_loadbv(BV_path, BV_file.name);
%     srate = EEG.srate;
%     % Stimtrack audio signal is included in the last channel
%     StimTrack = EEG.data (EEG.nbchan,:);
%     % detrending signal to obtain meaningful results with xcorr function
%     d_StimTrack = detrend(StimTrack);
%
%     %%%% BCI2000 data %%%%
%     %%%%%%%%%%%%%%%%%%%%%%
%     % Data was recorded in blocks. 4 blocks are expected.
%
%     % participant's path to data recorded w/ BCI2000
%     bci_path = char(A.BCI2000_folder(i));
%     cd(bci_path);
%     % getting the info of the blocks. recordings have '.dat' extension
%     BCI = dir('*.dat');
%
%     % Var to contain the each blocks onset on the continous BV recording
%     blockLatency = [];
%
%     for i1 = 1:length(BCI) % loop across (4) blocks
%
%         % loading data
%         EEG_BCI = pop_loadBCI2000(BCI(i1).name); %, {'PhaseInSequence','StimulusBegin','StimulusCode'});
%
%         % Saving BCI2000 data in eeglab format
%         pop_saveset( EEG_BCI, 'filename',['block_' num2str(i1) '_' BCI(i1).name(1:end-4) '.set'],'filepath',bci_path);
%
%         %%% Linking the BCI2000 recordings with the BV recording
%         % StimTrack channel is included in the before last channel
%         Fragment = EEG_BCI.data (EEG_BCI.nbchan-1,:);
%         % detrending signal to obtain meaningful results with xcorr function
%         d_Fragment = detrend(Fragment);
%
%         % Calculating the cross correlation between the audio channel in
%         % the BCI2000 recording and the BV recording
%         [xCorr,lags] = xcorr(d_StimTrack,d_Fragment);
%         % creating figure of the lagged correlation values
%         figure; plot(lags/srate,xCorr); grid; xlabel('Lags (s)');
%
%         % max corr (peak) indicates where signals align
%         [~,locs] = max(abs(xCorr));
%         % getting the exact lag
%         maxt = lags(locs);
%
%         % var with the audio in the block ???
%         Fs = 1/EEG_BCI.srate;
%         Time = 0:1/Fs:(length(StimTrack)-1)/Fs;
%         Trial = NaN(size(StimTrack));
%         Trial(maxt+1:maxt+length(Fragment)) = Fragment;
%
%         % Visual inspection to confirm it was detected correctly
%         figure('Name',['Block ' num2str(i1)]);
%         plot(Time,StimTrack,Time,Trial);
%         xlabel('Time (s)');
%
%         % Feedback about correctness for locating audio fragment
%         answer = questdlg('Is it correct?');
%         switch answer
%             case 'Yes'
%                 blockLatency = cat(2,blockLatency,locs); % Block latencies
%                 points = EEG_BCI.pnts;
%                 save([BCI(i1).name(1:end-4) '.mat'],'locs','points');
%         end
%         clear I answer Fs maxt Time Trial lags xCorr Fragment d_Fragment points;
%
%     end % loop blocks
%
%     blockLabels  = {'b1' 'b2' 'b3' 'b4'}; % Block labels
%     duration = num2cell([1 1 1 1]);
%     urevent  = num2cell([1 2 3 4]);
%
%     % eliminating the original structure
%     EEG.event = []; EEG.urevent = [];
%     % Build a minimal but valid EEG.event from scratch (Makoto)
%
%     % latency
%     latencyInCell = num2cell(blockLatency);
%     [EEG.event(1:length(blockLabels)).latency] = latencyInCell{:};
%     % type
%     [EEG.event(1:length(blockLabels)).type] = blockLabels{:};
%     % duration
%     [EEG.event(1:length(blockLabels)).duration] = duration{:};
%     % Also creating the urevent field
%     EEG.urevent = EEG.event;
%     % urevent
%     [EEG.event(1:length(blockLabels)).urevent] = urevent{:};
%
%     clear blockLatency latencyInCell blockLabels duration urevent StimTrack d_StimTrack BCI;
%     pop_saveset( EEG, 'filename',[BV_file.name(1:end-5) '.set'],'filepath',BV_path);
%
% end % loop participants

%%    Pre-processing 1 ICA (deprecated)
% clear;
% eeglab;
% clear ALLEEG EEG CURRENTSET ALLCOM CURRENTSTUDY globalvars LASTCOM PLUGINLIST STUDY;
% close all;
%
% % data at Surrey's OneDrive (pre-defined)
% data_path = 'C:\Users\ap0082\OneDrive - University of Surrey\Alejandro\Experiment3_closed-loop\';
% cd(data_path);
%
% % getting folders corresponding to each subject
% A = dir ('subj*');
%
% % loop across the participants (16)
% for i=5:length(A)
%
%     %%%% Brain Recorder data - continous %%%%
%     BV_path = [data_path A(i).name filesep 'Experiment\BrainRecorder\']; % path to data recorded w/ Brain Recorder
%     cd(BV_path);
%     BV_file = dir('*.set'); % Only one recording is expected.
%     EEG = pop_loadset('filename', BV_file.name,'filepath',BV_path); % Will produce an error if more than one.
%     EEG = pop_resample(EEG, 250);
%     EEG = pop_eegfiltnew(EEG, 'locutoff',2,'channels', ...
%         {'Fp1','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6','T7', ...
%         'C3','Cz','C4','T8','TP9','CP5','CP1','CP2','CP6','TP10','P7','P3', ...
%         'Pz','P4','P8','PO9','O1','Oz','O2','PO10'});
%     EEG = pop_runica( EEG,'chanind',[1:32] );
%     pop_saveset( EEG, 'filename',[BV_file.name(1:end-4) '_ica.set'],'filepath',BV_path);
%
% end % loop participants

%% Section 1: Data organisation
% % (Commented out to prevent code execution redundancy).
% % This section creates a table containing file information for all filenames.
% % The table is saved in the main directory where the data is located.
%
% % Data Directory
% % The data is located at D:\Tecnalia_data
% dataDir = 'D:\Tecnalia_data\';
% cd(dataDir);
%
% % Data Structure:
% % Within the data directory, individual participant folders are organized,
% % named according to the format: subj##_year-month-day_#########
% % Example: subj01_2021-11-25_131903001 and subj16_2022-02-18_152819001
% subjName = dir('subj*');
% % Each participant folder contains EEG data recorded using BCI2000,
% % divided into 4 recordings corresponding to blocks.
% % The EEG data is already in EEGLab format with .set and .fdt files.
% %
% % Block Naming:
% % The raw BCI2000 data was initially imported and saved using a deprecated
% % section of this script.
% % Blocks are named consistently with the prefix block_#_, followed by the folder name.
% % Example: block_4_subj16_2022-02-18_152819S001R10.set
% %
% % Initializing the table to store filename information across participants.
% sz = [16 9];
% varTypes = ["string","string","string","string","string","string","string","string","string"];
% varNames = ["participant", ...
%     "EEG_block_1","EEG_block_2","EEG_block_3","EEG_block_4", ...
%     "SavedData_block_1","SavedData_block_2","SavedData_block_3","SavedData_block_4"];
% DataTable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
%
% % Loop across participants and block to fulfill the table
% for subj = 1:length(subjName)
%
%     % change directory to the participant's folder
%     cd([dataDir subjName(subj).name]);
%
%     % Assigning the participant's codename
%     DataTable.participant(subj) = subjName(subj).name;
%
%     % Getting all filenames for the .mat data
%     saveddata = dir('SavedData_*');
%
%     % loop across the blocks
%     for i=1:4
%         temp = dir(['block_' num2str(i) '*']);
%         % Conditional to deal with possible differences in the eeglab format:
%         % (only .set) vs. (.set + .fdt)
%         if size(temp,1)==1
%             temp = temp(1).name;
%         elseif size(temp,1)==2
%             temp = temp(2).name;
%         end
%
%         % Assigning the filename of each EEG data block
%         eval(['DataTable.EEG_block_' num2str(i) '(subj) = temp;']);
%         clear temp;
%
% % Assigning the filename for each .mat data file.
% % The filenames include the exact recording hour in military format
% % (00:00 to 23:00). This order corresponds to the 'saveddata' variable order.
%         eval(['DataTable.SavedData_block_' num2str(i) '(subj) = saveddata(i).name;']);
%
%     end % loop across blocks
%     clear saveddata i;
% end % loop across subject
%
% save([dataDir 'Table_with_all_filenames.mat'],'DataTable');

%% Pre-processing step 1 ASR and ICA

% Data Directory
% The data directory is located at D:\Tecnalia_data
dataDir = 'D:\Tecnalia_data\';
cd(dataDir);

% loading var containing filenames and paths created on previous step
load Table_with_all_filenames;
% loading variable containing channel info: chanlocs & chaninfo
% (Created independently. Check you have the file.)
load electrode_layout_Tecnalia.mat;

for subj = 1:size(DataTable,1)

    % initialising eeglab for each participant
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    % Block 1
    EEG = pop_loadset('filename',char(DataTable.EEG_block_1(subj)),'filepath',[dataDir char(DataTable.participant(subj))]);
    EEG.chanlocs = chanlocs;
    EEG.chaninfo = chaninfo;
    % saving StimTrak data
    EEG1 = pop_select( EEG, 'channel',{'Stimtrack'});
    pop_saveset( EEG1, 'filename','Stimtrack_B1.set','filepath',[dataDir char(DataTable.participant(subj))]);
    clear EEG1;
    % removing StimTrak and Trigger channels
    EEG = pop_select( EEG, 'rmchannel',{'Stimtrack','Trig'});
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

    % Block 2
    EEG = pop_loadset('filename',char(DataTable.EEG_block_2(subj)),'filepath',[dataDir char(DataTable.participant(subj))]);
    EEG.chanlocs = chanlocs;
    EEG.chaninfo = chaninfo;
    % saving StimTrak data
    EEG2 = pop_select( EEG, 'channel',{'Stimtrack'});
    pop_saveset( EEG2, 'filename','Stimtrack_B2.set','filepath',[dataDir char(DataTable.participant(subj))]);
    clear EEG2;
    % removing StimTrak and Trigger channels
    EEG = pop_select( EEG, 'rmchannel',{'Stimtrack','Trig'});
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

    % Block 3
    EEG = pop_loadset('filename',char(DataTable.EEG_block_3(subj)),'filepath',[dataDir char(DataTable.participant(subj))]);
    EEG.chanlocs = chanlocs;
    EEG.chaninfo = chaninfo;
    % saving StimTrak data
    EEG3 = pop_select( EEG, 'channel',{'Stimtrack'});
    pop_saveset( EEG3, 'filename','Stimtrack_B3.set','filepath',[dataDir char(DataTable.participant(subj))]);
    clear EEG3;
    % removing StimTrak and Trigger channels
    EEG = pop_select( EEG, 'rmchannel',{'Stimtrack','Trig'});
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

    % Block 4
    EEG = pop_loadset('filename',char(DataTable.EEG_block_4(subj)),'filepath',[dataDir char(DataTable.participant(subj))]);
    EEG.chanlocs = chanlocs;
    EEG.chaninfo = chaninfo;
    % saving StimTrak data
    EEG4 = pop_select( EEG, 'channel',{'Stimtrack'});
    pop_saveset( EEG4, 'filename','Stimtrack_B4.set','filepath',[dataDir char(DataTable.participant(subj))]);
    clear EEG4;
    % removing StimTrak and Trigger channels
    EEG = pop_select( EEG, 'rmchannel',{'Stimtrack','Trig'});
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

    % merging the datasets
    EEG = pop_mergeset( ALLEEG, [1  2  3  4], 0);

    % data recorded at the CBU has a lot of line noise
    if subj > 13
        EEG = pop_eegfiltnew(EEG,'hicutoff',40,'plotfreqz',0);
    end

    % Cleaning data using ASR
    [EEG_clean,~,~,removed_channels] = clean_artifacts(EEG);
    pop_saveset( EEG_clean, 'filename','cleanedEEG.set','filepath',[dataDir char(DataTable.participant(subj))]);
    save([dataDir char(DataTable.participant(subj)) '\removed_channels.mat'],'removed_channels');
    clear removed_channels;
    vis_artifacts(EEG_clean,EEG);

    EEG_ica = pop_runica( EEG_clean );
    pop_saveset( EEG_ica, 'filename','ica_EEG.set','filepath',[dataDir char(DataTable.participant(subj))]);

    clear ALLEEG EEG CURRENTSET ALLCOM EEG_ica EEG_clean;
    close all;

end

%%    Pre-processing step 2 dipole fitting
clear;
eeglab;
clear ALLEEG EEG CURRENTSET ALLCOM CURRENTSTUDY globalvars LASTCOM PLUGINLIST STUDY;
close all;

dataDir = 'D:\Tecnalia_data\';
cd(dataDir);
load Table_with_all_filenames;

% Variables for Dipole Fitting (Adding Folder to the Path)
% Check if the dipfit plugin is located in a folder named 'dipfit'
% or a different folder, such as 'dipfit5.0'.
path_dipfit = what('dipfit');
path_dipfit = path_dipfit.path;
addpath(genpath(path_dipfit));
templateChannelFilePath       = [fileparts(which('standard_1005.elc')) filesep 'standard_1005.elc'];
hdmFilePath                   = [fileparts(which('standard_vol.mat')) filesep 'standard_vol.mat'];
MRIfile                       = [fileparts(which('standard_mri.mat')) filesep 'standard_mri.mat'];
% This was calculated using the eeglab GUI
coordinateTransformParameters = [-2.8846e-06 -9.4577e-06 8.3647e-07 1.5183e-08 2.5437e-07 -1.5708 1 1 1];

% loop across the participants (16)
for subj=1:size(DataTable,1)

    % loading ica data
    cd([dataDir char(DataTable.participant(subj))]);
    ica_file = dir('ica_*'); % (.set and .fdt files are expected)
    EEG = pop_loadset('filename', ica_file(2).name,'filepath',[dataDir char(DataTable.participant(subj))]);

    % automatic coordinate transformation (deprecated due to issues)
    %     [~,coordinateTransformParameters] = coregister(EEG.chanlocs, templateChannelFilePath, 'warp', 'auto', 'manual', 'off');

    EEG = pop_dipfit_settings( EEG, 'hdmfile',hdmFilePath,'mrifile',MRIfile, ...
        'chanfile',templateChannelFilePath,'coordformat','MNI', ...
        'coord_transform',coordinateTransformParameters ,'chansel',1:EEG.nbchan );

    EEG = pop_multifit(EEG, 1:EEG.nbchan,'threshold', 100, 'dipplot','off','plotopt',{'normlen' 'on'}); % 'dipplot','on',

    % labelling the components
    EEG = iclabel(EEG, 'default');

    % Inspecting the components
    pop_viewprops(EEG,0);
    pause; close;

    % The following threshold will select components if they fall within the
    % categories of eye or muscle with a specified level of confidence.
    threshold = [NaN NaN;0.48 1;0.9 1;0.8 1;0.8 1;0.8 1;0.5 1];
    EEG = pop_icflag(EEG, threshold);

    % Manually input the components that need to be rejected.
    EEG = pop_subcomp( EEG );

    % Post-process to update EEG.icaact.
    EEG.icaact = [];
    EEG = eeg_checkset(EEG, 'ica');

    close all;
    % sanity check
    formatSpec = 'Subject %s with %s channels and %s components\n';
    fprintf(formatSpec,EEG.filename,num2str(EEG.nbchan),num2str(size(EEG.icaact,1)));

    pop_saveset( EEG, 'filename','final_EEG.set','filepath',[dataDir char(DataTable.participant(subj))]);

end % end loop participants

%% Organising digit audio signals
% % (Commented out to prevent code execution redundancy).

% % Following comment is Not necessary for the current operations but it is
% % included here for organizational purposes:
% % 'debruijn_seq3.mat' holds the sequence of digits to be presented in the
% % experiment. This presentation order was generated using a DeBruijn algorithm.

% % Similarly, following code is a bonus in case it is needed
% % debruijn_seq = load('debruijn_seq3.mat','sequence');
% % debruijn_seq = debruijn_seq.sequence;
%
% % Spanish_digits = {'0_cero_es_es_1.mp3', '1_uno_es_es_1.mp3', '2_dos_es_es_1.mp3', ...
% %     '3_tres_es_es_1.mp3', '4_cuatro_es_es_1.mp3', '5_cinco_es_es_1.mp3', ...
% %     '6_seis_es_es_1.mp3', '7_siete_es_es_1.mp3', '8_ocho_es_es_1.mp3', ...
% %     '9_nueve_es_es_1.mp3'};
%
% % Italian_digits = {'0_cero_it_norm.wav','1_uno_it_norm.wav','2_dos_it_norm.wav', ...
% %     '3_tres_it_norm.wav','4_cuatro_it_norm.wav','5_cinco_it_norm.wav', ...
% %     '6_seis_it_norm.wav','7_siete_it_norm.wav','8_ocho_it_norm.wav', ...
% %     '9_nueve_it_norm.wav'};

% % Loading the 'audio_matrix.mat' variable, which contains the audio data
% % as utilized in the experiment, and subsequently resampling it to match
% % the EEG sampling rate of 1000Hz.
% cd('D:\TECNALIA_2021\Hyperscanning_the same as BCI2000 program(Closed Loop Tecnalia)\Codes\Audios\');
% spanish_audios = load('audio_matrix.mat');
% Fs_spanish = spanish_audios.Fs;
% spanish_audios = spanish_audios.all_audios;
% for i=1:10
% temp = spanish_audios{i};
% y = resample(temp,1000,22050);
% spanish_audios{i} = y;
% clear temp y;
% end
%
% % Loading the 'audio_matrix_italian.mat' variable, which contains the
% % audio data in italian that were used for two participants in the experiment:
% % subj04 and subj16
% % It is the Same as loading the var 'audio_matrix.mat' before.
% cd('D:\TECNALIA_2021\Hyperscanning_the same as BCI2000 program(Closed Loop Tecnalia)\Codes\Audios\ita\');
% italian_audios = load('audio_matrix_italian.mat');
% Fs_italian = italian_audios.Fs;
% italian_audios = italian_audios.all_audios;
% for i=1:10
% temp = italian_audios{i};
% y = resample(temp,1000,22050);
% italian_audios{i} = y;
% clear temp y;
% end
%
% % saving the audio data in the working directory
% save('D:\Tecnalia_data\audio_matrixes_1000Hz.mat','italian_audios','spanish_audios');

% % Copying the volume vector data used on each participant, from
% % TECNALIA_2021 to the Tecnalia_data folder.
% % This data will be used to insert the real audio data presented,
% % including the amplitude (volume) in the StimTrak files.
% data_origen = 'D:\TECNALIA_2021\';
% data_destination = 'D:\Tecnalia_data\';
% load('D:\Tecnalia_data\Table_with_all_filenames.mat');
%
% for i=1:size(DataTable,1)
%     cd([data_origen char(DataTable.participant(i)) '\Calibration']);
%     V = dir('Volume_*');
%     volume = V.name;
%     copyfile(volume, [data_destination char(DataTable.participant(i))]);
%     clear V volume;
% end
% clear data_origen data_destination;

%% Working with the audio data (volumes and digits)

% Processing the StimTrak recording for each block to generate an EEGlab-format signal
% that encapsulates the actual audio data at the specific volume levels it was presented.

clear;
close all;

% data path
dataDir = 'D:\Tecnalia_data\';
cd(dataDir);
load audio_matrixes_1000Hz.mat;

% working with one participant at a time. Selecting the folder.
selpath = uigetdir(dataDir,'Select participant folder');

% loading data
cd(selpath);

% loading the volume file which contains the four different amplitude values
% corresponding to eac participant: minimun audio, intermediate 1,
% intermediate 2 and maximun audio.
V = dir('Volume_*');
volume = V.name;
load(volume); % loads the var loud_vector
clear V volume dataDir;
close all;

% This code lists the SavedData files, each corresponding to one of the
% four blocks. These files contain various information, including data on
% the sequence of digits presented and the associated volume levels.
% SavedData = dir('SavedData_2021*');
SavedData = dir('SavedData_2022*'); % for the last 3 participants

for block = 1:4 % loop across blocks

    % loading the StimTrak
    EEG = pop_loadset('filename',['Stimtrack_B' num2str(block) '.set'],'filepath',selpath);

    % Removing all the unnecesary event markers (this time on StimTrak data)
    % First we keep only events type StimulusCode
    allEventTypes = {EEG.event.type}';
    strFindResult = strfind(allEventTypes, 'StimulusCode');
    probeIdx      = find(cellfun(@isempty, strFindResult));
    EEG.event(probeIdx) = [];
    % Second we keep only event with number 3 in the position field
    allEventTypes = {EEG.event.position}';
    allEventTypes = cell2mat(allEventTypes);
    probeIdx = find(allEventTypes~=4);
    EEG.event(probeIdx) = [];

    % updating the fields EEG.event.urevent and EEG.urevent
    for i= 1: size(EEG.event,2)
        EEG.event(i).urevent = i;
    end
    EEG.urevent = EEG.event;
    EEG.urevent = rmfield(EEG.urevent,'urevent');
    clear allEventTypes strFindResult probeIdx i;

    % we assume listed order corresponds to block order dur to the naming convention date-hour
    load(SavedData(block).name); % loads Vector_vol, Vector_digit, etc.

    % Conditional handling for cases where Vector_vol Is a Vector (Due to an Error)
    % The variable Vector_digit will be utilized to resize Vector_vol to the
    % correct dimensions. Vector_digit was previously resized in a manual step.

    if size(Vector_vol,1)==1
        % the error is because one element is missing and it was impossible
        % to create a matrix. We add a 999 at the end to resize the vector
        Vector_vol(end+1) = 999;
        [row,col] = size(Vector_digit);
        Vector_vol = reshape(Vector_vol,col,row); % reshape operates first in the rows
        Vector_vol = Vector_vol'; % transposing the matrix
        clear row col;
    end

    % We'll start by working with the 'Vector_vol' variable from the SavedData file.
    % 'Vector_vol' contains the volume levels presented to the participants.
    % In this section of the script, we will extensively debug 'Vector_vol' as
    % it plays a crucial role. One issue is that it appears to miss the first
    % volume value in the sequence and includes a last value that is meaningless.
    % Mathematically, it seems that each digit sequence is contained from rowN(end)
    % to rowN+1(end-1). We will address this issue below with a solution.

    % Original matrix size (for reshaping back)
    originalSize = size(Vector_vol);
    % Transpose the matrix and then convert it into a vector by concatenating rows
    Vector_vol = reshape(Vector_vol',[],1);
    % Adding a zero in the beggining
    Vector_vol = cat(1, 0, Vector_vol);
    % Removing the last number in the data which seems useless
    Vector_vol(end) = [];
    % Reshape the vector back to the original matrix size
    Vector_vol_new = reshape(Vector_vol, originalSize);

    % The second var within SavedData file to work with is 'Vector_digit'
    % Vector_digit contains the digits presented to the participants
    % Preparing the Vector of digits
    Vector_digit = reshape(Vector_digit',[],1);

    % preallocating var to savetiming info of the shifts
    latencies = nan(length(EEG.event),1);

    % Inspecting one by one the correct association between StimTrak signal
    % and the corresonding audio (digit/amplitude)
    for i=2:length(EEG.event) % loop across trials. Starting in 2 since we don't have correct loudness info for the first in the block
        t1 = EEG.event(i).latency; % onset latency of each trial
        t2 = t1 + 1600; % offset latency
        trial  = EEG.data(t1:t2); % data chunk of the StimTrak data in eeglab format
        digit = Vector_digit(i) + 1; % digit presented in the corresponding trial. Since it goes from 0-9, we must add 1.
        %         audio = spanish_audios{digit}; % the audio file corresponding to the digit
        audio = italian_audios{digit}; % for the cases 4 and 16, in italian
        % Calculating the cross correlation between the audio file and the Stimtrak signal
        % these two signals should be highly correlated in -500 lag and therefore showing a clear peak
        [r,lags] = xcorr(audio, trial);
        loudness_idx = Vector_vol(i); % index of the amplitude or loudness for this specific trial (1,2,3 or 4).
        loudness_value = loud_vector(loudness_idx); % loudness values corresponding to the index (e.g., 0.0658    0.3511    1.8738   10.0000)

        % Incorporating Loudness Value and Audio Digit Information into EEG.event.type
        % Format: Loud_x_dig_y, where 'x' represents loudness (1-4) and 'y' represents
        % the digit (0-9). The event type will be used for epoching purposes.
        EEG.event(i).type = ['Loud_' num2str(loudness_idx) '_dig_' num2str(digit-1)];
        % finding the exact lag with the larger correlation value and
        % therefore the better alignment
        new_r = zeros(size(r));
        % We restrict our search to this interval as a precaution for cases where
        % the audio is unclear in the StimTrak signal, resulting in very low correlation values.
        % The expected alignment should occur around lag -500 for reliable results.
        %         new_r(1001:1200) = r(1001:1200); % lags -600 to -400 in the cross correlation
        new_r(751:1100) = r(751:1100); % lags -850 to -500 in the cross correlation. This is for the last three cases.
        % finding index of larger r value
        larger_r = max(abs(new_r));
        idx = find(abs(new_r)==larger_r);
        shift = lags(idx);
        % printing on the screen to inspect the values
        disp(shift);
        latencies(i) = shift;
        % Scaling the audio data by multiplying it with the loudness value, reflecting
        % how the audio was presented to the participant.
        audio_with_loudness = audio*loudness_value;
        % getting the length of the audio
        audio_length = size(audio,1);
        % creating a new zeroed trial to include the real audio
        new_trial = zeros(size(trial));
        % including the audio exactly where it happened
        new_trial(abs(shift):audio_length+abs(shift)-1) = audio_with_loudness;
        % to solve a length issue on participant 16 (unclear cause)
        if size(new_trial,2) > size(trial,2)
            new_trial = new_trial(1,1:size(trial,2));
        end

        % Visualisation
        figure('Name', ['Trial ' num2str(i) ' Volume ' num2str(loudness_idx)],'NumberTitle','off');
        subplot(3,1,1);
        stem(lags,r); % showing cross correlation vales for the different lags
        subplot(3,1,2);
        plot(trial); % plotting the StimTrak signal
        subplot(3,1,3);
        plot(new_trial,'k');

        EEG.data(1,t1:t2) = new_trial;

    end % loop across trials

    pause;
    pop_saveset( EEG, 'filename',['StimTrack_corrected_B' num2str(block) '.set'],'filepath',selpath);
    save(['latencies_audios_B' num2str(block) '.mat'],'latencies');
    close all;

end % loop blocks

%% next step (commented after run once)

% % The following information lists the trials to be removed for each
% % participant and block. These trials are being removed due to our
% % inability to align the real audio data with the corresponding part of the
% % StimTrak channel. This alignment issue arose because the audio signal was
% % obscured by baseline noise within the StimTrak recording. Additionally,
% % the first trial of each block is being removed since it lacks
% % amplitude information for the audio. The identification of trials where
% % reliable alignment was not possible was determined through visual
% % inspection, as described in the previous section, and was annotated manually.

%
% subject_1 = struct();
% subject_1.block1 = [85,49,47,10];
% subject_1.block2 = [112,93,91,82,49,11,7,5];
% subject_1.block3 = [91,89,88,38,19,17,4,2];
% subject_1.block4 = [78,70,68,67,65,62,60,46,6,4];
%
% subject_2 = struct();
% subject_2.block1 = [136,89,83,64,29,23,22,20,18,16];
% subject_2.block2 = [133,129,127,78,74,67,56,49,48,45,20,19,18,17,16,15];
% subject_2.block3 = [125,121,112,110,103,97,73,72];
% subject_2.block4 = [140,136,135,120,114,97,74,71,64,61,60,54];
%
% subject_3 = struct();
% subject_3.block1 = [131,112,101,97,87,22,17,16,10,8,4];
% subject_3.block2 = [112,109,98,86,76,43,42,33,12];
% subject_3.block3 = [125,79,77,67,59,57,31,30,29];
% subject_3.block4 = [121,116,106,104,100,99,69,55,12,4];
%
% subject_4 = struct();
% subject_4.block1 = [78,74,73,72,71,70,49,47,43,40,39];
% subject_4.block2 = [128,125,122,113,106,96,62,59,36,20,19,17,9];
% subject_4.block3 = [82,78,77,75,72,63,61,60,59,56,54,52,33];
% subject_4.block4 = [119,115,99,98,49,43,41,36,5,4];
%
% subject_5 = struct();
% subject_5.block1 = [95,94,93,92,87,85,78,74,68,60,57,54,53,50,32];
% subject_5.block2 = [131,124,120,119,118,116,112,103,98,91,48,44,43,42,41,40,39,38,36,13,11];
% subject_5.block3 = [132,81,80,72,68,67,66,64,63,61,60,51,49,48];
% subject_5.block4 = [119,118,117,116,110,94,88,72,36,35,29,21,12,10,6,5,4,2];
%
% subject_6 = struct();
% subject_6.block1 = [66,53,2];
% subject_6.block2 = [21,20];
% subject_6.block3 = [65,61];
% subject_6.block4 = [126,59,44];
%
% subject_7 = struct();
% subject_7.block1 = [70,57,53,50,47,46,19];
% subject_7.block2 = [103,100,76,72,69,66,55,53,50,25,2];
% subject_7.block3 = [130,106,88,82,75,67,64,10];
% subject_7.block4 = [132,130,129,119,106,89,85,84,80,66,63,58,55,51,27];
%
% subject_8 = struct();
% subject_8.block1 = [];
% subject_8.block2 = [107,28];
% subject_8.block3 = [];
% subject_8.block4 = [120,77,60,39];
%
% subject_9 = struct();
% subject_9.block1 = [119,109,107,106,105,95,94,93,83,80,66,27];
% subject_9.block2 = [87,85,73,72,67,66,65,63,55,51,49,43,40];
% subject_9.block3 = [130,129,67,66,64,60,57,46,45,34,30,19,15];
% subject_9.block4 = [139,108,107,106,105,87,41,39,38,25,5,4];
%
% subject_10 = struct();
% subject_10.block1 = 28;
% subject_10.block2 = 36;
% subject_10.block3 = 91;
% subject_10.block4 = 29;
%
% subject_11 = struct();
% subject_11.block1 = [106,104,100,73,68,64,63,62,53,41,39];
% subject_11.block2 = [106,105,95,78,77,76,75,65,56,55,17,9,8,7,6,4,3];
% subject_11.block3 = [119,110,109,107,81,36,35,34,33,32,21,17,11,6,2];
% subject_11.block4 = [130,129,116,97,74,69,64,62,60,58,47,45,20,11,10];
%
% subject_12 = struct();
% subject_12.block1 = [42,41,40,27];
% subject_12.block2 = [76,6];
% subject_12.block3 = [34,35,24];
% subject_12.block4 = [100,69];
%
% subject_13 = struct();
% subject_13.block1 = [140,133,114,85,10];
% subject_13.block2 = [137,116,101,61,60,44,24];
% subject_13.block3 = [133,92,61,7];
% subject_13.block4 = [137,113,81,77,73,60,44,5];
%
% subject_14 = struct();
% subject_14.block1 = [126,124,115,114,113,112,95,91,89,87,68,63,56,52,34,33,30,29,17,15,5];
% subject_14.block2 = [128,127,121,114,109,82,78,76,75,61,57,46,42,31,20,15,14,13,10,4];
% subject_14.block3 = [130,112,104,101,91,78,72,30,29,27,17,13,9];
% subject_14.block4 = [122,113,112,97,81,79,40,37,36,25,23,22];
%
% subject_15 = struct();
% subject_15.block1 = [168,143,131,124,123,115,93,88,81,67,66,53,50,46,44,32,28,21,10];
% subject_15.block2 = [175,171,166,151,146,126,122,119,77,71,68,64,60,59,55,4539,34,20,18,17,15,12,8];
% subject_15.block3 = [175,165,164,160,156,136,130,117,114,113,104,95,86,76,74,67,57,52,50,42,41,16,12];
% subject_15.block4 = [170,165,161,155,154,141,139,130,125,119,106,104,77,70,68,50,48,34,21,13];
%
% subject_16 = struct();
% subject_16.block1 = [];
% subject_16.block2 = [];
% subject_16.block3 = [];
% subject_16.block4 = [];
%
% dataDir = 'D:\Tecnalia_data\';
% save([dataDir 'rejected_trials_per_block_allsubj.mat'], "subject_16", ...
%     "subject_15","subject_14","subject_13","subject_12","subject_11", ...
%     "subject_10","subject_9","subject_8","subject_7","subject_6", ...
%     "subject_5","subject_4","subject_3","subject_2","subject_1");

%% Transferring Information from StimTrak to EEG Files

clear;
close all;
dataDir = 'D:\Tecnalia_data\';
cd(dataDir);
load Table_with_all_filenames;

for subj=1:size(DataTable,1)

    selpath = [dataDir char(DataTable.participant(subj))];
    clear ALLEEG EEG CURRENTSET ALLCOM CURRENTSTUDY globalvars LASTCOM PLUGINLIST STUDY;
    cd(selpath);
    eeglab;

    % Working with the StimTrack_corrected files. They contain the real audio
    % data, and the  Stimulus Code information for each trial
    % (e.g., 'Loud_3_dig_4', 'Loud_1_dig_5', or 'Loud_2_dig_3').
    % We concatenate all StimulusCodes across blocks. Additionally,
    % Loading and storing data from the four blocks
    StimulusCode = [];
    for block=1:4
        EEG = pop_loadset('filename', ['StimTrack_corrected_B' num2str(block) '.set']);
        eventTypes = {EEG.event.type};
        StimulusCode = cat(2,StimulusCode, eventTypes);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    end % loop block
    clear eventTypes;

    % merging the datasets in one file
    EEG = pop_mergeset( ALLEEG, [1  2  3  4], 0);

    % For subj6, some data points were removed during the clean_artifacts step.
    % The issue was identified using the latency component of the sanity check.
    % To determine the missing data points, you can access them by evaluating
    % find(EEG.etc.clean_sample_mask == 0) in the 'cleanedEEG.set'.
    % In this case, the missing points range from 148921 to 149920.
    % Based on this information, we have implemented the following condition:

    if subj==6
        EEG = pop_select( EEG, 'rmpoint',[148921 149920] );
    end

    % removing the 'boundary' markers from the concatenated StimTrak
    % recordings
    allEventTypes = {EEG.event.type}';
    strFindResult = strfind(allEventTypes, 'boundary');
    probeIdx      = find(~cellfun(@isempty,strFindResult));
    EEG.event(probeIdx) = [];
    EEG.urevent(probeIdx) = [];

    % loading the processed EEG file for the participant
    final_EEG = pop_loadset('filename', 'final_EEG.set','filepath',[dataDir char(DataTable.participant(subj))]);
    % Removing all the unnecesary event markers
    % First we keep only events type StimulusCode
    allEventTypes = {final_EEG.event.type}';
    strFindResult = strfind(allEventTypes, 'StimulusCode');
    probeIdx      = find(cellfun(@isempty, strFindResult));
    final_EEG.event(probeIdx) = [];
    % Second we keep only event with number 3 in the position field
    allEventTypes = {final_EEG.event.position}';
    allEventTypes = cell2mat(allEventTypes);
    probeIdx = find(allEventTypes~=4);
    final_EEG.event(probeIdx) = [];

    % updating the fields EEG.event.urevent and EEG.urevent
    for i = 1: size(final_EEG.event,2)
        final_EEG.event(i).urevent = i;
    end
    final_EEG.urevent = final_EEG.event;
    final_EEG.urevent = rmfield(final_EEG.urevent,'urevent');
    clear allEventTypes strFindResult probeIdx i;

    % Sanity check about the amount of trials on the concatenated Stimtrak
    % files and in the EEG
    if length(StimulusCode) ~=  length(final_EEG.event)
        warningMessage = ['subject ' num2str(subj) ' has ' num2str(length(StimulusCode)) ' trials in the StimTrak but ' num2str(length(final_EEG.event)) ' in the EEG file'];
        warndlg(warningMessage, 'Warning');
    end

    % The StimulusCode information from the StimTrack_corrected files will
    % be transferred to the final_EEG files.
    % Assuming StimulusCode is a cell and final_EEG.event is a struct array
    % Convert StimulusCode cell array to a cell array of strings
    StimulusCode = cellfun(@char, StimulusCode, 'UniformOutput', false);
    % Assign StimulusCode to final_EEG.event.type
    [final_EEG.event.type] = StimulusCode{:};
    clear StimulusCode;

    % More sanity check, this time about the amount of trials on the
    % concatenated Stimtrak files and in the EEG
    if length(EEG.event) ~=  length(final_EEG.event)
        warningMessage = ['subject ' num2str(subj) ' has ' num2str(length(EEG.event)) ' trials in the StimTrak but ' num2str(length(final_EEG.event)) ' in the EEG file'];
        warndlg(warningMessage, 'Warning');
    end
    if ~isempty(find(cell2mat({EEG.event.latency}') - cell2mat({final_EEG.event.latency}')))
        warningMessage = ['subject ' num2str(subj) ' has different latencies for trials in the StimTrak than in the EEG file'];
        warndlg(warningMessage, 'Warning');
    end

    pop_saveset( EEG, 'filename','StimTrack_final.set','filepath',selpath);
    pop_saveset( final_EEG, 'filename','final_EEG4gcmi.set','filepath',selpath);

end % loop subject

%% Removing trials from the EEG and StimTrak
% This step involves further processing the StimTrak and EEG data by
% rejecting trials where perfect alignment of the audio signal was
% not possible. Information regarding rejected trials is stored in the
% 'rejected_trials_per_block_allsubj.mat' variable. Additionally, the
% speech envelope will be estimated.

clear;
close all;
dataDir = 'D:\Tecnalia_data\';
cd(dataDir);
load Table_with_all_filenames;
load rejected_trials_per_block_allsubj;

% Initializing a table to store information about number of trials for each
% condition (loud1, loud2, ...) and file type (EEG & StimTrak) across participants.
sz = [16 9];
varTypes = ["string","double","double","double","double","double","double","double","double"];
varNames = ["Participant", ...
    "EEG_loud1","StimTrak_loud1","EEG_loud2","StimTrak_loud2","EEG_loud3", ...
    "StimTrak_loud3","EEG_loud4","StimTrak_loud4"];
trials_DataTable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
clear varTypes varNames sz;

for subj=1:size(DataTable,1)
    selpath = [dataDir char(DataTable.participant(subj))];
    cd(selpath);
    eeglab;

    % Assigning the participant's codename
    trials_DataTable.Participant(subj) = char(DataTable.participant(subj));

    if subj<14
        SavedData = dir('SavedData_2021*');
    elseif subj>=14
        SavedData = dir('SavedData_2022*'); % for the last 3 participants
    end

    % Loading StimTrak files with the audios
    EEG = pop_loadset('filename', 'StimTrack_final.set');
    % calculating the speech envelope using the function from mTRF-Toolbox
    % https://github.com/mickcrosse/mTRF-Toolbox.git
    EEG.data = mTRFenvelope(EEG.data,EEG.srate,EEG.srate);
    EEG.data = EEG.data';

    bad_trials = eval(['subject_' (num2str(subj))]);

    % Initialising variables
    trials2rej = [];
    cumtrials = 0;

    % This code processes all the trials for rejection based on their block
    % structure and converts them for use in continuous recording. The trial
    % numbering is adjusted to consider the respective block they belong to.
    for block=1:4
        % We assume that the listed order corresponds to the block order
        % based on the naming convention of date and hour.
        load(SavedData(block).name); % loads Vector_vol, Vector_digit, etc.

        [rep,span]=size(Vector_digit);
        if span==0
            warningMessage = ['subject ' num2str(subj) ' has a mistake in block ' num2str(i) ];
            warndlg(warningMessage, 'Warning');
        end

        % calculating the total number of digits presented on a block
        trialblock = rep * span;

        % Getting the bad trials corresponding to each block
        %         temp = eval(['bad_trials.block' num2str(block)]); % old code
        temp = bad_trials. ("block" + num2str(block)); % improved version
        % Adding the first trial of each block for rejection (and sorting the order for clarity)
        temp = [1 sort(temp)];

        % for the first block trial order is the same
        if block==1
            trials2rej = temp;
        else
            cumtrials = cumtrials + trialblock;
            temp = temp + cumtrials;
            trials2rej = cat(2,trials2rej,temp);
        end
        clear temp trialblock rep span;
    end % block loop

    save('trials_to_reject_continous_format.mat','trials2rej');

    % creating a folder for saving the files to be directly loaded for the
    % gcmi analyses.
    mkdir('gcmi_folder');

    %%%% Epoching StimTrak data, rejecting trials and epoching based on loudness condition %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EEG = pop_epoch( EEG, { 'Loud_1_dig_0' 'Loud_1_dig_1' 'Loud_1_dig_2' 'Loud_1_dig_3' 'Loud_1_dig_4' 'Loud_1_dig_5' 'Loud_1_dig_6'...
        'Loud_1_dig_7' 'Loud_1_dig_8' 'Loud_1_dig_9' 'Loud_2_dig_0' 'Loud_2_dig_1' 'Loud_2_dig_2' 'Loud_2_dig_3' 'Loud_2_dig_4'...
        'Loud_2_dig_5' 'Loud_2_dig_6' 'Loud_2_dig_7' 'Loud_2_dig_8' 'Loud_2_dig_9' 'Loud_3_dig_0' 'Loud_3_dig_1' 'Loud_3_dig_2'...
        'Loud_3_dig_3' 'Loud_3_dig_4' 'Loud_3_dig_5' 'Loud_3_dig_6' 'Loud_3_dig_7' 'Loud_3_dig_8' 'Loud_3_dig_9' 'Loud_4_dig_0'...
        'Loud_4_dig_1' 'Loud_4_dig_2' 'Loud_4_dig_3' 'Loud_4_dig_4' 'Loud_4_dig_5' 'Loud_4_dig_6' 'Loud_4_dig_7' 'Loud_4_dig_8'...
        'Loud_4_dig_9' 'StimulusCode' }, [0    1.6], 'newname', 'Epoched', 'epochinfo', 'yes');
    EEG = pop_selectevent( EEG, 'omitepoch',trials2rej ,'deleteevents','off','deleteepochs','on','invertepochs','off');
    pop_saveset( EEG, 'filename','StimTrack_envelope.set','filepath',selpath);

    StimTrak_loud1 =  pop_epoch( EEG, { 'Loud_1_dig_0' 'Loud_1_dig_1' 'Loud_1_dig_2' 'Loud_1_dig_3' 'Loud_1_dig_4' 'Loud_1_dig_5' 'Loud_1_dig_6'...
        'Loud_1_dig_7' 'Loud_1_dig_8' 'Loud_1_dig_9' }, [0    1.6], 'newname', 'Epochs loud1', 'epochinfo', 'yes');
    pop_saveset( StimTrak_loud1, 'filename','StimTrak_loud1.set','filepath',[selpath filesep 'gcmi_folder']);
    trials_DataTable.StimTrak_loud1(subj) = StimTrak_loud1.trials;

    StimTrak_loud2 =  pop_epoch( EEG, { 'Loud_2_dig_0' 'Loud_2_dig_1' 'Loud_2_dig_2' 'Loud_2_dig_3' 'Loud_2_dig_4' 'Loud_2_dig_5' 'Loud_2_dig_6'...
        'Loud_2_dig_7' 'Loud_2_dig_8' 'Loud_2_dig_9' }, [0    1.6], 'newname', 'Epochs loud2', 'epochinfo', 'yes');
    pop_saveset( StimTrak_loud2, 'filename','StimTrak_loud2.set','filepath',[selpath filesep 'gcmi_folder']);
    trials_DataTable.StimTrak_loud2(subj) = StimTrak_loud2.trials;

    StimTrak_loud3 =  pop_epoch( EEG, { 'Loud_3_dig_0' 'Loud_3_dig_1' 'Loud_3_dig_2' 'Loud_3_dig_3' 'Loud_3_dig_4' 'Loud_3_dig_5' 'Loud_3_dig_6'...
        'Loud_3_dig_7' 'Loud_3_dig_8' 'Loud_3_dig_9' }, [0    1.6], 'newname', 'Epochs loud3', 'epochinfo', 'yes');
    pop_saveset( StimTrak_loud3, 'filename','StimTrak_loud3.set','filepath',[selpath filesep 'gcmi_folder']);
    trials_DataTable.StimTrak_loud3(subj) = StimTrak_loud3.trials;

    StimTrak_loud4 =  pop_epoch( EEG, { 'Loud_4_dig_0' 'Loud_4_dig_1' 'Loud_4_dig_2' 'Loud_4_dig_3' 'Loud_4_dig_4' 'Loud_4_dig_5' 'Loud_4_dig_6'...
        'Loud_4_dig_7' 'Loud_4_dig_8' 'Loud_4_dig_9' }, [0    1.6], 'newname', 'Epochs loud4', 'epochinfo', 'yes');
    pop_saveset( StimTrak_loud4, 'filename','StimTrak_loud4.set','filepath',[selpath filesep 'gcmi_folder']);
    trials_DataTable.StimTrak_loud4(subj) = StimTrak_loud4.trials;

    clear EEG StimTrak_loud*;

    %%%% Epoching EEG data, rejecting trials and epoching based on loudness condition %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EEG = pop_loadset('filename','final_EEG4gcmi.set');
    EEG = pop_epoch( EEG, { 'Loud_1_dig_0' 'Loud_1_dig_1' 'Loud_1_dig_2' 'Loud_1_dig_3' 'Loud_1_dig_4' 'Loud_1_dig_5' 'Loud_1_dig_6'...
        'Loud_1_dig_7' 'Loud_1_dig_8' 'Loud_1_dig_9' 'Loud_2_dig_0' 'Loud_2_dig_1' 'Loud_2_dig_2' 'Loud_2_dig_3' 'Loud_2_dig_4'...
        'Loud_2_dig_5' 'Loud_2_dig_6' 'Loud_2_dig_7' 'Loud_2_dig_8' 'Loud_2_dig_9' 'Loud_3_dig_0' 'Loud_3_dig_1' 'Loud_3_dig_2'...
        'Loud_3_dig_3' 'Loud_3_dig_4' 'Loud_3_dig_5' 'Loud_3_dig_6' 'Loud_3_dig_7' 'Loud_3_dig_8' 'Loud_3_dig_9' 'Loud_4_dig_0'...
        'Loud_4_dig_1' 'Loud_4_dig_2' 'Loud_4_dig_3' 'Loud_4_dig_4' 'Loud_4_dig_5' 'Loud_4_dig_6' 'Loud_4_dig_7' 'Loud_4_dig_8'...
        'Loud_4_dig_9' 'StimulusCode' }, [0    1.6], 'newname', 'Epoched', 'epochinfo', 'yes');
    EEG = pop_selectevent( EEG, 'omitepoch',trials2rej ,'deleteevents','off','deleteepochs','on','invertepochs','off');
    pop_saveset( EEG, 'filename','final_EEG4gcmi_epoched.set','filepath',selpath);

    EEG_loud1 =  pop_epoch( EEG, { 'Loud_1_dig_0' 'Loud_1_dig_1' 'Loud_1_dig_2' 'Loud_1_dig_3' 'Loud_1_dig_4' 'Loud_1_dig_5' 'Loud_1_dig_6'...
        'Loud_1_dig_7' 'Loud_1_dig_8' 'Loud_1_dig_9' }, [0    1.6], 'newname', 'Epochs loud1', 'epochinfo', 'yes');
    pop_saveset( EEG_loud1, 'filename','EEG_loud1.set','filepath',[selpath filesep 'gcmi_folder']);
    trials_DataTable.EEG_loud1(subj) = EEG_loud1.trials;

    EEG_loud2 =  pop_epoch( EEG, { 'Loud_2_dig_0' 'Loud_2_dig_1' 'Loud_2_dig_2' 'Loud_2_dig_3' 'Loud_2_dig_4' 'Loud_2_dig_5' 'Loud_2_dig_6'...
        'Loud_2_dig_7' 'Loud_2_dig_8' 'Loud_2_dig_9' }, [0    1.6], 'newname', 'Epochs loud2', 'epochinfo', 'yes');
    pop_saveset( EEG_loud2, 'filename','EEG_loud2.set','filepath',[selpath filesep 'gcmi_folder']);
    trials_DataTable.EEG_loud2(subj) = EEG_loud2.trials;

    EEG_loud3 =  pop_epoch( EEG, { 'Loud_3_dig_0' 'Loud_3_dig_1' 'Loud_3_dig_2' 'Loud_3_dig_3' 'Loud_3_dig_4' 'Loud_3_dig_5' 'Loud_3_dig_6'...
        'Loud_3_dig_7' 'Loud_3_dig_8' 'Loud_3_dig_9' }, [0    1.6], 'newname', 'Epochs loud3', 'epochinfo', 'yes');
    pop_saveset( EEG_loud3, 'filename','EEG_loud3.set','filepath',[selpath filesep 'gcmi_folder']);
    trials_DataTable.EEG_loud3(subj) = EEG_loud3.trials;

    EEG_loud4 =  pop_epoch( EEG, { 'Loud_4_dig_0' 'Loud_4_dig_1' 'Loud_4_dig_2' 'Loud_4_dig_3' 'Loud_4_dig_4' 'Loud_4_dig_5' 'Loud_4_dig_6'...
        'Loud_4_dig_7' 'Loud_4_dig_8' 'Loud_4_dig_9' }, [0    1.6], 'newname', 'Epochs loud4', 'epochinfo', 'yes');
    pop_saveset( EEG_loud4, 'filename','EEG_loud4.set','filepath',[selpath filesep 'gcmi_folder']);
    trials_DataTable.EEG_loud4(subj) = EEG_loud4.trials;

    clear badtrials block selpath Vector* EEG* trials2rej clear ALLEEG ...
        CURRENTSET ALLCOM CURRENTSTUDY globalvars LASTCOM PLUGINLIST STUDY;

end % subj loop

save([dataDir 'Table_number_trials_across_conditions.mat'],'trials_DataTable');
