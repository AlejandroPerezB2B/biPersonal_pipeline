function biPer_cleaning_v2(folder_name)
% biPer_cleaning() % Using EEGLAB functions, performs the semiautomatic
%                    removal of non-neural activity from pre-processed EEG signal.
%                    It runs twice inside the folder, each time for a
%                    different participant (xxx_A.set and xxx_B.set) 
%%INPUT
% folder_name     [string] Full path to the folder containing the two
%                 EEG recordings e.g. 'C:\Users\jojo\Downloads\CNT01'
%%OUTPUT
%                 EEGLAB struct (.set) containing the pruned EEG signal
%
%   The following EEGLAB plugins are needed:
%   1.- fitTwoDipoles (v0.01)
%   2.- Fieldtrip-lite20200922
%   3.- ICLabel1.2.5
%   4.- Viewprops1.5.4
%   5.- automagic (v2.3.8)

% Author: Alejandro Perez
% version 1: University of Toronto, UTSC, June 7th 2018.
% version 2: MRC-CBU, University of Cambridge, March 21, 2020.

foldersep = strfind(folder_name,filesep);
name      = folder_name(foldersep(end)+1:end);
clear foldersep;

for subj=[1 2] % loop across the dyad
    switch(subj)
        case {1}
            AorB = 'A';
        case {2}
            AorB = 'B';
    end
    
    EEG = pop_loadset('filename',[folder_name filesep name '_' AorB '_ICA_dipfit.set']);
    
    % Search for and estimate symmetrically constrained bilateral dipoles
    EEG = fitTwoDipoles(EEG, 'LRR', 35);
    pop_saveset( EEG, 'filename', [folder_name filesep name '_' AorB '_ICA_dipfit_two.set']);
    
% %if we are estimating more than one model with AMICA 
% outdir = [ folder_name filesep 'AMICA_' AorB];
% modout = loadmodout15(outdir);
% % load individual ICA model into EEG structure
% model_index = 1;
% EEG.icawinv = modout.A(:,:,model_index);
% EEG.icaweights = modout.W(:,:,model_index);
% EEG.icasphere = modout.S;
% EEG = eeg_checkset(EEG);
% pop_topoplot(EEG,0); % plot IC scalp maps of ICA model #1
% % compute model probability
% model_prob = 10 .^ modout.v; % modout.v (#models x #samples)
% figure; imagesc(model_prob(:,1:10*EEG.srate)); % model probability changes in first 10s

    % labelling the components
    EEG = iclabel(EEG, 'default');
    
    % Finding ICLabel scores with >70% brain
    brainIdx  = find(EEG.etc.ic_classification.ICLabel.classifications(:,1) >= 0.7);
    
    % Finding residual variance of the IC scalp maps <15%
    rvList    = [EEG.dipfit.model.rv];
    goodRvIdx = find(rvList < 0.15)';
    
    % Finding dipoles located inside brain
    dipoleXyz = zeros(length(EEG.dipfit.model),3);
    for icIdx = 1:length(EEG.dipfit.model)
        dipoleXyz(icIdx,:) = EEG.dipfit.model(icIdx).posxyz(1,:);
    end
    load(EEG.dipfit.hdmfile); % loading the var 'vol'.
    depth = ft_sourcedepth(dipoleXyz, vol);
    depthThreshold = 1;
    insideBrainIdx = find(depth<=depthThreshold);
    
    % Taking the three criteria.
    goodIcIdx = intersect(brainIdx, goodRvIdx);
    goodIcIdx = intersect(goodIcIdx, insideBrainIdx);
    
    % Performing IC rejection.
    EEG = pop_subcomp(EEG, goodIcIdx, 0, 1);
    
    % Post-process to update ICLabel data structure.
    EEG.etc.ic_classification.ICLabel.classifications = EEG.etc.ic_classification.ICLabel.classifications(goodIcIdx,:);
    
    % Post-process to update EEG.icaact.
    EEG.icaact = [];
    EEG = eeg_checkset(EEG, 'ica');
    
    % Inspecting the components
    pop_viewprops(EEG,0);
    pause; close;
    
    %     % For manually introducing those components to be rejected
    %     EEG = pop_subcomp( EEG );
    
    pop_saveset( EEG, 'filename', [name '_' AorB '_ICA_dipfit_two_pruned1.set'],'filepath', folder_name);
    close all;
    formatSpec = 'Subject %s with %s channels and %s components\n';
    fprintf(formatSpec,EEG.filename,num2str(EEG.nbchan),num2str(size(EEG.icaact,1))); % sanity check
    clear EEG;
end

end
