%% Path setting
eegpath = which('eeglab.m');
eegpath = eegpath(1:end-8);

%% load data
filepath = 'D:\ssvep dataset new\';
filename = 's05_hm_visual_oddball_2m.xdf';
EEG = pop_loadxdf([filepath,filename]);


%% preprocessing
% band pass
EEG = pop_eegfiltnew(EEG,1,50);
% look up channel location
% Rename channels label
chLabel_new = {'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3',...
               'FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1',...
               'P3','P5','P7','P9','PO7','PO3','O1','Iz','Oz','POz','Pz','CPz',...
               'Fpz','Fp2','AF8','AF4','AFz','Fz','F2','F4','F6','F8',...
               'FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8',...
               'TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2'};
EEG = pop_select(EEG,'nochannel',{'Trig1','EX3','EX4','EX5','EX6','EX7','EX8'});
EEG = pop_reref(EEG,{'EX1','EX2'});
[EEG.chanlocs.labels] = deal(chLabel_new{:});
% Add channel location    
EEG = pop_chanedit(EEG, 'lookup',[eegpath,'plugins/dipfit/standard_BEM/elec/standard_1005.elc']); % MNI
% re-center channel location
EEG = pop_chanedit(EEG, 'eval','chans = pop_chancenter( chans, [],[]);');
% ASR
thresFlatChannel = 5;
highPassBand = -1;
thresPoorCorrChannel = 0.7;
thresLineNoiseChannel = 4;
thresASR = 10;
thresWindow = -1;
EEG = clean_rawdata(EEG, thresFlatChannel, highPassBand, thresPoorCorrChannel, thresLineNoiseChannel, thresASR, thresWindow);
% ICA
EEG_ica = pop_runica(EEG_prep,'icatype','runica','extend',1);
% ICLabel and remove eye, muscle comp
EEG_ica = pop_iclabel(EEG_ica,'default');
EEG_ica = pop_icflag(EEG_ica, [NaN, NaN; 0.8, 1; 0.8, 1; 0.8, 1; 0.8, 1;0.8, 1;NaN, NaN;]);
EEG_ica = pop_subcomp(EEG_ica,find(EEG_ica.reject.gcompreject));
