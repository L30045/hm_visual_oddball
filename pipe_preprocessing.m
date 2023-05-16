%% preprocessing pipeline
% remove non-EEG channels
% bandpass 1-50 Hz
% reref to M1, M2
% bad channel removal
% ASR with cutoff parameter k = 10

savepath = '/data/projects/yuan/2021 HM_visual_oddball/dataset/preproc_data_ASR5/';
filepath = '/data/projects/yuan/2021 HM_visual_oddball/dataset/oddball/';
filename_list = {dir([loadpath,'*Oddball*.xdf']).name};

%% Without head movement
subj_list = [1,4:6,8:10];
parfor j = 1:length(subj_list)
    subj_i = subj_list(j);
    for i = 1:2
%         if subj_i == 8 && i==1
%             break
%         end
        cond_i = i;
        switch i
            case 1
                condname = 'Inner';
            case 2
                condname = 'Outer';
        end
        filename = sprintf('hm_visual_oddball_s%02d_cond%d.xdf',subj_i,cond_i);
        icaname = sprintf('s%02d_%s_resampled_250.set',subj_i, condname);

        %% load data
        EEG = pop_loadxdf([filepath,filename]);

        %% preprocessing
        % remove gyro channel
        rmCh_labels = {'GyroX','GyroY','GyroZ'};
        EEG = pop_select(EEG, 'nochannel', rmCh_labels);
        % band pass
        EEG = pop_eegfiltnew(EEG,1,50);
        % resample to 250Hz
        EEG = pop_resample(EEG,250);
        % rereference
        EEG = pop_reref(EEG,{'M1','M2'});
        % Add channel location    
        EEG = pop_chanedit(EEG, 'lookup','D:\Research\eeglab\plugins/dipfit/standard_BEM/elec/standard_1005.elc'); % MNI
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
        EEG_ica = pop_runica(EEG,'icatype','runica','extend',1);
        % ICLabel and remove eye, muscle comp
        EEG_ica = pop_iclabel(EEG_ica,'default');
        EEG_ica = pop_icflag(EEG_ica, [NaN, NaN; 0.8, 1; 0.8, 1; 0.8, 1; 0.8, 1;0.8, 1;NaN, NaN;]);
        EEG = pop_subcomp(EEG_ica,find(EEG_ica.reject.gcompreject));

        % savefile
        pop_saveset(EEG, [savepath,icaname]);
        % fprintf('Completed %s\n',filename{i});
    end
end
