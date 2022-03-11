%% preprocessing pipeline
% remove non-EEG channels
% bandpass 1-50 Hz
% reref to M1, M2
% bad channel removal
% ASR with cutoff parameter k = 10

addpath('../dataset/')
addpath(genpath('dependencies/'))
%% Without head movement
parfor subj_i = 8:10
    for i = 1:2
        if subj_i == 8 && i==1
            break
        end
        cond_i = i;
        filename = sprintf('hm_visual_oddball_s%02d_cond%d.xdf',subj_i,cond_i);
        icaname = sprintf('s%02d_cond%d_ica_k10.set',subj_i, cond_i);

        %% load data
        [~, EEG, ~, ~, ~, ~] = load_eyetracking_hm(filename);
        % re-center channel location
        EEG = pop_chanedit(EEG, 'eval','chans = pop_chancenter( chans, [],[]);');

        %% preprocessing
        EEG_prep = preproc_EEG_hm(EEG);
        rmCh = setdiff({EEG.chanlocs.labels},{EEG_prep.chanlocs.labels});
        EEG_ica = pop_runica(EEG_prep,'icatype','runica','extend',1);
        % ICLabel and remove eye, muscle comp
        EEG_ica = pop_iclabel(EEG_ica,'default');
        EEG_ica = pop_icflag(EEG_ica, [NaN, NaN; 0.8, 1; 0.8, 1; 0.8, 1; 0.8, 1;0.8, 1;NaN, NaN;]);
        EEG_ica = pop_subcomp(EEG_ica,find(EEG_ica.reject.gcompreject));

        pop_saveset(EEG_ica, ['../dataset/', icaname]);
    end
end
