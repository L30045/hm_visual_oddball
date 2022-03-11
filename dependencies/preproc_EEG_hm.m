function [EEG_prep, cali_5m] = preproc_EEG_hm(ori_EEG, varargin)
if ~isempty(varargin)
    ori_cali_EEG = varargin{1};
else
    ori_cali_EEG = [];
end
%% EEG processing
%% parameter setting
% band pass frequencies
f_bp = [1, 50];
% rereference channel
ref_Ch = {'M1','M2'};
% bad channel removal
thresFlatChannel = 5;
highPassBand = -1;
thresPoorCorrChannel = 0.7;
thresLineNoiseChannel = 4;
thresASR = 10;
thresWindow = -1;

%% preprocessing
% remove gyro channel
rmCh_labels = {'GyroX','GyroY','GyroZ'};
EEG = pop_select(ori_EEG, 'nochannel', rmCh_labels);
% centerized channels
EEG = pop_chanedit(EEG, 'eval','chans = pop_chancenter( chans, [],[]);');
% bandpass [1 50]
EEG = pop_eegfiltnew(EEG,f_bp(1),f_bp(2));
% reref
if isempty(ref_Ch)
    % remove mastoid channels if any
    EEG = pop_select(EEG, 'nochannel', {'M1','M2','A1','A2'});
    if ~isempty(ori_cali_EEG)
        cali_EEG = pop_chanedit(ori_cali_EEG, 'eval','chans = pop_chancenter( chans, [],[]);');
        cali_EEG = pop_eegfiltnew(cali_EEG,f_bp(1),f_bp(2));
        % remove bad channels and ASR
        EEG = clean_rawdata(EEG, thresFlatChannel, highPassBand, thresPoorCorrChannel, thresLineNoiseChannel, -1, thresWindow);
        cali_5m = pop_select(cali_EEG,'channel',{EEG.chanlocs.labels});
        EEG = clean_asr(EEG,thresASR,[],[],[],cali_5m);       
    else
        % remove bad channels and ASR
        EEG = clean_rawdata(EEG, thresFlatChannel, highPassBand, thresPoorCorrChannel, thresLineNoiseChannel, thresASR, thresWindow);
    end
    EEG  = pop_reref(EEG,ref_Ch);
else
    EEG = pop_reref(EEG,find(ismember({EEG.chanlocs.labels},ref_Ch)));
    if ~isempty(ori_cali_EEG)
        cali_EEG = pop_chanedit(ori_cali_EEG, 'eval','chans = pop_chancenter( chans, [],[]);');
        cali_EEG = pop_eegfiltnew(cali_EEG,f_bp(1),f_bp(2));
        cali_EEG = pop_reref(cali_EEG,find(ismember({EEG.chanlocs.labels},ref_Ch)));
        % remove bad channels and ASR
        EEG = clean_rawdata(EEG, thresFlatChannel, highPassBand, thresPoorCorrChannel, thresLineNoiseChannel, -1, thresWindow);
        cali_5m = pop_select(cali_EEG,'channel',{EEG.chanlocs.labels});       
        EEG = clean_asr(EEG,thresASR,[],[],[],cali_5m);
    else
        % remove bad channels and ASR
        EEG = clean_rawdata(EEG, thresFlatChannel, highPassBand, thresPoorCorrChannel, thresLineNoiseChannel, thresASR, thresWindow);
    end
end

EEG_prep = EEG;
disp('Done')

end