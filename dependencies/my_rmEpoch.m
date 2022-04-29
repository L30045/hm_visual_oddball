function [output,rm_idx] = my_rmEpoch(EEG)
EEG = pop_jointprob(EEG,1,1,5,5,0,0,0,[],0);
EEG = pop_eegthresh(EEG,1,1,-100,100,-0.5,0.998,2,0);
EEG = pop_rejtrend(EEG,1,1,750,50,0.3,2,0);
EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
%% reject data by mean amplitude
data = EEG.data;
% norm_data = max(abs(squeeze(mean(data,2))));
% norm_data = (norm_data-mean(norm_data))./std(norm_data);
% rm_epoch_idx = abs(norm_data)>=3;

% reject data by max amplitude
norm_data = max(squeeze(max(abs(data),[],2)));
norm_data = (norm_data-mean(norm_data))./std(norm_data);
rm_epoch_idx = abs(norm_data)>=2;
EEG.reject.rejglobal = EEG.reject.rejglobal|rm_epoch_idx;

%%
output = pop_rejepoch(EEG,EEG.reject.rejglobal,0);
rm_idx = EEG.reject.rejglobal;

end