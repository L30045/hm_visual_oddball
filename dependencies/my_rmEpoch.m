function output = my_rmEpoch(epoch_struct)
%% remove trial with too short/long GIP time
output = epoch_struct;
diff_gip_std = epoch_struct.event_time.diff_gip_std;
diff_gip_dev = epoch_struct.event_time.diff_gip_dev;
rm_idx_std = diff_gip_std<100 |diff_gip_std>1000;
rm_idx_dev = diff_gip_dev<100 |diff_gip_dev>1000;
output.event_time.diff_gip_std = diff_gip_std(~rm_idx_std);
output.event_time.diff_gip_dev = diff_gip_dev(~rm_idx_dev);

%% remove EEG trial
output.std_epoch = pop_rejepoch(epoch_struct.std_epoch,rm_idx_std,0);
output.dev_epoch = pop_rejepoch(epoch_struct.dev_epoch,rm_idx_dev,0);
output.gip_std = pop_rejepoch(epoch_struct.gip_std,rm_idx_std,0);
output.gip_dev = pop_rejepoch(epoch_struct.gip_dev,rm_idx_dev,0);

%%
% EEG = pop_jointprob(EEG,1,1,5,5,0,0,0,[],0);
% EEG = pop_eegthresh(EEG,1,1,-100,100,-0.5,0.998,2,0);
% EEG = pop_rejtrend(EEG,1,1,750,50,0.3,2,0);
% EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
% EEG.reject.rejglobal = EEG.reject.rejglobal|rm_epoch_idx;

%%
% output = pop_rejepoch(EEG,EEG.reject.rejglobal,0);
% rm_idx = EEG.reject.rejglobal;

end