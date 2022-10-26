function output = my_rmEpoch(epoch_struct, thres_time, varargin)
%% 
if ~isempty(varargin)
    tarCh = varargin{1};
    rm_thres = varargin{2};
    [rm_idx_stim, rm_idx_gip] = my_rmbase(stim_epoch, gip_epoch, event_time, tarCh, rm_thres);
end

%% remove trial with too short/long GIP time
output = epoch_struct;
diff_gip_std = epoch_struct.event_time.diff_gip_std;
diff_gip_dev = epoch_struct.event_time.diff_gip_dev;
rm_idx_std = diff_gip_std<thres_time(1) |diff_gip_std>thres_time(2);
rm_idx_dev = diff_gip_dev<thres_time(1) |diff_gip_dev>thres_time(2);
output.event_time.std_time = epoch_struct.event_time.std_time(~rm_idx_std);
output.event_time.dev_time = epoch_struct.event_time.dev_time(~rm_idx_dev);
output.event_time.diff_gip_std = diff_gip_std(~rm_idx_std);
output.event_time.diff_gip_dev = diff_gip_dev(~rm_idx_dev);
output.event_time.gipStd_time = epoch_struct.event_time.gipStd_time(~rm_idx_std);
output.event_time.gipDev_time = epoch_struct.event_time.gipDev_time(~rm_idx_dev);

% grab
% if ~isempty(output.event_time.grab_time)
%     output.event_time.grab_time = epoch_struct.event_time.grab_time(~rm_idx_std);
%     output.event_time.diff_stim_grab = epoch_struct.event_time.diff_stim_grab(~rm_idx_std);
%     output.grab_epoch = pop_rejepoch(epoch_struct.grab_epoch,rm_idx_std,0);
% end
if ~isempty(output.event_time.fixStd_time)
    output.event_time.fixStd_time = epoch_struct.event_time.fixStd_time(~rm_idx_std);
    output.event_time.fixDev_time = epoch_struct.event_time.fixDev_time(~rm_idx_dev);
    if ~isempty(output.fix_epoch)
        output.fix_std = pop_rejepoch(epoch_struct.fix_std,rm_idx_std,0);
        output.fix_dev = pop_rejepoch(epoch_struct.fix_dev,rm_idx_dev,0);
    end
end

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