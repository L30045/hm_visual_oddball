function processed_epoch = epoch_ssvep_internal(EEG, nbEEGchan, ev_name, len_epoch, rmbase_flag, rmtrial_flag)
%% epoching for ssvep data
% 1. epoch data
% 2. remove baseline for EEG data
% 3. remove bad trials
processed_epoch = pop_epoch(EEG,ev_name,len_epoch/1000, 'epochinfo', 'yes');

%% remove baseline
if rmbase_flag
    % save behavioral data
    behavi_data = processed_epoch.data(nbEEGchan+1:end,:,:);
    % remove baseline
    processed_epoch = pop_rmbase(processed_epoch,[max(len_epoch(1),processed_epoch.times(1)) 0],[]);
    % restore behavioral data
    processed_epoch.data(nbEEGchan+1:end,:,:) = behavi_data;
end

%% epoch auto rejection (UNDER CONSTRUCTION)
if rmtrial_flag
% std_epoch = pop_autorej(std_epoch,'electrodes',1:EEG_ica.nbchan-2,'nogui','on');
% dev_epoch = pop_autorej(dev_epoch,'electrodes',1:plt_EEG.nbchan-2,'nogui','on');
% grab_epoch = pop_autorej(grab_epoch,'electrodes',1:plt_EEG.nbchan-2,'nogui','on');
% fix_epoch = pop_autorej(fix_epoch,'electrodes',1:plt_EEG.nbchan-2,'nogui','on');
% fix_std = pop_autorej(fix_std,'electrodes',1:plt_EEG.nbchan-2,'nogui','on');
% fix_dev = pop_autorej(fix_dev,'electrodes',1:plt_EEG.nbchan-2,'nogui','on');
% ud_epoch = pop_autorej(ud_epoch,'electrodes',1:plt_EEG.nbchan-2,'nogui','on');
% lr_epoch = pop_autorej(lr_epoch,'electrodes',1:plt_EEG.nbchan-2,'nogui','on');
% gip_std = pop_autorej(gip_std,'electrodes',1:plt_EEG.nbchan-2,'nogui','on');
% gip_dev = pop_autorej(gip_dev,'electrodes',1:plt_EEG.nbchan-2,'nogui','on');
end

end