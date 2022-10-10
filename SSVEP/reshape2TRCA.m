function output = reshape2TRCA(epoch, tarCh)
%% Reshape data for TRCA
% Input: structure for epoched eeg 
% Output: TRCA format data
% [# of targets, # of channels, # of sampling points, # of blocks] = size(output);


% Parameter settings
% epoch time
time = epoch.up.times;
% sampling rate
fs = epoch.up.srate;
% target channels index
idx_ch = ismember({epoch.up.chanlocs.labels},tarCh);
% sample point for time equals to 0
len_zero_smpl = find(time==0,1); % sample
% SSVEP offset
len_delay_s = 0.13; % second
% SSVEP duration 
len_gaze_s = 0.5; % second
% Data length [samples]
len_gaze_smpl = round(len_gaze_s*fs);           
% Visual latency [samples]
len_delay_smpl = round(len_delay_s*fs);    
% segment data
segment_data = [len_delay_smpl+1:len_delay_smpl+len_gaze_smpl]+len_zero_smpl;

% reshape format ([right, up, left, down])
output = zeros(4,length(tarCh),length(segment_data),20);
output(1,:,:,:) = epoch.right.data(idx_ch,segment_data,:);
output(2,:,:,:) = epoch.up.data(idx_ch,segment_data,:);
output(3,:,:,:) = epoch.left.data(idx_ch,segment_data,:);
output(4,:,:,:) = epoch.down.data(idx_ch,segment_data,:);

end
