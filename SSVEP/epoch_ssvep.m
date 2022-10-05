function epoch_struct = epoch_ssvep(filename,rmbase_flag,rmtrial_flag)
%% This function adds event markers to the EEG structure and interprets behavioral data as pseudo-EEG signals.
% Input:
%       filename:
%           file to be processed (.xdf)
%       EEG:
%           EEG structure to be added markers and epoched. If EEG is not
%           given, the non-preprocessed EEG in filename will be used.
%           (NOTE: I run the preprocessing for EEG before epoching.)

% Output:
%       data_struct:
%           contain the raw behavioral data from eye tracker. See details
%           in >> help load_eyetracking_hm
%       fix_struct:
%           contain the results of fixation calculation. See details in >>
%           help cal_fix_pupil
%       t_c:
%           Time stamps of the eye gaze stream
%       t_std:
%           Time stamps of the standard trials after sychronizing
%           with t_c
%       t_dev:
%           Time stamps of the deviant trials after sychronizing
%           with t_c
%       epoch_struct:
%           EEG epochs time-lock to different events.
%           (1) std_epoch: standard trials time-lock to stimulus onset
%           (2) dev_epoch: deviant trials time-lock to stimulus onset
%           (3) grab_epoch: standard trials time-lock to response
%           (4) fix_epoch: fixation trials time-lock to fixation onset
%           (5) fix_std: standard trials time-lock to fixation onset
%           (6) fix_dev: deviant trials time-lock to fixation onset
%           (7) ud_epoch: up and down trials time-lock to stimulus onset
%           (8) lr_epoch: left and right trials time-lock to stimulus onset
%           (9) gip_std: standard trials time-lock to GIP onset
%           (10) gip_dev: deviant trials time-lock to GIP onset
%           (11) event_time: Event onset time and epoch length. The order
%                are the same as epoch_struct.


%% parameter setting
% epoch length
len_epoch = [-200 800]; % epoch length for stimulus-lock epoch
len_epoch_grab = [-500 1000]; % epoch length for GIP/Response/Fixation-lock epoch

% calculate fixation
len_blink = 100;
%   [thres_ang]: angular higher than this threshold will be marked as
%   saccade) (Default: 0.5 deg/s ref. Nystrom 2010.) (deg)
%   [thres_ang_v]: angular velocity higher than this threshold will be marked
%   as saccade) (Default: 30 deg/s ref. tobii, 130 deg/s ref. Eye tracking
%   2017.) (deg/s)
thres_ang = 1; 
thres_ang_v = 30;
% fix_selection = 'dispersion';
fix_selection = 'velocity';
max_fix_interval = 75;
max_fix_ang = 1;
min_fix_len = 200;

%% load data
[data_struct,~, ~, ~, ~, ~] = load_eyetracking_hm(filename);
EEG = pop_loadset([filename(1:end-4),'.set']);
nbchan = EEG.nbchan;

%% find out stimulus onset time
idx_ring1 = cellfun(@(x) ~isempty(regexp(x, 'Ring\ 1', 'ONCE')), {EEG.event.type});
idx_ring2 = cellfun(@(x) ~isempty(regexp(x, 'Ring\ 2', 'ONCE')), {EEG.event.type});
idx_blue = cellfun(@(x) ~isempty(regexp(x, 'blue_cube', 'ONCE')), {EEG.event.type});

%% calculate fixation
% fix_struct = cal_fix_pupil(test_data,data_struct.srate);
test_data = [data_struct.ori_eye_3D_pos;data_struct.eye_openess_idx];
fix_struct = cal_fix_pupil(test_data,data_struct.srate,'blink_length',len_blink,...
                           'thres_ang',thres_ang,'thres_ang_v',thres_ang_v,'fix_selection',fix_selection,...
                           'max_fix_interval',max_fix_interval,'max_fix_ang',max_fix_ang,'min_fix_len',min_fix_len);
eye_fix_idx = fix_struct.eye_fixation.eye_fix_idx;

%% get time stramps from eye tracker
t_c = fix_struct.time_stamps*1000; % change unit to ms
t_ring1 = ([EEG.event(idx_ring1).latency]-1)/EEG.srate*1000; % convert from sample point to time (ms) 
t_ring2 = ([EEG.event(idx_ring2).latency]-1)/EEG.srate*1000; % convert from sample point to time (ms)
t_blue_ori = ([EEG.event(idx_blue).latency]-1)/EEG.srate*1000; % convert from sample point to time (ms)

%% find up/down left/right event
up_idx = cellfun(@(x) ~isempty(regexp(x,'index\(.\)\:\ 1','ONCE')),{EEG.event.type});
down_idx = cellfun(@(x) ~isempty(regexp(x,'index\(.\)\:\ 3','ONCE')),{EEG.event.type});
left_idx = cellfun(@(x) ~isempty(regexp(x,'index\(.\)\:\ 2','ONCE')),{EEG.event.type});
right_idx = cellfun(@(x) ~isempty(regexp(x,'index\(.\)\:\ 0','ONCE')),{EEG.event.type});

ring1_up = up_idx & idx_ring1;
ring1_down = down_idx & idx_ring1;
ring1_left = left_idx & idx_ring1;
ring1_right = right_idx & idx_ring1;
ring2_up = up_idx & idx_ring2;
ring2_down = down_idx & idx_ring2;
ring2_left = left_idx & idx_ring2;
ring2_right = right_idx & idx_ring2;

ring1_ev_up = {EEG.event(ring1_up).type};
ring1_ev_down = {EEG.event(ring1_down).type};
ring1_ev_left = {EEG.event(ring1_left).type};
ring1_ev_right = {EEG.event(ring1_right).type};
ring2_ev_up = {EEG.event(ring2_up).type};
ring2_ev_down = {EEG.event(ring2_down).type};
ring2_ev_left = {EEG.event(ring2_left).type};
ring2_ev_right = {EEG.event(ring2_right).type};

%% find GIP onset time
gip_up_idx = cellfun(@(x) ~isempty(regexp(x,'Up','match','ONCE')),{EEG.event.type});
gip_down_idx = cellfun(@(x) ~isempty(regexp(x,'Bottom','match','ONCE')),{EEG.event.type});
gip_left_idx = cellfun(@(x) ~isempty(regexp(x,'Left','match','ONCE')),{EEG.event.type});
gip_right_idx = cellfun(@(x) ~isempty(regexp(x,'Right','match','ONCE')),{EEG.event.type});
ring1_stim_idx = [ring1_up;ring1_down;ring1_left;ring1_right]';
ring2_stim_idx = [ring2_up;ring2_down;ring2_left;ring2_right]';
gip_idx = [gip_up_idx ; gip_down_idx ; gip_left_idx ; gip_right_idx]';
gip_name = {'up','down','left','right'};

% find the first fixation after stimulus
% find out stimulus onset time
gipRing1_time = zeros(sum(ring1_up),4);
gipRing2_time = gipRing1_time;
stimRing1_time = gipRing1_time;
stimRing2_time = gipRing1_time;
for i = 1:4
    [EEG, stimTime, gipTime] = add_gip_evmarker(EEG, ring1_stim_idx(:,i), gip_idx(:,i), ['ring1_gip_',gip_name{i}]);
    stimRing1_time(:,i) = stimTime;
    gipRing1_time(:,i) = gipTime;
    [EEG, stimTime, gipTime] = add_gip_evmarker(EEG, ring2_stim_idx(:,i), gip_idx(:,i), ['ring2_gip_',gip_name{i}]);
    stimRing2_time(:,i) = stimTime;
    gipRing2_time(:,i) = gipTime;
end

%% add fixation event into EEG structure
% find fixation start and end
% fix_idx = [0,eye_fix_idx(2:end)-eye_fix_idx(1:end-1)];
% fprintf('Number of fixation: %d\n', sum(abs(fix_idx))/2);
% fix_duration = diff(fix_struct.time_stamps([find(fix_idx==1)' find(fix_idx==-1)']),[],2);
% fprintf('Median of fixation duration: %f (ms)\n', median(fix_duration*1000));
% fix_start = t_c(fix_idx==1); % ms
% fix_end = t_c(fix_idx==-1); % ms
% % adding fixation event
% for i = 1:length(fix_start)
%     len_event = length(EEG.event);
%     EEG.event(len_event+1).type = 'fix_start';
%     EEG.event(len_event+1).latency = fix_start(i)/1000*EEG.srate; % change to sample point
%     EEG.event(len_event+1).urevent = len_event+1;
% end
% for i = 1:length(fix_end)
%     len_event = length(EEG.event);
%     EEG.event(len_event+1).type = 'fix_end';
%     EEG.event(len_event+1).latency = fix_end(i)/1000*EEG.srate; % change to sample point
%     EEG.event(len_event+1).urevent = len_event+1;
% end
% 
% % find the first fixation after stimulus
% % find out stimulus onset time
% t_f_ring1 = zeros(1,length(t_ring1));
% for i = 1:length(t_ring1)
%     t_f = fix_start(find(fix_start>t_ring1(i),1));
%     if ~isempty(t_f)
%         t_f_ring1(i) = t_f;
%         len_event = length(EEG.event);
%         EEG.event(len_event+1).type = 'ring1_fix_start';
%         EEG.event(len_event+1).latency = t_f/1000*EEG.srate; % change to sample point
%         EEG.event(len_event+1).urevent = len_event+1;
%     else
%         t_f_ring1(i) = NaN;
%         len_event = length(EEG.event);
%         EEG.event(len_event+1).type = 'ring1_fix_missing';
%         EEG.event(len_event+1).latency = t_f/1000*EEG.srate; % change to sample point
%         EEG.event(len_event+1).urevent = len_event+1;
%     end
% end
% t_f_ring2 = zeros(1,length(t_ring2));
% for i = 1:length(t_ring2)
%     t_f = fix_start(find(fix_start>t_ring2(i),1));
%     if ~isempty(t_f)
%         t_f_ring2(i) = t_f;
%         len_event = length(EEG.event);
%         EEG.event(len_event+1).type = 'ring2_fix_start';
%         EEG.event(len_event+1).latency = t_f/1000*EEG.srate; % change to sample point
%         EEG.event(len_event+1).urevent = len_event+1;
%     else
%         t_f_ring2(i) = NaN;
%         len_event = length(EEG.event);
%         EEG.event(len_event+1).type = 'ring2_fix_missing';
%         EEG.event(len_event+1).latency = t_f/1000*EEG.srate; % change to sample point
%         EEG.event(len_event+1).urevent = len_event+1;
%     end
% end
% t_f_blue = zeros(1,length(t_blue_ori));
% for i = 1:length(t_blue_ori)
%     t_f = fix_start(find(fix_start>t_blue_ori(i),1));
%     if ~isempty(t_f)
%         t_f_blue(i) = t_f;
%         len_event = length(EEG.event);
%         EEG.event(len_event+1).type = 'blue_fix_start';
%         EEG.event(len_event+1).latency = t_f/1000*EEG.srate; % change to sample point
%         EEG.event(len_event+1).urevent = len_event+1;
%     else
%         t_f_blue(i) = NaN;
%         len_event = length(EEG.event);
%         EEG.event(len_event+1).type = 'blue_fix_missing';
%         EEG.event(len_event+1).latency = t_f/1000*EEG.srate; % change to sample point
%         EEG.event(len_event+1).urevent = len_event+1;
%     end
% end

%% Linear interpret eye gaze stream into EEG stream
% smoothing using moving average
mv_hd = smoothing_mv_avg(data_struct.ori_head_direct, true(1,size(data_struct.ori_head_direct,2)), 5);
mv_gip = smoothing_mv_avg(data_struct.ori_gip_3D_pos, true(1,size(data_struct.ori_gip_3D_pos,2)), 5);
mv_hloc = smoothing_mv_avg(data_struct.ori_head_3D_loc, true(1,size(data_struct.ori_head_3D_loc,2)), 5);
% get eye blink index
blk_idx = false(1,length(fix_struct.time_stamps));
dLose_idx = false(size(blk_idx));
if ~isempty(fix_struct.gap_detection.blink_idx)
    blk_idx([fix_struct.gap_detection.blink_idx{:}]) = true;
end
if ~isempty(fix_struct.gap_detection.dataLose_idx) 
    dLose_idx([fix_struct.gap_detection.dataLose_idx{:}]) = true;
end

%% Interpretation
% find out time point to interpt
[interp_p, rm_idx] = linear_interp_idx(t_c, EEG.times);
data_bus = [mv_hloc; mv_hd; mv_gip; blk_idx; dLose_idx;...
            fix_struct.reconstruct_eye_pos.test_data(1:3,:);...
            data_struct.ori_pupil_diameter];
name_list = {'HeadLoc_x','HeadLoc_y','HeadLoc_z',...
    'HeadDirect_x','HeadDirect_y','HeadDirect_z',...
    'GIP_x','GIP_y','GIP_z',...
    'blink_idx','dataLose_idx',...
    'EyeDirect_x','EyeDirect_y','EyeDirect_z'...
    'Pupil left', 'Pupil right'};
interp_data = const_interp_data(data_bus, interp_p, rm_idx);

%% add this information into EEG structure and create virtual channels
for n_i = 1:length(name_list)
    EEG = insert_data(EEG, interp_data(n_i,:), name_list{n_i});
end

%% epoching
ring1_epoch_up = epoch_ssvep_internal(EEG,nbchan,ring1_ev_up,len_epoch,rmbase_flag,rmtrial_flag);
ring1_epoch_down = epoch_ssvep_internal(EEG,nbchan,ring1_ev_down,len_epoch,rmbase_flag,rmtrial_flag);
ring1_epoch_left = epoch_ssvep_internal(EEG,nbchan,ring1_ev_left,len_epoch,rmbase_flag,rmtrial_flag);
ring1_epoch_right = epoch_ssvep_internal(EEG,nbchan,ring1_ev_right,len_epoch,rmbase_flag,rmtrial_flag);
ring1_gip_up = epoch_ssvep_internal(EEG,nbchan,{'ring1_gip_up_start'},len_epoch_grab,rmbase_flag,rmtrial_flag);
ring1_gip_down = epoch_ssvep_internal(EEG,nbchan,{'ring1_gip_down_start'},len_epoch_grab,rmbase_flag,rmtrial_flag);
ring1_gip_left = epoch_ssvep_internal(EEG,nbchan,{'ring1_gip_left_start'},len_epoch_grab,rmbase_flag,rmtrial_flag);
ring1_gip_right = epoch_ssvep_internal(EEG,nbchan,{'ring1_gip_right_start'},len_epoch_grab,rmbase_flag,rmtrial_flag);

ring2_epoch_up = epoch_ssvep_internal(EEG,nbchan,ring2_ev_up,len_epoch,rmbase_flag,rmtrial_flag);
ring2_epoch_down = epoch_ssvep_internal(EEG,nbchan,ring2_ev_down,len_epoch,rmbase_flag,rmtrial_flag);
ring2_epoch_left = epoch_ssvep_internal(EEG,nbchan,ring2_ev_left,len_epoch,rmbase_flag,rmtrial_flag);
ring2_epoch_right = epoch_ssvep_internal(EEG,nbchan,ring2_ev_right,len_epoch,rmbase_flag,rmtrial_flag);
ring2_gip_up = epoch_ssvep_internal(EEG,nbchan,{'ring2_gip_up_start'},len_epoch_grab,rmbase_flag,rmtrial_flag);
ring2_gip_down = epoch_ssvep_internal(EEG,nbchan,{'ring2_gip_down_start'},len_epoch_grab,rmbase_flag,rmtrial_flag);
ring2_gip_left = epoch_ssvep_internal(EEG,nbchan,{'ring2_gip_left_start'},len_epoch_grab,rmbase_flag,rmtrial_flag);
ring2_gip_right = epoch_ssvep_internal(EEG,nbchan,{'ring2_gip_right_start'},len_epoch_grab,rmbase_flag,rmtrial_flag);

ring1_stim = struct('up', ring1_epoch_up, 'down', ring1_epoch_down, 'left', ring1_epoch_left, 'right', ring1_epoch_right);
ring2_stim = struct('up', ring2_epoch_up, 'down', ring2_epoch_down, 'left', ring2_epoch_left, 'right', ring2_epoch_right);
ring1_gip = struct('up', ring1_gip_up, 'down', ring1_gip_down, 'left', ring1_gip_left, 'right', ring1_gip_right);
ring2_gip = struct('up', ring2_gip_up, 'down', ring2_gip_down, 'left', ring2_gip_left, 'right', ring2_gip_right);
ring1_epoch = struct('stim',ring1_stim,'gip',ring1_gip);
ring2_epoch = struct('stim',ring2_stim,'gip',ring2_gip);

%% fix related epoch
if false
    fix_epoch = pop_epoch(EEG,{'fix_start'},len_epoch_grab/1000, 'epochinfo', 'yes');
    fix_std = pop_epoch(EEG,{'circle_fix_start'},len_epoch_grab/1000, 'epochinfo', 'yes');
    fix_dev = pop_epoch(EEG,{'triangle_fix_start'},len_epoch_grab/1000, 'epochinfo', 'yes');
    fix_blue = pop_epoch(EEG,{'blue_fix_start'},len_epoch_grab/1000, 'epochinfo', 'yes');
    behavi_fixAll = fix_epoch.data(nbchan+1:end,:,:);
    behavi_fstd = fix_std.data(nbchan+1:end,:,:);
    behavi_fdev = fix_dev.data(nbchan+1:end,:,:);
    behavi_fblue = fix_blue.data(nbchan+1:end,:,:);
    fix_epoch = pop_rmbase(fix_epoch,[max(len_epoch_grab(1),fix_epoch.times(1)) 0],[]);
    fix_std = pop_rmbase(fix_std,[max(len_epoch_grab(1),fix_std.times(1)) 0],[]);
    fix_dev = pop_rmbase(fix_dev,[max(len_epoch_grab(1),fix_dev.times(1)) 0],[]);
    fix_blue = pop_rmbase(fix_blue,[max(len_epoch_grab(1),fix_blue.times(1)) 0],[]);
    fix_epoch.data(nbchan+1:end,:,:) = behavi_fixAll;
    fix_std.data(nbchan+1:end,:,:) = behavi_fstd;
    fix_dev.data(nbchan+1:end,:,:) = behavi_fdev;
    fix_blue.data(nbchan+1:end,:,:) = behavi_fblue;
end

%% record event time
diff_gip_stim_ring1 = gipRing1_time - stimRing1_time;
diff_gip_stim_ring2 = gipRing2_time - stimRing2_time;
% fixAll_time = fix_start; % msec
% fixRing1_time = t_f_ring1; % msec
% fixRing2_time = t_f_ring2; % msec
% fixBlue_time = t_f_blue; % msec
% event_time = struct('diff_gip_stim_ring1',diff_gip_stim_ring1,'diff_gip_stim_ring2',diff_gip_stim_ring2,...
%                     'fixAll_time',fixAll_time,'fixRing1_time',fixRing1_time,'fixRing2_time',fixRing2_time,'fixBlue_time',fixBlue_time,...
%                     'gipRing1_time',gipRing1_time,'gipRing2_time',gipRing2_time,...
%                     'ring1_up',ring1_up,'ring1_down',ring1_down,'ring1_left',ring1_left,'ring1_right',ring1_right,...
%                     'ring2_up',ring2_up,'ring2_down',ring2_down,'ring2_left',ring2_left,'ring2_right',ring2_right);
event_time = struct('diff_gip_stim_ring1',diff_gip_stim_ring1,'diff_gip_stim_ring2',diff_gip_stim_ring2,...
                    'ring1_up',ring1_up,'ring1_down',ring1_down,'ring1_left',ring1_left,'ring1_right',ring1_right,...
                    'ring2_up',ring2_up,'ring2_down',ring2_down,'ring2_left',ring2_left,'ring2_right',ring2_right);

                
%% output
epoch_struct = struct('ring1',ring1_epoch, 'ring2',ring2_epoch,'event_time',event_time);
                  
end