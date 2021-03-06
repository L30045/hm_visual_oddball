function [data_struct, fix_struct, t_c, t_std, t_dev, epoch_struct, EEG] = epoch_ez(filename, varargin)
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
%           Corrected time stamps of the eye gaze stream after sychronizing
%           with EEG stream
%       t_std:
%           Corrected time stamps of the standard trials after sychronizing
%           with t_c
%       t_dev:
%           Corrected time stamps of the deviant trials after sychronizing
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

%% parameter setting
% epoch length
len_epoch = [-200 1000]; % epoch length for stimulus-lock epoch
len_epoch_grab = [-500 1000]; % epoch length for GIP/Response/Fixation-lock epoch

% calculate fixation
len_blink = 100;
%   [thres_ang]: angular higher than this threshold will be marked as
%   saccade) (Default: 0.5 deg/s ref. Nystrom 2010.) (deg)
%   [thres_ang_v]: angular velocity higher than this threshold will be marked
%   as saccade) (Default: 30 deg/s ref. tobii, 130 deg/s ref. Eye tracking
%   2017.) (deg/s)
thres_ang = 0.5; 
thres_ang_v = 30;
fix_selection = 'dispersion';
% fix_selection = 'velocity';
max_fix_interval = 75;
max_fix_ang = 1;
min_fix_len = 150;

%% load data
[data_struct, EEG_ori, ~, s_eyeGaze, ~, s_EEG] = load_eyetracking_hm(filename);
if isempty(varargin)
    EEG = EEG_ori;
else
    EEG = varargin{1};
end

%% find out stimulus onset time
idx_dev = cellfun(@(x) ~isempty(regexp(x, 'Deviant', 'ONCE')), {EEG.event.type});
idx_std = cellfun(@(x) ~isempty(regexp(x, 'Standard', 'ONCE')), {EEG.event.type});
idx_tri = cellfun(@(x) ~isempty(regexp(x, 'Trigger', 'ONCE')), {EEG.event.type});
idx_tri(find(idx_tri,2,'last')) = false;
std_ev = {EEG.event(idx_std).type};
dev_ev = {EEG.event(idx_dev).type};
grab_ev = [EEG.event(idx_tri).latency];

%% calculate fixation
% fix_struct = cal_fix_pupil(test_data,data_struct.srate);
test_data = [data_struct.ori_eye_3D_pos;data_struct.eye_openess_idx];
fix_struct = cal_fix_pupil(test_data,data_struct.srate,'blink_length',len_blink,...
                           'thres_ang',thres_ang,'thres_ang_v',thres_ang_v,'fix_selection',fix_selection,...
                           'max_fix_interval',max_fix_interval,'max_fix_ang',max_fix_ang,'min_fix_len',min_fix_len);
eye_fix_idx = fix_struct.eye_fixation.eye_fix_idx;

%% synchronize eye gaze stream and EEG stream
t_offset = str2double(s_eyeGaze.info.first_timestamp) - str2double(s_EEG.info.first_timestamp); % sec
% correct time stamps of the eye gaze stream
t_c = (fix_struct.time_stamps - t_offset)*1000; % change unit to ms

% synchronize EEG stream to eye gaze stream
% sort event by pupil displacement/ angular velocity
% round up EEG stream time stamps with eye gaze stream time stamps
t_std_ori = ([EEG.event(idx_std).latency]-1)/EEG.srate*1000; % convert from sample point to time (ms) 
t_dev_ori = ([EEG.event(idx_dev).latency]-1)/EEG.srate*1000; % convert from sample point to time (ms)
t_std = zeros(size(t_std_ori));
t_dev = zeros(size(t_dev_ori));
for t_i = 1:length(t_std)
    t_std(t_i) = t_c(find(t_c>t_std_ori(t_i),1,'first'));
end
fprintf('RMSE of circle trials: %.f ms\n',sqrt(sum((t_std-t_std_ori).^2))/length(t_std))
fprintf('Max error: %.f ms\n',max(t_std-t_std_ori))
for t_i = 1:length(t_dev)
    t_dev(t_i) = t_c(find(t_c>t_dev_ori(t_i),1,'first'));
end
fprintf('RMSE of triangle trials: %.f ms\n',sqrt(sum((t_dev-t_dev_ori).^2))/length(t_dev))
fprintf('Max error: %.f ms\n',max(t_dev-t_dev_ori))

%% find up/down left/right event
%boolean array where trail starts and ends
bool_start = idx_dev + idx_std; % You can use the second cube as start points
bool_end = cellfun(@(x) ~isempty(regexp(x,'End Trial','match','ONCE')),{EEG.event.type});

%boolean arrays where the LEFT/RIGHT or UP/DOWN event markers occured (I assumed
%they were mutually trail exclusive i.e. both do not occur in the same trail)
condition_1 = cellfun(@(x) ~isempty(regexp(x,'(0.862, 2.101, 2.345)','match','ONCE')),{EEG.event.type}) + ...
    cellfun(@(x) ~isempty(regexp(x,'(0.862, 1.751, 2.345)','match','ONCE')),{EEG.event.type});
condition_2 = cellfun(@(x) ~isempty(regexp(x,'(0.687, 1.926, 2.345)','match','ONCE')),{EEG.event.type}) + ...
    cellfun(@(x) ~isempty(regexp(x,'(1.037, 1.926, 2.345)','match','ONCE')),{EEG.event.type});
idx_ud = find_trails(bool_start,bool_end,condition_1); %This is a boolean array denoting trails going Up or Down
idx_lr = find_trails(bool_start,bool_end,condition_2); %This is a boolean array denoting trails going Left or Right

ud_ev = [EEG.event(logical(idx_ud)).latency];
lr_ev = [EEG.event(logical(idx_lr)).latency];

%% find GIP onset time
% up_idx = cellfun(@(x) ~isempty(regexp(x,'(0.862, 2.101, 2.345)','match','ONCE')),{EEG_ica.event.type});
% bot_idx = cellfun(@(x) ~isempty(regexp(x,'(0.862, 1.751, 2.345)','match','ONCE')),{EEG_ica.event.type});
% left_idx = cellfun(@(x) ~isempty(regexp(x,'(0.687, 1.926, 2.345)','match','ONCE')),{EEG_ica.event.type});
% right_idx = cellfun(@(x) ~isempty(regexp(x,'(1.037, 1.926, 2.345)','match','ONCE')),{EEG_ica.event.type});
up_idx = cellfun(@(x) ~isempty(regexp(x,'Up','match','ONCE')),{EEG.event.type});
bot_idx = cellfun(@(x) ~isempty(regexp(x,'Bottom','match','ONCE')),{EEG.event.type});
left_idx = cellfun(@(x) ~isempty(regexp(x,'Left','match','ONCE')),{EEG.event.type});
right_idx = cellfun(@(x) ~isempty(regexp(x,'Right','match','ONCE')),{EEG.event.type});
gip_idx = find(up_idx | bot_idx | left_idx | right_idx);

% find the first fixation after stimulus
% find out stimulus onset time
f_std = find(idx_std);
f_dev = find(idx_dev);
for i = 1:length(f_std)
    t_f = gip_idx(find(gip_idx>f_std(i),1));
    if ~isempty(t_f)
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'circle_gip_start';
        EEG.event(len_event+1).latency = EEG.event(t_f).latency; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
    else
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'circle_gip_missing';
        EEG.event(len_event+1).latency = EEG.event(f_std(i)).latency; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
    end
end
for i = 1:length(f_dev)
    t_f = gip_idx(find(gip_idx>f_dev(i),1));
    if ~isempty(t_f)
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'triangle_gip_start';
        EEG.event(len_event+1).latency = EEG.event(t_f).latency; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
    else
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'triangle_gip_missing';
        EEG.event(len_event+1).latency = EEG.event(f_dev(i)).latency; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
    end
end

%% adding event marker to EEG structure
% adding grab marker
for i = 1:length(grab_ev)
    len_event = length(EEG.event);
    EEG.event(len_event+1).type = 'grab';
    EEG.event(len_event+1).latency = grab_ev(i); 
    EEG.event(len_event+1).urevent = len_event+1;
end
% adding up/down marker
for i = 1:length(ud_ev)
    len_event = length(EEG.event);
    EEG.event(len_event+1).type = 'upDown';
    EEG.event(len_event+1).latency = ud_ev(i); 
    EEG.event(len_event+1).urevent = len_event+1;
end
% adding left/right marker
for i = 1:length(lr_ev)
    len_event = length(EEG.event);
    EEG.event(len_event+1).type = 'leftRight';
    EEG.event(len_event+1).latency = lr_ev(i); 
    EEG.event(len_event+1).urevent = len_event+1;
end


%% add fixation event into EEG structure
% find fixation start and end
fix_idx = [0,eye_fix_idx(2:end)-eye_fix_idx(1:end-1)];
fprintf('Number of fixation: %d\n', sum(abs(fix_idx))/2);
fix_duration = diff(fix_struct.time_stamps([find(fix_idx==1)' find(fix_idx==-1)']),[],2);
% figure; histogram(fix_duration*1000,'binwidt',50,'normalization','probability')
fprintf('Median of fixation duration: %f (ms)\n', median(fix_duration*1000));
fix_start = t_c(fix_idx==1); % ms
fix_end = t_c(fix_idx==-1); % ms
% adding fixation event
for i = 1:length(fix_start)
    len_event = length(EEG.event);
    EEG.event(len_event+1).type = 'fix_start';
    EEG.event(len_event+1).latency = fix_start(i)/1000*EEG.srate; % change to sample point
    EEG.event(len_event+1).urevent = len_event+1;
end
for i = 1:length(fix_end)
    len_event = length(EEG.event);
    EEG.event(len_event+1).type = 'fix_end';
    EEG.event(len_event+1).latency = fix_end(i)/1000*EEG.srate; % change to sample point
    EEG.event(len_event+1).urevent = len_event+1;
end

% find the first fixation after stimulus
% find out stimulus onset time
for i = 1:length(t_std_ori)
    t_f = fix_start(find(fix_start>t_std_ori(i),1));
    len_event = length(EEG.event);
    EEG.event(len_event+1).type = 'circle_fix_start';
    EEG.event(len_event+1).latency = t_f/1000*EEG.srate; % change to sample point
    EEG.event(len_event+1).urevent = len_event+1;
end
for i = 1:length(t_dev_ori)
    t_f = fix_start(find(fix_start>t_dev_ori(i),1));
    len_event = length(EEG.event);
    EEG.event(len_event+1).type = 'triangle_fix_start';
    EEG.event(len_event+1).latency = t_f/1000*EEG.srate; % change to sample point
    EEG.event(len_event+1).urevent = len_event+1;
end

%% Linear interpret eye gaze stream into EEG stream
%% head rotation
velocity_smooth_win_len = 40;
srate = round(fix_struct.srate);
% smoothing using moving average
mv_hd = smoothing_mv_avg(data_struct.ori_head_direct, true(1,size(data_struct.ori_head_direct,2)), 5);
% calculate angular velocity
ang = nan(1,size(mv_hd,2));
v_ang = nan(1,size(mv_hd,2));
vel_win_len = round(0.001*velocity_smooth_win_len*srate/2);
for v_i = vel_win_len+1:size(mv_hd,2)-vel_win_len
    if ~isnan(mv_hd(1,v_i+[-vel_win_len,vel_win_len]))
        nomi = double(mv_hd(:,v_i-vel_win_len)'*mv_hd(:,v_i+vel_win_len));
        denomi = double((norm(mv_hd(:,v_i-vel_win_len))*norm(mv_hd(:,v_i+vel_win_len))));
        tolerance = 1 - cos(pi/180*0.1); % adding a 0.1 degree tolerance when calculating acos
        if nomi/denomi > 1 + tolerance
            error(sprintf('v_i = %d',v_i))
        end
        ang(v_i) = acos(nomi/denomi - tolerance)/pi*180;
        v_ang(v_i) = ang(v_i) / (0.001*velocity_smooth_win_len)/pi*180;
    end
end

%% Interpretation
% find out time point to interpt
[interp_p, rm_idx] = linear_interp_idx(t_c, EEG.times);
data_bus = [v_ang; fix_struct.eye_movement.right_ang_vel; data_struct.ori_pupil_diameter];
name_list = {'Head Rot.', 'Eye Rot.', 'Pupil left', 'Pupil right'};
interp_data = linear_interp_data(data_bus, interp_p, rm_idx);

%% add this information into EEG structure and create virtual channels
for n_i = 1:length(name_list)
    EEG = insert_data(EEG, interp_data(n_i,:), name_list{n_i});
end

%% epoching
std_epoch = pop_epoch(EEG,std_ev,len_epoch/1000, 'epochinfo', 'yes');
std_epoch = pop_rmbase(std_epoch,[len_epoch(1) 0],[]);
dev_epoch = pop_epoch(EEG,dev_ev,len_epoch/1000, 'epochinfo', 'yes');
dev_epoch = pop_rmbase(dev_epoch,[len_epoch(1) 0],[]);
ud_epoch = pop_epoch(EEG,dev_ev,len_epoch/1000, 'epochinfo', 'yes');
ud_epoch = pop_rmbase(ud_epoch,[len_epoch(1) 0],[]);
lr_epoch = pop_epoch(EEG,dev_ev,len_epoch/1000, 'epochinfo', 'yes');
lr_epoch = pop_rmbase(lr_epoch,[len_epoch(1) 0],[]);

std_time = [EEG.event(idx_std).latency];
grab_time = [EEG.event(idx_tri).latency];
miss_grab_count = length(std_time)-length(grab_time);
fprintf('Number of missing responses: %d\n', miss_grab_count);
diff_time = zeros(1,min(length(grab_time),length(std_time)));
k = 1; % counter for local_grab_time
while ~isempty(grab_time) && ~isempty(std_time)
    l_std = std_time(1);
    std_time(1) = [];
    if isempty(std_time) || std_time(1)>grab_time(1)
        diff_time(k) = grab_time(1)-l_std;
        k = k+1;
        while ~isempty(grab_time) && grab_time(1)<std_time(1)
            grab_time(1) = [];
        end
    end
end
diff_time = diff_time/EEG.srate;
fprintf('Mean reaction time: %d ms\n',round(median(diff_time*1000)));
% 631 ms for no headmovement
% 765 ms for headmovement

% grab
if any(ismember({EEG.event.type},'grab'))
    grab_epoch = pop_epoch(EEG,{'grab'},len_epoch_grab/1000, 'epochinfo', 'yes');
    grab_epoch = pop_rmbase(grab_epoch,[len_epoch_grab(1) 0],[]);
else
    grab_epoch = [];
end
fix_epoch = pop_epoch(EEG,{'fix_start'},len_epoch_grab/1000, 'epochinfo', 'yes');
fix_epoch = pop_rmbase(fix_epoch,[len_epoch_grab(1) 0],[]);
fix_std = pop_epoch(EEG,{'circle_fix_start'},len_epoch_grab/1000, 'epochinfo', 'yes');
fix_std = pop_rmbase(fix_std,[len_epoch_grab(1) 0],[]);
fix_dev = pop_epoch(EEG,{'triangle_fix_start'},len_epoch_grab/1000, 'epochinfo', 'yes');
fix_dev = pop_rmbase(fix_dev,[len_epoch_grab(1) 0],[]);
gip_std = pop_epoch(EEG,{'circle_gip_start'},len_epoch_grab/1000, 'epochinfo', 'yes');
gip_std = pop_rmbase(gip_std,[len_epoch_grab(1) 0],[]);
gip_dev = pop_epoch(EEG,{'triangle_gip_start'},len_epoch_grab/1000, 'epochinfo', 'yes');
gip_dev = pop_rmbase(gip_dev,[len_epoch_grab(1) 0],[]);


%% epoch auto rejection
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

%% output
epoch_struct = struct('std_epoch',std_epoch, 'dev_epoch',dev_epoch, 'grab_epoch',grab_epoch,...
                      'fix_epoch',fix_epoch, 'fix_std',fix_std, 'fix_dev',fix_dev, ...
                      'ud_epoch',ud_epoch, 'lr_epoch',lr_epoch, 'gip_std',gip_std, 'gip_dev',gip_dev);
end