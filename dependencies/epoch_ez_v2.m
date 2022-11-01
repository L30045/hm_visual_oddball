function [epoch_struct, EEG] = epoch_ez_v2(filename, varargin)
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
%           (11) event_time: Event onset time and epoch length. The order
%                are the same as epoch_struct.


%% parameter setting
% epoch length
len_epoch = [-200 2000]; % epoch length for stimulus-lock epoch
len_epoch_grab = [-1000 1000]; % epoch length for GIP/Response/Fixation-lock epoch

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
[data_struct, EEG_ori, ~, ~, ~, ~] = load_eyetracking_hm(filename);
if isempty(varargin)
    EEG = EEG_ori;
else
    EEG = varargin{1};
end
nbchan = EEG.nbchan;

% % Prevent event marker error 
condname = EEG.filename(5);
if EEG.filename(3)=='1'
    switch condname
        case 'I'
            reg_txt = 'Ring 0';
        case 'O'
            reg_txt = 'Ring 1';
    end
else
    reg_txt = 'Ring 0';
end

%% find out stimulus onset time
if ~strcmp(EEG.filename(1),'s')
    idx_dev = cellfun(@(x) ~isempty(regexp(x, '\(R\)', 'ONCE')), {EEG.event.type});
    idx_std = cellfun(@(x) ~isempty(regexp(x, '\(L\)', 'ONCE')), {EEG.event.type});
else
    idx_cond = cellfun(@(x) ~isempty(regexp(x, reg_txt, 'ONCE')), {EEG.event.type}); 
    idx_dev = cellfun(@(x) ~isempty(regexp(x, 'Deviant', 'ONCE')), {EEG.event.type});
    idx_std = cellfun(@(x) ~isempty(regexp(x, 'Standard', 'ONCE')), {EEG.event.type});
    idx_std = idx_std & idx_cond;
    idx_dev = idx_dev & idx_cond;
end
if sum(idx_std|idx_dev)~=100
    error('[Trial Number Incorrect]: Total number of trials is not 100.')
end
idx_tri = cellfun(@(x) ~isempty(regexp(x, 'Trigger', 'ONCE')), {EEG.event.type});
% idx_tri(find(idx_tri,2,'last')) = false;
idx_blue = cellfun(@(x) ~isempty(regexp(x, 'blue_cube', 'ONCE')), {EEG.event.type});
std_ev = {EEG.event(idx_std).type};
dev_ev = {EEG.event(idx_dev).type};
std_time = [EEG.event(idx_std).latency]/EEG.srate*1000; % msec
dev_time = [EEG.event(idx_dev).latency]/EEG.srate*1000; % msec

%% calculate fixation
% fix_struct = cal_fix_pupil(test_data,data_struct.srate);
test_data = [data_struct.ori_eye_3D_pos;data_struct.eye_openess_idx];
fix_struct = cal_fix_pupil(test_data,data_struct.srate,'blink_length',len_blink,...
                           'thres_ang',thres_ang,'thres_ang_v',thres_ang_v,'fix_selection',fix_selection,...
                           'max_fix_interval',max_fix_interval,'max_fix_ang',max_fix_ang,'min_fix_len',min_fix_len);
eye_fix_idx = fix_struct.eye_fixation.eye_fix_idx;

%% synchronize eye gaze stream and EEG stream
t_c = fix_struct.time_stamps*1000; % change unit to ms
t_std_ori = ([EEG.event(idx_std).latency]-1)/EEG.srate*1000; % convert from sample point to time (ms) 
t_dev_ori = ([EEG.event(idx_dev).latency]-1)/EEG.srate*1000; % convert from sample point to time (ms)
t_blue_ori = ([EEG.event(idx_blue).latency]-1)/EEG.srate*1000; % convert from sample point to time (ms)

%% find up/down left/right event
%boolean arrays where the LEFT/RIGHT or UP/DOWN event markers occured (I assumed
%they were mutually trail exclusive i.e. both do not occur in the same trail)
if ~strcmp(EEG.filename(1),'s')
    up_idx = cellfun(@(x) ~isempty(regexp(x,'index\(.\)\:\ 1','ONCE')),{EEG.event.type});
    down_idx = cellfun(@(x) ~isempty(regexp(x,'index\(.\)\:\ 3','ONCE')),{EEG.event.type});
    left_idx = cellfun(@(x) ~isempty(regexp(x,'index\(.\)\:\ 2','ONCE')),{EEG.event.type});
    right_idx = cellfun(@(x) ~isempty(regexp(x,'index\(.\)\:\ 0','ONCE')),{EEG.event.type});
    % find target loc
    key = 'Position: (';
    tmp_up = EEG.event(find(up_idx,1)).type;
    upLoc = str2double(split(tmp_up(regexp(tmp_up,key,'once')+length(key):end-1)));
    tmp_down = EEG.event(find(down_idx,1)).type;
    downLoc = str2double(split(tmp_down(regexp(tmp_down,key,'once')+length(key):end-1)));
    tmp_left = EEG.event(find(left_idx,1)).type;
    leftLoc = str2double(split(tmp_left(regexp(tmp_left,key,'once')+length(key):end-1)));
    tmp_right = EEG.event(find(right_idx,1)).type;
    rightLoc = str2double(split(tmp_right(regexp(tmp_right,key,'once')+length(key):end-1)));
else
    % gather 4 location
    tar_ev = unique({EEG.event(cellfun(@(x) ~isempty(regexp(x,reg_txt,'ONCE')),{EEG.event.type})).type});
    tar_ev = cellfun(@split ,unique(cellfun(@(x) x(regexp(x,'(')+1:end-1), tar_ev, 'uniformoutput',0)),'uniformoutput',0);
    num_loc = zeros(3,4);
    for i = 1:4
        tmp = tar_ev{i};
        num_loc(1,i) = str2double(tmp{1}(1:end-1));
        num_loc(2,i) = str2double(tmp{2}(1:end-1));
        num_loc(3,i) = str2double(tmp{3});
    end
    num_loc = sortrows(num_loc');

    upLoc = num_loc(3,:);
    downLoc = num_loc(2,:);
    leftLoc = num_loc(1,:);
    rightLoc = num_loc(4,:);
    up_idx = cellfun(@(x) ~isempty(regexp(x,sprintf('(%.4g, %.4g, %.4g)',upLoc),'match','ONCE')),{EEG.event.type});
    down_idx = cellfun(@(x) ~isempty(regexp(x,sprintf('(%.4g, %.4g, %.4g)',downLoc),'match','ONCE')),{EEG.event.type});
    left_idx = cellfun(@(x) ~isempty(regexp(x,sprintf('(%.4g, %.4g, %.4g)',leftLoc),'match','ONCE')),{EEG.event.type});
    right_idx = cellfun(@(x) ~isempty(regexp(x,sprintf('(%.4g, %.4g, %.4g)',rightLoc),'match','ONCE')),{EEG.event.type});
end

std_up = up_idx(idx_std);
std_down = down_idx(idx_std);
std_left = left_idx(idx_std);
std_right = right_idx(idx_std);
dev_up = up_idx(idx_dev);
dev_down = down_idx(idx_dev);
dev_left = left_idx(idx_dev);
dev_right = right_idx(idx_dev);

% check if all trials have been assigned a direcetion
if ~all(std_up|std_down|std_left|std_right)
    error('[Trial Information Missing]: Some trials do not have the direction information. (std)')
end
if ~all(dev_up|dev_down|dev_left|dev_right)
    error('[Trial Information Missing]: Some trials do not have the direction information. (dev)')
end


%% find GIP onset time
% find the first GIP after stimulus
% find out stimulus onset time
[gip_up_idx, gip_down_idx, gip_left_idx, gip_right_idx,miss_idx] = find_1st_gip(EEG,idx_std);
gip_1st_std = find(gip_up_idx|gip_down_idx|gip_left_idx|gip_right_idx);
miss_idx_std = find(miss_idx);
[gip_up_idx, gip_down_idx, gip_left_idx, gip_right_idx,miss_idx] = find_1st_gip(EEG,idx_dev);
gip_1st_dev = find(gip_up_idx|gip_down_idx|gip_left_idx|gip_right_idx);
miss_idx_dev = find(miss_idx);
f_std = find(idx_std);
f_dev = find(idx_dev);


gipStd_time = zeros(1,length(f_std));
count = 1;
for i = 1:length(f_std)
    if ~ ismember(f_std(i), miss_idx_std)
        t_f = gip_1st_std(count);
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'circle_gip_start';
        gipStd_time(i) = EEG.event(t_f).latency/EEG.srate*1000; % msec
        EEG.event(len_event+1).latency = EEG.event(t_f).latency; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
        count = count+1;
    else
        gipStd_time(i) = NaN;
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'circle_gip_missing';
        EEG.event(len_event+1).latency = EEG.event(f_std(i)).latency; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
    end
end


gipDev_time = zeros(1,length(f_dev));
count = 1;
for i = 1:length(f_dev)
    if ~ismember(f_dev(i), miss_idx_dev)
        t_f = gip_1st_dev(count);
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'triangle_gip_start';
        gipDev_time(i) = EEG.event(t_f).latency/EEG.srate*1000; % msec
        EEG.event(len_event+1).latency = EEG.event(t_f).latency; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
        count = count+1;
    else
        gipDev_time(i) = NaN;
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'triangle_gip_missing';
        EEG.event(len_event+1).latency = EEG.event(f_dev(i)).latency; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
    end
end

%% adding event marker to EEG structure
% adding grab marker
grab_1st = find_1st_grab(EEG,find(idx_tri),f_std); % event index
grab_1st_lat = nan(size(grab_1st));
grab_time = nan(size(grab_1st));
for i = 1:length(grab_1st)
    if ~isnan(grab_1st)
        grab_1st_lat(i) = EEG.event(grab_1st(i)).latency;
        grab_time(i) = grab_1st_lat(i)/EEG.srate*1000;
    end
end
for i = 1:length(grab_1st)
    if ~isnan(grab_1st(i))
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'grab';
        EEG.event(len_event+1).latency = grab_1st_lat(i); 
        EEG.event(len_event+1).urevent = len_event+1;
    else
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'grab missing';
        EEG.event(len_event+1).latency = f_std(i); 
        EEG.event(len_event+1).urevent = len_event+1;
    end
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
t_f_std = find_1st_fix(fix_start,t_std_ori); % fixation time, nan if missing
t_f_dev = find_1st_fix(fix_start,t_dev_ori);

for i = 1:length(t_f_std)
    if ~isnan(t_f_std(i))
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'circle_fix_start';
        EEG.event(len_event+1).latency = t_f_std(i)/1000*EEG.srate; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
    else
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'circle_fix_missing';
        EEG.event(len_event+1).latency = t_std_ori(i)/1000*EEG.srate; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
    end
end

for i = 1:length(t_f_dev)
    if ~isnan(t_f_dev(i))
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'triangle_fix_start';
        EEG.event(len_event+1).latency = t_f_dev(i)/1000*EEG.srate; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
    else
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'triangle_fix_missing';
        EEG.event(len_event+1).latency = t_dev_ori(i)/1000*EEG.srate; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
    end
end

t_f_blue = find_1st_fix(fix_start,t_blue_ori);
for i = 1:length(t_f_blue)
    if ~isnan(t_f_blue(i))
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'blue_fix_start';
        EEG.event(len_event+1).latency = t_f_blue(i)/1000*EEG.srate; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
    else
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = 'blue_fix_missing';
        EEG.event(len_event+1).latency = t_blue_ori(i)/1000*EEG.srate; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
    end
end

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
std_epoch = pop_epoch(EEG,std_ev,len_epoch/1000, 'epochinfo', 'yes');
dev_epoch = pop_epoch(EEG,dev_ev,len_epoch/1000, 'epochinfo', 'yes');

% grab
if any(ismember({EEG.event.type},'grab'))
    grab_epoch = pop_epoch(EEG,{'grab'},len_epoch_grab/1000, 'epochinfo', 'yes');
else
    grab_epoch = [];
end

gip_std = pop_epoch(EEG,{'circle_gip_start'},len_epoch_grab/1000, 'epochinfo', 'yes');
gip_dev = pop_epoch(EEG,{'triangle_gip_start'},len_epoch_grab/1000, 'epochinfo', 'yes');

%% remove baseline
% save behavioral data
behavi_std = std_epoch.data(nbchan+1:end,:,:);
behavi_dev = dev_epoch.data(nbchan+1:end,:,:);
behavi_gstd = gip_std.data(nbchan+1:end,:,:);
behavi_gdev = gip_dev.data(nbchan+1:end,:,:);
% remove baseline
std_base = mean(std_epoch.data(:,1:find(std_epoch.times==0,1),:),2);
gip_count = 1;
for t_i = 1:size(std_epoch.data,3)
    std_epoch.data(:,:,t_i) = std_epoch.data(:,:,t_i) - std_base(:,:,t_i);
    if ~isnan(gipStd_time(t_i))
        gip_std.data(:,:,gip_count) = gip_std.data(:,:,gip_count) - std_base(:,:,t_i);
        gip_count = gip_count+1;
    end
end
dev_base = mean(dev_epoch.data(:,1:find(dev_epoch.times==0,1),:),2);
gip_count = 1;
for t_i = 1:size(dev_epoch.data,3)
    dev_epoch.data(:,:,t_i) = dev_epoch.data(:,:,t_i) - dev_base(:,:,t_i);
    if ~isnan(gipDev_time(t_i))
        gip_dev.data(:,:,gip_count) = gip_dev.data(:,:,gip_count) - dev_base(:,:,t_i);
        gip_count = gip_count+1;
    end
end    
% restore behavioral data
std_epoch.data(nbchan+1:end,:,:) = behavi_std;
dev_epoch.data(nbchan+1:end,:,:) = behavi_dev;
gip_std.data(nbchan+1:end,:,:) = behavi_gstd;
gip_dev.data(nbchan+1:end,:,:) = behavi_gdev;
% some data missing grab event markers
if ~isempty(grab_epoch)
    behavi_grab = grab_epoch.data(nbchan+1:end,:,:);
    % remove baseline
    grab_count = 1;
    for t_i = 1:size(grab_epoch.data,3)
        if ~isnan(grab_time(t_i))
            grab_epoch.data(:,:,t_i) = grab_epoch.data(:,:,t_i) - std_base(:,:,t_i);
            grab_count = grab_count+1;
        end
    end
    grab_epoch.data(nbchan+1:end,:,:) = behavi_grab;
end

%% fix related epoch
if all(ismember({'circle_fix_start','triangle_fix_start'},{EEG.event.type}))
    fix_epoch = pop_epoch(EEG,{'fix_start'},len_epoch_grab/1000, 'epochinfo', 'yes');
    fix_std = pop_epoch(EEG,{'circle_fix_start'},len_epoch_grab/1000, 'epochinfo', 'yes');
    fix_dev = pop_epoch(EEG,{'triangle_fix_start'},len_epoch_grab/1000, 'epochinfo', 'yes');
    fix_blue = pop_epoch(EEG,{'blue_fix_start'},len_epoch_grab/1000, 'epochinfo', 'yes');
    behavi_fixAll = fix_epoch.data(nbchan+1:end,:,:);
    behavi_fstd = fix_std.data(nbchan+1:end,:,:);
    behavi_fdev = fix_dev.data(nbchan+1:end,:,:);
    behavi_fblue = fix_blue.data(nbchan+1:end,:,:);
    % remove baseline
    fix_epoch = pop_rmbase(fix_epoch,[max(len_epoch_grab(1),fix_epoch.times(1)) 0],[]);
    fix_blue = pop_rmbase(fix_blue,[max(len_epoch_grab(1),fix_blue.times(1)) 0],[]);
    % remove baseline
    std_base = mean(std_epoch.data(:,1:find(std_epoch.times==0,1),:),2);
    fix_count = 1;
    for t_i = 1:size(std_epoch.data,3)
        if ~isnan(t_f_std(t_i))
            fix_std.data(:,:,fix_count) = fix_std.data(:,:,fix_count) - std_base(:,:,t_i);
            fix_count = fix_count+1;
        end
    end
    dev_base = mean(dev_epoch.data(:,1:find(dev_epoch.times==0,1),:),2);
    fix_count = 1;
    for t_i = 1:size(dev_epoch.data,3)
        if ~isnan(t_f_dev(t_i))
            fix_dev.data(:,:,fix_count) = fix_dev.data(:,:,fix_count) - dev_base(:,:,t_i);
            fix_count = fix_count+1;
        end
    end    
    fix_epoch.data(nbchan+1:end,:,:) = behavi_fixAll;
    fix_std.data(nbchan+1:end,:,:) = behavi_fstd;
    fix_dev.data(nbchan+1:end,:,:) = behavi_fdev;
    fix_blue.data(nbchan+1:end,:,:) = behavi_fblue;
else
    fix_epoch = [];
    fix_std = [];
    fix_dev = [];
    fix_blue = [];
end

%% record event time
diff_gip_std = gipStd_time - std_time;
diff_gip_dev = gipDev_time - dev_time;
diff_grab_std = grab_time - std_time;
fixAll_time = fix_start; % msec
fixStd_time = t_f_std; % msec
fixDev_time = t_f_dev; % msec
fixBlue_time = t_f_blue; % msec
event_time = struct('std_time',std_time,'dev_time',dev_time,'grab_time',grab_time,...
                    'diff_stim_grab',diff_grab_std,'diff_gip_std',diff_gip_std,'diff_gip_dev',diff_gip_dev,...
                    'fixAll_time',fixAll_time,'fixStd_time',fixStd_time,'fixDev_time',fixDev_time,'fixBlue_time',fixBlue_time,...
                    'gipStd_time',gipStd_time,'gipDev_time',gipDev_time,...
                    'std_up',std_up,'std_down',std_down,'std_left',std_left,'std_right',std_right,...
                    'dev_up',dev_up,'dev_down',dev_down,'dev_left',dev_left,'dev_right',dev_right,...
                    'upLoc',upLoc,'downLoc',downLoc,'leftLoc',leftLoc,'rightLoc',rightLoc);

%% output
epoch_struct = struct('std_epoch',std_epoch, 'dev_epoch',dev_epoch, 'grab_epoch',grab_epoch,...
                      'fix_epoch',fix_epoch, 'fix_std',fix_std, 'fix_dev',fix_dev, 'fix_blue',fix_blue,...
                      'gip_std',gip_std, 'gip_dev',gip_dev,'event_time',event_time);
                  
end