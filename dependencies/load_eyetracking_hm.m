function [data_struct, EEG, s_eyeMarker, s_eyeGaze, s_event, s_EEG] = load_eyetracking_hm(filename)
%% load eye tracking data
% load head rotation stream
streams = load_xdf(filename);

%% extract streams
s_eyeMarker = streams{cellfun(@(x) strcmp(x.info.name,'ProEyeMarker'), streams)};
s_eyeGaze = streams{cellfun(@(x) strcmp(x.info.name,'ProEyeGaze'), streams)};
s_event = streams{cellfun(@(x) strcmp(x.info.name,'EventMarker'), streams)};
s_EEG = streams{cellfun(@(x) strcmp(x.info.name,'EEG'), streams)};
EEG = pop_loadxdf(filename,'streamname','EEG');

%% ======= data stream in eye stream (39 Chs) ======== 05/31/2021
% // 0 - 2d coordinate of left eye
% // 2 - 2d coordinate of right eye
% // 4 - 3d direction of left eye
% // 7 - 3d direction of right eye
% // 10 - 3d position of combined hit spot
% // 13 - 3d position of head
% // 16 - 3d forward direction of head
% // 19 - 3d velocity of head
% // 22 - 3d angular velocity of head
% // 25 - left eye openness
% // 26 - right eye openness
% // 27 - 3d position of chest IMU
% // 30 - 3d forward direction of chest IMU
% // 33 - 3d velocity of chest IMU
% // 36 - heatmap value (right now, the heat map is raw values)d
% // 37 - left eye pupil diameter
% // 38 - right eye pupil diameter


%% extract data from eye Gaze stream
seg_range = s_eyeGaze.segments(1).index_range;
ori_eye_2D_pos = s_eyeGaze.time_series(1:4,seg_range(1):seg_range(2)); % left_xy, right_xy
ori_eye_3D_pos = s_eyeGaze.time_series(5:10,seg_range(1):seg_range(2)); % left_xyz, right_xyz
ori_gip_3D_pos = s_eyeGaze.time_series(11:13,seg_range(1):seg_range(2));
ori_head_3D_loc = s_eyeGaze.time_series(14:16,seg_range(1):seg_range(2));
ori_head_direct = s_eyeGaze.time_series(17:19,seg_range(1):seg_range(2));
ori_head_vel = s_eyeGaze.time_series(20:22,seg_range(1):seg_range(2));
ori_head_rot = s_eyeGaze.time_series(23:25,seg_range(1):seg_range(2));
eye_openess_idx = s_eyeGaze.time_series(26:27,seg_range(1):seg_range(2));
ori_pupil_diameter = s_eyeGaze.time_series(38:39,seg_range(1):seg_range(2)); % left, right
eyeMarker_name = s_eyeMarker.time_series;
pt_eg = s_eyeGaze.time_stamps(seg_range(1):seg_range(2))-s_eyeGaze.time_stamps(seg_range(1));
pt_em = s_eyeMarker.time_stamps-s_eyeGaze.time_stamps(seg_range(1));
pt_em(pt_em > max(pt_eg)) = [];
% round up eye marker stream time stamps with eye gaze stream time stamps
for t_i = 1:length(pt_em)
    pt_em(t_i) = pt_eg(find(pt_eg>pt_em(t_i),1,'first'));
end
pt_event = s_event.time_stamps-s_eyeGaze.time_stamps(seg_range(1));
pt_event(pt_event > max(pt_eg)) = [];
% round up grab stream time stamps with eye gaze stream time stamps
for t_i = 1:length(pt_event)
    pt_event(t_i) = pt_eg(find(pt_eg>pt_event(t_i),1,'first'));
end
% sampling rate
srate = 1/mean(diff(pt_eg));

%% output struct
data_struct = struct('ori_eye_2D_pos',ori_eye_2D_pos,'ori_eye_3D_pos',ori_eye_3D_pos,...
    'ori_gip_3D_pos',ori_gip_3D_pos,'ori_head_3D_loc',ori_head_3D_loc,'ori_head_direct',ori_head_direct,...
    'ori_head_vel',ori_head_vel,'ori_head_rot',ori_head_rot,'eye_openess_idx',eye_openess_idx,...
    'ori_pupil_diameter',ori_pupil_diameter,'eyeMarker_name',{eyeMarker_name},...
    'pt_eg',pt_eg,'pt_em',pt_em,'pt_event',pt_event,'srate',srate);

end