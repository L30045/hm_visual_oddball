%% analysis pipeline
%% Pilot study for head-movement-included visual oddball
c_path = pwd;
eegpath = which('eeglab');
cd(eegpath(1:regexp(eegpath,'eeglab.m')-1))
eeglab
cd(c_path)
filepath = '../dataset/';
addpath(genpath('dependencies/'))


%% epoch EZ
subj_i = 4;
EEG_noHm = pop_loadset([filepath,sprintf('s%02d_cond1_ica_k10.set',subj_i)]);
EEG_Hm = pop_loadset([filepath,sprintf('s%02d_cond2_ica_k10.set',subj_i)]);
[data_struct_noHm, fix_struct_noHm, t_c_noHm, t_std_noHm, t_dev_noHm, epoch_struct_noHm, EEG_noHm]...
= epoch_ez([filepath,sprintf('hm_visual_oddball_s%02d_cond1.xdf',subj_i)], EEG_noHm);
[data_struct_Hm, fix_struct_Hm, t_c_Hm, t_std_Hm, t_dev_Hm, epoch_struct_Hm, EEG_Hm]...
= epoch_ez([filepath,sprintf('hm_visual_oddball_s%02d_cond2.xdf',subj_i)], EEG_Hm);

%%
dev_name_noHm = { 'Ring 0; Trial 100; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 101; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 102; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 103; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 104; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 105; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 106; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 107; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 108; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 109; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 110; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 111; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 112; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 113; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 114; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 115; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 116; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 117; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 118; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 119; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 11; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 13; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 14; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 15; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 17; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 20; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 22; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 23; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 25; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 26; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 27; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 28; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 2; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 30; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 31; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 32; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 33; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 35; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 36; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 37; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 39; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 3; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 40; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 42; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 43; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 45; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 46; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 47; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 48; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 49; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 4; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 50; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 51; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 52; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 53; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 55; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 57; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 58; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 59; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 5; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 60; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 62; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 63; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 64; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 65; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 69; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 6; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 70; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 71; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 72; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 73; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 75; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 76; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 77; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 78; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 79; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 7; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 81; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 82; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 83; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 85; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 86; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 88; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 89; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 8; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 90; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 91; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 93; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 94; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 95; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 96; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 97; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 98; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 99; Deviant cube: (0.862, 1.751, 2.345)'};
std_name_noHm = { 'Ring 0; Trial 0; Standard cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 10; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 12; Standard cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 16; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 18; Standard cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 19; Standard cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 21; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 24; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 29; Standard cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 34; Standard cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 38; Standard cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 41; Standard cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 44; Standard cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 54; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 56; Standard cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 61; Standard cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 66; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 67; Standard cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 68; Standard cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 74; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 80; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 84; Standard cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 87; Standard cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 92; Standard cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 9; Standard cube: (0.862, 1.751, 2.345)'};
dev_name_Hm = { 'Ring 1; Trial 102; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 104; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 105; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 106; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 107; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 108; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 109; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 10; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 110; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 111; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 112; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 113; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 116; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 117; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 118; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 119; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 12; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 13; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 15; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 17; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 18; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 19; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 20; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 21; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 22; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 24; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 25; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 26; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 27; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 28; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 29; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 30; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 31; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 32; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 33; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 34; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 35; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 36; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 37; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 38; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 42; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 45; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 46; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 47; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 48; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 49; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 51; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 53; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 59; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 5; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 60; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 61; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 63; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 64; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 65; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 66; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 68; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 6; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 72; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 73; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 74; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 75; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 76; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 77; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 78; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 79; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 7; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 80; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 81; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 82; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 84; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 85; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 86; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 87; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 88; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 89; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 8; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 90; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 91; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 92; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 93; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 95; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 96; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 97; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 98; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 99; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 9; Deviant cube: (0.862, 0.7712, 2.345)'};
std_name_Hm = { 'Ring 1; Trial 101; Standard cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 114; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 115; Standard cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 14; Standard cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 23; Standard cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 39; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 41; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 43; Standard cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 44; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 52; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 54; Standard cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 55; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 56; Standard cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 57; Standard cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 58; Standard cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 62; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 67; Standard cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 69; Standard cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 70; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 71; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 83; Standard cube: (0.862, 3.081, 2.345)'};

smooth_idx = 10;
% plt_EEG = pop_mergeset(epoch_struct_Hm.gip_dev, epoch_struct_Hm.gip_std);
plt_EEG = epoch_struct_noHm.dev_epoch;
tar_Ch = 'Cz';
ch_idx = find(ismember({plt_EEG.chanlocs.labels},tar_Ch));
figure; pop_erpimage(plt_EEG,1, [ch_idx],[[]],tar_Ch,smooth_idx,1,...
    {'triangle_gip_start'},[],'latency' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [ch_idx] plt_EEG.chanlocs plt_EEG.chaninfo } );

%%
% ssvep has a delay of 140ms after fixation onset so epoch from 140ms to 1140ms, no baseline
% required
% 1. PSD (FFT)
% 2. get the time course across trials
% 3. CCA by canoncorr
% 4. distinguish up/down and left/right

%% plot fixation onset time
% extract fixation epoch
plt_struct = fix_struct_noHm;
ev_duration = [-0.2 1]; % sec
srate = round(plt_struct.srate);
plt_t = [fliplr(0:-1/srate:ev_duration(1)),1/srate:1/srate:ev_duration(2)]*1000;
plt_t(1) = [];
% visualize
ev_ep = cal_epoch_ev(plt_struct.eye_fixation.eye_fix_idx, t_c_noHm, t_std_noHm, ev_duration, srate);
plt_ev = ev_ep';
[nr,nc] = size(plt_ev);
figure;
hold on
h = pcolor([plt_ev, nan(nr,1);nan(1,nc+1)]);
xline(find(plt_t==0,1,'first'),'linewidth',3)
set(h,'edgecolor','none')
set(gca,'xTick',1:5:nc)
set(gca,'xTickLabel',round(plt_t(1:5:end)))
ylabel('Trial')
xlabel('Time')
ylim([1 size(plt_ev,1)])
set(gca,'fontsize',10)
title('Fixation index - Circle')

%% visualize pupil activities (right pupil)
plt_noHm = cal_epoch_ev(fix_struct_noHm.eye_movement.right_ang_vel, t_c_noHm, [t_std_noHm, t_dev_noHm], ev_duration, srate);
plt_Hm = cal_epoch_ev(fix_struct_Hm.eye_movement.right_ang_vel, t_c_Hm, [t_std_Hm, t_dev_Hm], ev_duration, srate);
% plt_angv_std = cal_epoch_ev(fix_struct.eye_movement.right_ang_vel, t_c, t_std, ev_duration, srate);
% plt_angv_dev = cal_epoch_ev(fix_struct.eye_movement.right_ang_vel, t_c, t_dev, ev_duration, srate);

figure
hold on
grid on
% fs = shadedErrorBar(plt_t, plt_ang_std_noHm', {@nanmean, @nanstd},'lineprops','c-');
% fs.mainLine.LineWidth = 3;
% fs.mainLine.DisplayName = 'Circle (noHm)';
% fd = shadedErrorBar(plt_t, plt_ang_dev_noHm', {@nanmean, @nanstd},'lineprops','m-');
% fd.mainLine.LineWidth = 3;
% fd.mainLine.DisplayName = 'Triangle (noHm)';
fs = shadedErrorBar(plt_t, plt_noHm', {@nanmean, @nanstd},'lineprops','b-');
fs.mainLine.LineWidth = 3;
fs.mainLine.DisplayName = 'no Hm';
% fd = shadedErrorBar(plt_t, plt_Hm', {@nanmean, @nanstd},'lineprops','r-');
% fd.mainLine.LineWidth = 3;
% fd.mainLine.DisplayName = 'Hm';

ylabel('Degree/ sec')
xlabel('Time (ms)')
legend(findobj(gca,'-regexp', 'DisplayName','[^'']'))
set(gca,'fontsize',20)
title('Angular speed')

%% plot head rotation
velocity_smooth_win_len = 40;
% smoothing using moving average
mv_hd = smoothing_mv_avg(data_struct_Hm.ori_head_direct, true(1,size(data_struct_Hm.ori_head_direct,2)), 5);
% mv_hd = smoothing_mv_avg(data_struct_noHm.ori_head_direct, true(1,size(data_struct_noHm.ori_head_direct,2)), 5);
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
% check if angle are all real number
% if ~isreal(ang_Hm)
%     error('[Calculate angular velocity]: dot product of 2 head directions is greater than 1.')
% end

plt_std = cal_epoch_ev(v_ang, t_c_Hm, t_std_Hm, ev_duration, srate);
plt_dev = cal_epoch_ev(v_ang, t_c_Hm, t_dev_Hm, ev_duration, srate);
norm_hr = [plt_std,plt_dev];
norm_hr = (norm_hr-min(norm_hr,[],'all'))./...
    (max(norm_hr,[],'all')-min(norm_hr,[],'all'));

norm_er = plt_Hm;
norm_er = (norm_er-min(norm_er,[],'all'))./...
    (max(norm_er,[],'all')-min(norm_er,[],'all'));


figure
hold on
grid on
fs = shadedErrorBar(plt_t, norm_hr', {@nanmean, @nanstd},'lineprops','b-');
fs.mainLine.LineWidth = 3;
fs.mainLine.DisplayName = 'Head rot.';
fd = shadedErrorBar(plt_t, norm_er', {@nanmean, @nanstd},'lineprops','r-');
fd.mainLine.LineWidth = 3;
fd.mainLine.DisplayName = 'Pupil rot.';
% fs = shadedErrorBar(plt_t, plt_ang_std', {@nanmean, @nanstd},'lineprops','c-.');
% fs.mainLine.LineWidth = 3;
% fs.mainLine.DisplayName = 'Circle (Hm)';
% fd = shadedErrorBar(plt_t, plt_ang_dev', {@nanmean, @nanstd},'lineprops','m-.');
% fd.mainLine.LineWidth = 3;
% fd.mainLine.DisplayName = 'Triangle (Hm)';

ylabel('Normalized degree/ sec')
xlabel('Time (ms)')
legend(findobj(gca,'-regexp', 'DisplayName','[^'']'))
set(gca,'fontsize',20)
title('Rotation timing comparison')

%%
% ERPImage of pupil data (stimulus-locked, response-locked, fixation-locked)
% plot GIP-locked event
















