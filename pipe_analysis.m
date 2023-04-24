%% analysis pipeline
%% Pilot study for head-movement-included visual oddball
c_path = pwd;
eegpath = which('eeglab');
cd(eegpath(1:regexp(eegpath,'eeglab.m')-1))
eeglab
cd(c_path)
filepath = '/home/yuan/Documents/2021 HM_visual_oddball/dataset/';
savepath = [filepath,'new epoch/']; % Correct the synchronization issue in epoch_ez and create new epoch.
savepath = 'C:\Users\Yuan\OneDrive\Desktop\visualOddball-NicoleXinDataCollection/new epoch/';
filepath = 'C:\Users\Yuan\OneDrive\Desktop\visualOddball-NicoleXinDataCollection\preproc_data/';
addpath(genpath('dependencies/'))


%% epoch EZ
% subj_list = [1,4:6,8:10];
% % subj_i = 4;
% for subj_i = subj_list
%     EEG_noHm = pop_loadset([filepath,sprintf('s%02d_cond1_ica_k10.set',subj_i)]);
%     EEG_Hm = pop_loadset([filepath,sprintf('s%02d_cond2_ica_k10.set',subj_i)]);
%     [data_struct_noHm, fix_struct_noHm, t_c_noHm, t_std_noHm, t_dev_noHm, epoch_struct_noHm, EEG_noHm]...
%     = epoch_ez([filepath,sprintf('hm_visual_oddball_s%02d_cond1.xdf',subj_i)], EEG_noHm);
%     [data_struct_Hm, fix_struct_Hm, t_c_Hm, t_std_Hm, t_dev_Hm, epoch_struct_Hm, EEG_Hm]...
%     = epoch_ez([filepath,sprintf('hm_visual_oddball_s%02d_cond2.xdf',subj_i)], EEG_Hm);
%     save([savepath,sprintf('s%02d_epoch.mat',subj_i)],'-v7.3','epoch_struct_noHm','epoch_struct_Hm');
% end
% disp('Done')

filename = '2004_Oddball_1163';
loadpath = ['C:\Users\Yuan\OneDrive\Desktop\visualOddball-NicoleXinDataCollection\sub',filename(end-3:end),'\'];
EEG_noHm = pop_loadset([filepath,sprintf('%s_Inner.set',filename)]);
EEG_Hm = pop_loadset([filepath,sprintf('%s_Outer.set',filename)]);
[~, ~, ~, ~, ~, epoch_struct_noHm, ~]...
= epoch_ez([loadpath,sprintf('%s_Inner.xdf',filename)], EEG_noHm);
[~, ~, ~, ~, ~, epoch_struct_Hm, ~]...
= epoch_ez([loadpath,sprintf('%s_Outer.xdf',filename)], EEG_Hm);
save([savepath,'new_s03_epoch.mat'],'-v7.3','epoch_struct_noHm','epoch_struct_Hm');
disp('Done')

%% calculate head rotation velocity and eye rotation velocity
cond_i = 1;
switch cond_i
    case 1
        epoch_struct = epoch_struct_noHm;
        
    case 2
        epoch_struct = epoch_struct_Hm;
end
upLoc = epoch_struct.event_time.upLoc;
downLoc = epoch_struct.event_time.downLoc;
leftLoc = epoch_struct.event_time.leftLoc;
rightLoc = epoch_struct.event_time.rightLoc;
tar_lib = [upLoc;downLoc;leftLoc;rightLoc];
ev_list = fieldnames(epoch_struct);
ev_list([4:6,9]) = [];
ev_direct = {[epoch_struct.event_time.std_up;epoch_struct.event_time.std_down;...
              epoch_struct.event_time.std_left; epoch_struct.event_time.std_right];...
             [epoch_struct.event_time.dev_up;epoch_struct.event_time.dev_down;...
              epoch_struct.event_time.dev_left; epoch_struct.event_time.dev_right]};
dir_idx = [1,2,1,1,2];
    
e_i = 5;
switch e_i
    case 1
        tname = 'Stim Lock (Circle)';
        tlock = 'Stim';
    case 2
        tname = 'Stim Lock (Triangle)';
        tlock = 'Stim';
    case 3
        tname = 'Response Lock (Circle)';
        tlock = 'Response';
    case 4
        tname = 'GIP Lock (Circle)';
        tlock = 'GIP';
    case 5
        tname = 'GIP Lock (Triangle)';
        tlock = 'GIP';
end
        
tar_epoch = epoch_struct.(ev_list{e_i});
nbchan = find(ismember({tar_epoch.chanlocs.labels},'HeadLoc_x'));
tar_direct = ev_direct{dir_idx(e_i)};
if e_i == 4
    s_t = isnan(epoch_struct.event_time.diff_gip_std);
elseif e_i == 5
    s_t = isnan(epoch_struct.event_time.diff_gip_dev);
else
    s_t = false(1,size(tar_direct,2));
end
tar_direct(:,s_t) = [];
headLoc = tar_epoch.data(nbchan:nbchan+2,:,:);
headDirect = tar_epoch.data(nbchan+3:nbchan+5,:,:);
GIP = tar_epoch.data(nbchan+6:nbchan+8,:,:);
%     blink_idx = tar_epoch.data(nbchan+9,:,:);
%     dataLose_idx = tar_epoch.data(nbchan+10,:,:);
eyeDirect = tar_epoch.data(nbchan+11:nbchan+13,:,:);
% calculate GIP distance
dist_gip2box = dist2Box(GIP, headLoc, tar_lib, tar_direct);

% calculate rotation
[headAng, headRot, headAngDiff, headAngCumsum] = cal_rot(headDirect, tar_direct, tar_epoch.srate);
[eyeAng, eyeRot, eyeAngDiff, eyeAngCumsum] = cal_rot(eyeDirect, tar_direct, tar_epoch.srate);
[gipAng, gipRot, gipAngDiff, gipAngCumsum] = cal_rot(GIP, tar_direct, tar_epoch.srate);


%% sanity check
shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(length(subj_list)))};
figure;
% normalized distance
plt_dist = (dist_gip2box-min(dist_gip2box,[],1))./(max(dist_gip2box,[],1)+min(dist_gip2box,[],1));
scale = 1; %max(mean(eyeAngDiff,2));
plt_dist = plt_dist * scale;
shadedErrorBar(tar_epoch.times,eyeAngDiff', shaded_method, 'lineprops',{'b-','DisplayName','EyeAngDiff','linewidth',3});
grid on; hold on;
shadedErrorBar(tar_epoch.times,headAngDiff', shaded_method, 'lineprops',{'r-','DisplayName','HeadAngDiff','linewidth',3})
shadedErrorBar(tar_epoch.times,gipAngDiff', shaded_method, 'lineprops',{'k-','DisplayName','GIPAngDiff','linewidth',3})
shadedErrorBar(tar_epoch.times,plt_dist', shaded_method, 'lineprops',{'g-','DisplayName','dist2box','linewidth',3})
xline(0,'k--','linewidth',3,'DisplayName',tlock);
title(tname);
set(gca,'fontsize',20)
xlabel('Time (ms)')
ylabel('Angle (deg)')
legend(findobj(gca,'-regexp','DisplayName', '[^'']'));

figure;
plt_dist = (dist_gip2box-min(dist_gip2box,[],1))./(max(dist_gip2box,[],1)+min(dist_gip2box,[],1));
scale = max(mean(eyeAng,2,'omitnan'));
plt_dist = plt_dist * scale;
shadedErrorBar(tar_epoch.times,eyeAng', shaded_method, 'lineprops',{'b-','DisplayName','EyeAng','linewidth',3});
grid on; hold on;
shadedErrorBar(tar_epoch.times,headAng', shaded_method, 'lineprops',{'r-','DisplayName','HeadAng','linewidth',3})
shadedErrorBar(tar_epoch.times,gipAng', shaded_method, 'lineprops',{'k-','DisplayName','GIPAng','linewidth',3})
shadedErrorBar(tar_epoch.times,plt_dist', shaded_method, 'lineprops',{'g-','DisplayName','dist2box','linewidth',3})
xline(0,'k--','linewidth',3,'DisplayName',tlock);
title(tname);
set(gca,'fontsize',20)
xlabel('Time (ms)')
ylabel('Angle (deg)')
legend(findobj(gca,'-regexp','DisplayName', '[^'']'));


%%
dev_name_noHm = { 'Ring 0; Trial 100; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 101; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 102; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 103; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 104; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 105; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 106; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 107; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 108; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 109; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 110; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 111; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 112; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 113; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 114; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 115; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 116; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 117; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 118; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 119; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 11; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 13; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 14; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 15; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 17; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 20; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 22; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 23; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 25; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 26; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 27; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 28; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 2; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 30; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 31; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 32; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 33; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 35; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 36; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 37; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 39; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 3; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 40; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 42; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 43; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 45; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 46; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 47; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 48; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 49; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 4; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 50; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 51; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 52; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 53; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 55; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 57; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 58; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 59; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 5; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 60; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 62; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 63; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 64; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 65; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 69; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 6; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 70; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 71; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 72; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 73; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 75; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 76; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 77; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 78; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 79; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 7; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 81; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 82; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 83; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 85; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 86; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 88; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 89; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 8; Deviant cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 90; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 91; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 93; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 94; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 95; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 96; Deviant cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 97; Deviant cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 98; Deviant cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 99; Deviant cube: (0.862, 1.751, 2.345)'};
std_name_noHm = { 'Ring 0; Trial 0; Standard cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 10; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 12; Standard cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 16; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 18; Standard cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 19; Standard cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 21; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 24; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 29; Standard cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 34; Standard cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 38; Standard cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 41; Standard cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 44; Standard cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 54; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 56; Standard cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 61; Standard cube: (0.687, 1.926, 2.345)' 'Ring 0; Trial 66; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 67; Standard cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 68; Standard cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 74; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 80; Standard cube: (0.862, 1.751, 2.345)' 'Ring 0; Trial 84; Standard cube: (0.862, 2.101, 2.345)' 'Ring 0; Trial 87; Standard cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 92; Standard cube: (1.037, 1.926, 2.345)' 'Ring 0; Trial 9; Standard cube: (0.862, 1.751, 2.345)'};
dev_name_Hm = { 'Ring 1; Trial 102; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 104; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 105; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 106; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 107; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 108; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 109; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 10; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 110; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 111; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 112; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 113; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 116; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 117; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 118; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 119; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 12; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 13; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 15; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 17; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 18; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 19; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 20; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 21; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 22; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 24; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 25; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 26; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 27; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 28; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 29; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 30; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 31; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 32; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 33; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 34; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 35; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 36; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 37; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 38; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 42; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 45; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 46; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 47; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 48; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 49; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 51; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 53; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 59; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 5; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 60; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 61; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 63; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 64; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 65; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 66; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 68; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 6; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 72; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 73; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 74; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 75; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 76; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 77; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 78; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 79; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 7; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 80; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 81; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 82; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 84; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 85; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 86; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 87; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 88; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 89; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 8; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 90; Deviant cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 91; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 92; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 93; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 95; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 96; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 97; Deviant cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 98; Deviant cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 99; Deviant cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 9; Deviant cube: (0.862, 0.7712, 2.345)'};
std_name_Hm = { 'Ring 1; Trial 101; Standard cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 114; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 115; Standard cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 14; Standard cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 23; Standard cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 39; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 41; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 43; Standard cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 44; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 52; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 54; Standard cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 55; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 56; Standard cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 57; Standard cube: (-0.2928, 1.926, 2.345)' 'Ring 1; Trial 58; Standard cube: (2.017, 1.926, 2.345)' 'Ring 1; Trial 62; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 67; Standard cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 69; Standard cube: (0.862, 0.7712, 2.345)' 'Ring 1; Trial 70; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 71; Standard cube: (0.862, 3.081, 2.345)' 'Ring 1; Trial 83; Standard cube: (0.862, 3.081, 2.345)'};

smooth_idx = 10;
% plt_EEG = pop_mergeset(epoch_struct_Hm.gip_dev, epoch_struct_Hm.gip_std);
plt_EEG = epoch_struct_noHm.std_epoch;
tar_Ch = 'Cz';
ch_idx = find(ismember({plt_EEG.chanlocs.labels},tar_Ch));
figure; pop_erpimage(plt_EEG,1, [ch_idx],[[]],tar_Ch,smooth_idx,1,...
    {'circle_gip_start'},[],'latency' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [ch_idx] plt_EEG.chanlocs plt_EEG.chaninfo } );

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
srate = EEG_noHm.srate;
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

%% compare standard and deviant event
subj_i = 1;
load(sprintf('s%02d_epoch.mat',subj_i));
cond_name = 'Hm';
eval(sprintf('cond_struct = epoch_struct_%s;',cond_name));
% get mean gip time
tri_epoch = cond_struct.dev_epoch;
cir_epoch = cond_struct.std_epoch;
tri_lat = nan(1,tri_epoch.trials);
cir_lat = nan(1,cir_epoch.trials);
for e_i = 1:length(tri_epoch.epoch)
    if any(ismember(tri_epoch.epoch(e_i).eventtype, 'triangle_gip_start'))
        tri_lat(e_i) = tri_epoch.epoch(e_i).eventlatency{find(ismember(tri_epoch.epoch(e_i).eventtype, 'triangle_gip_start'),1)};
    end
end
for e_i = 1:length(cir_epoch.epoch)
    if any(ismember(cir_epoch.epoch(e_i).eventtype, 'circle_gip_start'))
        cir_lat(e_i) = cir_epoch.epoch(e_i).eventlatency{find(ismember(cir_epoch.epoch(e_i).eventtype, 'circle_gip_start'),1)};
    end
end
figure
histogram(tri_lat,'binwidth',50)
title(sprintf('Triangle (%d)',sum(~isnan(tri_lat))))
figure
histogram(cir_lat,'binwidth',50)
title(sprintf('Circle (%d)',sum(~isnan(cir_lat))))
switch cond_name
    case 'Hm'
        tri_lat_Hm = tri_lat;
        cir_lat_Hm = cir_lat;
        save([filepath,sprintf('s%02d_epoch.mat',subj_i)],'tri_lat_Hm','cir_lat_Hm','-append')
    case 'noHm'
        tri_lat_noHm = tri_lat;
        cir_lat_noHm = cir_lat;
        save([filepath,sprintf('s%02d_epoch.mat',subj_i)],'tri_lat_noHm','cir_lat_noHm','-append')
end

%%
ev_name = 'gip';
tar_Ch = 'Cz';
switch ev_name
    case 'stim'
        tri_epoch = cond_struct.dev_epoch;
        cir_epoch = cond_struct.std_epoch;
        plt_t = tri_epoch.times(end);
        ch_idx = find(ismember({tri_epoch.chanlocs.labels},tar_Ch));
        plt_idx = find(tri_epoch.times >= plt_t, 1)-1;
        plt_t = tri_epoch.times(1:plt_idx);
        tri_data = tri_epoch.data(ch_idx,1:plt_idx,:);
        cir_data = cir_epoch.data(ch_idx,1:plt_idx,:);
        lat_name = 'GIP';
    case 'gip'
        tri_epoch = cond_struct.gip_dev;
        cir_epoch = cond_struct.gip_std;
        plt_t = tri_epoch.times(end);
        ch_idx = find(ismember({tri_epoch.chanlocs.labels},tar_Ch));
        plt_idx = find(tri_epoch.times >= plt_t, 1)-1;
        plt_t = tri_epoch.times(1:plt_idx);
        tri_data = tri_epoch.data(ch_idx,1:plt_idx,:);
        cir_data = cir_epoch.data(ch_idx,1:plt_idx,:);
        lat_name = 'Stim';
%         tri_lat = -tri_lat;
%         cir_lat = -cir_lat;
end
% eye and head rotation
tri_eye = squeeze(tri_epoch.data(end-2,1:end-1,:));
tri_head = squeeze(tri_epoch.data(end-3,1:end-1,:));
cir_eye = squeeze(cir_epoch.data(end-2,1:end-1,:));
cir_head = squeeze(cir_epoch.data(end-3,1:end-1,:));
% normalize
sf = 50;
ntri_eye = sf*((tri_eye-min(tri_eye,[],'all'))./...
    (max(tri_eye,[],'all')-min(tri_eye,[],'all')));
ntri_head = sf*((tri_head-min(tri_head,[],'all'))./...
    (max(tri_head,[],'all')-min(tri_head,[],'all')));
ncir_eye = sf*((cir_eye-min(cir_eye,[],'all'))./...
    (max(cir_eye,[],'all')-min(cir_eye,[],'all')));
ncir_head = sf*((cir_head-min(cir_head,[],'all'))./...
    (max(cir_head,[],'all')-min(cir_head,[],'all')));
eye_rot = [tri_eye,cir_eye];
head_rot = [tri_head, cir_head];
neye_rot = sf*(eye_rot-min(eye_rot,[],'all'))./...
    (max(eye_rot,[],'all')-min(eye_rot,[],'all'));
nhead_rot = sf*(head_rot-min(head_rot,[],'all'))./...
    (max(head_rot,[],'all')-min(head_rot,[],'all'));
plt_eye = nanmean(neye_rot,2);
plt_eye = plt_eye - min(plt_eye);
plt_head = nanmean(nhead_rot,2);
plt_head = plt_head - min(plt_head);
%
figure
% plot(plt_t, nanmean(tri_data,3),'b-','DisplayName',sprintf('Triangle (%d)',sum(~all(isnan(tri_data),2))),'linewidth',3)
ht = shadedErrorBar(plt_t, squeeze(tri_data)', {@nanmean, @nanstd},'lineprops',...
    {'color','b','linewidth',3,'DisplayName',sprintf('Triangle (%d)',sum(~all(isnan(tri_data),2)))});
ht.patch.FaceAlpha = 0.3;
grid on
hold on
% plot(plt_t, nanmean(cir_data,3),'r-','DisplayName',sprintf('Circle (%d)',sum(~all(isnan(cir_data),2))),'linewidth',3)
hc = shadedErrorBar(plt_t, squeeze(cir_data)', {@nanmean, @nanstd},'lineprops',...
    {'color','r','linewidth',3,'DisplayName',sprintf('Circle (%d)',sum(~all(isnan(cir_data),2)))});
hc.patch.FaceAlpha = 0.3;
plot(plt_t,plt_eye,'m--','DisplayName','eye rot.','linewidth',3)
plot(plt_t,plt_head,'g--','DisplayName','head rot.','linewidth',3)
xline(0,'k-','DisplayName',ev_name,'linewidth',3)
xline(nanmedian(tri_lat),'b--','DisplayName',lat_name,'linewidth',3)
xline(nanmedian(cir_lat),'r--','DisplayName',lat_name,'linewidth',3)
legend(findobj(gca,'-regexp','DisplayName', '[^'']'),'location','northwest')
set(gca,'fontsize',30)
set(gca,'xtick',plt_t(1):100:plt_t(end))
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
title(sprintf('%s lock (%s)', ev_name, tar_Ch))


