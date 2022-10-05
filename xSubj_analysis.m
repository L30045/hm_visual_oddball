%% cross subjects erp
% load epoch information
% s05_epoch.mat seems to have strong artifacts in Hm_gip_Cz condition
% subj_list = [1, 4:6, 8:10];
filepath = '/home/yuan/Documents/2021 HM_visual_oddball/dataset/new epoch/';
subj_list = {dir([filepath, '*_epoch*']).name};
epoch_lib = cell(2,length(subj_list));
for i = 1:length(subj_list)
    load([filepath,subj_list{i}]);
    epoch_lib{1,i} = epoch_struct_noHm;
    epoch_lib{2,i} = epoch_struct_Hm;
end
% savepath = filepath;
% save([savepath,'epoch_lib_20221003.mat'],'-v7.3','epoch_lib');
disp('Done')

%%
filepath = '/home/yuan/Documents/2021 HM_visual_oddball/dataset/new epoch/';
savepath = filepath;
load([savepath,'epoch_lib_20221003.mat'],'epoch_lib');

%%
cond_name = 'noHm';
ev_name = 'gip';
tar_Ch = 'POz';
filepath = '/home/yuan/Documents/2021 HM_visual_oddball/dataset/new epoch/';
subj_list = {dir([filepath, '*_epoch*']).name};
% subj_list = [1,4:6,8:10];


len_stim = epoch_lib{1}.std_epoch.pnts;
len_gip = epoch_lib{1}.gip_std.pnts;

switch ev_name
    case 'stim'
        tri_lib = zeros(length(subj_list),len_stim);
        cir_lib = zeros(length(subj_list),len_stim);
        eye_lib = zeros(length(subj_list),len_stim);
        head_lib = zeros(length(subj_list),len_stim);
    case 'gip'
        tri_lib = zeros(length(subj_list),len_gip);
        cir_lib = zeros(length(subj_list),len_gip);
        eye_lib = zeros(length(subj_list),len_gip);
        head_lib = zeros(length(subj_list),len_gip);
    case 'fix'
        tri_lib = zeros(length(subj_list),len_gip);
        cir_lib = zeros(length(subj_list),len_gip);
        eye_lib = zeros(length(subj_list),len_gip);
        head_lib = zeros(length(subj_list),len_gip);
end

cir_stack = [];
tri_stack = [];
cir_time = [];
tri_time = [];

% for i = 1:length(subj_list)
for i = 1:8
    switch cond_name
        case 'noHm'
            cond_struct = epoch_lib{1,i};
        case 'Hm'
            cond_struct = epoch_lib{2,i};
    end
    switch ev_name
        case 'stim'
            tri_epoch = cond_struct.dev_epoch;
            cir_epoch = cond_struct.std_epoch;
        case 'gip'
            tri_epoch = cond_struct.gip_dev;
            cir_epoch = cond_struct.gip_std;
        case 'fix'
            tri_epoch = cond_struct.fix_dev;
            cir_epoch = cond_struct.fix_std;
            
    end
    ch_idx = find(ismember({tri_epoch.chanlocs.labels},tar_Ch));
    plt_t = tri_epoch.times;
    tri_lib(i,:) = mean(squeeze(tri_epoch.data(ch_idx,:,:)),2);
    cir_lib(i,:) = mean(squeeze(cir_epoch.data(ch_idx,:,:)),2);
    eye_lib(i,:) = mean([squeeze(cir_epoch.data(end-2,:,:)),squeeze(tri_epoch.data(end-2,:,:))],2,'omitnan');
    head_lib(i,:) = mean([squeeze(cir_epoch.data(end-3,:,:)),squeeze(tri_epoch.data(end-2,:,:))],2,'omitnan' );
    cir_stack = [cir_stack;squeeze(cir_epoch.data(ch_idx,:,:))'];
    tri_stack = [tri_stack;squeeze(tri_epoch.data(ch_idx,:,:))'];
    cir_time = [cir_time, cond_struct.event_time.diff_gip_std];
    tri_time = [tri_time, cond_struct.event_time.diff_gip_dev];
    
end

%
sf = 5;
neye_lib = sf*((eye_lib-min(eye_lib,[],'all'))./(max(eye_lib,[],'all')-min(eye_lib,[],'all')));
neye_lib = neye_lib - min(neye_lib,[],'all');
nhead_lib = sf*((head_lib-min(head_lib,[],'all'))./(max(head_lib,[],'all')-min(head_lib,[],'all')));
nhead_lib = nhead_lib - min(nhead_lib,[],'all');
%
shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(length(subj_list)))};
% shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan'))};
figure
% >> Triangle
ht = shadedErrorBar(plt_t, tri_lib, shaded_method,'lineprops',...
    {'color','b','linewidth',3,'DisplayName','Standard'});
ht.patch.FaceAlpha = 0.1;
grid on
hold on
% >> Circle
hc = shadedErrorBar(plt_t, cir_lib, shaded_method,'lineprops',...
    {'color','r','linewidth',3,'DisplayName','Deviant'});
hc.patch.FaceAlpha = 0.1;

%
% cond_struct = epoch_lib{2,2};
% cir_epoch = cond_struct.gip_std;
% tri_epoch = cond_struct.gip_dev;
% plt_cir = squeeze(cir_epoch.data(ch_idx,:,:));
% plt_tri = squeeze(tri_epoch.data(ch_idx,:,:));
% plt_t = tri_epoch.times;
% shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(length(subj_list)))};
% % shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan'))};
% figure
% % >> Triangle
% ht = shadedErrorBar(plt_t, plt_tri', shaded_method,'lineprops',...
%     {'color','b','linewidth',3,'DisplayName','Standard'});
% ht.patch.FaceAlpha = 0.1;
% grid on
% hold on
% % >> Circle
% hc = shadedErrorBar(plt_t, plt_cir', shaded_method,'lineprops',...
%     {'color','r','linewidth',3,'DisplayName','Deviant'});
% hc.patch.FaceAlpha = 0.1;

% >> Behavioral plots
% >> shadedErrorBar 
% he = shadedErrorBar(plt_t, neye_lib, {@nanmean, @nanstd},...
%     {'color','m','linewidth',1,'DisplayName','Eye rot.'});
% he.patch.FaceAlpha = 0.3;
% hh = shadedErrorBar(plt_t, nhead_lib, {@nanmean, @nanstd},...
%     {'color','g','linewidth',1,'DisplayName','Head rot.'});
% hh.patch.FaceAlpha = 0.3;
% >> Normalized plots
% plot(plt_t, mean(neye_lib,1,'omitnan'),'color','m','linewidth',1,'DisplayName','Eye rot.');
% plot(plt_t, mean(nhead_lib,1,'omitnan'),'color','g','linewidth',1,'DisplayName','Head rot.');
% >> Raw data
% plot(plt_t, mean(eye_lib,1,'omitnan'),'color','m','linewidth',1,'DisplayName','Eye rot.');
% plot(plt_t, mean(head_lib,1,'omitnan'),'color','g','linewidth',1,'DisplayName','Head rot.');

xline(0,'k-','DisplayName',ev_name,'linewidth',3)
legend(findobj(gca,'-regexp','DisplayName', '[^'']'),'location','northwest')
set(gca,'fontsize',30)
set(gcf,'color',[1 1 1])
set(gca,'xtick',round(plt_t(1):100:plt_t(end)))
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
title(sprintf('%s lock - %s (%s)', ev_name, cond_name, tar_Ch))

% sanity check on ERPImage
figure
[sort_time,sort_idx] = sort(cir_time);
plt_data = cir_stack(sort_idx,:);
% plt_data = cir_lib(1:8,:);
subplot(2,1,1)
h = pcolor([plt_data,nan(size(plt_data,1),1);nan(1,size(plt_data,2)+1)]);
set(h,'edgecolor','none')
hold on
[~,v_pos] = min(abs(plt_t));
vline(v_pos,'k--','linewidth',50);
colorbar
caxis([min([cir_lib(:);tri_lib(:)]),max([cir_lib(:);tri_lib(:)])])
set(gca,'xtick',1:50:size(plt_data,2))
set(gca,'xticklabel',round(plt_t(1:50:end)))
set(gca,'fontsize',20)
title(sprintf('Cir, %s lock - %s (%s)', ev_name, cond_name, tar_Ch))
subplot(2,1,2)
[~,sort_idx] = sort(tri_time);
plt_data = tri_stack(sort_idx,:);
% plt_data = tri_lib(1:8,:);
h = pcolor([plt_data,nan(size(plt_data,1),1);nan(1,size(plt_data,2)+1)]);
set(h,'edgecolor','none')
hold on 
[~,v_pos] = min(abs(plt_t));
vline(v_pos,'k--','linewidth',50);
colorbar
caxis([min([cir_lib(:);tri_lib(:)]),max([cir_lib(:);tri_lib(:)])])
xlabel('Time (ms)')
set(gca,'xtick',1:50:size(plt_data,2))
set(gca,'xticklabel',round(plt_t(1:50:end)))
title(sprintf('Tri, %s lock - %s (%s)', ev_name, cond_name, tar_Ch))
set(gca,'fontsize',20)
set(gcf,'color',[1 1 1])

%% Laterialized Readinese Potential (LRP)
cond_name = 'Hm';
C3_lib = zeros(length(subj_list),len_gip);
C4_lib = zeros(length(subj_list),len_gip);
        

for i = 1:length(subj_list)
    switch cond_name
        case 'noHm'
            cond_struct = epoch_lib{1,i};
        case 'Hm'
            cond_struct = epoch_lib{2,i};
    end
    cir_epoch = cond_struct.grab_epoch;
         
    C3_idx = find(ismember({tri_epoch.chanlocs.labels},'C3'));
    C4_idx = find(ismember({tri_epoch.chanlocs.labels},'C4'));
    
    plt_t = tri_epoch.times;
    C3_lib(i,:) = mean(squeeze(cir_epoch.data(C3_idx,:,:)),2);
    C4_lib(i,:) = mean(squeeze(cir_epoch.data(C4_idx,:,:)),2);
end

%
shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(length(subj_list)))};
% shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan'))};
figure
% >> Circle
hc = shadedErrorBar(plt_t, C3_lib, shaded_method,'lineprops',...
    {'color','r','linewidth',3,'DisplayName','C3'});
hc.patch.FaceAlpha = 0.1;
grid on
hold on
hc = shadedErrorBar(plt_t, C4_lib, shaded_method,'lineprops',...
    {'color','b','linewidth',3,'DisplayName','C4'});
hc.patch.FaceAlpha = 0.1;
xline(0,'k-','DisplayName','Response','linewidth',3)
legend(findobj(gca,'-regexp','DisplayName', '[^'']'),'location','northwest')
set(gca,'fontsize',30)
set(gcf,'color',[1 1 1])
set(gca,'xtick',plt_t(1):100:plt_t(end))
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
title('Response lock')



%% calculate head rotation velocity and eye rotation velocity
ang_lib = cell(4, length(subj_list), 2); % time lock * subj * condition
angDiff_lib = cell(4, length(subj_list), 2); % time lock * subj * condition
dist_lib = cell(4, length(subj_list), 2); % time lock * subj * condition

for subj_i = 1:length(subj_list)
    for cond_i = 1:2    
        epoch_struct = epoch_lib{cond_i,subj_i};
        upLoc = epoch_struct.event_time.upLoc;
        downLoc = epoch_struct.event_time.downLoc;
        leftLoc = epoch_struct.event_time.leftLoc;
        rightLoc = epoch_struct.event_time.rightLoc;
        tar_lib = [upLoc;downLoc;leftLoc;rightLoc];
%         ev_list = fieldnames(epoch_struct);
        ev_list = {'std_epoch','dev_epoch','gip_std','gip_dev'};
        ev_direct = {[epoch_struct.event_time.std_up;epoch_struct.event_time.std_down;...
                      epoch_struct.event_time.std_left; epoch_struct.event_time.std_right];...
                     [epoch_struct.event_time.dev_up;epoch_struct.event_time.dev_down;...
                      epoch_struct.event_time.dev_left; epoch_struct.event_time.dev_right]};
        dir_idx = [1,2,1,2];

        for e_i = 1:4
            tar_epoch = epoch_struct.(ev_list{e_i});
            nbchan = find(ismember({tar_epoch.chanlocs.labels},'HeadLoc_x'));
            tar_direct = ev_direct{dir_idx(e_i)};
            if e_i == 3
                s_t = isnan(epoch_struct.event_time.diff_gip_std);
            elseif e_i == 4
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
            
            % save data to lib
            dist_lib{e_i,subj_i,cond_i} = dist_gip2box;
            ang_lib{e_i,subj_i,cond_i} = {headAng,eyeAng,gipAng};
            angDiff_lib{e_i,subj_i,cond_i} = {headAngDiff,eyeAngDiff,gipAngDiff};
        end
    end
end
disp('Done')
com = 'Dimension: [std_epoch,dev_epoch,gip_std,gip_dev] * subject * condition. Ang_lib cell: [head, eye, gip]';
% save('behav_lib_v2.mat','ang_lib','angDiff_lib','dist_lib','com');


%% plot behavior for each condition
cond_i = 1;
ev_name = 'stim';

switch ev_name
    case 'stim'
        ev_idx = [1,2];
        plt_t = epoch_struct.dev_epoch.times;
    case 'gip'
        ev_idx = [3,4];
        plt_t = epoch_struct.gip_dev.times;
end

switch cond_i
    case 1
        cond_name = 'noHm';
    case 2
        cond_name = 'Hm';
end

for e_i = 2
    switch e_i
        case 1
            trial_name = 'DEV';
        case 2
            trial_name = 'STD';
    end
    plt_dist_ori = cell2mat(cellfun(@(x) double(mean(x,2,'omitnan')),...
        dist_lib(ev_idx(e_i),:,cond_i),'uniformoutput',0));
    tmp_ang = ang_lib(ev_idx(e_i),:,cond_i);
    tmp_angDiff = angDiff_lib(ev_idx(e_i),:,cond_i);
    % angle
    plt_h_a = cell2mat(cellfun(@(x) mean(x{1},2,'omitnan'),...
        tmp_ang, 'uniformoutput',0));
    plt_e_a = cell2mat(cellfun(@(x) mean(x{2},2,'omitnan'),...
        tmp_ang, 'uniformoutput',0));
    plt_g_a = cell2mat(cellfun(@(x) mean(x{3},2,'omitnan'),...
        tmp_ang, 'uniformoutput',0));
    % angle diff
    plt_h_ad = cell2mat(cellfun(@(x) mean(x{1},2,'omitnan'),...
        tmp_angDiff, 'uniformoutput',0));
    plt_e_ad = cell2mat(cellfun(@(x) mean(x{2},2,'omitnan'),...
        tmp_angDiff, 'uniformoutput',0));
    plt_g_ad = cell2mat(cellfun(@(x) mean(x{3},2,'omitnan'),...
        tmp_angDiff, 'uniformoutput',0));
    % angle speed
    plt_h_as = plt_h_ad*500;
    plt_e_as = plt_e_ad*500;
    plt_g_as = plt_g_ad*500;

    shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(length(subj_list)))};
    my_norm = @(x)((x-min(x,[],1))./(max(x,[],1)-min(x,[],1)));
    
    % angle difference
%     figure;
%     % normalized distance
%     plt_dist = my_norm(plt_dist_ori);
%     scale = max(mean(plt_e_ad,2,'omitnan'));
%     plt_dist = plt_dist * scale;
%     shadedErrorBar(plt_t,plt_e_ad', shaded_method, 'lineprops',{'b-','DisplayName','EyeAngDiff','linewidth',3});
%     grid on; hold on;
%     shadedErrorBar(plt_t,plt_h_ad', shaded_method, 'lineprops',{'r-','DisplayName','HeadAngDiff','linewidth',3})
% %     shadedErrorBar(plt_t,plt_g_ad', shaded_method, 'lineprops',{'k-','DisplayName','GIPAngDiff','linewidth',3})
%     shadedErrorBar(plt_t,plt_dist', shaded_method, 'lineprops',{'g-','DisplayName','dist2box','linewidth',3})
%     xline(0,'k--','linewidth',3,'DisplayName',ev_name);
%     title(sprintf('%s lock - %s (%s)',ev_name, trial_name, cond_name));
%     set(gca,'fontsize',20)
%     set(gcf,'color','w')
%     xlabel('Time (ms)')
%     ylabel('Angle (deg)')
%     legend(findobj(gca,'-regexp','DisplayName', '[^'']'));
    
    % angle 
    figure;
    plt_dist = my_norm(plt_dist_ori);
    scale = max(mean(plt_e_a,2,'omitnan'));
    plt_dist = plt_dist * scale;
    shadedErrorBar(plt_t,plt_e_a', shaded_method, 'lineprops',{'b-','DisplayName','EyeAng','linewidth',3});
    grid on; hold on;
    shadedErrorBar(plt_t,plt_h_a', shaded_method, 'lineprops',{'r-','DisplayName','HeadAng','linewidth',3})
%     shadedErrorBar(plt_t,plt_g_a', shaded_method, 'lineprops',{'k-','DisplayName','GIPAng','linewidth',3})
    shadedErrorBar(plt_t,plt_dist', shaded_method, 'lineprops',{'g-','DisplayName','dist2box','linewidth',3})
    xline(0,'k--','linewidth',3,'DisplayName',ev_name);
    title(sprintf('%s lock - %s (%s)',ev_name, trial_name, cond_name));
    set(gca,'fontsize',20)
    set(gcf,'color','w')
    xlabel('Time (ms)')
    ylabel('Angle (deg)')
    legend(findobj(gca,'-regexp','DisplayName', '[^'']'));
    
    % angle speed
    figure;
    plt_dist = my_norm(plt_dist_ori);
    scale = max(mean(plt_e_as,2,'omitnan'));
    plt_dist = plt_dist * scale;
    shadedErrorBar(plt_t,plt_e_as', shaded_method, 'lineprops',{'b-','DisplayName','EyeAng','linewidth',3});
    grid on; hold on;
    shadedErrorBar(plt_t,plt_h_as', shaded_method, 'lineprops',{'r-','DisplayName','HeadAng','linewidth',3})
%     shadedErrorBar(plt_t,plt_g_a', shaded_method, 'lineprops',{'k-','DisplayName','GIPAng','linewidth',3})
    shadedErrorBar(plt_t,plt_dist', shaded_method, 'lineprops',{'g-','DisplayName','dist2box','linewidth',3})
    xline(0,'k--','linewidth',3,'DisplayName',ev_name);
    title(sprintf('%s lock - %s (%s)',ev_name, trial_name, cond_name));
    set(gca,'fontsize',20)
    set(gcf,'color','w')
    xlabel('Time (ms)')
    ylabel('Angular Speed (deg/sec)')
    legend(findobj(gca,'-regexp','DisplayName', '[^'']'));

end

%% compare std dev behavior
cond_i = 2;
ev_name = 'stim';

switch ev_name
    case 'stim'
        ev_idx = [1,2];
        plt_t = epoch_struct.std_epoch.times;
    case 'gip'
        ev_idx = [3,4];
        plt_t = epoch_struct.gip_std.times;
end

switch cond_i
    case 1
        cond_name = 'noHm';
    case 2
        cond_name = 'Hm';
end

plt_dist_ori = cell(1,2);
plt_h_ad = cell(1,2);
plt_e_ad = cell(1,2);

for e_i = 1:2
    switch e_i
        case 1
            trial_name = 'DEV';
        case 2
            trial_name = 'STD';
    end
    plt_dist_ori{e_i} = cell2mat(cellfun(@(x) double(mean(x,2,'omitnan')),...
        dist_lib(ev_idx(e_i),:,cond_i),'uniformoutput',0));
    tmp_angDiff = angDiff_lib(ev_idx(e_i),:,cond_i);
    % angle diff
    plt_h_ad{e_i} = cell2mat(cellfun(@(x) mean(x{1},2,'omitnan'),...
        tmp_angDiff, 'uniformoutput',0));
    plt_e_ad{e_i} = cell2mat(cellfun(@(x) mean(x{2},2,'omitnan'),...
        tmp_angDiff, 'uniformoutput',0));
end

shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(length(subj_list)))};
my_norm = @(x)((x-min(x,[],1))./(max(x,[],1)-min(x,[],1)));
%
figure;
% normalized distance
h_std = plt_h_ad{1};
h_dev = plt_h_ad{2};
e_std = plt_e_ad{1};
e_dev = plt_e_ad{2};
shadedErrorBar(plt_t,h_std', shaded_method, 'lineprops',{'b--','DisplayName','DEV (head)','linewidth',3});
grid on; hold on;
shadedErrorBar(plt_t,h_dev', shaded_method, 'lineprops',{'r--','DisplayName','STD (head)','linewidth',3});
shadedErrorBar(plt_t,e_std', shaded_method, 'lineprops',{'b-','DisplayName','DEV (eye)','linewidth',3});
shadedErrorBar(plt_t,e_dev', shaded_method, 'lineprops',{'r-','DisplayName','STD (eye)','linewidth',3});
xline(0,'k--','linewidth',3,'DisplayName',ev_name);
title(sprintf('%s lock - %s (%s)',ev_name, trial_name, cond_name));
set(gca,'fontsize',20)
xlabel('Time (ms)')
ylabel('Angle (deg)')
legend(findobj(gca,'-regexp','DisplayName', '[^'']'));



%% plot response time
% Need to fix epoch_ez for extracting diff_stim_grab

cond_i = 2;
switch cond_i
    case 1
        cond_name = 'w/o Head Moving';
    case 2
        cond_name = 'w/ Head Moving';
end
time_thres = 2000; %msec
diff_gip_std = cell2mat(cellfun(@(x) x.event_time.diff_gip_std, epoch_lib(cond_i,:),'uniformoutput',0));
diff_gip_dev = cell2mat(cellfun(@(x) x.event_time.diff_gip_dev, epoch_lib(cond_i,:),'uniformoutput',0));
diff_stim_grab = cell2mat(cellfun(@(x) x.event_time.diff_stim_grab, epoch_lib(cond_i,:),'uniformoutput',0)); % grab time - stim time (sec)
diff_stim_grab = diff_stim_grab*1000;
diff_grab_gip = diff_stim_grab-diff_gip_std;

diff_gip_std(diff_gip_std>time_thres) = [];
diff_gip_std(isnan(diff_gip_std)) = [];
diff_gip_dev(diff_gip_dev>time_thres)= [];
diff_gip_dev(isnan(diff_gip_dev)) = [];
diff_stim_grab(diff_stim_grab>time_thres) = [];
diff_stim_grab(isnan(diff_stim_grab)) = [];
diff_grab_gip(diff_grab_gip>time_thres) = [];
diff_grab_gip(isnan(diff_grab_gip)) = [];


figure
[~,~,p] = statcond({diff_gip_std,diff_gip_dev},'method','bootstrap','naccu',1000);
hs = histogram(diff_gip_std,'binwidth',50,'Normalization','probability','DisplayName','Deviant',...
    'DisplayStyle','Stairs','linewidth',3,'linestyle','-','edgecolor','r');
hold on; grid on
hd = histogram(diff_gip_dev,'binwidth',50,'Normalization','probability','DisplayName','Standard',...
    'DisplayStyle','Stairs','linewidth',3,'linestyle','-','edgecolor','b');
plot(nan,'DisplayName',sprintf('p = %g',p));
legend(findobj(gca,'-regexp','DisplayName', '[^'']'));
xlabel('latency (ms)')
ylabel('Probability')
set(gca,'fontsize',20)
title(['GIP latency ',cond_name])
set(gca,'fontsize',20)
set(gcf,'color','w')

figure
hs = histogram(diff_grab_gip,'binwidth',50,'Normalization','probability','DisplayName','Response');

figure
hs = histogram(diff_gip_std,'binwidth',50,'Normalization','probability','DisplayName','GIP',...
    'DisplayStyle','Stairs','linewidth',3,'linestyle','-','edgecolor','r');
hold on; grid on
hd = histogram(diff_stim_grab,'binwidth',50,'Normalization','probability','DisplayName','Response',...
    'DisplayStyle','Stairs','linewidth',3,'linestyle','-','edgecolor','b');
legend(findobj(gca,'-regexp','DisplayName', '[^'']'));
xlabel('Latency(ms)')
ylabel('Probability')
set(gca,'fontsize',20)
title(['GIP and Response ',cond_name])
set(gca,'fontsize',20)
set(gcf,'color','w')

%% ERPImage
tar_Ch = {'Cz','POz'};
cond_i = 2;
ev_name = 'gip';
all_event = unique({epoch_lib{1}.std_epoch.event.type});


switch ev_name
    case 'stim'
        std_name = 'std_epoch';
        dev_name = 'dev_epoch';
        sort_std = 'circle_gip_start';
        sort_dev = 'triangle_gip_start';
    case 'gip'
        std_name = 'gip_std';
        dev_name = 'gip_dev';
end

std_merge = eeg_emptyset();
dev_merge = eeg_emptyset();

subj_list = [1,4:6,8:10];
tar_list = [1,4,6,8,9,10];
for subj_i = 1:length(tar_list)
    tmp_std = epoch_lib{cond_i,subj_list==tar_list(subj_i)}.(std_name);
    tmp_dev = epoch_lib{cond_i,subj_list==tar_list(subj_i)}.(dev_name);
    ch_idx = ismember({tmp_std.chanlocs.labels},tar_Ch);
    tmp_std = pop_select(tmp_std,'channel',find(ch_idx));
    tmp_dev = pop_select(tmp_dev,'channel',find(ch_idx));
    if isempty(std_merge.data)
        std_merge = tmp_std;
        dev_merge = tmp_dev;
    else
        std_merge = pop_mergeset(std_merge, tmp_std);
        dev_merge = pop_mergeset(dev_merge, tmp_dev);
    end
end


[std_clean,rm_idx_std] = my_rmEpoch(std_merge);
[dev_clean,rm_idx_dev] = my_rmEpoch(dev_merge);
fprintf('STD remove rate = %2.1f%% (%d/ %d)\n',sum(rm_idx_std)/std_merge.trials*100,sum(rm_idx_std),std_merge.trials);
fprintf('DEV remove rate = %2.1f%% (%d/ %d)\n',sum(rm_idx_dev)/dev_merge.trials*100,sum(rm_idx_dev),dev_merge.trials);

%%
smooth = 10;
figure;
pop_erpimage(std_clean,1, [1],[[]],'Cz',smooth,1,{ 'circle_gip_start'},[],'latency' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [1] std_clean.chanlocs std_clean.chaninfo } );
figure
pop_erpimage(dev_clean,1, [1],[[]],'Cz',smooth,1,{ 'triangle_gip_start'},[],'latency' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [1] dev_clean.chanlocs dev_clean.chaninfo } );

%% Sanity check. Number of trials in epoch_lib Ring 0 and 1 both appear in the file
nb_std_epoch_cond1 = cellfun(@(x) size(x.std_epoch.data,3), epoch_lib(1,:));
nb_std_epoch_cond2 = cellfun(@(x) size(x.std_epoch.data,3), epoch_lib(2,:));
nb_dev_epoch_cond1 = cellfun(@(x) size(x.dev_epoch.data,3), epoch_lib(1,:));
nb_dev_epoch_cond2 = cellfun(@(x) size(x.dev_epoch.data,3), epoch_lib(2,:));
nb_gip_std_cond1 = cellfun(@(x) size(x.gip_std.data,3), epoch_lib(1,:));
nb_gip_std_cond2 = cellfun(@(x) size(x.gip_std.data,3), epoch_lib(2,:));
nb_gip_dev_cond1 = cellfun(@(x) size(x.gip_dev.data,3), epoch_lib(1,:));
nb_gip_dev_cond2 = cellfun(@(x) size(x.gip_dev.data,3), epoch_lib(2,:));

%% check event markers

%% TO DO LIST
% hm_visual_oddball
% plot the shadedErrorBar for the ERP
% sort the GIP event by stim onset
% check the event marker Ring 1 and Ring 0 appear together problem
% Regression artifact cleaning
% disable baseline removal
% 
% SSVEP
% Calculate PSD after 500ms to avoid 10Hz created by stimulus-onset ERP 
% Epoch by GIP onset
% disable baseline removal





