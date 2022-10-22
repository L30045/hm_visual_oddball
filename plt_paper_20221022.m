%% plot for paper
filepath = 'D:\Research\';
load([filepath,'behav_lib_20221022.mat']);
load([filepath,'epoch_lib_rmPreStim.mat']);
% filepath = '//hoarding/yuan/Documents/2021 HM_visual_oddball/dataset/new epoch/';
% subj_list = {dir([filepath, 'rmPreStim*']).name};
subj_list = 1:14;


%% plot behavior
cond_i = 2;
ev_name = 'gip';

switch ev_name
    case 'stim'
        ev_idx = [1,2];
        plt_t = epoch_lib{1}.dev_epoch.times;
    case 'gip'
        ev_idx = [3,4];
        plt_t = epoch_lib{1}.gip_dev.times;
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
            trial_name = 'CIR';
        case 2
            trial_name = 'TRI';
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

%     shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(length(subj_list)))};
    shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan'))};
    my_norm = @(x)((x-min(x,[],1))./(max(x,[],1)-min(x,[],1)));
    
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

%% Merge epoch_lib to perform ERPImage using EEGLAB function
tmp_ch = cellfun(@(x) {x.std_epoch.chanlocs.labels}, epoch_lib(:),'uniformoutput',0);
[~, tmp_idx] = max(cellfun(@length, tmp_ch));
common_ch = tmp_ch{tmp_idx};
for ch_i = 1:length(epoch_lib(:))
    common_ch = intersect(common_ch, tmp_ch{ch_i});
end

% tarCh = {'POz'};
tarCh = common_ch;
merge_stim_cir = {eeg_emptyset(),eeg_emptyset()};
merge_stim_tri = {eeg_emptyset(),eeg_emptyset()};
merge_gip_cir = {eeg_emptyset(),eeg_emptyset()};
merge_gip_tri = {eeg_emptyset(),eeg_emptyset()};

for subj_i = 1:size(epoch_lib,2)
    for cond_i = 1:2
        tmp_epoch = my_rmEpoch(epoch_lib{cond_i,subj_i});
        tmp_eeg_1 = pop_select(tmp_epoch.std_epoch,'channel',tarCh);
        tmp_eeg_2 = pop_select(tmp_epoch.dev_epoch,'channel',tarCh);
        tmp_eeg_3 = pop_select(tmp_epoch.gip_std,'channel',tarCh);
        tmp_eeg_4 = pop_select(tmp_epoch.gip_dev,'channel',tarCh);
        if merge_stim_cir{cond_i}.nbchan==0
            merge_stim_cir{cond_i} = tmp_eeg_1;
            merge_stim_tri{cond_i} = tmp_eeg_2;
            merge_gip_cir{cond_i} = tmp_eeg_3;
            merge_gip_tri{cond_i} = tmp_eeg_4;
        else
            merge_stim_cir{cond_i} = pop_mergeset(tmp_eeg_1,merge_stim_cir{cond_i});
            merge_stim_tri{cond_i} = pop_mergeset(tmp_eeg_2,merge_stim_tri{cond_i});
            merge_gip_cir{cond_i} = pop_mergeset(tmp_eeg_3,merge_gip_cir{cond_i});
            merge_gip_tri{cond_i} = pop_mergeset(tmp_eeg_4,merge_gip_tri{cond_i});
        end
    end
end
    
%% plot ERPImage
tarCh = 'Cz';
lock_name = 'gip';
ev_name = 'cir';
cond_name = 'Hm';
thres_amp = 100;
switch cond_name
    case 'noHm'
        cond_i = 1;
    case 'Hm'
        cond_i = 2;
end

eval(sprintf('plt_EEG=merge_%s_%s{cond_i};',lock_name,ev_name));
plt_EEG = pop_select(plt_EEG,'channel',{tarCh});
% [plt_EEG,rm_idx] = pop_autorej(plt_EEG,'threshold',50,'nogui','on');
[plt_EEG, rm_idx] = pop_eegthresh(plt_EEG,1,1,-thres_amp,thres_amp,plt_EEG.xmin,plt_EEG.xmax,0,0);
plt_EEG = pop_rejepoch(plt_EEG, rm_idx, 0);
idx_cir_1 = cellfun(@(x) ~isempty(regexp(x,'Standard','once')),{plt_EEG.event.type});
idx_cir_2 = cellfun(@(x) ~isempty(regexp(x,'\(L\)','once')),{plt_EEG.event.type});
idx_cir = idx_cir_1 | idx_cir_2;
idx_tri_1 = cellfun(@(x) ~isempty(regexp(x,'Deviant','once')),{plt_EEG.event.type});
idx_tri_2 = cellfun(@(x) ~isempty(regexp(x,'\(R\)','once')),{plt_EEG.event.type});
idx_tri = idx_tri_1 | idx_tri_2;
switch lock_name
    case 'stim'
        if strcmp(ev_name,'cir')
            sort_name = {'circle_gip_start'};
        else
            sort_name = {'triangle_gip_start'};
        end
    case 'gip'
        if strcmp(ev_name,'cir')
            sort_name = {plt_EEG.event(idx_cir).type};
        else
            sort_name = {plt_EEG.event(idx_tri).type};
        end
end
smoothing = 10;
figure; pop_erpimage(plt_EEG,1, [1],[[]],tarCh,smoothing,1,sort_name,[],'latency','yerplabel','\muV','erp','on','cbar','on','topo', { [1] plt_EEG.chanlocs plt_EEG.chaninfo } );

%% plot ERP
tarCh = 'Cz';
ev_name = 'gip';
cond_name = 'Hm';
switch cond_name
    case 'noHm'
        cond_i = 1;
    case 'Hm'
        cond_i = 2;
end
eval(sprintf('plt_EEG1 = merge_%s_cir{cond_i};',ev_name));
eval(sprintf('plt_EEG2 = merge_%s_tri{cond_i};',ev_name));
plt_EEG1 = pop_select(plt_EEG1,'channel',{tarCh});
plt_EEG2 = pop_select(plt_EEG2,'channel',{tarCh});
% plt_EEG1 = pop_eegfiltnew(plt_EEG1,5,20);
% plt_EEG2 = pop_eegfiltnew(plt_EEG2,5,20);
plt_EEG1 = pop_autorej(plt_EEG1,'threshold',100,'nogui','on');
plt_EEG2 = pop_autorej(plt_EEG2,'threshold',100,'nogui','on');
% shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(length(subj_list)))};
% shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan'))};
shaded_method = {@(x)(mean(x,'omitnan')),@(x)([quantile(x,0.8)-mean(x,'omitnan');mean(x,'omitnan')-quantile(x,0.2)])};
figure
plt_t = plt_EEG1.times;
% >> Triangle
ht = shadedErrorBar(plt_t, squeeze(plt_EEG2.data)', shaded_method,'lineprops',...
    {'color','b','linewidth',3,'DisplayName','Standard'});
ht.patch.FaceAlpha = 0.1;
grid on
hold on
% >> Circle
hc = shadedErrorBar(plt_t, squeeze(plt_EEG1.data)', shaded_method,'lineprops',...
    {'color','r','linewidth',3,'DisplayName','Deviant'});
hc.patch.FaceAlpha = 0.1;
xline(0,'k-','DisplayName',ev_name,'linewidth',3)
legend(findobj(gca,'-regexp','DisplayName', '[^'']'),'location','northwest')
set(gca,'fontsize',30)
set(gcf,'color',[1 1 1])
set(gca,'xtick',round(plt_t(1):100:plt_t(end)))
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
title(sprintf('%s lock - %s (%s)', ev_name, cond_name, tarCh))