%% plot for paper
filepath = 'D:\Research\';
% load([filepath,'behav_lib_20221022.mat']);
load([filepath,'epoch_lib_rmPreStim_new.mat']);
% filepath = '//hoarding/yuan/Documents/2021 HM_visual_oddball/dataset/new epoch/';
% subj_list = {dir([filepath, 'rmPreStim*']).name};
subj_list = 1:14;

% print_multi([saveFigPath,'v_fig1_',tName],{'pdf'})

%% plot the mean of mean ERP across subject
plt_epoch_lib = epoch_lib;
thres_time = [100 1000];
rm_thres = 30;
nb_subj = size(plt_epoch_lib,2);
% shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(nb_subj))};
shaded_method = {@(x)(median(x,'omitnan')),@(x)([quantile(x,0.8)-median(x,'omitnan');median(x,'omitnan')-quantile(x,0.2)])};
savepath = 'D:\Research\oddball_fig\xSubj\mean ERP\';

subplot_lib = cell(1,6);
count = 1;
for cond_name = {'noHm','Hm'}
    for ev_name = {'stim','gip','fix'}
        for tarCh = {'CPz'}
%             loc_path = sprintf('%s%s/',savepath,tarCh{:});
%             loc_path = 'C:\Users\Yuan\OneDrive\Desktop\graduate!\fig\vr\';
%             if ~exist(loc_path,'dir')
%                 mkdir(loc_path)
%             end
            [fig,plt_t,cir_lib,tri_lib] = plt_erp_meanXsubj(plt_epoch_lib, cond_name{:}, ev_name{:}, tarCh{:}, thres_time, rm_thres,shaded_method, false);
            % for printing PDF
%             set(gca,'units','centimeters')
%             pos = get(gca,'Position');
%             ti = get(gca,'TightInset');
%             set(gcf, 'PaperUnits','centimeters');
%             set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%             set(gcf, 'PaperPositionMode', 'manual');
%             set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% %             saveas(fig, sprintf('%serp_%s_%s_%s.png',loc_path,cond_name{:},ev_name{:},tarCh{:}));
% %             saveas(fig, sprintf('%serp_%s_%s_%s.pdf',loc_path,cond_name{:},ev_name{:},tarCh{:}));
%             close(fig)
            subplot_lib{count} = {plt_t, cir_lib, tri_lib};
            count = count+1;
        end
    end
end

%% plot above figures in one big subplot
plt_lib = reshape(reshape(subplot_lib,3,2)',[],1); % for 3 by 2 subplot
shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(nb_subj))};
% shaded_method = {@(x)(median(x,'omitnan')),@(x)([quantile(x,0.8)-median(x,'omitnan');median(x,'omitnan')-quantile(x,0.2)])};
% plt_lib = subplot_lib; % for 2 by 3 subplot
ev_name = {'stim','stim','gip','gip','fix','fix'};
fig = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
plt_time_idx = [-500 1000];
ax_lib = cell(1,6);
for p_i = 1:6
    ax_lib{p_i} = subplot(3,2,p_i);
    plt_t = plt_lib{p_i}{1};
    plt_idx = plt_t >= plt_time_idx(1) & plt_t<=plt_time_idx(2);
    plt_t = plt_t(plt_idx);
    cir_lib = plt_lib{p_i}{2}(:,plt_idx);
    tri_lib = plt_lib{p_i}{3}(:,plt_idx);
    
    [~, p] = ttest(tri_lib,cir_lib);
%     [corrected_p, h] = bonf_holm(p);
    h = p <= 0.05;
    
    ht = shadedErrorBar(plt_t, tri_lib, shaded_method,'lineprops',...
        {'color','b','linewidth',3,'DisplayName','Standard'});
    ht.patch.FaceAlpha = 0.1;
    grid on
    hold on
    % >> Circle
    hc = shadedErrorBar(plt_t, cir_lib, shaded_method,'lineprops',...
        {'color','r','linewidth',3,'DisplayName','Deviant'});
    hc.patch.FaceAlpha = 0.1;
%     xline(0,'k-','DisplayName',sprintf('%s onset',ev_name{p_i}),'linewidth',3)
    xline(0,'k-','DisplayName','Event onset','linewidth',3)
%     plot(plt_t(h2), mean(cir_lib(:,h2)), 'gx','linewidth',2,'markersize',15);
%     plot(plt_t(h), mean(cir_lib(:,h)), 'kx','linewidth',2,'markersize',15);
    set(gca,'fontsize',15)
    set(gcf,'color',[1 1 1])
    set(gca,'xtick',round(plt_t(1):100:plt_t(end)))
    if p_i == 3
        ylabel('Amplitude (\muV)')
    end
    if p_i == 5||p_i==6
        xlabel('Time (ms)')
    end
    if p_i == 2
        legend(findobj(gca,'-regexp','DisplayName', '[^'']'),'location','northeast')
    end
    
end
linkaxes([ax_lib{:}],'y')

fig.Units = 'centimeters';
fig.PaperUnits = 'centimeters';
fig.PaperSize = fig.Position(3:4);
savespath = 'D:\Research\oddball_fig\xSubj\mean ERP\';
% saveas(fig, sprintf('%smeanERP_O2.png',savepath));
% saveas(fig, sprintf('%smeanERP_O2.pdf',savepath));

%% anova p300
ev_name = {'stim','stim','gip','gip','fix','fix'};
cond_name = {'noHm','Hm'};
sig_name = {'noSig','Sig'};
for p_i = 6
    ax_lib{p_i} = subplot(3,2,p_i);
    plt_t = plt_lib{p_i}{1};
    plt_idx = plt_t >= plt_time_idx(1) & plt_t<=plt_time_idx(2);
    plt_t = plt_t(plt_idx);
    cir_lib = plt_lib{p_i}{2}(:,plt_idx);
    tri_lib = plt_lib{p_i}{3}(:,plt_idx);
    
    % extract amplitude around 300 ms
    p3_idx = plt_t>=250 & plt_t<=350;
    tri_p3 = mean(tri_lib(:,p3_idx),2);
    cir_p3 = mean(cir_lib(:,p3_idx),2);
    p = anova1([cir_p3,tri_p3]);
    
    title(sprintf('%s - %s (%s)',ev_name{p_i},cond_name{mod(p_i-1,2)+1}),sig_name{(p<=0.05)+1})
    
end


%% calculate behavior
behav_lib_ori = cell(size(epoch_lib(:,1:8)));
for subj_i = 1:size(behav_lib_ori,2)
    for cond_i = 1:2
        behav_lib_ori{cond_i,subj_i} = cal_behav(epoch_lib{cond_i,subj_i});
    end
end
disp('Done')

%% plot behavior
behav_lib = behav_lib_ori;
plt_epoch_lib = epoch_lib(:,1:8);
% avoid outlier near the edges
rm_edge = round(100*0.001*plt_epoch_lib{1}.std_epoch.srate); % sample points
% [fix_subj_idx, grab_subj_idx] = find_if_device(epoch_lib);
% behav_lib = behav_lib_ori(:,sum_fix_subj_idx==2);
% savepath = 'D:\Research\oddball_fig\xSubj\behav\';
savepath = 'C:\Users\Yuan\OneDrive\Desktop\graduate!\fig\vr\';

ang_lib = cell(2,6);
v_ang_lib = cell(2,6);

for ev_name_loop = {'std','dev'}
    count = 1;
    for cond_name_loop = {'noHm','Hm'}
        for lock_name_loop = {'stim','gip','fix'}
            cond_name = cond_name_loop{:};
            lock_name = lock_name_loop{:};
            ev_name = ev_name_loop{:};
            switch cond_name
                case 'noHm'
                    cond_i = 1;
                case 'Hm'
                    cond_i = 2;
            end
            switch ev_name
                case 'std'
                    trial_name = 'CIR';
                    l_i = 1;
                case 'dev'
                    trial_name = 'TRI';
                    l_i = 2;
            end
            [fix_subj_idx, grab_subj_idx] = find_if_device(plt_epoch_lib);
            fix_subj_idx = sum(fix_subj_idx)==2;
            grab_subj_idx = sum(grab_subj_idx)==2;
            switch lock_name
                case 'stim'
                    epoch_name = sprintf('%s_epoch',ev_name);
                    plt_t = plt_epoch_lib{1}.dev_epoch.times;
                    include_subj = true(1,size(behav_lib,2));
                case 'gip'
                    epoch_name = sprintf('gip_%s',ev_name);
                    plt_t = plt_epoch_lib{1}.gip_dev.times;
                    include_subj = true(1,size(behav_lib,2));
                case 'fix'
                    epoch_name = sprintf('fix_%s',ev_name);
                    plt_t = plt_epoch_lib{1}.gip_dev.times;
                    include_subj = fix_subj_idx;
            end
            % head 
            plt_h_a = cell2mat(cellfun(@(x) mean(x.(sprintf('%s',epoch_name)).headAng,2,'omitnan'),...
                              behav_lib(cond_i,include_subj),'uniformoutput',0));
            plt_h_ad = cell2mat(cellfun(@(x) mean(x.(sprintf('%s',epoch_name)).headRot,2,'omitnan'),...
                              behav_lib(cond_i,include_subj),'uniformoutput',0));
            % eye
            plt_e_a = cell2mat(cellfun(@(x) mean(x.(sprintf('%s',epoch_name)).eyeAng,2,'omitnan'),...
                              behav_lib(cond_i,include_subj),'uniformoutput',0));
            plt_e_ad = cell2mat(cellfun(@(x) mean(x.(sprintf('%s',epoch_name)).eyeRot,2,'omitnan'),...
                              behav_lib(cond_i,include_subj),'uniformoutput',0));
            % gip
            plt_g_a = cell2mat(cellfun(@(x) mean(x.(sprintf('%s',epoch_name)).gipAng,2,'omitnan'),...
                              behav_lib(cond_i,include_subj),'uniformoutput',0));
            plt_g_ad = cell2mat(cellfun(@(x) mean(x.(sprintf('%s',epoch_name)).gipRot,2,'omitnan'),...
                              behav_lib(cond_i,include_subj),'uniformoutput',0));
            % dist
            plt_dist_ori = cell2mat(cellfun(@(x) mean(x.(sprintf('%s',epoch_name)).dist,2,'omitnan'),...
                              behav_lib(cond_i,include_subj),'uniformoutput',0));
            % remove edge
            plt_t = plt_t(rm_edge:end-rm_edge);
            plt_h_a = plt_h_a(rm_edge:end-rm_edge,:);
            plt_h_ad = plt_h_ad(rm_edge:end-rm_edge,:);
            plt_e_a = plt_e_a(rm_edge:end-rm_edge,:);
            plt_e_ad = plt_e_ad(rm_edge:end-rm_edge,:);
            plt_g_a = plt_g_a(rm_edge:end-rm_edge,:);
            plt_g_ad = plt_g_ad(rm_edge:end-rm_edge,:);
            plt_dist_ori = plt_dist_ori(rm_edge:end-rm_edge,:);
            shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(size(plt_epoch_lib,2)))};
%             shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan'))};
            my_norm = @(x)((x-min(x,[],1))./(max(x,[],1)-min(x,[],1)));

            % angle 
%             fig = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
            plt_dist = my_norm(plt_dist_ori);
            scale = max(mean(plt_e_a,2,'omitnan'));
            plt_dist = plt_dist * scale;
%             shadedErrorBar(plt_t,plt_e_a', shaded_method, 'lineprops',{'b-','DisplayName','EyeAng','linewidth',3});
%             grid on; hold on;
%             shadedErrorBar(plt_t,plt_h_a', shaded_method, 'lineprops',{'r-','DisplayName','HeadAng','linewidth',3})
%             shadedErrorBar(plt_t,plt_g_a', shaded_method, 'lineprops',{'k-','DisplayName','GIPAng','linewidth',3})
%             shadedErrorBar(plt_t,plt_dist', shaded_method, 'lineprops',{'g-','DisplayName','dist2box','linewidth',3})
%             xline(0,'k--','linewidth',3,'DisplayName',sprintf('%s onset',lock_name));
% %             title(sprintf('%s lock - %s (%s)',lock_name, trial_name, cond_name));
%             set(gca,'fontsize',20)
%             set(gcf,'color','w')
%             xlabel('Time (ms)')
%             ylabel('Angle (deg)')
%             legend(findobj(gca,'-regexp','DisplayName', '[^'']'));
%             % for printing PDF
%             set(gca,'units','centimeters')
%             pos = get(gca,'Position');
%             ti = get(gca,'TightInset');
%             set(gcf, 'PaperUnits','centimeters');
%             set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%             set(gcf, 'PaperPositionMode', 'manual');
%             set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% %             saveas(fig,sprintf('%s%s_%s_%s_ang.png',savepath,cond_name,lock_name,ev_name));
%             saveas(fig,sprintf('%s%s_%s_%s_ang.pdf',savepath,cond_name,lock_name,ev_name));
%             close(fig)
            ang_lib{l_i,count} = {plt_t,plt_h_a,plt_e_a,plt_g_a,plt_dist};
%             % angle speed
%             fig = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
            plt_dist = my_norm(plt_dist_ori);
            scale = max(mean(plt_e_ad,2,'omitnan'));
            plt_dist = plt_dist * scale;
%             shadedErrorBar(plt_t,plt_e_ad', shaded_method, 'lineprops',{'b-','DisplayName','EyeAng','linewidth',3});
%             grid on; hold on;
%             shadedErrorBar(plt_t,plt_h_ad', shaded_method, 'lineprops',{'r-','DisplayName','HeadAng','linewidth',3})
%             shadedErrorBar(plt_t,plt_g_ad', shaded_method, 'lineprops',{'k-','DisplayName','GIPAng','linewidth',3})
%             shadedErrorBar(plt_t,plt_dist', shaded_method, 'lineprops',{'g-','DisplayName','dist2box','linewidth',3})
%             xline(0,'k--','linewidth',3,'DisplayName',sprintf('%s onset',lock_name));
% %             title(sprintf('%s lock - %s (%s)',lock_name, trial_name, cond_name));
%             set(gca,'fontsize',20)
%             set(gcf,'color','w')
%             xlabel('Time (ms)')
%             ylabel('Angular Speed (deg/sec)')
%             legend(findobj(gca,'-regexp','DisplayName', '[^'']'));
%             % for printing PDF
%             set(gca,'units','centimeters')
%             pos = get(gca,'Position');
%             ti = get(gca,'TightInset');
%             set(gcf, 'PaperUnits','centimeters');
%             set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%             set(gcf, 'PaperPositionMode', 'manual');
%             set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% %             saveas(fig,sprintf('%s%s_%s_%s_vang.png',savepath,cond_name,lock_name,ev_name));
%             saveas(fig,sprintf('%s%s_%s_%s_vang.pdf',savepath,cond_name,lock_name,ev_name));
%             close(fig)
            v_ang_lib{l_i,count} = {plt_t,plt_h_ad,plt_e_ad,plt_g_ad,plt_dist};
            count=count+1;
        end
    end
end

%% plot above figure into a big subplot
% ylabel_name = 'Angular (deg/s)';
ylabel_name = 'Angle (deg)';
plt_lib_std = reshape(reshape(ang_lib(1,:),3,2)',[],1); % for 3 by 2 subplot
plt_lib_dev = reshape(reshape(ang_lib(2,:),3,2)',[],1); % for 3 by 2 subplot
% plt_lib = subplot_lib; % for 2 by 3 subplot
ev_name = {'stim','stim','gip','gip','fix','fix'};
fig = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
ax_lib = cell(1,6);
plt_time_idx = [-500 1000];
for p_i = 1:6
    ax_lib{p_i} = subplot(3,2,p_i);
    plt_t = plt_lib_std{p_i}{1};
    plt_idx = plt_t >= plt_time_idx(1)&plt_t<=plt_time_idx(2);
    plt_t = plt_t(plt_idx);
    plt_h_a = plt_lib_std{p_i}{2}(plt_idx,:);
    plt_e_a = plt_lib_std{p_i}{3}(plt_idx,:);
    plt_g_a = plt_lib_std{p_i}{4}(plt_idx,:);
    plt_dist = plt_lib_std{p_i}{5}(plt_idx,:);
    
    shadedErrorBar(plt_t,plt_e_a', shaded_method, 'lineprops',{'b-','DisplayName','Eye','linewidth',3});
    grid on; hold on;
    shadedErrorBar(plt_t,plt_h_a', shaded_method, 'lineprops',{'r-','DisplayName','Head','linewidth',3})
%     shadedErrorBar(plt_t,plt_g_a', shaded_method, 'lineprops',{'k-','DisplayName','Gaze','linewidth',3})
%     shadedErrorBar(plt_t,plt_dist', shaded_method, 'lineprops',{'g-','DisplayName','Dist.','linewidth',3})
    xline(0,'k-','linewidth',3,'DisplayName','Event onset');
    set(gca,'fontsize',20)
    set(gcf,'color','w')
    if p_i == 3
        ylabel(ylabel_name)
    end
    if p_i==2
        legend(findobj(gca,'-regexp','DisplayName', '[^'']'));    
    end
    if ~ismember(p_i,[1,2])
        set(gca,'xtick',plt_time_idx(1):100:plt_time_idx(2))
    else
        set(gca,'xtick',plt_t(2):100:plt_t(end))
    end
    if ismember(p_i,[5,6])
        xlabel('Time (ms)')
    end
end
% linkaxes([ax_lib{:}],'y')

fig.Units = 'centimeters';
fig.PaperUnits = 'centimeters';
fig.PaperSize = fig.Position(3:4);
savepath = 'D:\Research\oddball_fig\xSubj\behav\';
% saveas(fig, sprintf('%sxSubj_vangle_dev.png',savepath));
% saveas(fig, sprintf('%sxSubj_vangle_dev.pdf',savepath));

%% compare behavior one by one
ylabel_name = 'Angular Velocity (deg/s)';
% ylabel_name = 'Angle (deg)';
plt_lib_std = reshape(reshape(ang_lib(1,:),3,2)',[],1); % for 3 by 2 subplot
plt_lib_dev = reshape(reshape(ang_lib(2,:),3,2)',[],1); % for 3 by 2 subplot
% plt_lib = subplot_lib; % for 2 by 3 subplot
ev_name = {'stim','stim','gip','gip','fix','fix'};
fig = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for p_i = 1:6
    ax_lib{p_i} = subplot(3,2,p_i);
    plt_t = plt_lib_std{p_i}{1};
    plt_idx = plt_t >= plt_time_idx(1)&plt_t<=plt_time_idx(2);
    plt_t = plt_t(plt_idx);
    plt_h_a_s = plt_lib_std{p_i}{2}(plt_idx,:);
    plt_e_a_s = plt_lib_std{p_i}{3}(plt_idx,:);
    plt_h_a_d = plt_lib_dev{p_i}{2}(plt_idx,:);
    plt_e_a_d = plt_lib_dev{p_i}{3}(plt_idx,:);
    % gaze
    plt_g_a_s = plt_lib_std{p_i}{4}(plt_idx,:);
    plt_g_a_d = plt_lib_dev{p_i}{4}(plt_idx,:);
    
    
    [~,p_h] = ttest(plt_h_a_s',plt_h_a_d');
    [~,p_e] = ttest(plt_e_a_s',plt_e_a_d');
    [~, h_h]=bonf_holm(p_h);
    [~, h_e]=bonf_holm(p_e);
    [~,p_g] = ttest(plt_g_a_s',plt_g_a_d');
    [~,h_g]=bonf_holm(p_g);
    
    
    
%     shadedErrorBar(plt_t,plt_e_a_s', shaded_method, 'lineprops',{'b-','DisplayName','Eye (Deviant)','linewidth',3});
%     grid on; hold on;
%     shadedErrorBar(plt_t,plt_e_a_d', shaded_method, 'lineprops',{'-','color',[0.5 0.5 1],'DisplayName','Eye (Standard)','linewidth',3});
%     shadedErrorBar(plt_t,plt_h_a_s', shaded_method, 'lineprops',{'r-','DisplayName','Head (Deviant)','linewidth',3});
%     shadedErrorBar(plt_t,plt_h_a_d', shaded_method, 'lineprops',{'-','color',[1 0.5 0.5],'DisplayName','Head (Deviant)','linewidth',3});
%     % significant
%     plot(plt_t(h_h), mean(plt_h_a_s(h_h,:),2), 'rx', 'linewidth',3,'markersize',10)
%     plot(plt_t(h_e), mean(plt_e_a_s(h_e,:),2), 'bo', 'linewidth',3,'markersize',10)
    shadedErrorBar(plt_t,plt_g_a_s', shaded_method, 'lineprops',{'b-','DisplayName','Deviant','linewidth',3});
    grid on; hold on;
    shadedErrorBar(plt_t,plt_g_a_d', shaded_method, 'lineprops',{'r-','DisplayName','Standard','linewidth',3});
    plot(plt_t(h_g), mean(plt_g_a_s(h_g,:),2), 'rx', 'linewidth',3,'markersize',10)
    xline(0,'k-','linewidth',3,'DisplayName','Event onset');
    set(gca,'fontsize',20)
    set(gcf,'color','w')
    if p_i == 3
        ylabel(ylabel_name)
    end
    if p_i==2
        legend(findobj(gca,'-regexp','DisplayName', '[^'']'));    
    end
    if ~ismember(p_i,[1,2])
        set(gca,'xtick',plt_time_idx(1):100:plt_time_idx(2))
    else
        set(gca,'xtick',plt_t(2):100:plt_t(end))
    end
    if ismember(p_i,[5,6])
        xlabel('Time (ms)')
    end
end






%% Merge epoch_lib to perform ERPImage using EEGLAB function
[fix_subj_idx, grab_subj_idx] = find_if_device(epoch_lib);
preserve_idx = sum(fix_subj_idx,1)==2;
merged_lib = merge_epoch_lib(epoch_lib(:,preserve_idx));
% merged_lib_new = merge_epoch_lib(epoch_lib(:,1:8));
% merged_lib_old = merge_epoch_lib(epoch_lib(:,9:end));

%% plot ERPImage
thres_amp = 30;
savepath = 'D:\Research\oddball_fig\xSubj\Biosemi only\ERPImage\';
for tarCh = {'CPz'}
    loc_path = sprintf('%s%s/',savepath,tarCh{:});
    if ~exist(loc_path,'dir')
        mkdir(loc_path)
    end
    for lock_name = {'stim','gip','fix'}
        for ev_name = {'cir','tri'}
            for cond_name = {'noHm','Hm'}
                switch cond_name{:}
                    case 'noHm'
                        cond_i = 1;
                    case 'Hm'
                        cond_i = 2;
                end
                switch lock_name{:}
                    case 'stim'
                        if strcmp(ev_name{:},'cir')
                            plt_EEG = merged_lib{cond_i}.std_epoch;
                        else
                            plt_EEG = merged_lib{cond_i}.dev_epoch;
                        end
                    case 'gip'
                        if strcmp(ev_name{:},'cir')
                            plt_EEG = merged_lib{cond_i}.gip_std;
                        else
                            plt_EEG = merged_lib{cond_i}.gip_dev;
                        end
                    case 'fix'
                        if strcmp(ev_name{:},'cir')
                            plt_EEG = merged_lib{cond_i}.fix_std;
                        else
                            plt_EEG = merged_lib{cond_i}.fix_dev;
                        end
                    case 'grab'
                        plt_EEG = merged_lib{cond_i}.grab_epoch;
                end

            plt_EEG = pop_select(plt_EEG,'channel',tarCh);
            % [plt_EEG,rm_idx] = pop_autorej(plt_EEG,'threshold',50,'nogui','on');
%             [plt_EEG, rm_idx] = pop_eegthresh(plt_EEG,1,1,-thres_amp,thres_amp,plt_EEG.xmin,plt_EEG.xmax,0,0);
%             plt_EEG = pop_rejepoch(plt_EEG, rm_idx, 0);
            idx_cir_1 = cellfun(@(x) ~isempty(regexp(x,'Standard','once')),{plt_EEG.event.type});
            idx_cir_2 = cellfun(@(x) ~isempty(regexp(x,'\(L\)','once')),{plt_EEG.event.type});
            idx_cir = idx_cir_1 | idx_cir_2;
            idx_tri_1 = cellfun(@(x) ~isempty(regexp(x,'Deviant','once')),{plt_EEG.event.type});
            idx_tri_2 = cellfun(@(x) ~isempty(regexp(x,'\(R\)','once')),{plt_EEG.event.type});
            idx_tri = idx_tri_1 | idx_tri_2;
            switch lock_name{:}
                case 'stim'
                    if strcmp(ev_name{:},'cir')
                        sort_name = {'circle_gip_start'};
                    else
                        sort_name = {'triangle_gip_start'};
                    end
                case 'gip'
                    if strcmp(ev_name{:},'cir')
                        sort_name = {plt_EEG.event(idx_cir).type};
                    else
                        sort_name = {plt_EEG.event(idx_tri).type};
                    end
                case 'fix'
                    if strcmp(ev_name{:},'cir')
                        sort_name = {'circle_gip_start'};
                    else
                        sort_name = {'triangle_gip_start'};
                    end
                case 'grab'
            %         sort_name = {'circle_gip_start'};
                    sort_name = {plt_EEG.event(idx_cir).type};
            end
            smoothing = 10;
            fig = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
            pop_erpimage(plt_EEG,1, [1],[[]],tarCh{:},smoothing,1,sort_name,[],'latency','yerplabel','\muV','erp','on','cbar','on','topo', { [1] plt_EEG.chanlocs plt_EEG.chaninfo } );
            % for printing PDF
            fig.Units = 'centimeters';
            fig.PaperUnits = 'centimeters';
            fig.PaperSize = fig.Position(3:4);
            saveas(fig, sprintf('%s%s_%s_%s_%s.png',loc_path,cond_name{:},lock_name{:},ev_name{:},tarCh{:}));
            close(fig)
            end
        end
    end
end

%% plot ERP
tarCh = 'CPz';
lock_name = 'stim';
cond_name = 'noHm';
switch cond_name
    case 'noHm'
        cond_i = 1;
    case 'Hm'
        cond_i = 2;
end
switch lock_name
    case 'stim'
        plt_EEG_cir = merged_lib{cond_i}.std_epoch;
        plt_EEG_tri = merged_lib{cond_i}.dev_epoch;
    case 'gip'
        plt_EEG_cir = merged_lib{cond_i}.gip_std;
        plt_EEG_tri = merged_lib{cond_i}.gip_dev;
    case 'fix'
        plt_EEG_cir = merged_lib{cond_i}.fix_std;
        plt_EEG_tri = merged_lib{cond_i}.fix_dev;
end
% shaded_method = {@(x)(mean(x,'omitnan')),@(x)([quantile(x,0.8)-mean(x,'omitnan');mean(x,'omitnan')-quantile(x,0.2)])};
% shaded_method = {@(x)(mean(x,'omitnan')),@(x)(std(x,'omitnan'))};
shaded_method = {@(x)(mean(x,'omitnan')),@(x)(std(x,'omitnan')/sqrt(size(x,1)))};
fig = plt_erp(plt_EEG_cir,plt_EEG_tri,tarCh,lock_name,shaded_method);

%% ttest between condition
p_lib = cell(2,3);
tarCh = 'Cz';
ch_idx = ismember({merged_lib{1}.dev_epoch.chanlocs.labels},tarCh);
for cond_i = 1:2
    for ev_i = 1:3
        switch ev_i
            case 1
                tri_data = squeeze(merged_lib{cond_i}.dev_epoch.data(ch_idx,:,:))';
                cir_data = squeeze(merged_lib{cond_i}.std_epoch.data(ch_idx,:,:))';
                plt_t = merged_lib{cond_i}.dev_epoch.times;
            case 2
                tri_data = squeeze(merged_lib{cond_i}.gip_dev.data(ch_idx,:,:))';
                cir_data = squeeze(merged_lib{cond_i}.gip_std.data(ch_idx,:,:))';
                plt_t = merged_lib{cond_i}.gip_dev.times;
            case 3
                tri_data = squeeze(merged_lib{cond_i}.fix_dev.data(ch_idx,:,:))';
                cir_data = squeeze(merged_lib{cond_i}.fix_std.data(ch_idx,:,:))';
                plt_t = merged_lib{cond_i}.fix_dev.times;
        end
        [~,p] = ttest2(tri_data,cir_data);
        p_lib{cond_i,ev_i} = p;
        figure
        plot(plt_t,mean(tri_data,1),'b-','linewidth',3);
        hold on
        grid on
        plot(plt_t,mean(cir_data,1),'r-','linewidth',3);
        plot(plt_t(p<=0.05), mean(cir_data(:,p<=0.05),1),'kx','linewidth',3,'markersize',15)
        
    end
end


%% plot event timing
% savepath = 'D:\Research\oddball_fig\xSubj\event latency\';
savepath = 'C:\Users\Yuan\OneDrive\Desktop\graduate!\fig\vr\';
thres_time = [100 1000];
for cond_name = {'noHm','Hm'}
    for lock_name = {'gip','fix'}
        switch cond_name{:}
            case 'noHm'
                cond_i = 1;
            case 'Hm'
                cond_i = 2;
        end
        diff_gip_std = merged_lib{cond_i}.event_time.gipStd_time - merged_lib{cond_i}.event_time.std_time;
        diff_gip_dev = merged_lib{cond_i}.event_time.gipDev_time - merged_lib{cond_i}.event_time.dev_time;
        diff_fix_std = merged_lib{cond_i}.event_time.fixStd_time - merged_lib{cond_i}.event_time.std_time;
        diff_fix_dev = merged_lib{cond_i}.event_time.fixDev_time - merged_lib{cond_i}.event_time.dev_time;
        diff_gip_std = diff_gip_std(diff_gip_std>=thres_time(1)&diff_gip_std<=thres_time(2));
        diff_gip_dev = diff_gip_dev(diff_gip_dev>=thres_time(1)&diff_gip_dev<=thres_time(2));
        diff_fix_std = diff_fix_std(diff_fix_std>=thres_time(1)&diff_fix_std<=thres_time(2));
        diff_fix_dev = diff_fix_dev(diff_fix_dev>=thres_time(1)&diff_fix_dev<=thres_time(2));
        if strcmp(lock_name{:},'gip')
            plt_std = diff_gip_std;
            plt_dev = diff_gip_dev;
        else
            plt_std = diff_fix_std;
            plt_dev = diff_fix_dev;
        end
        plt_std = plt_std(~isnan(plt_std));
        plt_dev = plt_dev(~isnan(plt_dev));

        fig = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
        [~,~,p] = statcond({plt_std,plt_dev},'method','bootstrap','naccu',1000);
        hs = histogram(plt_std,'binwidth',50,'Normalization','probability','DisplayName',sprintf('Deviant (med.=%dms)',round(median(plt_std))),...
            'DisplayStyle','Stairs','linewidth',3,'linestyle','-','edgecolor','r');
        hold on; grid on
        hd = histogram(plt_dev,'binwidth',50,'Normalization','probability','DisplayName',sprintf('Standard (med.=%dms)',round(median(plt_dev))),...
            'DisplayStyle','Stairs','linewidth',3,'linestyle','-','edgecolor','b');
        plot(nan,'DisplayName',sprintf('p = %g',p));
        legend(findobj(gca,'-regexp','DisplayName', '[^'']'));
        xlabel('latency (ms)')
        ylabel('Probability')
        set(gca,'fontsize',20)
%         title(sprintf('%s latency %s',lock_name{:},cond_name{:}))
        set(gca,'fontsize',20)
        set(gcf,'color','w')
        % for printing PDF
        set(gca,'units','centimeters')
        pos = get(gca,'Position');
        ti = get(gca,'TightInset');
        set(gcf, 'PaperUnits','centimeters');
        set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%         saveas(fig, sprintf('%s%s_%s_latency.png',savepath,cond_name{:},lock_name{:}));
%         saveas(fig, sprintf('%s%s_%s_latency.pdf',savepath,cond_name{:},lock_name{:}));
%         close(fig)
        
    end
end

%% investigate grab timing for merged condition
savepath = 'D:\Research\oddball_fig\xSubj\event latency\';
plt_func = @mean;
plt_lib = merged_lib;
fig = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
ax_lib = cell(1,2);
for cond_i = 1:2
    switch cond_i
        case 1
            cond_name = 'noHm';
            t_name = 'Eye-shifting';
        case 2
            cond_name = 'Hm';
            t_name = 'Freely-moving';
    end

    output = find_ev_time(plt_lib{cond_i});

    plt_grab_gip = output.cir.diff_grab_stim - output.cir.diff_gip_stim;
    plt_grab_fix = output.cir.diff_grab_stim - output.cir.diff_fix_stim;
    plt_fix_gip = output.cir.diff_fix_stim - output.cir.diff_gip_stim;
    
    ax1 = subplot(3,2,cond_i);
    [f,xi]=ksdensity(plt_fix_gip);
%     hs = histogram(plt_fix_gip,'binwidth',100,'Normalization','probability','DisplayName',sprintf('Fix-GIP (mean.=%dms)',round(plt_func(plt_fix_gip))),...
%                 'DisplayStyle','Stairs','linewidth',3,'linestyle','-','edgecolor','m');
    plot(xi,f,'m-','DisplayName',sprintf('mean = %dms',round(plt_func(plt_fix_gip))),'linewidth',3);
%     title(ax1, sprintf('%s',t_name));
    legend(findobj(gca,'-regexp','DisplayName', '[^'']'),'location','northwest');
    set(gca,'fontsize',20)
%     ylabel('Probability')
    grid on
    hold on
    ax2 = subplot(3,2,cond_i*2+(cond_i==1));
    grid on
    hold on
    [f,xi] = ksdensity(plt_grab_gip);
%     hs = histogram(plt_grab_gip,'binwidth',100,'Normalization','probability','DisplayName',sprintf('Grab-GIP (mean.=%dms)',round(plt_func(plt_grab_gip))),...
%                 'DisplayStyle','Stairs','linewidth',3,'linestyle','-','edgecolor','b');
    plot(xi, f, 'b-', 'DisplayName',sprintf('mean = %dms',round(plt_func(plt_grab_gip))),'linewidth',3);
    legend(findobj(gca,'-regexp','DisplayName', '[^'']'),'location','northwest');
    set(gca,'fontsize',20)
    if cond_i==1
        ylabel('Probability')
    end
    ax3 = subplot(3,2,cond_i*3+2*(cond_i==1));
    grid on
    hold on
    
    legend(findobj(gca,'-regexp','DisplayName', '[^'']'),'location','northwest');
    xlabel('latency (ms)')
%     ylabel('Probability')
    set(gca,'fontsize',20)
    set(gcf,'color','w')
    [f,xi] = ksdensity(plt_grab_fix);
%     hs = histogram(plt_grab_fix,'binwidth',100,'Normalization','probability','DisplayName',sprintf('Grab-Fix (mean.=%dms)',round(plt_func(plt_grab_fix))),...
%                 'DisplayStyle','Stairs','linewidth',3,'linestyle','-','edgecolor','r');
    plot(xi,f,'r-','DisplayName',sprintf('mean = %dms',round(plt_func(plt_grab_fix))),'linewidth',3);
    ax_lib{cond_i} = [ax1,ax2,ax3];
end

linkaxes([ax_lib{:}],'x')
fig.Units = 'centimeters';
fig.PaperUnits = 'centimeters';
fig.PaperSize = fig.Position(3:4);
% saveas(fig, sprintf('%sevent_latency_allDir_cir.png',savepath));
% saveas(fig, sprintf('%sevent_latency_allDir_cir.pdf',savepath));
% close(fig)

%% scatter plot among fixation, gip, and response
savepath = 'D:\Research\oddball_fig\xSubj\event latency\';
plt_func = @mean;
plt_lib = merged_lib;
fig = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
ax_lib = cell(1,2);
for cond_i = 1:2
    switch cond_i
        case 1
            cond_name = 'noHm';
            t_name = 'Eye-shifting';
        case 2
            cond_name = 'Hm';
            t_name = 'Freely-moving';
    end

    output = find_ev_time(plt_lib{cond_i});

    plt_grab_gip = output.cir.diff_grab_stim - output.cir.diff_gip_stim;
    plt_grab_fix = output.cir.diff_grab_stim - output.cir.diff_fix_stim;
    plt_fix_gip = output.cir.diff_fix_stim - output.cir.diff_gip_stim;
    
    ax1 = subplot(1,2,cond_i);
    x = plt_fix_gip;
    y = plt_grab_fix;
    [corr_R, corr_p] = corrcoef(x,y);
    plot(x, y,'bo','DisplayName',sprintf('Corr. = %.3f, p val. = %.3f', corr_R(1,2),corr_p(1,2)),'linewidth',3,'markersize',15);
    set(gca,'fontsize',20)
    xlabel('GIP-Fix')
    ylabel('Response-Fix')
    title('Latency: GIP-Fix vs Response-Fix')
    hold on
    grid on
    pfit = polyfit(x,y,1);
    f = polyval(pfit,x);
    r2 = 1 - sum((y-f).^2)/sum((y-mean(y)).^2);
    plot(x,f,'k-','DisplayName',sprintf('R^2 = %.3f',r2),'linewidth',3);
    legend(findobj(gca,'-regexp','DisplayName', '[^'']'),'location','northeast');
    pbaspect([1 1 1])
end

fig.Units = 'centimeters';
fig.PaperUnits = 'centimeters';
fig.PaperSize = fig.Position(3:4);
% saveas(fig, sprintf('%sevent_latency_allDir_cir.png',savepath));
% saveas(fig, sprintf('%sevent_latency_allDir_cir.pdf',savepath));
% close(fig)

%% investigate event onset timing for merged condition
plt_t_lib = cell(2,2);
plt_merged_lib = merged_lib_old;
for cond_i = 1:2
    for ev_i = 1:2
        if ev_i== 1
            ev_name = 'circle_fix_start';
            stim_epoch = plt_merged_lib{cond_i}.std_epoch;
        else
            ev_name = 'triangle_fix_start';
            stim_epoch = plt_merged_lib{cond_i}.dev_epoch;
        end
        % get event latency
        stim_idx = cellfun(@(x) find(cellfun(@(y) ~isempty(regexp(y,'Ring','ONCE')), x)), {stim_epoch.epoch.eventtype}, 'uniformoutput',0);
        gip_idx = cellfun(@(x,x_stim) find(cellfun(@(y) ismember(y,{'Up','Bottom','Left','Right'}),x(x_stim:end)),1)+x_stim-1, {stim_epoch.epoch.eventtype},stim_idx,'uniformoutput',0);
        fix_idx = cellfun(@(x,x_stim) find(cellfun(@(y) strcmp(y,ev_name),x(x_stim:end)),1)+x_stim-1, {stim_epoch.epoch.eventtype},stim_idx,'uniformoutput',0);
        rm_idx = cellfun(@isempty, gip_idx)|cellfun(@isempty,fix_idx);
        if strcmp(ev_name, 'circle_fix_start')
            grab_idx = cellfun(@(x) find(cellfun(@(y) strcmp(y,'grab'),x),1), {stim_epoch.epoch.eventtype},'uniformoutput',0);
            rm_idx = rm_idx|cellfun(@isempty, grab_idx);
            grab_idx = grab_idx(~rm_idx);
            grab_t = zeros(1,length(grab_idx));
        else
            grab_t = [];
        end
        gip_idx = gip_idx(~rm_idx);
        fix_idx = fix_idx(~rm_idx);
        stim_idx = stim_idx(~rm_idx);
        stim_epoch = pop_rejepoch(stim_epoch,rm_idx,0);
        % find time

        gip_t =  zeros(1,length(gip_idx));
        fix_t =  zeros(1,length(gip_idx));
        for i = 1:length(gip_idx)
            if strcmp(ev_name, 'circle_fix_start')
                grab_t(i) = stim_epoch.epoch(i).eventlatency{grab_idx{i}}; % ms
            end
            gip_t(i) = stim_epoch.epoch(i).eventlatency{gip_idx{i}}; % ms
            fix_t(i) = stim_epoch.epoch(i).eventlatency{fix_idx{i}}; % ms
        end
        plt_t_lib{cond_i,ev_i} = {gip_t,fix_t,grab_t};
    end
end

savepath = 'D:\Research\oddball_fig\xSubj\event latency\';
plt_func = @mean;
fig = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
ax_lib = cell(2,2);
for cond_i = 1:2
    switch cond_i
        case 1
            cond_name = 'noHm';
            t_name = 'Eye-shifting';
        case 2
            cond_name = 'Hm';
            t_name = 'Head-turning';
    end
    % std
    gip_t_std = plt_t_lib{cond_i,1}{1};
    fix_t_std = plt_t_lib{cond_i,1}{2};
    grab_t_std = plt_t_lib{cond_i,1}{3};
    [f_gip_std, xi_gip_std] = ksdensity(gip_t_std);
    [f_fix_std, xi_fix_std] = ksdensity(fix_t_std);
    [f_grab_std, xi_grab_std] = ksdensity(grab_t_std);
    % dev
    gip_t_dev= plt_t_lib{cond_i,2}{1};
    fix_t_dev = plt_t_lib{cond_i,2}{2};
    [f_gip_dev, xi_gip_dev] = ksdensity(gip_t_dev);
    [f_fix_dev, xi_fix_dev] = ksdensity(fix_t_dev);
    % max
    [~, mi] = max(f_gip_std);
    max_gip_std = xi_gip_std(mi);
         [~, mi] = max(f_fix_std);
    max_fix_std = xi_fix_std(mi);
    [~, mi] = max(f_grab_std);
    max_grab_std = xi_grab_std(mi);
    [~, mi] = max(f_gip_dev);
    max_gip_dev = xi_gip_dev(mi);
    [~, mi] = max(f_fix_dev);
    max_fix_dev = xi_fix_dev(mi);

    ax = subplot(2,1,cond_i);
%     title(ax1, sprintf('%s',t_name));
    plot(xi_grab_std, f_grab_std, 'c-', 'DisplayName',sprintf('Response (mean = %dms)',round(max_grab_std)),'linewidth',3);
    hold on
    grid on
    plot(xi_fix_std, f_fix_std, 'b-', 'DisplayName',sprintf('Fix (mean = %dms)',round(max_fix_std)),'linewidth',3);
    plot(xi_gip_std, f_gip_std, 'k-', 'DisplayName',sprintf('GIP (mean = %dms)',round(max_gip_std)),'linewidth',3);
%     plot(xi_gip_dev, f_gip_dev, 'b--', 'DisplayName',sprintf('GIP Standard (mean = %dms)',round(max_gip_dev)),'linewidth',3);
%     plot(xi_fix_dev, f_fix_dev, 'r--', 'DisplayName',sprintf('Fix Standard (mean = %dms)',round(max_fix_dev)),'linewidth',3);
    
    
    legend(findobj(gca,'-regexp','DisplayName', '[^'']'),'location','northeast');
    set(gca,'fontsize',30)
    ylabel('Probability')

    ax_lib{cond_i} = ax;
end
set(gcf,'color','w')
xlabel('Time (ms)')
linkaxes([ax_lib{:}],'x')
fig.Units = 'centimeters';
fig.PaperUnits = 'centimeters';
fig.PaperSize = fig.Position(3:4);
% saveas(fig, sprintf('%sevent_onset.png',savepath));
% saveas(fig, sprintf('%sevent_onset.pdf',savepath));

%% compare behavior
% merged
ev_name = 'stim';
cond_i = 2;
plt_t = epoch_lib{1}.std_epoch.times;

plt_behav = struct('h_a_cir',[],'h_a_dev',[],'e_a_cir',[],'e_a_dev',[],'g_a_cir',[],'g_a_dev',[],...
                   'h_ad_cir',[],'h_ad_dev',[],'e_ad_cir',[],'e_ad_dev',[],'g_ad_cir',[],'g_ad_dev',[]);
for i = 1:length(behav_lib(cond_i,:))
    switch ev_name
        case 'stim'
            e_name_std = 'std_epoch';
            e_name_dev = 'dev_epoch';
    end
            
    h_a_cir = behav_lib{cond_i,i}.(e_name_std).headAng;
    h_ad_cir = behav_lib{cond_i,i}.(e_name_std).headRot;
    e_a_cir = behav_lib{cond_i,i}.(e_name_std).eyeAng;
    e_ad_cir = behav_lib{cond_i,i}.(e_name_std).eyeRot;
    g_a_cir = behav_lib{cond_i,i}.(e_name_std).gipAng;
    g_ad_cir = behav_lib{cond_i,i}.(e_name_std).gipRot;
    
    h_a_dev = behav_lib{cond_i,i}.(e_name_dev).headAng;
    h_ad_dev = behav_lib{cond_i,i}.(e_name_dev).headRot;
    e_a_dev = behav_lib{cond_i,i}.(e_name_dev).eyeAng;
    e_ad_dev = behav_lib{cond_i,i}.(e_name_dev).eyeRot;
    g_a_dev = behav_lib{cond_i,i}.(e_name_dev).gipAng;
    g_ad_dev = behav_lib{cond_i,i}.(e_name_dev).gipRot;
    
    plt_behav.h_a_cir = [plt_behav.h_a_cir, h_a_cir];
    plt_behav.h_ad_cir = [plt_behav.h_ad_cir, h_ad_cir];
    plt_behav.e_a_cir = [plt_behav.e_a_cir, e_a_cir];
    plt_behav.e_ad_cir = [plt_behav.e_ad_cir, e_ad_cir];
    plt_behav.g_a_cir = [plt_behav.g_a_cir, g_a_cir];
    plt_behav.g_ad_cir = [plt_behav.g_ad_cir, g_ad_cir];
    
    plt_behav.h_a_dev = [plt_behav.h_a_dev, h_a_dev];
    plt_behav.h_ad_dev = [plt_behav.h_ad_dev, h_ad_dev];
    plt_behav.e_a_dev = [plt_behav.e_a_dev, e_a_dev];
    plt_behav.e_ad_dev = [plt_behav.e_ad_dev, e_ad_dev];
    plt_behav.g_a_dev = [plt_behav.g_a_dev, g_a_dev];
    plt_behav.g_ad_dev = [plt_behav.g_ad_dev, g_ad_dev];
    
end


%%
plt_1 = 'e_ad_dev';
plt_2 = 'e_ad_cir';
figure
h = ttest2(plt_behav.(plt_1)',plt_behav.(plt_2)');
h(isnan(h)) = 0;
h = logical(h);
plot(plt_t, mean(plt_behav.(plt_1),2,'omitnan'), 'b-', 'linewidth',3);
hold on
grid on
plot(plt_t, mean(plt_behav.(plt_2),2,'omitnan'), 'r-', 'linewidth',3);
plot(plt_t(h), mean(plt_behav.(plt_2)(h,:),2,'omitnan'),'kx','linewidth',3,'markersize',15);
