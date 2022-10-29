%% plot for paper
filepath = 'D:\Research\';
% load([filepath,'behav_lib_20221022.mat']);
load([filepath,'epoch_lib_rmPreStim.mat']);
% filepath = '//hoarding/yuan/Documents/2021 HM_visual_oddball/dataset/new epoch/';
% subj_list = {dir([filepath, 'rmPreStim*']).name};
subj_list = 1:14;

% print_multi([saveFigPath,'v_fig1_',tName],{'pdf'})

%% plot the mean of mean ERP across subject
thres_time = [100 1000];
rm_thres = 30;
nb_subj = size(epoch_lib,2);
shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(nb_subj))};
% shaded_method = {@(x)(mean(x,'omitnan')),@(x)([quantile(x,0.8)-mean(x,'omitnan');mean(x,'omitnan')-quantile(x,0.2)])};
savepath = 'D:\Research\oddball_fig\xSubj\mean ERP\';
for cond_name = {'noHm','Hm'}
    for ev_name = {'stim','gip','fix'}
        for tarCh = {'Cz'}
%             loc_path = sprintf('%s%s/',savepath,tarCh{:});
            loc_path = 'C:\Users\Yuan\OneDrive\Desktop\graduate!\fig\vr\';
            if ~exist(loc_path,'dir')
                mkdir(loc_path)
            end
            fig = plt_erp_meanXsubj(epoch_lib, cond_name{:}, ev_name{:}, tarCh{:}, thres_time, rm_thres,shaded_method);
            % for printing PDF
            set(gca,'units','centimeters')
            pos = get(gca,'Position');
            ti = get(gca,'TightInset');
            set(gcf, 'PaperUnits','centimeters');
            set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%             saveas(fig, sprintf('%serp_%s_%s_%s.png',loc_path,cond_name{:},ev_name{:},tarCh{:}));
            saveas(fig, sprintf('%serp_%s_%s_%s.pdf',loc_path,cond_name{:},ev_name{:},tarCh{:}));
            close(fig)
        end
    end
end

%% calculate behavior
behav_lib = cell(size(epoch_lib));
for subj_i = 1:size(behav_lib,2)
    for cond_i = 1:2
        behav_lib{cond_i,subj_i} = cal_behav(epoch_lib{cond_i,subj_i});
    end
end
        

%% plot behavior
% avoid outlier near the edges
rm_edge = round(100*0.001*epoch_lib{1}.std_epoch.srate); % sample points
% savepath = 'D:\Research\oddball_fig\xSubj\behav\';
savepath = 'C:\Users\Yuan\OneDrive\Desktop\graduate!\fig\vr\';

for cond_name_loop = {'noHm','Hm'}
    for lock_name_loop = {'stim','gip','fix'}
        for ev_name_loop = {'dev'}
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
                case 'dev'
                    trial_name = 'TRI';
            end
            [fix_subj_idx, grab_subj_idx] = find_if_device(epoch_lib);
            fix_subj_idx = sum(fix_subj_idx)==2;
            grab_subj_idx = sum(grab_subj_idx)==2;
            switch lock_name
                case 'stim'
                    epoch_name = sprintf('%s_epoch',ev_name);
                    plt_t = epoch_lib{1}.dev_epoch.times;
                    include_subj = true(1,size(behav_lib,2));
                case 'gip'
                    epoch_name = sprintf('gip_%s',ev_name);
                    plt_t = epoch_lib{1}.gip_dev.times;
                    include_subj = true(1,size(behav_lib,2));
                case 'fix'
                    epoch_name = sprintf('fix_%s',ev_name);
                    plt_t = epoch_lib{1}.gip_dev.times;
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
            %     shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(length(subj_list)))};
            shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan'))};
            my_norm = @(x)((x-min(x,[],1))./(max(x,[],1)-min(x,[],1)));

            % angle 
            fig = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
            plt_dist = my_norm(plt_dist_ori);
            scale = max(mean(plt_e_a,2,'omitnan'));
            plt_dist = plt_dist * scale;
            shadedErrorBar(plt_t,plt_e_a', shaded_method, 'lineprops',{'b-','DisplayName','EyeAng','linewidth',3});
            grid on; hold on;
            shadedErrorBar(plt_t,plt_h_a', shaded_method, 'lineprops',{'r-','DisplayName','HeadAng','linewidth',3})
            shadedErrorBar(plt_t,plt_g_a', shaded_method, 'lineprops',{'k-','DisplayName','GIPAng','linewidth',3})
            shadedErrorBar(plt_t,plt_dist', shaded_method, 'lineprops',{'g-','DisplayName','dist2box','linewidth',3})
            xline(0,'k--','linewidth',3,'DisplayName',sprintf('%s onset',lock_name));
%             title(sprintf('%s lock - %s (%s)',lock_name, trial_name, cond_name));
            set(gca,'fontsize',20)
            set(gcf,'color','w')
            xlabel('Time (ms)')
            ylabel('Angle (deg)')
            legend(findobj(gca,'-regexp','DisplayName', '[^'']'));
            % for printing PDF
            set(gca,'units','centimeters')
            pos = get(gca,'Position');
            ti = get(gca,'TightInset');
            set(gcf, 'PaperUnits','centimeters');
            set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%             saveas(fig,sprintf('%s%s_%s_%s_ang.png',savepath,cond_name,lock_name,ev_name));
            saveas(fig,sprintf('%s%s_%s_%s_ang.pdf',savepath,cond_name,lock_name,ev_name));
            close(fig)

            % angle speed
            fig = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
            plt_dist = my_norm(plt_dist_ori);
            scale = max(mean(plt_e_ad,2,'omitnan'));
            plt_dist = plt_dist * scale;
            shadedErrorBar(plt_t,plt_e_ad', shaded_method, 'lineprops',{'b-','DisplayName','EyeAng','linewidth',3});
            grid on; hold on;
            shadedErrorBar(plt_t,plt_h_ad', shaded_method, 'lineprops',{'r-','DisplayName','HeadAng','linewidth',3})
            shadedErrorBar(plt_t,plt_g_ad', shaded_method, 'lineprops',{'k-','DisplayName','GIPAng','linewidth',3})
            shadedErrorBar(plt_t,plt_dist', shaded_method, 'lineprops',{'g-','DisplayName','dist2box','linewidth',3})
            xline(0,'k--','linewidth',3,'DisplayName',sprintf('%s onset',lock_name));
%             title(sprintf('%s lock - %s (%s)',lock_name, trial_name, cond_name));
            set(gca,'fontsize',20)
            set(gcf,'color','w')
            xlabel('Time (ms)')
            ylabel('Angular Speed (deg/sec)')
            legend(findobj(gca,'-regexp','DisplayName', '[^'']'));
            % for printing PDF
            set(gca,'units','centimeters')
            pos = get(gca,'Position');
            ti = get(gca,'TightInset');
            set(gcf, 'PaperUnits','centimeters');
            set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%             saveas(fig,sprintf('%s%s_%s_%s_vang.png',savepath,cond_name,lock_name,ev_name));
            saveas(fig,sprintf('%s%s_%s_%s_vang.pdf',savepath,cond_name,lock_name,ev_name));
            close(fig)
        end
    end
end


%% Merge epoch_lib to perform ERPImage using EEGLAB function
output = merge_epoch_lib(epoch_lib);

%% plot ERPImage
thres_amp = 30;
savepath = 'D:\Research\oddball_fig\xSubj\ERPImage\';
for tarCh = {'Cz','CPz','POz','O2'}
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
                            plt_EEG = output{cond_i}.std_epoch;
                        else
                            plt_EEG = output{cond_i}.dev_epoch;
                        end
                    case 'gip'
                        if strcmp(ev_name{:},'cir')
                            plt_EEG = output{cond_i}.gip_std;
                        else
                            plt_EEG = output{cond_i}.gip_dev;
                        end
                    case 'fix'
                        if strcmp(ev_name{:},'cir')
                            plt_EEG = output{cond_i}.fix_std;
                        else
                            plt_EEG = output{cond_i}.fix_dev;
                        end
                    case 'grab'
                        plt_EEG = output{cond_i}.grab_epoch;
                end

            plt_EEG = pop_select(plt_EEG,'channel',tarCh);
            % [plt_EEG,rm_idx] = pop_autorej(plt_EEG,'threshold',50,'nogui','on');
            [plt_EEG, rm_idx] = pop_eegthresh(plt_EEG,1,1,-thres_amp,thres_amp,plt_EEG.xmin,plt_EEG.xmax,0,0);
            plt_EEG = pop_rejepoch(plt_EEG, rm_idx, 0);
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
            set(gca,'units','centimeters')
            pos = get(gca,'Position');
            ti = get(gca,'TightInset');
            set(gcf, 'PaperUnits','centimeters');
            set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
            saveas(fig, sprintf('%s%s_%s_%s_%s.png',loc_path,cond_name{:},lock_name{:},ev_name{:},tarCh{:}));
            close(fig)
            end
        end
    end
end

%% plot ERP
tarCh = 'Cz';
lock_name = 'fix';
cond_name = 'Hm';
switch cond_name
    case 'noHm'
        cond_i = 1;
    case 'Hm'
        cond_i = 2;
end
switch lock_name
    case 'stim'
        plt_EEG_cir = output{cond_i}.std_epoch;
        plt_EEG_tri = output{cond_i}.dev_epoch;
    case 'gip'
        plt_EEG_cir = output{cond_i}.gip_std;
        plt_EEG_tri = output{cond_i}.gip_dev;
    case 'fix'
        plt_EEG_cir = output{cond_i}.fix_std;
        plt_EEG_tri = output{cond_i}.fix_dev;
end
shaded_method = {@(x)(mean(x,'omitnan')),@(x)([quantile(x,0.8)-mean(x,'omitnan');mean(x,'omitnan')-quantile(x,0.2)])};
% shaded_method = {@(x)(mean(x,'omitnan')),@(x)(std(x,'omitnan'))};
fig = plt_erp(plt_EEG_cir,plt_EEG_tri,tarCh,lock_name,shaded_method);

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
        diff_gip_std = output{cond_i}.event_time.gipStd_time - output{cond_i}.event_time.std_time;
        diff_gip_dev = output{cond_i}.event_time.gipDev_time - output{cond_i}.event_time.dev_time;
        diff_fix_std = output{cond_i}.event_time.fixStd_time - output{cond_i}.event_time.std_time;
        diff_fix_dev = output{cond_i}.event_time.fixDev_time - output{cond_i}.event_time.dev_time;
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
        saveas(fig, sprintf('%s%s_%s_latency.pdf',savepath,cond_name{:},lock_name{:}));
        close(fig)
        
    end
end



%% Single subject analysis
saveFigPath = 'D:\Research\oddball_fig\';
if ~exist(saveFigPath,'dir')
    mkdir(saveFigPath)
end
tarCh = {'Cz','CPz','POz','O2'};
rm_thres = 30;
filepath = '//hoarding/yuan/Documents/2021 HM_visual_oddball/dataset/new epoch/';
subj_list = {dir([filepath, 'rmPreStim*']).name};


for subj_i = 1:size(epoch_lib,2)
    subj_savepath = sprintf('%s%s/',saveFigPath,subj_list{subj_i}(1:end-4));
    if ~exist(subj_savepath,'dir')
        mkdir(subj_savepath)
    end
    for cond_i = 1:2
        switch cond_i
            case 1
                cond_savepath = sprintf('%snoHm/',subj_savepath);
                cond_name = 'noHm';
            case 2
                cond_savepath = sprintf('%sHm/',subj_savepath);
                cond_name = 'Hm';
        end
        if ~exist(cond_savepath,'dir')
            mkdir(cond_savepath)
        end 
        % get data
        tmp_epoch = my_rmEpoch(epoch_lib{cond_i,subj_i});
        std_epoch = tmp_epoch.std_epoch;
        dev_epoch = tmp_epoch.dev_epoch;
        gip_std = tmp_epoch.gip_std;
        gip_dev= tmp_epoch.gip_dev;
        
        for ch_i = 1:length(tarCh)
            ch_savepath = sprintf('%s%s/',cond_savepath,tarCh{ch_i});
            if ~exist(ch_savepath,'dir')
                mkdir(ch_savepath)
            end 
            %=================
            % remove trial with bad baseline
            [rm_idx_stim, rm_idx_gip] = my_rmbase(std_epoch, gip_std, tmp_epoch.event_time.gipStd_time,...
                                                  tarCh{ch_i}, rm_thres);
            std_epoch = pop_rejepoch(std_epoch,rm_idx_stim,0);
            gip_std = pop_rejepoch(gip_std,rm_idx_gip,0);
            [rm_idx_stim, rm_idx_gip] = my_rmbase(dev_epoch, gip_dev, tmp_epoch.event_time.gipDev_time,...
                                                  tarCh{ch_i}, rm_thres);
            dev_epoch = pop_rejepoch(dev_epoch,rm_idx_stim,0);
            gip_dev = pop_rejepoch(gip_dev,rm_idx_gip,0);
            %=================
            % plot ERP
            % stim lock
            shaded_method = {@(x)(mean(x,'omitnan')),@(x)([quantile(x,0.8)-mean(x,'omitnan');mean(x,'omitnan')-quantile(x,0.2)])};
            fig = plt_erp(std_epoch,dev_epoch,tarCh{ch_i},'stim',shaded_method);
            saveas(fig, sprintf('%serp_%s_stim_%s.png',ch_savepath,cond_name, tarCh{ch_i}));
            close(fig)
            % GIP lock
            shaded_method = {@(x)(mean(x,'omitnan')),@(x)([quantile(x,0.8)-mean(x,'omitnan');mean(x,'omitnan')-quantile(x,0.2)])};
            fig = plt_erp(gip_std,gip_dev,tarCh{ch_i},'gip',shaded_method);
            saveas(fig, sprintf('%serp_%s_gip_%s.png',ch_savepath,cond_name, tarCh{ch_i}));
            close(fig)
            %=================
            % plot ERSP
            % stim lock, circle trial
            smoothing = 1;
            fig = plt_ersp(std_epoch,tarCh{ch_i},'stim','cir',smoothing);
            saveas(fig, sprintf('%sersp_%s_stim_cir_%s.png',ch_savepath,cond_name, tarCh{ch_i}));
            close(fig)
            % stim lock, triangle trial
            fig = plt_ersp(dev_epoch,tarCh{ch_i},'stim','tri',smoothing);
            saveas(fig, sprintf('%sersp_%s_stim_tri_%s.png',ch_savepath,cond_name, tarCh{ch_i}));
            close(fig)
            % gip lock, circle trial
            fig = plt_ersp(gip_std,tarCh{ch_i},'gip','cir',smoothing);
            saveas(fig, sprintf('%sersp_%s_gip_cir_%s.png',ch_savepath,cond_name, tarCh{ch_i}));
            close(fig)
            % gip lock, triangle trial
            fig = plt_ersp(gip_dev,tarCh{ch_i},'gip','tri',smoothing);
            saveas(fig, sprintf('%sersp_%s_gip_tri_%s.png',ch_savepath,cond_name, tarCh{ch_i}));
            close(fig)
            %=================
            % plot ERP diff
            % stim lock
            fig = plt_erp_diff(std_epoch,dev_epoch,tarCh{ch_i},'stim');
            saveas(fig, sprintf('%serpDiff_%s_stim_%s.png',ch_savepath,cond_name, tarCh{ch_i}));
            close(fig)
            % GIP lock
            fig = plt_erp_diff(gip_std,gip_dev,tarCh{ch_i},'gip');
            saveas(fig, sprintf('%serpDiff_%s_gip_%s.png',ch_savepath,cond_name, tarCh{ch_i}));
            close(fig)
        end
    end
end
            
            
    
    
