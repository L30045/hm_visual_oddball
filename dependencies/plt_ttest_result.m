% plot t-test result for the selected channels
%% plot for paper
filepath = 'D:\Research\';
% load([filepath,'behav_lib_20221022.mat']);
load([filepath,'epoch_lib_rmPreStim_new.mat']);
disp('Done')

%% remove those don't have fixation
filepath = '/data/projects/yuan/2021 HM_visual_oddball/dataset/20230508/';
subj_list = cellfun(@(x) x(end-7:end-4),{dir([filepath, 'rmPreStim*']).name},'uniformoutput',0);
[fix_subj_idx, grab_subj_idx] = find_if_device(epoch_lib);
% include subject without grab event (all subjects are included)
preserve_idx = sum(fix_subj_idx,1)==2;
% & sum(grab_subj_idx,1)==2;
% remove 1998 2029
rm_subj = {'1998','2029','2059','2085'};
rm_noFRP = {'2057','2030','2050'};
preserve_idx(ismember(subj_list,rm_subj)) = 0;
% preserve_subj = {'1145','1151','1988','1997','2000','2004','2022','2031','2050','2056'};
% preserve_idx(~ismember(subj_list,preserve_subj)) = 0;
subj_list = subj_list(preserve_idx);


%% find cross channels
% Rename channels label
cross_channels = {'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3',...
               'FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1',...
               'P3','P5','P7','P9','PO7','PO3','O1','Iz','Oz','POz','Pz','CPz',...
               'Fpz','Fp2','AF8','AF4','AFz','Fz','F2','F4','F6','F8',...
               'FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8',...
               'TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2'};
tmp_epoch = epoch_lib(:);
for epoch_i = 1:length(tmp_epoch)
    tmp_ch = {tmp_epoch{epoch_i}.std_epoch.chanlocs.labels};
    cross_channels = intersect(cross_channels,tmp_ch);
end

% subselect channels
select_ch = {'Fz','FCz','Cz','CPz','Pz','POz','Oz','O1','O2','CP3','CP4','C3','C4'};
% select_ch = {'Fz','FCz','Cz','CPz','Pz','POz','Oz','O1','O2','P3','P4','CP3','CP4',...
%     'C3','C4'};
% select_ch = {'Fz','FCz','Cz','CPz','Pz','POz','Oz','O1','O2','PO4','PO3','P3','P7','P4','P8','CP3','CP4',...
%     'C3','C4','F5','F4'};
% select_ch = {'Cz'};

%% plot the mean of mean ERP across subject
plt_epoch_lib = epoch_lib(:,preserve_idx);
thres_time = [100 1000];
rm_thres = 30;
nb_subj = size(plt_epoch_lib,2);
shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(nb_subj))};
% shaded_method = {@(x)(median(x,'omitnan')),@(x)([quantile(x,0.8)-median(x,'omitnan');median(x,'omitnan')-quantile(x,0.2)])};
savepath = 'D:\Research\oddball_fig\xSubj\mean ERP\';

ch_lib = cell(1,length(select_ch));
parfor tarch_i = 1:length(select_ch)
    subplot_lib = cell(1,6);
    count = 1;
    for cond_i = 1:2
        switch cond_i
            case 1
                cond_name = 'noHm';
            case 2
                cond_name = 'Hm';
        end
        for ev_i = 1:3
            switch ev_i
                case 1
                    ev_name = 'stim';
                case 2
                    ev_name = 'gip';
                case 3
                    ev_name = 'fix';
            end
                   
            [fig,plt_t,cir_lib,tri_lib] = plt_erp_meanXsubj(plt_epoch_lib, cond_name, ev_name, select_ch{tarch_i}, thres_time, rm_thres,shaded_method, false);
            subplot_lib{count} = {plt_t, cir_lib, tri_lib};
            count = count+1;
        end
    end
    ch_lib{tarch_i} = subplot_lib;
end
disp('Done')

%% plot channel ERP
figure
tar_ch = 'Cz';
select_subj = true(1,size(plt_epoch_lib,2));
% select_subj = 1:13;
shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(nb_subj))};
% shaded_method = {@(x)(median(x,'omitnan')),@(x)([quantile(x,0.8)-median(x,'omitnan');median(x,'omitnan')-quantile(x,0.2)])};
tar_ch_idx = ismember(select_ch, tar_ch);
plt_lib = reshape(reshape(ch_lib{tar_ch_idx},3,2)',[],1); % for 3 by 2 subplot
ax_lib = cell(1,6);
t_name = reshape(repmat({'Stim','GIP','Fix'},2,1),1,[]);
for p_i = 1:6
    ax_lib{p_i} = subplot(3,2,p_i);
    plt_t = plt_lib{p_i}{1};
    plt_cir = plt_lib{p_i}{2}(select_subj,:);
    plt_tri = plt_lib{p_i}{3}(select_subj,:);
    %======== shaded
    ht = shadedErrorBar(plt_t, plt_tri, shaded_method,'lineprops',...
            {'color','b','linewidth',3,'DisplayName','Standard'});
    ht.patch.FaceAlpha = 0.1;
    hold on
    grid on
    hc = shadedErrorBar(plt_t, plt_cir, shaded_method,'lineprops',...
            {'color','r','linewidth',3,'DisplayName','Standard'});
    hc.patch.FaceAlpha = 0.1;
    xlabel('Time (ms)')
    ylabel('Amplitude (\muV)')
    xline(0,'k','linewidth',3)
    set(gcf,'color','w')
    set(gca,'xtick',plt_t(1:25:end));
    set(gca,'xticklabels',plt_t(1:25:length(plt_t)));

    %======= stat test
    plt_idx = plt_t >= plt_time_idx(1) & plt_t<=plt_time_idx(2);
    plot(plt_t(plt_idx),ttest_lib{p_i},'kx','linewidth',3);
    %======= pcolor
%     [nr, nc] = size(plt_cir);
%     ht = pcolor([plt_cir,nan(nr,1);nan(1,nc+1)]);
%     set(ht,'edgecolor','none')
%     colorbar
    
    
    title(t_name{p_i})
end
linkaxes([ax_lib{:}],'y')

%% ttest
ttest_lib = cell(1,6);
pval_lib = cell(1,6);
plt_t_lib = cell(1,6);
plt_time_idx = [-500 1000];
test_time_idx = [-100 1000];
is_ttest = true;
is_correct = true;
% sub_ch = {'Cz','CPz','Pz','POz','Oz'};
sub_ch = {'Cz'};
sub_ch_idx = ismember(select_ch,sub_ch);
sub_ch_lib = ch_lib(sub_ch_idx);

for p_i = 1:6
    tmp_ttest_result = [];
    tmp_p_val = [];
    for t_i = 1:length(sub_ch_lib)
        plt_lib = reshape(reshape(sub_ch_lib{t_i},3,2)',[],1); % for 3 by 2 subplot
        plt_t = plt_lib{p_i}{1};
        h_tmp = zeros(1,length(plt_t));
        p_tmp = nan(1,length(plt_t));
        plt_idx = plt_t >= plt_time_idx(1) & plt_t<=plt_time_idx(2);
        test_idx = plt_t >= test_time_idx(1) & plt_t<=test_time_idx(2);
        plt_t = plt_t(plt_idx);
        cir_lib = plt_lib{p_i}{2}(:,test_idx);
        tri_lib = plt_lib{p_i}{3}(:,test_idx);

        if is_ttest
            [~, p] = ttest(tri_lib,cir_lib);
        %     [corrected_p, h] = bonf_holm(p);
            fdrType = 'parametric';
        else
            p = rowfun(@ranksum,table(cir_lib',tri_lib'));
            p = p{:,:}';
            fdrType = 'nonparametric';
        end
        
        if is_correct
            [p_fdr, h] = fdr(p, 0.05, fdrType);
%             [p_fdr, h] = fdr(p, 0.05);
        else
            h = p <= 0.05;
        end  
        
        h_tmp(test_idx) = h;
        p_tmp(test_idx) = p;
        tmp_ttest_result = [tmp_ttest_result;h_tmp(plt_idx)];
        tmp_p_val = [tmp_p_val;p_tmp(plt_idx)];
    end
    ttest_lib{p_i} = tmp_ttest_result;
    pval_lib{p_i} = tmp_p_val;
    plt_t_lib{p_i} = plt_t;
end
disp('Done')

%% extract mean amplitude
t_noHm_stim = [50, 100; 300, 350; 700, 800];
t_noHm_gip = [-100, 0; 25, 50; 250, 325; 450, 500; 650, 750; 875, 925];
t_noHm_fix = [-100, 0; 300, 375; 450, 500; 650, 750; 875, 925];
t_Hm_stim = [50, 100; 200, 225; 300, 350; 800, 900];
t_Hm_gip = [-200, -100; 90,150; 450, 500; 700, 800];
t_Hm_fix = [-200, -100; 350, 450; 450, 525; 950, 1000];

t_window = {t_noHm_stim, t_Hm_stim, t_noHm_gip, t_Hm_gip, t_noHm_fix, t_Hm_fix};
amp_cir = cell(6,1);
amp_tri = cell(size(amp_cir));

for p_i = 1:6
    t_win = t_window{p_i};
    plt_t = ch_lib{1}{p_i}{1};
    tmp_cir = [];
    tmp_tri = [];
    for t_i = 1:size(t_win,1)
        t_idx = plt_t>=t_win(t_i,1) & plt_t<=t_win(t_i,2);
        tmp_cir = [tmp_cir, cellfun(@(x) mean(x{p_i}{2}(:,t_idx),2), ch_lib, 'uniformoutput',0)'];
        tmp_tri = [tmp_tri, cellfun(@(x) mean(x{p_i}{3}(:,t_idx),2), ch_lib, 'uniformoutput',0)'];
    end
    amp_cir{p_i} = tmp_cir;
    amp_tri{p_i} = tmp_tri;
end
    
% save('D:\Research\meanAmp.mat','t_window','amp_cir','amp_tri','select_ch');

%% visualize each subplot
ax_lib = cell(1,6);
t_name = reshape(repmat({'Stim','GIP','Fix'},2,1),1,[]);
for p_i = 1:6
    ax_lib{p_i} = subplot(3,2,p_i);
    plt_result = ttest_lib{p_i};
%     plt_result = pval_lib{p_i};
%     plt_result(plt_result>0.05) = nan;
    plt_t = plt_t_lib{p_i};
    h = pcolor([plt_result,nan(size(plt_result,1),1);nan(1,size(plt_result,2)+1)]);
    hold on
    xlabel('Time (ms)')
    xline(find(plt_t==0),'k','linewidth',3)
    xline(find(plt_t==test_time_idx(1)),'r--','linewidth',2)
    xline(find(plt_t==test_time_idx(2)),'r--','linewidth',2)
    set(h,'edgecolor','none')
    set(gca,'xtick',1:25:length(plt_t));
    set(gca,'xticklabels',plt_t(1:25:length(plt_t)));
    set(gca,'ytick',1.5:size(plt_result,1)+0.5);
    set(gca,'yticklabels',sub_ch);
%     title(t_name{p_i})
end
% colorbar
% caxis([0 0.05])

%% merge EEG for EEGLAB function
plt_epoch_lib = epoch_lib(:,preserve_idx);
nb_subj = size(plt_epoch_lib,2);
merged_lib = merge_epoch_lib(plt_epoch_lib,select_ch);
disp('Done')
ori_merged_lib = merged_lib;

%% remove trials with latency between stimulus and events in range of 400 ms to 100 ms
thres_lat = [50, 400; 100, 600];
for cond_i = 1:2
    event_time = merged_lib{cond_i}.event_time;
    diff_gip_std = event_time.gipStd_time - event_time.std_time;
    diff_gip_dev = event_time.gipDev_time - event_time.dev_time;
    diff_fix_std = event_time.fixStd_time - event_time.std_time;
    diff_fix_dev = event_time.fixDev_time - event_time.dev_time;
    
    rm_std_idx = diff_gip_std < thres_lat(cond_i,1) |...
                 diff_gip_std > thres_lat(cond_i,2);
%                  diff_fix_std < thres_lat(cond_i,1)
%                  diff_fix_std > thres_lat(cond_i,2);
    merged_lib{cond_i}.std_epoch = pop_rejepoch(merged_lib{cond_i}.std_epoch, rm_std_idx,0);
    merged_lib{cond_i}.gip_std = pop_rejepoch(merged_lib{cond_i}.gip_std, rm_std_idx(~isnan(diff_gip_std)),0);
    merged_lib{cond_i}.fix_std = pop_rejepoch(merged_lib{cond_i}.fix_std, rm_std_idx(~isnan(diff_fix_std)),0);
    rm_dev_idx = diff_gip_dev < thres_lat(cond_i,1) |...
                 diff_gip_dev > thres_lat(cond_i,2);
%                  diff_fix_dev < thres_lat(cond_i,1)
%                  diff_fix_dev > thres_lat(cond_i,2);
    merged_lib{cond_i}.dev_epoch = pop_rejepoch(merged_lib{cond_i}.dev_epoch, rm_dev_idx,0);
    merged_lib{cond_i}.gip_dev = pop_rejepoch(merged_lib{cond_i}.gip_dev, rm_dev_idx(~isnan(diff_gip_dev)),0);
    merged_lib{cond_i}.fix_dev = pop_rejepoch(merged_lib{cond_i}.fix_dev, rm_dev_idx(~isnan(diff_fix_dev)),0);
end

%% topo
eeg_array = [];
cond_i = 2;
tmp_cir = pop_select(merged_lib{cond_i}.std_epoch,'channel',select_ch);
tmp_cir.setname = 'stim circle';
tmp_tri = pop_select(merged_lib{cond_i}.dev_epoch,'channel',select_ch);
tmp_tri.setname = 'stim triangle';
eeg_array = [eeg_array, tmp_cir, tmp_tri];
tmp_cir = pop_select(merged_lib{cond_i}.gip_std,'channel',select_ch);
tmp_cir.setname = 'gip circle';
tmp_tri = pop_select(merged_lib{cond_i}.gip_dev,'channel',select_ch);
tmp_tri.setname = 'gip triangle';
eeg_array = [eeg_array, tmp_cir, tmp_tri];
tmp_cir = pop_select(merged_lib{cond_i}.fix_std,'channel',select_ch);
tmp_cir.setname = 'fix circle';
tmp_tri = pop_select(merged_lib{cond_i}.fix_dev,'channel',select_ch);
tmp_tri.setname = 'fix triangle';
eeg_array = [eeg_array, tmp_cir, tmp_tri];
tmp_grab = pop_select(merged_lib{cond_i}.grab_epoch,'channel',select_ch);
tmp_grab.setname = 'grab';
eeg_array = [eeg_array, tmp_grab];

%% compare ERP
plt_idx = 3;
pop_comperp(eeg_array, 1, plt_idx*2-1,plt_idx*2,'addavg','on','addstd','off','subavg','on','diffavg','on','diffstd','off','alpha',0.05,'tplotopt',{'ydir',1});

%% plot ERP
tarCh = 'Cz';
lock_name = 'stim';
cond_name = 'Hm';
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
% shaded_method = {@(x)(median(x,'omitnan')),@(x)([quantile(x,0.8)-median(x,'omitnan');median(x,'omitnan')-quantile(x,0.2)])};
shaded_method = {@(x)(mean(x,'omitnan')),@(x)(std(x,'omitnan'))};
% shaded_method = {@(x)(mean(x,'omitnan')),@(x)(std(x,'omitnan')/sqrt(size(x,1)))};
fig = plt_erp(plt_EEG_cir,plt_EEG_tri,tarCh,lock_name,shaded_method);



%% plot ERP topo
% further remove trial based on high variance
plt_EEG = eeg_array(5);
var_thres = 0.95;
var_dist = reshape(mean(var(plt_EEG.data,[],2),1),[],1);
plt_EEG = pop_rejepoch(plt_EEG,var_dist> quantile(var_dist,var_thres),0);
pop_topoplot(plt_EEG, 1, [-400:50:600],'Hm fix circle',[5 5] ,0,'electrodes','on');

%% plot difference
plt_idx = 1;
plt_method = @mean;
if plt_idx == 1
    plt_t_idx = [-200:50:800];
else
    plt_t_idx = [-400:50:600];
end
cir_data = eeg_array(plt_idx*2-1).data;
tri_data = eeg_array(plt_idx*2).data;

% further remove trial based on high variance
var_thres = 0.95;
var_dist_cir = reshape(mean(var(cir_data,[],2),1),[],1);
var_dist_tri = reshape(mean(var(tri_data,[],2),1),[],1);
cir_data(:,:,var_dist_cir > quantile(var_dist_cir,var_thres)) = [];
tri_data(:,:,var_dist_tri > quantile(var_dist_tri,var_thres)) = [];

tmp_eeg = eeg_array(plt_idx*2-1); 
tmp_eeg.trials = 1;
tmp_eeg.data = plt_method(cir_data,3) - plt_method(tri_data,3);
pop_topoplot(tmp_eeg,1,plt_t_idx ,'Hm fix diff',[5 5] ,0,'electrodes','on');


%% ERPImage
smooth = 10;
var_thres = 0.95;
clim = [-10 10];
plt_EEG = eeg_array(5);
tarCh = 'Cz';
idx_tarCh = find(ismember({plt_EEG.chanlocs.labels},tarCh));
var_dist = reshape(mean(var(plt_EEG.data(idx_tarCh,:,:),[],2),1),[],1);
rm_trial = var_dist > quantile(var_dist,var_thres);
plt_EEG = pop_rejepoch(plt_EEG,rm_trial,0);
sort_ev = {plt_EEG.event(cellfun(@(x) ~isempty(regexp(x,'Ring 0','once')),{plt_EEG.event.type})).type};
% sort_ev = {plt_EEG.event(cellfun(@(x) ~isempty(regexp(x,'Right Trigger is holded','once')),{plt_EEG.event.type})).type};
% sort_ev = {plt_EEG.event(cellfun(@(x) ~isempty(regexp(x,'circle_gip_start','once')),{plt_EEG.event.type})).type};
figure; pop_erpimage(plt_EEG,1,[idx_tarCh],[[]],tarCh,smooth,1,sort_ev,[],...
    'latency' ,'yerplabel','\muV','erp','on','cbar','on','topo',...
    { [idx_tarCh] plt_EEG.chanlocs plt_EEG.chaninfo },'caxis',clim);

%% ERP comparison
tarCh = 'CPz';
method = @mean;
plt_EEG = eeg_array(5);
plt_t = plt_EEG.times;
idx_tarCh = find(ismember({plt_EEG.chanlocs.labels},tarCh));
var_dist = reshape(mean(var(plt_EEG.data(idx_tarCh,:,:),[],2),1),[],1);
rm_trial = var_dist > quantile(var_dist,var_thres);
plt_EEG = pop_rejepoch(plt_EEG,rm_trial,0);
fix_erp_cir = squeeze(method(plt_EEG.data(idx_tarCh,:,:),3));
plt_EEG = eeg_array(6);
idx_tarCh = find(ismember({plt_EEG.chanlocs.labels},tarCh));
var_dist = reshape(mean(var(plt_EEG.data(idx_tarCh,:,:),[],2),1),[],1);
rm_trial = var_dist > quantile(var_dist,var_thres);
plt_EEG = pop_rejepoch(plt_EEG,rm_trial,0);
fix_erp_tri = squeeze(method(plt_EEG.data(idx_tarCh,:,:),3));
fix_erp_diff = fix_erp_cir-fix_erp_tri;
%
plt_EEG = eeg_array(6);
idx_tarCh = find(ismember({plt_EEG.chanlocs.labels},tarCh));
var_dist = reshape(mean(var(plt_EEG.data(idx_tarCh,:,:),[],2),1),[],1);
rm_trial = var_dist > quantile(var_dist,var_thres);
plt_EEG = pop_rejepoch(plt_EEG,rm_trial,0);
gip_erp_cir = squeeze(method(plt_EEG.data(idx_tarCh,:,:),3));
plt_EEG = eeg_array(4);
idx_tarCh = find(ismember({plt_EEG.chanlocs.labels},tarCh));
var_dist = reshape(mean(var(plt_EEG.data(idx_tarCh,:,:),[],2),1),[],1);
rm_trial = var_dist > quantile(var_dist,var_thres);
plt_EEG = pop_rejepoch(plt_EEG,rm_trial,0);
gip_erp_tri = squeeze(method(plt_EEG.data(idx_tarCh,:,:),3));
gip_erp_diff = gip_erp_cir-gip_erp_tri;
%
figure
plot(plt_t,gip_erp_diff,'b')
hold on
plot(plt_t,fix_erp_diff,'r')
grid on

%% ERP with reaction time


%% ICA on fixation-locked trials
plt_EEG = eeg_array(5);
var_thres = 0.95;
var_dist = reshape(mean(var(plt_EEG.data,[],2),1),[],1);
plt_EEG = pop_rejepoch(plt_EEG,var_dist> quantile(var_dist,var_thres),0);

plt_EEG_ica = pop_runica(plt_EEG, 'icatype','runica','extended',1);
plt_EEG_ica = pop_iclabel(plt_EEG_ica,'default');

%% Check number of trial
trial_name = fieldnames(merged_lib{1}.nb_trial);
nbtrial_lib = {2,6};
nbdir_lib = {2,2};

for cond_i = 1:2
    nb_trial = merged_lib{cond_i}.nb_trial;
    tmp = zeros(6,nb_subj);
    for t_i = 1:length(trial_name)
        tmp(t_i,:) = [nb_trial.(trial_name{t_i})];
    end
    nbtrial_lib{cond_i} = tmp;
    
    nbdir_lib{cond_i,1} =cellfun(@(x) size(x,2), merged_lib{cond_i}.dir_trial.dir_std);
    nbdir_lib{cond_i,2} =cellfun(@(x) size(x,2), merged_lib{cond_i}.dir_trial.dir_dev);
    
end


%%  Cross subjects result with merge_lib
tarCh = 'O1';
var_thres = 0.95;
trial_name = fieldnames(merged_lib{1}.nb_trial);
epoch_name = {'std_epoch','dev_epoch',trial_name{3:end}};
plt_lib = cell(2,length(trial_name));
plt_t_lib = cell(2,length(trial_name));
for cond_i = 1:2
    nb_trial = merged_lib{cond_i}.nb_trial;
    for t_i = 1:length(trial_name)
        loc_t = [0, cumsum(nb_trial.(trial_name{t_i}))];
        tmp_epoch = merged_lib{cond_i}.(epoch_name{t_i});
        idx_tarCh = ismember({tmp_epoch.chanlocs.labels},tarCh);
        plt_t_lib{cond_i,t_i} = tmp_epoch.times;
        tmp_data = zeros(length(loc_t)-1,tmp_epoch.pnts);
        for c_i = 1:length(loc_t)-1
            loc_tmp = tmp_epoch.data(idx_tarCh,:,loc_t(c_i)+1:loc_t(c_i+1));
            var_dist = reshape(mean(var(loc_tmp,[],2),1),[],1);
            rm_trial = var_dist > quantile(var_dist,var_thres);
            tmp_data(c_i,:) = squeeze(mean(loc_tmp(:,:,~rm_trial),3));
        end
        plt_lib{cond_i,t_i} = tmp_data;
    end
end
    
figure
shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(nb_subj))};
% shaded_method = {@(x)(median(x,'omitnan')),@(x)([quantile(x,0.8)-median(x,'omitnan');median(x,'omitnan')-quantile(x,0.2)])};
ax_lib = cell(1,6);
t_name = reshape(repmat({'Stim','GIP','Fix'},2,1),1,[]);
plt_t_lib = reshape(plt_t_lib,[],1);
plt_lib = [plt_lib(1,1:2);plt_lib(2,1:2);plt_lib(1,3:4);plt_lib(2,3:4);plt_lib(1,5:6);plt_lib(2,5:6)];

for p_i = 1:6
    ax_lib{p_i} = subplot(3,2,p_i);
    plt_t = plt_t_lib{p_i*2};
    plt_cir = plt_lib{p_i,1};
    plt_tri = plt_lib{p_i,2};
    ht = shadedErrorBar(plt_t, plt_tri, shaded_method,'lineprops',...
            {'color','b','linewidth',3,'DisplayName','Standard'});
    ht.patch.FaceAlpha = 0.1;
    hold on
    grid on
    hc = shadedErrorBar(plt_t, plt_cir, shaded_method,'lineprops',...
            {'color','r','linewidth',3,'DisplayName','Standard'});
    hc.patch.FaceAlpha = 0.1;
    xlabel('Time (ms)')
    ylabel('Amplitude (\muV)')
    xline(0,'k','linewidth',3)
    set(gcf,'color','w')
    set(gca,'xtick',plt_t(1:25:end));
%     set(gca,'xticklabels',plt_t(1:25:length(plt_t)));
%     set(gca,'ytick',1.5:size(plt_result,1)+0.5);
%     set(gca,'yticklabels',select_ch);
    title(t_name{p_i})
end
linkaxes([ax_lib{:}],'y')

%% subject-wise ERPImage
plt_idx = 5;
is_std = 1;
plt_obj = plt_lib{plt_idx, is_std};
plt_t = plt_t_lib{plt_idx*2};
[nr, nc] = size(plt_obj);
figure
h = pcolor([plt_obj,nan(nr,1);nan(1,nc+1)]);
set(h,'edgecolor','none')
hold on
xline(find(plt_t==0),'k--','linewidth',3)
set(gca,'xtick',1:25:length(plt_t))
set(gca,'xticklabels',plt_t(1:25:end))
set(gca,'ytick',(1:(nc-1))+0.5)
set(gca,'yticklabels',1:(nc-1))
xlabel('Time (ms)')
ylabel('Subject ID')
set(gca,'fontsize',20)
colorbar
colormap('jet')
title('FIX - Hm - DEV')


%% separate into different directions
cond_i = 2;
ev_name = 'fix';
tarCh = 'O2';

% right, up, left, down
up_idx_cir = cell2mat(cellfun(@(x) x(2,:), merged_lib{cond_i}.dir_trial.dir_std,'uniformoutput',0));
right_idx_cir = cell2mat(cellfun(@(x) x(1,:), merged_lib{cond_i}.dir_trial.dir_std,'uniformoutput',0));
left_idx_cir = cell2mat(cellfun(@(x) x(3,:), merged_lib{cond_i}.dir_trial.dir_std,'uniformoutput',0));
down_idx_cir = cell2mat(cellfun(@(x) x(4,:), merged_lib{cond_i}.dir_trial.dir_std,'uniformoutput',0));
up_idx_tri = cell2mat(cellfun(@(x) x(2,:), merged_lib{cond_i}.dir_trial.dir_dev,'uniformoutput',0));
right_idx_tri = cell2mat(cellfun(@(x) x(1,:), merged_lib{cond_i}.dir_trial.dir_dev,'uniformoutput',0));
left_idx_tri = cell2mat(cellfun(@(x) x(3,:), merged_lib{cond_i}.dir_trial.dir_dev,'uniformoutput',0));
down_idx_tri = cell2mat(cellfun(@(x) x(4,:), merged_lib{cond_i}.dir_trial.dir_dev,'uniformoutput',0));

switch ev_name
    case 'stim'
        plt_cir = merged_lib{cond_i}.std_epoch;
        plt_tri = merged_lib{cond_i}.dev_epoch;
        plt_t = merged_lib{cond_i}.std_epoch.times;
        d_r_c = right_idx_cir;
        d_u_c = up_idx_cir;
        d_l_c = left_idx_cir;
        d_d_c = down_idx_cir;
        d_r_t = right_idx_tri;
        d_u_t = up_idx_tri;
        d_l_t = left_idx_tri;
        d_d_t = down_idx_tri;
    case 'gip'
        plt_cir = merged_lib{cond_i}.gip_std;
        plt_tri = merged_lib{cond_i}.gip_dev;
        plt_t = merged_lib{cond_i}.gip_std.times;
        d_r_c = right_idx_cir(~isnan(merged_lib{cond_i}.event_time.gipStd_time));
        d_u_c = up_idx_cir(~isnan(merged_lib{cond_i}.event_time.gipStd_time));
        d_l_c = left_idx_cir(~isnan(merged_lib{cond_i}.event_time.gipStd_time));
        d_d_c = down_idx_cir(~isnan(merged_lib{cond_i}.event_time.gipStd_time));
        d_r_t = right_idx_tri(~isnan(merged_lib{cond_i}.event_time.gipDev_time));
        d_u_t = up_idx_tri(~isnan(merged_lib{cond_i}.event_time.gipDev_time));
        d_l_t = left_idx_tri(~isnan(merged_lib{cond_i}.event_time.gipDev_time));
        d_d_t = down_idx_tri(~isnan(merged_lib{cond_i}.event_time.gipDev_time));
    case 'fix'
        plt_cir = merged_lib{cond_i}.fix_std;
        plt_tri = merged_lib{cond_i}.fix_dev;
        plt_t = merged_lib{cond_i}.fix_std.times;
        d_r_c = right_idx_cir(~isnan(merged_lib{cond_i}.event_time.fixStd_time));
        d_u_c = up_idx_cir(~isnan(merged_lib{cond_i}.event_time.fixStd_time));
        d_l_c = left_idx_cir(~isnan(merged_lib{cond_i}.event_time.fixStd_time));
        d_d_c = down_idx_cir(~isnan(merged_lib{cond_i}.event_time.fixStd_time));
        d_r_t = right_idx_tri(~isnan(merged_lib{cond_i}.event_time.fixDev_time));
        d_u_t = up_idx_tri(~isnan(merged_lib{cond_i}.event_time.fixDev_time));
        d_l_t = left_idx_tri(~isnan(merged_lib{cond_i}.event_time.fixDev_time));
        d_d_t = down_idx_tri(~isnan(merged_lib{cond_i}.event_time.fixDev_time));
end

% gather trials
idx_tarCh = ismember({merged_lib{cond_i}.std_epoch.chanlocs.labels},tarCh);
plt_dir_cir = cell(1,4);
plt_dir_cir{1} = squeeze(plt_cir.data(idx_tarCh,:,d_r_c))';
plt_dir_cir{2} = squeeze(plt_cir.data(idx_tarCh,:,d_u_c))';
plt_dir_cir{3} = squeeze(plt_cir.data(idx_tarCh,:,d_l_c))';
plt_dir_cir{4} = squeeze(plt_cir.data(idx_tarCh,:,d_d_c))';
plt_dir_tri = cell(1,4);
plt_dir_tri{1} = squeeze(plt_tri.data(idx_tarCh,:,d_r_t))';
plt_dir_tri{2} = squeeze(plt_tri.data(idx_tarCh,:,d_u_t))';
plt_dir_tri{3} = squeeze(plt_tri.data(idx_tarCh,:,d_l_t))';
plt_dir_tri{4} = squeeze(plt_tri.data(idx_tarCh,:,d_d_t))';

figure
shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(nb_subj))};
% shaded_method = {@(x)(median(x,'omitnan')),@(x)([quantile(x,0.8)-median(x,'omitnan');median(x,'omitnan')-quantile(x,0.2)])};
ax_lib = cell(1,4);
dir_name = {'right','up','left','down'};
for p_i = 1:4
    ax_lib{p_i} = subplot(2,2,p_i);
    plt_cir = plt_dir_cir{p_i};
    plt_tri = plt_dir_tri{p_i};
    ht = shadedErrorBar(plt_t, plt_tri, shaded_method,'lineprops',...
            {'color','b','linewidth',3,'DisplayName','Standard'});
    ht.patch.FaceAlpha = 0.1;
    hold on
    grid on
    hc = shadedErrorBar(plt_t, plt_cir, shaded_method,'lineprops',...
            {'color','r','linewidth',3,'DisplayName','Standard'});
    hc.patch.FaceAlpha = 0.1;
    xlabel('Time (ms)')
    ylabel('Amplitude (\muV)')
    xline(0,'k','linewidth',3)
    set(gcf,'color','w')
    set(gca,'xtick',plt_t(1:25:end));
%     set(gca,'xticklabels',plt_t(1:25:length(plt_t)));
%     set(gca,'ytick',1.5:size(plt_result,1)+0.5);
%     set(gca,'yticklabels',select_ch);
    title([ev_name,'-',dir_name{p_i}],'fontsize',20)
end
linkaxes([ax_lib{:}],'y')





%% Check Fixation related potential on Oz
tarCh = 'Oz';
smooth = 10;
savepath = '/data/projects/yuan/2021 HM_visual_oddball/figures/20230517/';
for cond_i = 1:2
    switch cond_i
        case 1
            condname = 'noHm';
        case 2
            condname = 'Hm';
    end
    for subj_i = 1:size(plt_epoch_lib,2)
        plt_EEG = plt_epoch_lib{cond_i,subj_i}.fix_epoch;
        idx_tarCh = find(ismember({plt_EEG.chanlocs.labels},tarCh));
        figure; 
        pop_erpimage(plt_EEG,1,idx_tarCh,[[]],tarCh,smooth,1,{},[],'' ,...
                    'yerplabel','\muV','erp','on','cbar','on','topo',...
                    {[idx_tarCh] plt_EEG.chanlocs plt_EEG.chaninfo } );
        saveas(gcf,[savepath, sprintf('FRP_%s_%s_%s.png',tarCh,condname,subj_list{subj_i})]);
        close(gcf)
    end
end

%% alignment
% aligned by FRP
rng('default');
cond_i = 1;
ev_name = 'gip';
tarCh = 'O1';
time_range = [-200,300];
switch cond_i
    case 1
        cond_name = 'noHm';
    case 2
        cond_name = 'Hm';
end
switch ev_name
    case 'stim'
        plt_cir = pop_select(merged_lib{cond_i}.std_epoch,'channel',{tarCh});
        plt_tri = pop_select(merged_lib{cond_i}.dev_epoch,'channel',{tarCh});        
        tname = 'Stim';
    case 'gip'
        plt_cir = pop_select(merged_lib{cond_i}.gip_std,'channel',{tarCh});
        plt_tri = pop_select(merged_lib{cond_i}.gip_dev,'channel',{tarCh});
        tname = 'GIP';
%         peak_polar = true(size(t_interest,1),1);
    case 'fix'
        plt_cir = pop_select(merged_lib{cond_i}.fix_std,'channel',{tarCh});
        plt_tri = pop_select(merged_lib{cond_i}.fix_dev,'channel',{tarCh});
        tname = 'Fix';
end
ref_frp = pop_select(merged_lib{cond_i}.fix_epoch,'channel',{tarCh});
plt_t = plt_cir.times;
plt_idx = plt_t>-500;
plt_cir = squeeze(plt_cir.data)';
plt_tri = squeeze(plt_tri.data)';
ref_frp = squeeze(ref_frp.data)';



% round(size(ref_frp,1)/10)
[align_epoch_cir, amp_lib_cir, lat_lib_cir, idx_rm_epoch_cir, reg_mat_cir] = align_erp(plt_cir, plt_t,...
    't_interest',time_range,'mvavg_winlen',10,'ref_epoch',ref_frp);
[align_epoch_tri, amp_lib_tri, lat_lib_tri, idx_rm_epoch_tri, reg_mat_tri] = align_erp(plt_tri, plt_t,...
    't_interest',time_range,'mvavg_winlen',10,'ref_epoch',ref_frp);
% visualize regressors
% sanity check on regression
[align_epoch_ref, amp_lib_ref, lat_lib_ref, idx_rm_epoch_ref, reg_mat_ref] = align_erp(ref_frp, plt_t,...
    't_interest',time_range,'mvavg_winlen',10,'ref_epoch',ref_frp);
figure
plot(plt_t(plt_idx), reg_mat_cir(plt_idx,:), 'linewidth',3, 'DisplayName','Regressors')
hold on
grid on 
xline(time_range(1),'k--','linewidth',3)
xline(time_range(2),'k--','linewidth',3)
plot(plt_t(plt_idx),mean(align_epoch_ref(:,plt_idx),1),'r','linewidth',3,'DisplayName',...
    sprintf('Aligned Ref (%d/%d)',sum(~isnan(lat_lib_ref)),length(lat_lib_ref)))
plot(plt_t(plt_idx),mean(ref_frp(:,plt_idx),1),'b','linewidth',3,'DisplayName','Avg. Ref')
legend

% visualize aligned ERP
shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan'))};
% shaded_method = {@(x)(median(x,'omitnan')),@(x)([quantile(x,0.8)-median(x,'omitnan');median(x,'omitnan')-quantile(x,0.2)])};
figure
% ht = shadedErrorBar(plt_t, align_epoch_tri, shaded_method,'lineprops',...
%             {'color','b','linewidth',3,'DisplayName','Standard'});
% ht.patch.FaceAlpha = 0.1;
plot(plt_t(plt_idx), mean(align_epoch_tri(:,plt_idx),1),'b','linewidth',3,'DisplayName',...
    sprintf('Standard (%d/%d)',sum(~isnan(lat_lib_tri)),length(lat_lib_tri)))
hold on
grid on
% hc = shadedErrorBar(plt_t, align_epoch_cir, shaded_method,'lineprops',...
%             {'color','r','linewidth',3,'DisplayName','Deviant'});
% hc.patch.FaceAlpha = 0.1;
plot(plt_t(plt_idx), mean(align_epoch_cir(:,plt_idx),1),'r','linewidth',3,'DisplayName',...
    sprintf('Deviant (%d/%d)',sum(~isnan(lat_lib_cir)),length(lat_lib_cir)))
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
xline(0,'k','linewidth',3)
set(gcf,'color','w')
set(gca,'xtick',plt_t(1:25:end));
legend(findobj(gca,'-regexp','DisplayName', '[^'']'));
set(gca,'fontsize',20)
title([ev_name,'-',cond_name],'fontsize',20)
ylim([-5,5])






