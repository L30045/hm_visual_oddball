% plot t-test result for the selected channels
%% plot for paper
filepath = 'D:\Research\';
% load([filepath,'behav_lib_20221022.mat']);
load([filepath,'epoch_lib_rmPreStim_new.mat']);
disp('Done')

%% remove those don't have fixation
[fix_subj_idx, grab_subj_idx] = find_if_device(epoch_lib);
preserve_idx = sum(fix_subj_idx,1)==2 & sum(grab_subj_idx,1)==2;
% remove subject 9
preserve_idx([5,9]) = 0;

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
select_ch = {'Fz','FCz','Cz','CPz','Pz','POz','Oz','O1','O2','P3','P4','CP3','CP4',...
    'C3','C4'};
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
for tarch_i = 1:length(select_ch)
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
select_subj = [1:3,5,6,8:13];
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
test_time_idx = [200 600];
is_ttest = false;
is_correct = true;
sub_ch = {'Cz','CPz','Pz','POz','Oz'};
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
tarCh = 'CPz';
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
plt_EEG = eeg_array(6);
tarCh = 'CPz';
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
tarCh = 'Cz';
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

for cond_i = 1:2
    nb_trial = merged_lib{cond_i}.nb_trial;
    tmp = zeros(6,nb_subj);
    for t_i = 1:length(trial_name)
        tmp(t_i,:) = [nb_trial.(trial_name{t_i})];
    end
    nbtrial_lib{cond_i} = tmp;
end


%%  Cross subjects result with merge_lib
tarCh = 'CPz';
var_thres = 0.95;
trial_name = fieldnames(merged_lib{1}.nb_trial);
epoch_name = fieldnames(merged_lib{1});
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
plt_lib = reshape(plt_lib,[],1);

for p_i = 1:6
    ax_lib{p_i} = subplot(3,2,p_i);
    plt_t = plt_t_lib{p_i*2};
    plt_cir = plt_lib{(p_i-1)*2+1};
    plt_tri = plt_lib{p_i*2};
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


%% separate into different directions






