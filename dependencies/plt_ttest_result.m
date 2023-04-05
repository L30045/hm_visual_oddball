% plot t-test result for the selected channels
%% plot for paper
filepath = 'D:\Research\';
% load([filepath,'behav_lib_20221022.mat']);
load([filepath,'epoch_lib_rmPreStim_new.mat']);
disp('Done')

%% remove those don't have fixation
[fix_subj_idx, grab_subj_idx] = find_if_device(epoch_lib);
preserve_idx = sum(fix_subj_idx,1)==2 & sum(grab_subj_idx,1)==2;

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
% shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(nb_subj))};
shaded_method = {@(x)(median(x,'omitnan')),@(x)([quantile(x,0.8)-median(x,'omitnan');median(x,'omitnan')-quantile(x,0.2)])};
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
tar_ch = 'CPz';
% shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(nb_subj))};
shaded_method = {@(x)(median(x,'omitnan')),@(x)([quantile(x,0.8)-median(x,'omitnan');median(x,'omitnan')-quantile(x,0.2)])};
tar_ch_idx = ismember(select_ch, tar_ch);
plt_lib = reshape(reshape(ch_lib{tar_ch_idx},3,2)',[],1); % for 3 by 2 subplot
ax_lib = cell(1,6);
t_name = reshape(repmat({'Stim','GIP','Fix'},2,1),1,[]);
for p_i = 1:6
    ax_lib{p_i} = subplot(3,2,p_i);
    plt_t = plt_lib{p_i}{1};
    plt_cir = plt_lib{p_i}{2};
    plt_tri = plt_lib{p_i}{3};
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
merged_lib = merge_epoch_lib(plt_epoch_lib);
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

%% compare ERP
plt_idx = 3;
pop_comperp(eeg_array, 1, plt_idx*2-1,plt_idx*2,'addavg','on','addstd','off','subavg','on','diffavg','on','diffstd','off','alpha',0.05,'tplotopt',{'ydir',1});

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
plt_EEG = eeg_array(4);
var_thres = 0.95;
var_dist = reshape(mean(var(plt_EEG.data,[],2),1),[],1);
plt_EEG = pop_rejepoch(plt_EEG,var_dist> quantile(var_dist,var_thres),0);
pop_topoplot(plt_EEG, 1, [-400:50:600],'Hm fix circle',[5 5] ,0,'electrodes','on');

%% plot difference
plt_idx = 1;
plt_method = @median;
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
smooth = 5;
var_thres = 0.95;
clim = [-10 10];
plt_EEG = eeg_array(5);
tarCh = 'Cz';
idx_tarCh = find(ismember({plt_EEG.chanlocs.labels},tarCh));
var_dist = reshape(mean(var(plt_EEG.data(idx_tarCh,:,:),[],2),1),[],1);
rm_trial = var_dist > quantile(var_dist,var_thres);
plt_EEG = pop_rejepoch(plt_EEG,rm_trial,0);
sort_ev = {plt_EEG.event(cellfun(@(x) ~isempty(regexp(x,'Ring 0','once')),{plt_EEG.event.type})).type};
figure; pop_erpimage(plt_EEG,1,[idx_tarCh],[[]],tarCh,smooth,1,sort_ev,[],...
    'latency' ,'yerplabel','\muV','erp','on','cbar','on','topo',...
    { [idx_tarCh] plt_EEG.chanlocs plt_EEG.chaninfo },'caxis',clim);
% figure; pop_erpimage(EEG,1, [15],[[]],'Cz',smooth,1,{ 'Ring 0; Trial 0; Cube index(L): 1; Frequency: 9; Position: (0.9208, 1.629, 2.345)' 'Ring 0; Trial 0; Cube index(L): 3; Frequency: 11; Position: (0.9091, 1.314, 2.345)' 'Ring 0; Trial 10; Cube index(L): 3; Frequency: 11; Position: (0.8541, 1.151, 2.345)' 'Ring 0; Trial 11; Cube index(L): 0; Frequency: 8; Position: (1.041, 1.444, 2.345)' 'Ring 0; Trial 12; Cube index(L): 1; Frequency: 9; Position: (0.8659, 1.619, 2.345)' 'Ring 0; Trial 13; Cube index(L): 0; Frequency: 8; Position: (1.233, 1.435, 2.345)' 'Ring 0; Trial 13; Cube index(L): 0; Position: (0.9164, 1.414, 2.345)' 'Ring 0; Trial 13; Cube index(L): 2; Frequency: 10; Position: (0.7458, 1.454, 2.345)' 'Ring 0; Trial 15; Cube index(L): 1; Position: (0.7414, 1.589, 2.345)' 'Ring 0; Trial 15; Cube index(L): 2; Frequency: 10; Position: (1.013, 1.483, 2.345)' 'Ring 0; Trial 15; Cube index(L): 3; Frequency: 11; Position: (0.8659, 1.269, 2.345)' 'Ring 0; Trial 16; Cube index(L): 1; Frequency: 9; Position: (0.9091, 1.664, 2.345)' 'Ring 0; Trial 16; Cube index(L): 1; Position: (0.7414, 1.589, 2.345)' 'Ring 0; Trial 16; Cube index(L): 3; Frequency: 11; Position: (0.9208, 1.279, 2.345)' 'Ring 0; Trial 17; Cube index(L): 1; Frequency: 9; Position: (0.9208, 1.629, 2.345)' 'Ring 0; Trial 17; Cube index(L): 3; Position: (0.9253, 1.35, 2.345)' 'Ring 0; Trial 18; Cube index(L): 3; Frequency: 11; Position: (0.9208, 1.279, 2.345)' 'Ring 0; Trial 19; Cube index(L): 2; Frequency: 10; Position: (1.013, 1.483, 2.345)' 'Ring 0; Trial 19; Cube index(L): 3; Frequency: 11; Position: (0.8541, 1.151, 2.345)' 'Ring 0; Trial 19; Cube index(L): 3; Position: (0.9253, 1.35, 2.345)' 'Ring 0; Trial 20; Cube index(L): 1; Frequency: 9; Position: (0.9091, 1.664, 2.345)' 'Ring 0; Trial 20; Cube index(L): 3; Frequency: 11; Position: (0.8659, 1.269, 2.345)' 'Ring 0; Trial 21; Cube index(L): 2; Frequency: 10; Position: (0.7341, 1.489, 2.345)' 'Ring 0; Trial 21; Cube index(L): 3; Position: (0.7414, 1.239, 2.345)' 'Ring 0; Trial 22; Cube index(L): 3; Frequency: 11; Position: (0.9208, 1.279, 2.345)' 'Ring 0; Trial 24; Cube index(L): 2; Frequency: 10; Position: (0.6909, 1.444, 2.345)' 'Ring 0; Trial 24; Cube index(L): 3; Frequency: 11; Position: (1.188, 1.308, 2.345)' 'Ring 0; Trial 25; Cube index(L): 1; Frequency: 9; Position: (0.8659, 1.619, 2.345)' 'Ring 0; Trial 25; Cube index(L): 2; Frequency: 10; Position: (0.6791, 1.326, 2.345)' 'Ring 0; Trial 25; Cube index(L): 3; Position: (0.7414, 1.239, 2.345)' 'Ring 0; Trial 26; Cube index(L): 2; Position: (0.5664, 1.414, 2.345)' 'Ring 0; Trial 27; Cube index(L): 0; Frequency: 8; Position: (1.096, 1.454, 2.345)' 'Ring 0; Trial 27; Cube index(L): 1; Position: (0.7414, 1.589, 2.345)' 'Ring 0; Trial 28; Cube index(L): 3; Frequency: 11; Position: (0.9208, 1.279, 2.345)' 'Ring 0; Trial 2; Cube index(L): 0; Position: (0.9164, 1.414, 2.345)' 'Ring 0; Trial 30; Cube index(L): 3; Frequency: 11; Position: (0.9208, 1.279, 2.345)' 'Ring 0; Trial 31; Cube index(L): 1; Frequency: 9; Position: (1.058, 1.61, 2.345)' 'Ring 0; Trial 31; Cube index(L): 1; Position: (0.9253, 1.7, 2.345)' 'Ring 0; Trial 32; Cube index(L): 2; Position: (0.5664, 1.414, 2.345)' 'Ring 0; Trial 35; Cube index(L): 3; Frequency: 11; Position: (0.8541, 1.151, 2.345)' 'Ring 0; Trial 35; Cube index(L): 3; Frequency: 11; Position: (0.9208, 1.279, 2.345)' 'Ring 0; Trial 36; Cube index(L): 0; Position: (0.9164, 1.414, 2.345)' 'Ring 0; Trial 36; Cube index(L): 1; Frequency: 9; Position: (1.188, 1.658, 2.345)' 'Ring 0; Trial 36; Cube index(L): 3; Frequency: 11; Position: (0.9208, 1.279, 2.345)' 'Ring 0; Trial 37; Cube index(L): 2; Frequency: 10; Position: (0.7458, 1.454, 2.345)' 'Ring 0; Trial 38; Cube index(L): 2; Position: (0.5664, 1.414, 2.345)' 'Ring 0; Trial 39; Cube index(L): 1; Frequency: 9; Position: (0.9208, 1.629, 2.345)' 'Ring 0; Trial 40; Cube index(L): 3; Frequency: 11; Position: (1.188, 1.308, 2.345)' 'Ring 0; Trial 43; Cube index(L): 2; Frequency: 10; Position: (1.013, 1.483, 2.345)' 'Ring 0; Trial 44; Cube index(L): 1; Position: (0.9253, 1.7, 2.345)' 'Ring 0; Trial 44; Cube index(L): 2; Frequency: 10; Position: (0.7458, 1.454, 2.345)' 'Ring 0; Trial 45; Cube index(L): 2; Frequency: 10; Position: (1.013, 1.483, 2.345)' 'Ring 0; Trial 46; Cube index(L): 2; Frequency: 10; Position: (0.883, 1.435, 2.345)' 'Ring 0; Trial 47; Cube index(L): 1; Frequency: 9; Position: (1.058, 1.61, 2.345)' 'Ring 0; Trial 49; Cube index(L): 0; Frequency: 8; Position: (1.363, 1.483, 2.345)' 'Ring 0; Trial 49; Cube index(L): 3; Position: (0.7414, 1.239, 2.345)' 'Ring 0; Trial 49; Cube index(L): 3; Position: (0.9253, 1.35, 2.345)' 'Ring 0; Trial 4; Cube index(L): 1; Frequency: 9; Position: (0.8541, 1.501, 2.345)' 'Ring 0; Trial 4; Cube index(L): 2; Position: (0.5664, 1.414, 2.345)' 'Ring 0; Trial 4; Cube index(L): 3; Frequency: 11; Position: (0.9208, 1.279, 2.345)' 'Ring 0; Trial 50; Cube index(L): 1; Frequency: 9; Position: (0.8659, 1.619, 2.345)' 'Ring 0; Trial 51; Cube index(L): 0; Frequency: 8; Position: (1.029, 1.326, 2.345)' 'Ring 0; Trial 51; Cube index(L): 0; Frequency: 8; Position: (1.096, 1.454, 2.345)' 'Ring 0; Trial 51; Cube index(L): 3; Frequency: 11; Position: (1.058, 1.26, 2.345)' 'Ring 0; Trial 52; Cube index(L): 0; Position: (1.1, 1.525, 2.345)' 'Ring 0; Trial 53; Cube index(L): 0; Position: (1.1, 1.525, 2.345)' 'Ring 0; Trial 53; Cube index(L): 2; Frequency: 10; Position: (0.883, 1.435, 2.345)' 'Ring 0; Trial 53; Cube index(L): 3; Frequency: 11; Position: (0.9091, 1.314, 2.345)' 'Ring 0; Trial 55; Cube index(L): 0; Frequency: 8; Position: (1.029, 1.326, 2.345)' 'Ring 0; Trial 55; Cube index(L): 1; Frequency: 9; Position: (1.188, 1.658, 2.345)' 'Ring 0; Trial 55; Cube index(L): 3; Frequency: 11; Position: (0.9091, 1.314, 2.345)' 'Ring 0; Trial 56; Cube index(L): 2; Frequency: 10; Position: (0.7458, 1.454, 2.345)' 'Ring 0; Trial 57; Cube index(L): 1; Frequency: 9; Position: (0.8541, 1.501, 2.345)' 'Ring 0; Trial 57; Cube index(L): 1; Frequency: 9; Position: (0.9091, 1.664, 2.345)' 'Ring 0; Trial 57; Cube index(L): 2; Position: (0.7503, 1.525, 2.345)' 'Ring 0; Trial 57; Cube index(L): 3; Frequency: 11; Position: (1.188, 1.308, 2.345)' 'Ring 0; Trial 58; Cube index(L): 3; Frequency: 11; Position: (1.058, 1.26, 2.345)' 'Ring 0; Trial 59; Cube index(L): 2; Position: (0.7503, 1.525, 2.345)' 'Ring 0; Trial 5; Cube index(L): 0; Position: (1.1, 1.525, 2.345)' 'Ring 0; Trial 5; Cube index(L): 1; Frequency: 9; Position: (0.8541, 1.501, 2.345)' 'Ring 0; Trial 61; Cube index(L): 2; Frequency: 10; Position: (0.6791, 1.326, 2.345)' 'Ring 0; Trial 62; Cube index(L): 0; Frequency: 8; Position: (1.041, 1.444, 2.345)' 'Ring 0; Trial 63; Cube index(L): 0; Position: (1.1, 1.525, 2.345)' 'Ring 0; Trial 63; Cube index(L): 2; Frequency: 10; Position: (0.6791, 1.326, 2.345)' 'Ring 0; Trial 64; Cube index(L): 1; Frequency: 9; Position: (1.188, 1.658, 2.345)' 'Ring 0; Trial 65; Cube index(L): 2; Frequency: 10; Position: (0.6909, 1.444, 2.345)' 'Ring 0; Trial 65; Cube index(L): 2; Frequency: 10; Position: (1.013, 1.483, 2.345)' 'Ring 0; Trial 66; Cube index(L): 2; Position: (0.5664, 1.414, 2.345)' 'Ring 0; Trial 67; Cube index(L): 0; Position: (0.9164, 1.414, 2.345)' 'Ring 0; Trial 67; Cube index(L): 1; Frequency: 9; Position: (0.8541, 1.501, 2.345)' 'Ring 0; Trial 67; Cube index(L): 3; Frequency: 11; Position: (0.8659, 1.269, 2.345)' 'Ring 0; Trial 68; Cube index(L): 0; Frequency: 8; Position: (1.096, 1.454, 2.345)' 'Ring 0; Trial 68; Cube index(L): 2; Frequency: 10; Position: (1.013, 1.483, 2.345)' 'Ring 0; Trial 69; Cube index(L): 0; Frequency: 8; Position: (1.084, 1.489, 2.345)' 'Ring 0; Trial 69; Cube index(L): 2; Frequency: 10; Position: (0.6909, 1.444, 2.345)' 'Ring 0; Trial 6; Cube index(L): 0; Position: (1.1, 1.525, 2.345)' 'Ring 0; Trial 70; Cube index(L): 0; Frequency: 8; Position: (1.363, 1.483, 2.345)' 'Ring 0; Trial 70; Cube index(L): 1; Frequency: 9; Position: (0.8541, 1.501, 2.345)' 'Ring 0; Trial 70; Cube index(L): 3; Frequency: 11; Position: (0.8659, 1.269, 2.345)' 'Ring 0; Trial 71; Cube index(L): 1; Frequency: 9; Position: (0.8541, 1.501, 2.345)' 'Ring 0; Trial 73; Cube index(L): 0; Frequency: 8; Position: (1.084, 1.489, 2.345)' 'Ring 0; Trial 74; Cube index(L): 2; Frequency: 10; Position: (0.883, 1.435, 2.345)' 'Ring 0; Trial 74; Cube index(L): 3; Frequency: 11; Position: (0.9208, 1.279, 2.345)' 'Ring 0; Trial 75; Cube index(L): 3; Frequency: 11; Position: (0.9208, 1.279, 2.345)' 'Ring 0; Trial 76; Cube index(L): 0; Frequency: 8; Position: (1.084, 1.489, 2.345)' 'Ring 0; Trial 76; Cube index(L): 1; Position: (0.7414, 1.589, 2.345)' 'Ring 0; Trial 76; Cube index(L): 2; Frequency: 10; Position: (0.7458, 1.454, 2.345)' 'Ring 0; Trial 76; Cube index(L): 2; Frequency: 10; Position: (0.883, 1.435, 2.345)' 'Ring 0; Trial 77; Cube index(L): 0; Frequency: 8; Position: (1.363, 1.483, 2.345)' 'Ring 0; Trial 78; Cube index(L): 0; Frequency: 8; Position: (1.041, 1.444, 2.345)' 'Ring 0; Trial 78; Cube index(L): 2; Frequency: 10; Position: (1.013, 1.483, 2.345)' 'Ring 0; Trial 79; Cube index(L): 0; Position: (0.9164, 1.414, 2.345)' 'Ring 0; Trial 79; Cube index(L): 3; Frequency: 11; Position: (0.9091, 1.314, 2.345)' 'Ring 0; Trial 7; Cube index(L): 1; Frequency: 9; Position: (1.058, 1.61, 2.345)' 'Ring 0; Trial 7; Cube index(L): 2; Frequency: 10; Position: (1.013, 1.483, 2.345)' 'Ring 0; Trial 80; Cube index(L): 2; Position: (0.5664, 1.414, 2.345)' 'Ring 0; Trial 82; Cube index(L): 1; Frequency: 9; Position: (0.9208, 1.629, 2.345)' 'Ring 0; Trial 83; Cube index(L): 1; Frequency: 9; Position: (1.058, 1.61, 2.345)' 'Ring 0; Trial 83; Cube index(L): 2; Frequency: 10; Position: (0.6791, 1.326, 2.345)' 'Ring 0; Trial 83; Cube index(L): 3; Position: (0.7414, 1.239, 2.345)' 'Ring 0; Trial 84; Cube index(L): 3; Position: (0.7414, 1.239, 2.345)' 'Ring 0; Trial 85; Cube index(L): 2; Frequency: 10; Position: (0.6909, 1.444, 2.345)' 'Ring 0; Trial 87; Cube index(L): 0; Frequency: 8; Position: (1.096, 1.454, 2.345)' 'Ring 0; Trial 87; Cube index(L): 1; Position: (0.9253, 1.7, 2.345)' 'Ring 0; Trial 87; Cube index(L): 3; Position: (0.7414, 1.239, 2.345)' 'Ring 0; Trial 89; Cube index(L): 3; Frequency: 11; Position: (0.9091, 1.314, 2.345)' 'Ring 0; Trial 89; Cube index(L): 3; Position: (0.9253, 1.35, 2.345)' 'Ring 0; Trial 8; Cube index(L): 1; Frequency: 9; Position: (0.9208, 1.629, 2.345)' 'Ring 0; Trial 8; Cube index(L): 1; Position: (0.9253, 1.7, 2.345)' 'Ring 0; Trial 92; Cube index(L): 1; Frequency: 9; Position: (0.9091, 1.664, 2.345)' 'Ring 0; Trial 93; Cube index(L): 1; Position: (0.9253, 1.7, 2.345)' 'Ring 0; Trial 94; Cube index(L): 2; Frequency: 10; Position: (0.6909, 1.444, 2.345)' 'Ring 0; Trial 95; Cube index(L): 3; Position: (0.7414, 1.239, 2.345)' 'Ring 0; Trial 96; Cube index(L): 0; Frequency: 8; Position: (1.233, 1.435, 2.345)' 'Ring 0; Trial 97; Cube index(L): 3; Frequency: 11; Position: (1.188, 1.308, 2.345)' 'Ring 0; Trial 99; Cube index(L): 2; Frequency: 10; Position: (0.6909, 1.444, 2.345)'},[],'latency' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [15] EEG.chanlocs EEG.chaninfo },'caxis',[-30 30] );
