function [fig,varargout] = plt_erp_meanXsubj_v2(merge_struct, cond_name, ev_name, tarCh, thres_time, rm_thres, shaded_method)
% this version used merged EEG struct instead of individual epoch struct
nb_subj = size(merge_struct,2);
len_stim = size(merge_struct{1}.std_epoch.data,2);
len_gip = size(merge_struct{1}.gip_std.data,2);
[fix_subj_idx, grab_subj_idx] = find_if_device(merge_struct);

switch ev_name
    case 'stim'
        tri_lib = zeros(nb_subj,len_stim);
        cir_lib = zeros(nb_subj,len_stim);
    case 'gip'
        tri_lib = zeros(nb_subj,len_gip);
        cir_lib = zeros(nb_subj,len_gip);
    case 'fix'
        merge_struct = merge_struct(:,fix_subj_idx(1,:) & fix_subj_idx(2,:));
        nb_subj = size(merge_struct,2);
        tri_lib = zeros(nb_subj,len_gip);
        cir_lib = zeros(nb_subj,len_gip);
    case 'grab'
        merge_struct = merge_struct(:,grab_subj_idx(1,:) & grab_subj_idx(2,:));
        nb_subj = size(merge_struct,2);
        cir_lib = zeros(nb_subj,len_gip);
end

for i = 1:nb_subj
    switch cond_name
        case 'noHm'
            cond_struct = merge_struct{1,i};
        case 'Hm'
            cond_struct = merge_struct{2,i};
    end
    % remove trials by GIP latency
    cond_struct = my_rmEpoch(cond_struct,thres_time);
    [rm_idx_stim_std, rm_idx_gip_std] = my_rmbase(cond_struct.std_epoch, cond_struct.gip_std, cond_struct.event_time.gipStd_time, tarCh, rm_thres);
    [rm_idx_stim_dev, rm_idx_gip_dev] = my_rmbase(cond_struct.dev_epoch, cond_struct.gip_dev, cond_struct.event_time.gipDev_time, tarCh, rm_thres);
    if strcmp(ev_name,'fix')
        [~, rm_idx_fix_std] = my_rmbase(cond_struct.std_epoch, cond_struct.fix_std, cond_struct.event_time.fixStd_time, tarCh, rm_thres);
        [~, rm_idx_fix_dev] = my_rmbase(cond_struct.dev_epoch, cond_struct.fix_dev, cond_struct.event_time.fixDev_time, tarCh, rm_thres);
        cond_struct.fix_std = pop_rejepoch(cond_struct.fix_std,rm_idx_fix_std,0);
        cond_struct.fix_dev = pop_rejepoch(cond_struct.fix_dev,rm_idx_fix_dev,0);
    elseif strcmp(ev_name, 'grab')
        [rm_idx_stim_std, rm_idx_grab] = my_rmbase(cond_struct.std_epoch, cond_struct.grab_epoch, cond_struct.event_time.diff_stim_grab, tarCh, rm_thres);        
        cond_struct.grab_epoch = pop_rejepoch(cond_struct.grab_epoch,rm_idx_grab,0);
    end
    % remove from epoch struct
    cond_struct.std_epoch = pop_rejepoch(cond_struct.std_epoch,rm_idx_stim_std,0);
    cond_struct.dev_epoch = pop_rejepoch(cond_struct.dev_epoch,rm_idx_stim_dev,0);
    cond_struct.gip_std = pop_rejepoch(cond_struct.gip_std,rm_idx_gip_std,0);
    cond_struct.gip_dev = pop_rejepoch(cond_struct.gip_dev,rm_idx_gip_dev,0);
    
    
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
        case 'grab'
            tri_epoch = [];
            cir_epoch = cond_struct.grab_epoch;     
    end
    % remove bad trials
    ch_idx = find(ismember({cir_epoch.chanlocs.labels},tarCh));
    plt_t = cir_epoch.times;
    cir_lib(i,:) = mean(cir_epoch.data(ch_idx,:,:),3);    
    if ~isempty(tri_epoch)
        tri_lib(i,:) = mean(tri_epoch.data(ch_idx,:,:),3);
    end
end

% shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan')/sqrt(nb_subj))};
% shaded_method = {@(x)(mean(x,'omitnan')), @(x)(std(x,'omitnan'))};
fig = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
% >> Triangle
if ~isempty(tri_epoch)
    ht = shadedErrorBar(plt_t, tri_lib, shaded_method,'lineprops',...
        {'color','b','linewidth',3,'DisplayName','Standard'});
    ht.patch.FaceAlpha = 0.1;
end
grid on
hold on
% >> Circle
hc = shadedErrorBar(plt_t, cir_lib, shaded_method,'lineprops',...
    {'color','r','linewidth',3,'DisplayName','Deviant'});
hc.patch.FaceAlpha = 0.1;
xline(0,'k-','DisplayName',sprintf('%s onset',ev_name),'linewidth',3)
legend(findobj(gca,'-regexp','DisplayName', '[^'']'),'location','northwest')
set(gca,'fontsize',30)
set(gcf,'color',[1 1 1])
set(gca,'xtick',round(plt_t(1):100:plt_t(end)))
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
% title(sprintf('%s lock - %s (%s)', ev_name, cond_name, tarCh))

varargout{1} = plt_t;
varargout{2} = cir_lib;
varargout{3} = tri_lib;

end