function [rm_idx_stim, rm_idx_gip] = my_rmbase(stim_epoch, gip_epoch, event_time, tarCh, rm_thres)
%% remove trial with noisy baseline
ch_idx = find(ismember({stim_epoch.chanlocs.labels},tarCh));
stim_epoch = pop_select(stim_epoch,'channel',ch_idx);
gip_epoch = pop_select(gip_epoch,'channel',ch_idx);
base = squeeze(stim_epoch.data(:,1:find(stim_epoch.times==0,1),:));
rm_idx_stim = max(base,[],1)>rm_thres | min(base,[],1)<-rm_thres;

rm_idx_gip = false(1,size(gip_epoch.data,3));
% gip
gip_count = 1;
for t_i = 1:length(rm_idx_stim)
    if ~isnan(event_time(t_i))
        if rm_idx_stim(t_i)
            rm_idx_gip(gip_count) = true;
        end
        gip_count = gip_count+1;
    end
end