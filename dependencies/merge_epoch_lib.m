function output = merge_epoch_lib(epoch_lib)
%% Merge epoch_lib to perform ERPImage using EEGLAB function
ev_t = struct('std_time',[],'dev_time',[],'grab_time',[],'gipStd_time',[],'gipDev_time',[],...
              'fixStd_time',[],'fixDev_time',[],'fixBlue_time',[]);
merge_noHm = struct('std_epoch',[],'dev_epoch',[],'grab_epoch',[],...
                'fix_epoch',[],'fix_std',[],'fix_dev',[],'fix_blue',[],...
                'gip_std',[],'gip_dev',[],'event_time',ev_t,...
                'nb_trial',struct('stim_std',[],'stim_dev',[],'gip_std',[],'gip_dev',[],'fix_std',[],'fix_dev',[]));
merge_Hm = merge_noHm;
output = {merge_noHm, merge_Hm};
            
tmp_ch = cellfun(@(x) {x.std_epoch.chanlocs.labels}, epoch_lib(:),'uniformoutput',0);
[~, tmp_idx] = max(cellfun(@length, tmp_ch));
common_ch = tmp_ch{tmp_idx};
for ch_i = 1:length(epoch_lib(:))
    common_ch = intersect(common_ch, tmp_ch{ch_i});
end

thres_time = [100, 1000]; %ms
rm_thres = 30;
tarCh = common_ch;

for subj_i = 1:size(epoch_lib,2)
    for cond_i = 1:2
        tmp_epoch = my_rmEpoch(epoch_lib{cond_i,subj_i}, thres_time);
        std_epoch = pop_select(tmp_epoch.std_epoch,'channel',tarCh);
        dev_epoch = pop_select(tmp_epoch.dev_epoch,'channel',tarCh);
        gip_std = pop_select(tmp_epoch.gip_std,'channel',tarCh);
        gip_dev = pop_select(tmp_epoch.gip_dev,'channel',tarCh);
        [rm_idx_stim_std, rm_idx_gip_std] = my_rmbase(std_epoch, gip_std, tmp_epoch.event_time.gipStd_time, 'Cz', rm_thres);
        [rm_idx_stim_dev, rm_idx_gip_dev] = my_rmbase(dev_epoch, gip_dev, tmp_epoch.event_time.gipDev_time, 'Cz', rm_thres);
        std_epoch = pop_rejepoch(std_epoch,rm_idx_stim_std,0);
        gip_std = pop_rejepoch(gip_std,rm_idx_gip_std,0);
        dev_epoch = pop_rejepoch(dev_epoch,rm_idx_stim_dev,0);
        gip_dev = pop_rejepoch(gip_dev,rm_idx_gip_dev,0);
        output{cond_i}.event_time.std_time = [output{cond_i}.event_time.std_time, tmp_epoch.event_time.std_time(~rm_idx_stim_std)];
        output{cond_i}.event_time.dev_time = [output{cond_i}.event_time.dev_time, tmp_epoch.event_time.dev_time(~rm_idx_stim_dev)];
        output{cond_i}.event_time.gipStd_time = [output{cond_i}.event_time.gipStd_time, tmp_epoch.event_time.gipStd_time(~rm_idx_gip_std)];
        output{cond_i}.event_time.gipDev_time = [output{cond_i}.event_time.gipDev_time, tmp_epoch.event_time.gipDev_time(~rm_idx_gip_dev)];
        if isempty(output{cond_i}.std_epoch)
            output{cond_i}.std_epoch = std_epoch;
            output{cond_i}.dev_epoch = dev_epoch;
            output{cond_i}.gip_std = gip_std;
            output{cond_i}.gip_dev = gip_dev;
        else
            output{cond_i}.std_epoch = pop_mergeset(std_epoch,output{cond_i}.std_epoch);
            output{cond_i}.dev_epoch = pop_mergeset(dev_epoch,output{cond_i}.dev_epoch);
            output{cond_i}.gip_std = pop_mergeset(gip_std,output{cond_i}.gip_std);
            output{cond_i}.gip_dev = pop_mergeset(gip_dev,output{cond_i}.gip_dev);
        end
        % grab
        if ~isempty(tmp_epoch.grab_epoch)
            grab_epoch = pop_select(tmp_epoch.grab_epoch,'channel',tarCh);
            [~, rm_idx_grab] = my_rmbase(std_epoch, grab_epoch, tmp_epoch.event_time.diff_stim_grab, 'Cz', rm_thres);
            grab_epoch = pop_rejepoch(grab_epoch,rm_idx_grab,0);
            output{cond_i}.event_time.grab_time = [output{cond_i}.event_time.grab_time, tmp_epoch.event_time.grab_time(~rm_idx_grab)];
            if isempty(output{cond_i}.grab_epoch)
                output{cond_i}.grab_epoch = grab_epoch;
            else
                output{cond_i}.grab_epoch = pop_mergeset(grab_epoch,output{cond_i}.grab_epoch);
            end 
        else
            output{cond_i}.event_time.grab_time = [output{cond_i}.event_time.grab_time, nan(1,length(tmp_epoch.event_time.std_time))];
        end
        % fix
        if ~isempty(tmp_epoch.fix_std)
            fix_epoch = pop_select(tmp_epoch.fix_epoch,'channel',tarCh);
            fix_std = pop_select(tmp_epoch.fix_std,'channel',tarCh);
            fix_dev = pop_select(tmp_epoch.fix_dev,'channel',tarCh);
            fix_blue = pop_select(tmp_epoch.fix_blue,'channel',tarCh);
            [~, rm_idx_fix_std] = my_rmbase(std_epoch, fix_std, tmp_epoch.event_time.fixStd_time, 'Cz', rm_thres);
            [~, rm_idx_fix_dev] = my_rmbase(dev_epoch, fix_dev, tmp_epoch.event_time.fixDev_time, 'Cz', rm_thres);
            fix_std = pop_rejepoch(fix_std,rm_idx_fix_std,0);
            fix_dev = pop_rejepoch(fix_dev,rm_idx_fix_dev,0);
            output{cond_i}.event_time.fixStd_time = [output{cond_i}.event_time.fixStd_time, tmp_epoch.event_time.fixStd_time(~rm_idx_fix_std)];
            output{cond_i}.event_time.fixDev_time = [output{cond_i}.event_time.fixDev_time, tmp_epoch.event_time.fixDev_time(~rm_idx_fix_dev)];
            output{cond_i}.event_time.fixBlue_time = [output{cond_i}.event_time.fixBlue_time, tmp_epoch.event_time.fixBlue_time];
            if isempty(output{cond_i}.fix_epoch)
                output{cond_i}.fix_epoch = fix_epoch;
                output{cond_i}.fix_std = fix_std;
                output{cond_i}.fix_dev = fix_dev;
                output{cond_i}.fix_blue = fix_blue;
            else
                output{cond_i}.fix_epoch = pop_mergeset(fix_epoch,output{cond_i}.fix_epoch);
                output{cond_i}.fix_std = pop_mergeset(fix_std,output{cond_i}.fix_std);
                output{cond_i}.fix_dev = pop_mergeset(fix_dev,output{cond_i}.fix_dev);
                output{cond_i}.fix_blue = pop_mergeset(fix_blue,output{cond_i}.fix_blue);
            end
        else
            output{cond_i}.event_time.fixStd_time = [output{cond_i}.event_time.fixStd_time, nan(1,length(tmp_epoch.event_time.std_time))];
            output{cond_i}.event_time.fixDev_time = [output{cond_i}.event_time.fixDev_time, nan(1,length(tmp_epoch.event_time.dev_time))];
        end
        % record number of trials
        output{cond_i}.nb_trial.stim_std(subj_i) = size(std_epoch.data,3);
        output{cond_i}.nb_trial.stim_dev(subj_i) = size(dev_epoch.data,3);
        output{cond_i}.nb_trial.gip_std(subj_i) = size(gip_std.data,3);
        output{cond_i}.nb_trial.gip_dev(subj_i) = size(gip_dev.data,3);
        output{cond_i}.nb_trial.fix_std(subj_i) = size(fix_std.data,3);
        output{cond_i}.nb_trial.fix_dev(subj_i) = size(fix_dev.data,3);
    end
end


end