function output = merge_epoch_lib(epoch_lib)
%% Merge epoch_lib to perform ERPImage using EEGLAB function
ev_t = struct('std_time',[],'dev_time',[],'grab_time',[],'gipStd_time',[],'gipDev_time',[],...
              'fixStd_time',[],'fixDev_time',[],'fixBlue_time',[]);
merge_noHm = struct('std_epoch',[],'dev_epoch',[],'grab_epoch',[],...
                'fix_epoch',[],'fix_std',[],'fix_dev',[],'fix_blue',[],...
                'gip_std',[],'gip_dev',[],'event_time',ev_t);
merge_Hm = merge_noHm;
output = {merge_noHm, merge_Hm};
            
tmp_ch = cellfun(@(x) {x.std_epoch.chanlocs.labels}, epoch_lib(:),'uniformoutput',0);
[~, tmp_idx] = max(cellfun(@length, tmp_ch));
common_ch = tmp_ch{tmp_idx};
for ch_i = 1:length(epoch_lib(:))
    common_ch = intersect(common_ch, tmp_ch{ch_i});
end

thres_time = [100, 1000]; %ms
tarCh = common_ch;

for subj_i = 1:size(epoch_lib,2)
    for cond_i = 1:2
        tmp_epoch = my_rmEpoch(epoch_lib{cond_i,subj_i}, thres_time);
        output{cond_i}.event_time.std_time = [output{cond_i}.event_time.std_time, tmp_epoch.event_time.std_time];
        output{cond_i}.event_time.dev_time = [output{cond_i}.event_time.dev_time, tmp_epoch.event_time.dev_time];
        output{cond_i}.event_time.gipStd_time = [output{cond_i}.event_time.gipStd_time, tmp_epoch.event_time.gipStd_time];
        output{cond_i}.event_time.gipDev_time = [output{cond_i}.event_time.gipDev_time, tmp_epoch.event_time.gipDev_time];
        std_epoch = pop_select(tmp_epoch.std_epoch,'channel',tarCh);
        dev_epoch = pop_select(tmp_epoch.dev_epoch,'channel',tarCh);
        gip_std = pop_select(tmp_epoch.gip_std,'channel',tarCh);
        gip_dev = pop_select(tmp_epoch.gip_dev,'channel',tarCh);
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
            output{cond_i}.event_time.grab_time = [output{cond_i}.event_time.grab_time, tmp_epoch.event_time.grab_time];
            grab_epoch = pop_select(tmp_epoch.grab_epoch,'channel',tarCh);
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
            output{cond_i}.event_time.fixStd_time = [output{cond_i}.event_time.fixStd_time, tmp_epoch.event_time.fixStd_time];
            output{cond_i}.event_time.fixDev_time = [output{cond_i}.event_time.fixDev_time, tmp_epoch.event_time.fixDev_time];
            output{cond_i}.event_time.fixBlue_time = [output{cond_i}.event_time.fixBlue_time, tmp_epoch.event_time.fixBlue_time];
            fix_epoch = pop_select(tmp_epoch.fix_epoch,'channel',tarCh);
            fix_std = pop_select(tmp_epoch.fix_std,'channel',tarCh);
            fix_dev = pop_select(tmp_epoch.fix_dev,'channel',tarCh);
            fix_blue = pop_select(tmp_epoch.fix_blue,'channel',tarCh);
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
    end
end









end