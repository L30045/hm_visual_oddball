function output = merge_epoch_lib(epoch_lib,selected_ch)
%% Merge epoch_lib to perform ERPImage using EEGLAB function
ev_t = struct('std_time',[],'dev_time',[],'grab_time',[],'gipStd_time',[],'gipDev_time',[],...
              'fixStd_time',[],'fixDev_time',[],'fixBlue_time',[]);
merge_noHm = struct('std_epoch',[],'dev_epoch',[],'grab_epoch',[],...
                'fix_epoch',[],'fix_std',[],'fix_dev',[],'fix_blue',[],...
                'gip_std',[],'gip_dev',[],'event_time',ev_t,...
                'nb_trial',struct('stim_std',[],'stim_dev',[],'gip_std',[],'gip_dev',[],'fix_std',[],'fix_dev',[]),...
                'dir_trial',struct('dir_std',{cell(1,size(epoch_lib,2))},'dir_dev',{cell(1,size(epoch_lib,2))}));
merge_Hm = merge_noHm;
output = {merge_noHm, merge_Hm};
            
tarCh = selected_ch;

for subj_i = 1:size(epoch_lib,2)
    for cond_i = 1:2
        tmp_epoch = epoch_lib{cond_i,subj_i};
        % collect direction
        tmp_dir_std = [tmp_epoch.event_time.std_right;...
                   tmp_epoch.event_time.std_up;...
                   tmp_epoch.event_time.std_left;...
                   tmp_epoch.event_time.std_down];
        tmp_dir_dev = [tmp_epoch.event_time.dev_right;...
                   tmp_epoch.event_time.dev_up;...
                   tmp_epoch.event_time.dev_left;...
                   tmp_epoch.event_time.dev_down];
        std_epoch = pop_select(tmp_epoch.std_epoch,'channel',tarCh);
        dev_epoch = pop_select(tmp_epoch.dev_epoch,'channel',tarCh);
        % check if trial number didn't match
        if size(tmp_dir_std,2)~=std_epoch.trials
            error(fprintf('%d, %d',cond_i, subj_i));
        end
        gip_std = pop_select(tmp_epoch.gip_std,'channel',tarCh);
        gip_dev = pop_select(tmp_epoch.gip_dev,'channel',tarCh);
        output{cond_i}.event_time.std_time = [output{cond_i}.event_time.std_time, tmp_epoch.event_time.std_time];
        output{cond_i}.event_time.dev_time = [output{cond_i}.event_time.dev_time, tmp_epoch.event_time.dev_time];
        output{cond_i}.event_time.gipStd_time = [output{cond_i}.event_time.gipStd_time, tmp_epoch.event_time.gipStd_time];
        output{cond_i}.event_time.gipDev_time = [output{cond_i}.event_time.gipDev_time, tmp_epoch.event_time.gipDev_time];
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
            output{cond_i}.event_time.grab_time = [output{cond_i}.event_time.grab_time, tmp_epoch.event_time.grab_time];
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
            output{cond_i}.event_time.fixStd_time = [output{cond_i}.event_time.fixStd_time, tmp_epoch.event_time.fixStd_time];
            output{cond_i}.event_time.fixDev_time = [output{cond_i}.event_time.fixDev_time, tmp_epoch.event_time.fixDev_time];
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
        % record direction of trials
        output{cond_i}.dir_trial.dir_std{subj_i} = tmp_dir_std;
        output{cond_i}.dir_trial.dir_dev{subj_i} = tmp_dir_dev;
    end
end


end