function merge_dir_lib = separate_direct(epoch_lib)
%% separate direction
epoch_struct = epoch_lib{1};
output = cell(size(epoch_lib));

tmp_ch = cellfun(@(x) {x.std_epoch.chanlocs.labels}, epoch_lib(:),'uniformoutput',0);
[~, tmp_idx] = max(cellfun(@length, tmp_ch));
common_ch = tmp_ch{tmp_idx};
for ch_i = 1:length(epoch_lib(:))
    common_ch = intersect(common_ch, tmp_ch{ch_i});
end

thres_time = [100, 1000]; %ms
tarCh = common_ch;

%% for each epoch structure
for subj_i = 1:size(epoch_lib,2)
    for cond_i = 1:2
        tmp_epoch = epoch_lib{cond_i,subj_i};
        tmp_output = struct('up',epoch_struct,'down',epoch_struct,'left',epoch_struct,'right',epoch_struct);
        std_up_idx = tmp_epoch.event_time.std_up;
        std_down_idx = tmp_epoch.event_time.std_down;
        std_left_idx = tmp_epoch.event_time.std_left;
        std_right_idx = tmp_epoch.event_time.std_right;
        dev_up_idx = tmp_epoch.event_time.dev_up;
        dev_down_idx = tmp_epoch.event_time.dev_down;
        dev_left_idx = tmp_epoch.event_time.dev_left;
        dev_right_idx = tmp_epoch.event_time.dev_right;
        
        %std
        if any(std_up_idx)
            tmp_output.up.std_epoch = pop_rejepoch(tmp_epoch.std_epoch,~std_up_idx,0);
            tmp_output.up.gip_std= pop_rejepoch(tmp_epoch.gip_std,~std_up_idx,0);
        else
            tmp_output.up.std_epoch = [];
        end
        if any(std_down_idx)
            tmp_output.down.std_epoch = pop_rejepoch(tmp_epoch.std_epoch,~std_down_idx,0);
            tmp_output.down.gip_std= pop_rejepoch(tmp_epoch.gip_std,~std_down_idx,0);
        else
            tmp_output.down.std_epoch = [];
        end
        if any(std_left_idx)
            tmp_output.left.std_epoch = pop_rejepoch(tmp_epoch.std_epoch,~std_left_idx,0);
            tmp_output.left.gip_std= pop_rejepoch(tmp_epoch.gip_std,~std_left_idx,0);
        else
            tmp_output.left.std_epoch = [];
        end
        if any(std_right_idx)
            tmp_output.right.std_epoch = pop_rejepoch(tmp_epoch.std_epoch,~std_right_idx,0);
            tmp_output.right.gip_std= pop_rejepoch(tmp_epoch.gip_std,~std_right_idx,0);
        else
            tmp_output.right.std_epoch = [];
        end
        
        % dev
        if any(dev_up_idx)
            tmp_output.up.dev_epoch = pop_rejepoch(tmp_epoch.dev_epoch,~dev_up_idx,0);
            tmp_output.up.gip_dev = pop_rejepoch(tmp_epoch.gip_dev,~dev_up_idx,0);
        else
            tmp_output.up.dev_epoch = [];
        end
        if any(dev_down_idx)
            tmp_output.down.dev_epoch = pop_rejepoch(tmp_epoch.dev_epoch,~dev_down_idx,0);
            tmp_output.down.gip_dev= pop_rejepoch(tmp_epoch.gip_dev,~dev_down_idx,0);
        else
            tmp_output.down.dev_epoch =[];
        end
        if any(dev_left_idx)
            tmp_output.left.dev_epoch = pop_rejepoch(tmp_epoch.dev_epoch,~dev_left_idx,0);
            tmp_output.left.gip_dev= pop_rejepoch(tmp_epoch.gip_dev,~dev_left_idx,0);
        else
            tmp_output.left.dev_epoch =[];
        end
        if any(dev_right_idx)
            tmp_output.right.dev_epoch = pop_rejepoch(tmp_epoch.dev_epoch,~dev_right_idx,0);
            tmp_output.right.gip_dev= pop_rejepoch(tmp_epoch.gip_dev,~dev_right_idx,0);
        else
            tmp_output.right.dev_epoch =[];
        end
        
        % grab
        if ~isempty(tmp_epoch.grab_epoch)
            %std
            if any(std_up_idx)
                tmp_output.up.grab_epoch = pop_rejepoch(tmp_epoch.grab_epoch,~std_up_idx,0);
            end
            if any(std_down_idx)
                tmp_output.down.grab_epoch = pop_rejepoch(tmp_epoch.grab_epoch,~std_down_idx,0);
            end
            if any(std_left_idx)
                tmp_output.left.grab_epoch = pop_rejepoch(tmp_epoch.grab_epoch,~std_left_idx,0);
            end
            if any(std_right_idx)
                tmp_output.right.grab_epoch = pop_rejepoch(tmp_epoch.grab_epoch,~std_right_idx,0);
            end
        end
        % fix
        if ~isempty(tmp_epoch.fix_std)
            %std
            if any(std_up_idx)
                tmp_output.up.fix_std= pop_rejepoch(tmp_epoch.fix_std,~std_up_idx,0);
            end
            if any(std_down_idx)
                tmp_output.down.fix_std= pop_rejepoch(tmp_epoch.fix_std,~std_down_idx,0);
            end
            if any(std_left_idx)
                tmp_output.left.fix_std= pop_rejepoch(tmp_epoch.fix_std,~std_left_idx,0);
            end
            if any(std_right_idx)
                tmp_output.right.fix_std= pop_rejepoch(tmp_epoch.fix_std,~std_right_idx,0);
            end

            % dev
            if any(dev_up_idx)
                tmp_output.up.fix_dev = pop_rejepoch(tmp_epoch.fix_dev,~dev_up_idx,0);
            end
            if any(dev_down_idx)
                tmp_output.down.fix_dev= pop_rejepoch(tmp_epoch.fix_dev,~dev_down_idx,0);
            end
            if any(dev_left_idx)
                tmp_output.left.fix_dev= pop_rejepoch(tmp_epoch.fix_dev,~dev_left_idx,0);
            end
            if any(dev_right_idx)
                tmp_output.right.fix_dev= pop_rejepoch(tmp_epoch.fix_dev,~dev_right_idx,0);
            end
        end
        output{cond_i,subj_i} = tmp_output;
    end
end

%% merge each direction
merge_dir_lib = cell(4,1);
thres_time = [100 1000];
for d_i = 1:4
    switch d_i
        case 1
            d_name = 'up';
        case 2
            d_name = 'down';
        case 3
            d_name = 'left';
        case 4
            d_name = 'right';
    end
    noHm_epoch = cellfun(@(x) x.(d_name), output(1,:),'uniformoutput',0);
    Hm_epoch = cellfun(@(x) x.(d_name), output(2,:),'uniformoutput',0);
    dir_epoch = [noHm_epoch;Hm_epoch];
    
    merge_noHm = struct('std_epoch',[],'dev_epoch',[],'grab_epoch',[],...
                'fix_epoch',[],'fix_std',[],'fix_dev',[],...
                'gip_std',[],'gip_dev',[]);
    merge_Hm = merge_noHm;
    merge_dir = {merge_noHm, merge_Hm};

    for subj_i = 1:size(dir_epoch,2)
        for cond_i = 1:2
%             tmp_epoch = my_rmEpoch(dir_epoch{cond_i,subj_i}, thres_time);
            tmp_epoch = dir_epoch{cond_i,subj_i};
            if ~isempty(tmp_epoch.std_epoch)
                std_epoch = pop_select(tmp_epoch.std_epoch,'channel',tarCh);
                gip_std = pop_select(tmp_epoch.gip_std,'channel',tarCh);
            end
            if ~isempty(tmp_epoch.dev_epoch)
                dev_epoch = pop_select(tmp_epoch.dev_epoch,'channel',tarCh);
                gip_dev = pop_select(tmp_epoch.gip_dev,'channel',tarCh);
            end
            if isempty(merge_dir{cond_i}.std_epoch)
                merge_dir{cond_i}.std_epoch = std_epoch;
                merge_dir{cond_i}.dev_epoch = dev_epoch;
                merge_dir{cond_i}.gip_std = gip_std;
                merge_dir{cond_i}.gip_dev = gip_dev;
            else
                merge_dir{cond_i}.std_epoch = pop_mergeset(std_epoch,merge_dir{cond_i}.std_epoch);
                merge_dir{cond_i}.dev_epoch = pop_mergeset(dev_epoch,merge_dir{cond_i}.dev_epoch);
                merge_dir{cond_i}.gip_std = pop_mergeset(gip_std,merge_dir{cond_i}.gip_std);
                merge_dir{cond_i}.gip_dev = pop_mergeset(gip_dev,merge_dir{cond_i}.gip_dev);
            end
            % grab
            if ~isempty(tmp_epoch.grab_epoch)
                grab_epoch = pop_select(tmp_epoch.grab_epoch,'channel',tarCh);
                if isempty(merge_dir{cond_i}.grab_epoch)
                    merge_dir{cond_i}.grab_epoch = grab_epoch;
                else
                    merge_dir{cond_i}.grab_epoch = pop_mergeset(grab_epoch,merge_dir{cond_i}.grab_epoch);
                end 
            end
            % fix
            if ~isempty(tmp_epoch.fix_std)
                fix_epoch = pop_select(tmp_epoch.fix_epoch,'channel',tarCh);
                fix_std = pop_select(tmp_epoch.fix_std,'channel',tarCh);
                fix_dev = pop_select(tmp_epoch.fix_dev,'channel',tarCh);
                fix_blue = pop_select(tmp_epoch.fix_blue,'channel',tarCh);
                if isempty(merge_dir{cond_i}.fix_epoch)
                    merge_dir{cond_i}.fix_epoch = fix_epoch;
                    merge_dir{cond_i}.fix_std = fix_std;
                    merge_dir{cond_i}.fix_dev = fix_dev;
                    merge_dir{cond_i}.fix_blue = fix_blue;
                else
                    merge_dir{cond_i}.fix_epoch = pop_mergeset(fix_epoch,merge_dir{cond_i}.fix_epoch);
                    merge_dir{cond_i}.fix_std = pop_mergeset(fix_std,merge_dir{cond_i}.fix_std);
                    merge_dir{cond_i}.fix_dev = pop_mergeset(fix_dev,merge_dir{cond_i}.fix_dev);
                    merge_dir{cond_i}.fix_blue = pop_mergeset(fix_blue,merge_dir{cond_i}.fix_blue);
                end
            end
        end
    end
    merge_dir_lib{d_i} = merge_dir;
end







end