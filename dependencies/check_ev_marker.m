%% check event markers one by one
filepath = 'D:\Research\oddball_epoch\';
subj_list = {dir([filepath, 'rmPreStim*']).name};
savepath = filepath;

for i = 1:length(subj_list)
    fprintf('Current Subject: %d / %d\n',i, length(subj_list));
    load([filepath,subj_list{i}]);
    % check gip 
    for j = 1:length(epoch_struct_noHm.gip_std.epoch)
        tmp = epoch_struct_noHm.gip_std.epoch(j).eventtype;
        if any(cellfun(@(x) ~isempty(regexp(x,'circle','once')), tmp)) && ...
           any(cellfun(@(x) ~isempty(regexp(x,'triangle','once')), tmp))
            fprintf('Error: Subj %s, noHm, GIP, std, Epoch %d\n', subj_list{i}(end-7:end-4),j);
        end
    end
    for j = 1:length(epoch_struct_noHm.gip_dev.epoch)
        tmp = epoch_struct_noHm.gip_dev.epoch(j).eventtype;
        if any(cellfun(@(x) ~isempty(regexp(x,'circle','once')), tmp)) && ...
           any(cellfun(@(x) ~isempty(regexp(x,'triangle','once')), tmp))
            fprintf('Error: Subj %s, noHm, GIP, dev, Epoch %d\n', subj_list{i}(end-7:end-4),j);
        end
    end
    % check fix
    for j = 1:length(epoch_struct_noHm.fix_std.epoch)
        tmp = epoch_struct_noHm.fix_std.epoch(j).eventtype;
        if any(cellfun(@(x) ~isempty(regexp(x,'circle','once')), tmp)) && ...
           any(cellfun(@(x) ~isempty(regexp(x,'triangle','once')), tmp))
            fprintf('Error: Subj %s, noHm, FIX, std, Epoch %d\n', subj_list{i}(end-7:end-4),j);
        end
    end
    for j = 1:length(epoch_struct_noHm.fix_dev.epoch)
        tmp = epoch_struct_noHm.fix_dev.epoch(j).eventtype;
        if any(cellfun(@(x) ~isempty(regexp(x,'circle','once')), tmp)) && ...
           any(cellfun(@(x) ~isempty(regexp(x,'triangle','once')), tmp))
            fprintf('Error: Subj %s, noHm, FIX, dev, Epoch %d\n', subj_list{i}(end-7:end-4),j);
        end
    end
    % Hm
    % check gip 
    for j = 1:length(epoch_struct_Hm.gip_std.epoch)
        tmp = epoch_struct_Hm.gip_std.epoch(j).eventtype;
        if any(cellfun(@(x) ~isempty(regexp(x,'circle','once')), tmp)) && ...
           any(cellfun(@(x) ~isempty(regexp(x,'triangle','once')), tmp))
            fprintf('Error: Subj %s, Hm, GIP, std, Epoch %d\n', subj_list{i}(end-7:end-4),j);
        end
    end
    for j = 1:length(epoch_struct_Hm.gip_dev.epoch)
        tmp = epoch_struct_Hm.gip_dev.epoch(j).eventtype;
        if any(cellfun(@(x) ~isempty(regexp(x,'circle','once')), tmp)) && ...
           any(cellfun(@(x) ~isempty(regexp(x,'triangle','once')), tmp))
            fprintf('Error: Subj %s, Hm, GIP, dev, Epoch %d\n', subj_list{i}(end-7:end-4),j);
        end
    end
    % check fix
    for j = 1:length(epoch_struct_Hm.fix_std.epoch)
        tmp = epoch_struct_Hm.fix_std.epoch(j).eventtype;
        if any(cellfun(@(x) ~isempty(regexp(x,'circle','once')), tmp)) && ...
           any(cellfun(@(x) ~isempty(regexp(x,'triangle','once')), tmp))
            fprintf('Error: Subj %s, Hm, FIX, std, Epoch %d\n', subj_list{i}(end-7:end-4),j);
        end
    end
    for j = 1:length(epoch_struct_Hm.fix_dev.epoch)
        tmp = epoch_struct_Hm.fix_dev.epoch(j).eventtype;
        if any(cellfun(@(x) ~isempty(regexp(x,'circle','once')), tmp)) && ...
           any(cellfun(@(x) ~isempty(regexp(x,'triangle','once')), tmp))
            fprintf('Error: Subj %s, Hm, FIX, dev, Epoch %d\n', subj_list{i}(end-7:end-4),j);
        end
    end
end