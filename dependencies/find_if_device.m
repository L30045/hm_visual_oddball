function [fix_subj_idx, grab_subj_idx] = find_if_device(epoch_lib)
fix_subj_idx = true(size(epoch_lib));
grab_subj_idx = true(size(epoch_lib));


for subj_i = 1:size(epoch_lib,2)
    for cond_i = 1:2
        tmp_epoch = epoch_lib{cond_i,subj_i};
        if isempty(tmp_epoch.fix_std)||isempty(tmp_epoch.fix_dev)
            fix_subj_idx(cond_i,subj_i) = false;
        end
        if isempty(tmp_epoch.grab_epoch)
            grab_subj_idx(cond_i,subj_i) = false;
        end
    end
end

end