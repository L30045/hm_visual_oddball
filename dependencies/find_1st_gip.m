function [gip_up_idx, gip_down_idx, gip_left_idx, gip_right_idx] = find_1st_gip(EEG,tar_ev_idx)
% this function find the first gip report after stimulus onset
%% find GIP onset time
gip_up= cellfun(@(x) ~isempty(regexp(x,'Up','match','ONCE')),{EEG.event.type});
gip_down= cellfun(@(x) ~isempty(regexp(x,'Bottom','match','ONCE')),{EEG.event.type});
gip_left= cellfun(@(x) ~isempty(regexp(x,'Left','match','ONCE')),{EEG.event.type});
gip_right= cellfun(@(x) ~isempty(regexp(x,'Right','match','ONCE')),{EEG.event.type});
gip_idx = find(gip_up|gip_down|gip_left|gip_right);

% find the first gip after stimulus
% find out stimulus onset time
f_std = find(tar_ev_idx);
gip_1st = false(size(tar_ev_idx));

for i = 1:length(f_std)
    t_f = gip_idx(find(gip_idx>f_std(i),1));
    if ~isempty(t_f)
        if i<length(f_std) && t_f < f_std(i+1)
            gip_1st(t_f) = true;
        end
    end
end

gip_up_idx = gip_1st & gip_up;
gip_down_idx = gip_1st & gip_down;
gip_left_idx = gip_1st & gip_left;
gip_right_idx = gip_1st & gip_right;


end