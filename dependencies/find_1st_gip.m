function [gip_up_idx, gip_down_idx, gip_left_idx, gip_right_idx,miss_idx] = find_1st_gip(EEG,tar_ev_idx)
% this function find the first gip report after stimulus onset
%% find GIP onset time
gip_up= cellfun(@(x) strcmp(x,'Up'),{EEG.event.type});
gip_down= cellfun(@(x) strcmp(x,'Bottom'),{EEG.event.type});
gip_left= cellfun(@(x) strcmp(x,'Left'),{EEG.event.type});
gip_right= cellfun(@(x) strcmp(x,'Right'),{EEG.event.type});
gip_idx = find(gip_up|gip_down|gip_left|gip_right);

% find the first gip after stimulus
% find out stimulus onset time
f_std = find(tar_ev_idx);
gip_1st = false(size(tar_ev_idx));
miss_idx = false(size(tar_ev_idx));

for i = 1:length(f_std)-1
    t_f = gip_idx(find(gip_idx>f_std(i),1));
    if ~isempty(t_f)
        if t_f < f_std(i+1)
            gip_1st(t_f) = true;
        else
            miss_idx(f_std(i)) = true;
        end
    else
        miss_idx(f_std(i)) = true;
    end
end
% last event
t_f = gip_idx(find(gip_idx>f_std(end),1));
if ~isempty(t_f)
    if (EEG.event(t_f).latency - EEG.event(f_std(end)).latency)/EEG.srate <= 2
        gip_1st(t_f) = true;
    else
        miss_idx(f_std(i)) = true;
    end
else
    miss_idx(f_std(i)) = true;
end

gip_up_idx = gip_1st & gip_up;
gip_down_idx = gip_1st & gip_down;
gip_left_idx = gip_1st & gip_left;
gip_right_idx = gip_1st & gip_right;


end