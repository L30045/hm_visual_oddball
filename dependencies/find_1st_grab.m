function grab_1st = find_1st_grab(EEG,grab_start,tar_t)
% both input is in the unit of index.
% this function find the first gip report after stimulus onset
%% find GIP onset time
grab_1st = nan(size(tar_t));

for i = 1:length(tar_t)-1
    t_f = grab_start(find(grab_start>tar_t(i),1));
    if ~isempty(t_f)
        if t_f < tar_t(i+1)
            grab_1st(i) = t_f;
        end
    end
end
% last event
t_f = grab_start(find(grab_start>tar_t(end),1));
if ~isempty(t_f)
    if (EEG.event(t_f).latency - EEG.event(tar_t(end)).latency)/EEG.srate <= 2
        grab_1st(end) = t_f;
    end
end


end