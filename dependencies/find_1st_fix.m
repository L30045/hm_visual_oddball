function fix_1st = find_1st_fix(fix_start,tar_t)
% both input is in the unit of ms.
% this function find the first gip report after stimulus onset
%% find GIP onset time
fix_1st = nan(size(tar_t));

for i = 1:length(tar_t)-1
    t_f = fix_start(find(fix_start>tar_t(i),1));
    if ~isempty(t_f)
        if t_f < tar_t(i+1)
            fix_1st(i) = t_f;
        end
    end
end
% last event
t_f = fix_start(find(fix_start>tar_t(end),1));
if ~isempty(t_f)
    if t_f-tar_t(end)<= 2000
        fix_1st(end) = t_f;
    end
end


end