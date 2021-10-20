function [interp_p, rm_idx] = linear_interp_idx(input_time, target_time)
% all time scales are in ms
rm_idx = [];
interp_p = false(size(target_time));
for t_i = 1:length(input_time)
    if input_time(t_i) >= 0
        % find the interpretation point
        tmp = find(target_time >= input_time(t_i),1);
        if isempty(tmp)
            rm_idx = [rm_idx, t_i];
        else
            if abs(target_time(tmp-1)-input_time(t_i))...
               < abs(target_time(tmp)-input_time(t_i))
               tmp = tmp-1;
            end
            interp_p(tmp) = true;
        end
    else
        rm_idx = [rm_idx, t_i];
    end
end

end