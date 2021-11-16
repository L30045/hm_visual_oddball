function [interp_p, rm_idx] = linear_interp_idx(input_time, target_time)
%% This function labels the corresponding time points of input_time on the scale of target_time.
% Input:
%       input_time:  time scale to be interpretted.
%       target_time: time scale to be projected on.
% Output
%       interp_p: labels of time points to interpret. (on the scale of
%       target_time)
%       rm_idx: index of the input_time time points which cannot be
%       interpretted.
% e.g.
%   input_time = [-1, 1.1, 3.5, 7.9, 10.5];
%   target_time = [1,2,3,4,5,6,7,8,9,10];
%   interp_p =    [1,0,0,1,0,0,0,1,0,0];
%   rm_idx = [1, 5];
%
% Note: all time scales are in ms

%%
% record input_time time points which cannot be interpretted.
rm_idx = zeros(size(input_time));
rm_i = 1;
% label time points to interpret on target_time
interp_p = false(size(target_time));
tmp = 1;
% loop over input_time to find the corresponding time points on target_time
for t_i = 1:length(input_time)
    % ignore input_time before target_time start to record
    if input_time(t_i) >= target_time(1)
        % find the interpretation point
        tmp = find(target_time(tmp+1:end) >= input_time(t_i),1)+tmp;
        if isempty(tmp)
            rm_idx(rm_i) = t_i;
            rm_i = rm_i+1;
        else
            % compare current time point and previous time point on
            % target_time and choose the closer time point to interpret.
            if abs(target_time(tmp-1)-input_time(t_i))...
               < abs(target_time(tmp)-input_time(t_i))
               tmp = tmp-1;
            end
            interp_p(tmp) = true;
        end
    else
        rm_idx(rm_i) = t_i;
        rm_i = rm_i+1;
    end
end
% remove unused rm_idx space
rm_idx(rm_i:end) = [];

end