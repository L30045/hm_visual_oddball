function ev_idx = find_near_ev(ev_start_t,tar_t)
% both input is in the unit of ms.
% Input:
%   ev_start: Start time of events to find (N by 1 array). If (N by 2) array is
%   given, treat first column as start time and second column as end time.
%   tar_t: target time points to find the event. If tar_t is nan, skip.
%   [find_method] 1st event after. Find the start time when target time
%   points are covered by the event duration.
% Output:
%   ev_idx


%% Data parsing
t_threshold = 1000; % If the duration between target time points and
% event found is larger than this threshold, discard the event found. (ms)
if size(ev_start_t,1)==1 || size(ev_start_t,2)==1
    ev_start_t = reshape(ev_start_t,[],1);
    ev_end_t = ev_start_t;
else
    ev_end_t = ev_start_t(:,2);
    ev_start_t = ev_start_t(:,1);
end
% intialize output
ev_idx = nan(size(tar_t));

for i = 1:length(tar_t)
    % find the first event onset after target
    f_idx = find(ev_start_t>tar_t(i),1);
    % find the first event offset after target
    b_idx = find(ev_end_t>tar_t(i),1);
    % check if event onset exist
    if ~isempty(f_idx)
        % check if event offset happens earlier than event onset
        if b_idx <= f_idx
            % replace f_idx by b_idx
            f_idx = b_idx;
        end
        % check if the duration exceed threshold.
        if abs(ev_start_t(f_idx)-tar_t(i))<=t_threshold
            ev_idx(i) = ev_start_t(f_idx);
        end
    end
end

end