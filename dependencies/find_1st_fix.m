function ev_idx = find_1st_fix(fix_t,gip_t, stim_t)
% all inputs are in the unit of ms.
% Input:
%   fix_start: Start time of events to find (N by 1 array). If (N by 2) array is
%   given, treat first column as start time and second column as end time.
%   tar_t: target time points to find the event. If tar_t is nan, skip.
%   [find_method] 1st event after. Find the start time when target time
%   points are covered by the event duration.
% Output:
%   ev_idx


%% Data parsing
t_threshold = 1000; % If the duration between target time points and
% event found is larger than this threshold, discard the event found. (ms)
fix_start = fix_t(:,1);
fix_end = fix_t(:,2);

% intialize output
ev_idx = nan(size(gip_t));

for i = 1:length(gip_t)
    % find the first event onset after target
    f_idx = find(fix_start>gip_t(i),1);
    % find the first event offset after target
    b_idx = find(fix_end>gip_t(i),1);
    % check if event onset exist
    if ~isempty(f_idx) || ~isempty(b_idx)
        % check if event offset happens earlier than event onset
        if isempty(f_idx) || b_idx <= f_idx
            % replace f_idx by b_idx
            f_idx = b_idx;
        end
        % check if the duration exceed threshold.
        if abs(fix_start(f_idx)-gip_t(i))<=t_threshold
            % check if fix happens earlier than stim
            if fix_start(f_idx) > stim_t(i)
                ev_idx(i) = fix_start(f_idx);
            end
        end
    end
end



end