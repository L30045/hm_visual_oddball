function output = find_trails(bool_start,bool_end,condition)
%FIND_TRIALS Function to find if trails contain a specific event
%   Trails are defined by a start and end index when specific event markers
%   appear in the EEG struct (EEG.events.type). To see if a specific event
%   occcured in these trails, a condition array must be provided. The
%   output is a boolean array of all start indices that had specified
%   event marker(s).

%Defining indices when trail starts and ends
idx_start = find(bool_start);
idx_end = find(bool_end);

%Output of function
output = bool_start;

%Iterating through start and end indices to see if event markers are in
%trail.
for i = 1:length(idx_start)
    % find the first end trial marker after bool_start
    loc_end = idx_end(find(idx_end>idx_start(i),1));
    if sum(condition(idx_start(i):loc_end)) > 0
        output(idx_start(i)) = 1;
    else
        output(idx_start(i)) = 0;
    end
end

