function [up_idx, down_idx, left_idx, right_idx, tar_ev_idx] = find_box_loc(EEG, reg_txt)

%% find up/down left/right event
%boolean arrays where the LEFT/RIGHT or UP/DOWN event markers occured (I assumed
%they were mutually trail exclusive i.e. both do not occur in the same trail)

% gather 4 location
tar_ev = unique({EEG.event(cellfun(@(x) ~isempty(regexp(x,reg_txt,'ONCE')),{EEG.event.type})).type});
tar_ev_idx = cellfun(@(x) ~isempty(regexp(x,reg_txt,'ONCE')),{EEG.event.type});
tar_ev = cellfun(@split ,unique(cellfun(@(x) x(regexp(x,'(')+1:end-1), tar_ev, 'uniformoutput',0)),'uniformoutput',0);
num_loc = zeros(3,4);
for i = 1:4
    tmp = tar_ev{i};
    num_loc(1,i) = str2double(tmp{1}(1:end-1));
    num_loc(2,i) = str2double(tmp{2}(1:end-1));
    num_loc(3,i) = str2double(tmp{3});
end
num_loc = sortrows(num_loc');

upLoc = num_loc(3,:);
downLoc = num_loc(2,:);
leftLoc = num_loc(1,:);
rightLoc = num_loc(4,:);

up_idx = cellfun(@(x) ~isempty(regexp(x,sprintf('(%.4g, %.4g, %.4g)',upLoc),'match','ONCE')),{EEG.event.type});
down_idx = cellfun(@(x) ~isempty(regexp(x,sprintf('(%.4g, %.4g, %.4g)',downLoc),'match','ONCE')),{EEG.event.type});
left_idx = cellfun(@(x) ~isempty(regexp(x,sprintf('(%.4g, %.4g, %.4g)',leftLoc),'match','ONCE')),{EEG.event.type});
right_idx = cellfun(@(x) ~isempty(regexp(x,sprintf('(%.4g, %.4g, %.4g)',rightLoc),'match','ONCE')),{EEG.event.type});

end