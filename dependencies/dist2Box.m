%   This function calculates the distance between target box and GIP vector
%   Inputs: GIP     -   in euclidean points
%           headLoc -   in euclidean points
%           tar_lib -   possible locations for the target box in euclidean points.
%                       rows 1 to 4 are for positions up/down/left/right
%           dir_lib -   trialwise location of the target box
%                       columns = trial index
%                       rows 1 to 4 indicate location (up/down/left/right)
%   Output: dist    -   distance in sample points x trials

function dist = dist2Box(GIP, headLoc, tar_lib, dir_lib)

    % find the correct target location for each trial
    for r = 1:4
        for i = find(dir_lib(r,:)==1)
           tar_loc(i,:) =  tar_lib(r,:);
        end
    end

    % get distance by cory's function distPt2Line - loop over trials
    for i = 1:size(GIP,3)
        P = tar_loc(i,:);
        Q0 = headLoc(:,:,i)';
        Q1 = GIP(:,:,i)';
        dist(:,i) = distPt2Ln(P,Q0,Q1,'ray');
    end

end
