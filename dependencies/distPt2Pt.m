% distPt2Pt - calculates the distance between two sets of points
%
% Usage:
%   >> distPt2Pt(pts1,pts2); %calculates the euclidan distance between the
%       points in pts1 and the corresponding points in pts2, faster
%       equivelant of diag(pdist2(pts1,pts2);
%
% Inputs:
%   >> pts1, pts2 - Nx2 or Nx3 matricies of N points in 2D/3D Eucldian
%       space. N is the same for both matricies, unless N=1 for only one
%       matrix, in which case it will measure the distance between that
%       point and each of the other points in N
%
% Outputs:
%   >> distances - an Nx1 vector of distances between each corresponding
%       point in pts1 and pts2 (unless N=1 for only 1 input)
%
%   Author: Cory
%   Created: 20170428
function [distances] = distPt2Pt(pts1,pts2)
    if (isvector(pts1) || isvector(pts2)) && any(size(pts1)~=size(pts2)) %compare 1 pt to many
        if size(pts1,1)>size(pts2,1)
            distances = pdist2(pts1,pts2);
        else
            distances = pdist2(pts2,pts1);
        end
    else %compare each pt to each pt
        distances = sqrt(sum((pts1-pts2).^2,2));
    end
    
end