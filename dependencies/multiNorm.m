% multiNorm() - a function to transform a set of vectors into unit vectors
%
% Usage:
%   >> multiNorm(vec); %finds the unit vectors corresponding to the
%       vectors in vec
%
% Inputs:
%   >> vec - a 3xN or Nx3 set of N vectors in 3-space (x,y,z) OR a 2xN or
%       Nx2 set of N vectors in 2-space (x,y)
%
% Outputs:
%   >> normalized - a Nx3 or Nx2 set of N unit vectors (||v||=1) in the 
%       direction of the corresponding input vectors (and in R-space).
%
%   NOTE: If one of the dimensions of vec is 3 then, multiNorm will treat
%       it as vectors in R3, so 2x3 and 3x2 matricies are 2 vectors of 3
%       components each
%
%   Author: Cory
%   Modified: 20160929 - added 2D possibilities
function normalized = multiNorm(vec)
    dim = find(size(vec)==3);
    z = 3;
    if isempty(dim)
        dim = find(size(vec)==2);
        z = 2;
    end
    if dim~=2
        vec = vec';
    end
    
    norms = sqrt(sum(vec.^2,2));
    normalized = vec./(repmat(norms,1,z));
end