%distPt2Ln() - Finds the shortest distance between points and lines/rays/seg
%
% Usage:
%   >> distPt2Ln(P,Q0,Q1); %finds the shortest Euclidian distance, in R2
%       or R3, between a point (or set of points) and a line (or set of
%       lines), and the closest points on those lines
%   >> distPt2Ln(P,Q0,Q1,MODE); %same as above but the distance is now to
%       rays or line segments
%
% Dependencies:
%   >> multiNorm(seg); %finds unit vectors
%   >> distPt2Pt;
%
% Inputs: 
%   >> P - a set of Mx3 (in R3) or Mx2 (in R2) Euclidan points from which
%       a distance to lines is to found
%   >> Q0 - a set of Nx3 (in R3) or Nx2 (in R2) points which along with the
%       points in Q1 define N lines.
%   >> Q1 - a set of Nx3 (in R3) or Nx2 (in R2) points which along with the
%       points in Q0 define N lines.
%   >> MODE - string entries to switch between lines, rays, segments
%       - 'line' - distance is found between the points and infinite lines,
%           this mode is the DEFAULT
%       - 'ray' - distance is found between the points and rays starting
%           at points Q0 and going through Q1
%       - 'seg' - distance is found between the points and line segments
%           with end points Q0 and Q1
%
% Outputs:
%   >> dst - an Jx1 vector of distances from J points to J lines
%   >> Qc - a Jx3 (in R3) or Jx2 (in R2) set of J points on the lines
%       defined by Q1 and Q2 that are the closest to the points in P;
%       dst = min(||Qc-P||)
%
%
%   NOTE: If there is one point in P (M=1), then the distance from P to
%       every line in Q0,Q1 will be returned (J=N). If there is one line
%       (N=1) and M points then the distances from those points to that
%       line will be calculated (J=M). If the number of points matches the
%       number of lines (M=N) then the pairwise distance between each point
%       and the corresponding line (and the corresponding Qc) will be
%       calculated.
%
%   Author: Cory
%   Created: 20160929
%   Modified: 20170428 - replaced diag(dist(X,Y')) with distPt2Pt(X,Y)

function [dst,Qc] = distPt2Ln(P,Q0,Q1,MODE)
    if size(P,1)~=size(Q0,1)
        if size(P,1)==1
           P = repmat(P,size(Q0,1),1);
        end
        if size(Q0,1)==1
           Q0 = repmat(Q0,size(P,1),1);
           Q1 = repmat(Q1,size(P,1),1);
        end
    end

    N = multiNorm(Q1-Q0);
    proj = dot(Q0-P,N,2);
    Qc = Q0 - repmat(proj,1,size(N,2)).*N;
    dst = distPt2Pt(P,Qc); %dst = diag(dist(P,Qc'));
    
    if exist('MODE','var')
        switch lower(MODE)
            case 'line'
                %do nothing
            case {'ray','rays'}
                dq0 = distPt2Pt(P,Q0); %dq0 = diag(dist(P,Q0'));
                Qc(proj>=0,:) = Q0(proj>=0,:);
                dst(proj>=0) = dq0(proj>=0);
            case {'seg','segs','segment','segments'}
                dq0 = distPt2Pt(P,Q0); %dq0 = diag(dist(P,Q0'));
                Qc(proj>=0,:) = Q0(proj>=0,:);
                dst(proj>=0) = dq0(proj>=0);
                
                sproj = dot(Q1-P,N,2);
                dq1 = distPt2Pt(P,Q1); %dq1 = diag(dist(P,Q1'));
                Qc(sproj<=0,:) = Q1(sproj<=0,:);
                dst(sproj<=0) = dq1(sproj<=0);
            otherwise
                %do nothing
        end
        %do nothing different if there is no mode
    end
end