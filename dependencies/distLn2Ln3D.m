% distLn2Ln3D() - finds the closest points between 2 sets of lines in 3D
%       space
%
% Usage:
%   >> distLn2Ln3D(P0,P1,Q0,Q1); %finds the closest points and distances
%       in 3D space between each of the lines defined by the points in 
%       P0 & P1 and and the corresponding lines in Q0 Q1
%   >> distLn2Ln3D(P0,P1,Q0,Q1,optionP,optionQ); %finds the same parameters
%       but for rays and line segments as specified by option1 and option2
%
% Inputs:
%   >> P0 - an Nx3 array containing N points (of form [x,y,z]) on each 
%       of N lines and with the corresponding points P1 to define lines P
%   >> P1 - a Nx3 array containing N points (of form [x,y,z]) on each of N
%       lines and with the corresponding points P0 define to define lines P
%   >> Q0 - an Nx3 array containing N points (of form [x,y,z]) on each 
%       of N lines and with the corresponding points Q1 to define lines Q,
%       different than lines P.
%   >> Q1 - a Nx3 array containing N points (of form [x,y,z]) on each of N
%       lines and with the corresponding points Q0 define to define lines
%       Q, different than lines P.
%   >> optionP,optionQ - string which specifies whehter P & Q are lines, 
%       rays, or line segments. optionP & optionQ need not be the same.
%       optionP corresponds to P, optionQ corresponds to Q
%       - 'line' - distance is found between lines
%       - 'ray' - distance is found between rays starting at points P0 or
%           Q0, and going through P1 or Q1
%       - 'seg' - distance is found between line segments with end points 
%           P0 and P1 or Q0 and Q1
%
% Outputs:
%   >> dist - a Nx1 vec containing the shortest distances between each of 
%       the 2 paired sets of input lines or rays
%   >> Pc - a Nx3 array of 3D points corresponding to the point on the
%       first line P nearest the corresponding line Q
%   >> Qc - a Nx3 array of 3D points corresponding to the point on the
%       second line Q nearest the corresponding P. The distance
%       between each Pc and Qc is dist.
%
%   Author: Cory
%   Created: 20160920
%   Modified: 20160930 - fixed ray/seg support, changed input vars
%   Modified: 20170606 - replaced diag(m*n') with more efficient dot(m,n,2)
function [dst,Pc,Qc] = distLn2Ln3D(P0,P1,Q0,Q1,option1,option2)
    
    %find relative vectors
    u = P1-P0;
    v = Q1-Q0;
    w0 = P0-Q0;
    
    %find linear equation parameters
    a = dot(u,u,2);
    b = dot(u,v,2);
    c = dot(v,v,2);
    d = dot(u,w0,2);
    e = dot(v,w0,2);    
%     a = diag(u*u');
%     b = diag(u*v');
%     c = diag(v*v');
%     d = diag(u*w0');
%     e = diag(v*w0');
    
    %find closest point in parametric form P(s) and Q(t)
    sc = (b.*e - c.*d)./(a.*c - b.^2);
    tc = (a.*e - b.*d)./(a.*c - b.^2);
    
    %sc2 = sc;
    tc2 = tc;
    %modify sc & tc based on line, ray, or seg
    %P(s) is ray or seg
    if exist('option1','var')
        switch lower(option1)
            case 'line'
                %do nothing different
            case {'ray','rays'}
                
                if any(sc<=0)
                    sc(sc<=0) = 0;
                    tc(sc<=0) = d(sc<=0) ./ b(sc<=0);
                end
            case {'seg','segs','segment','segments'}
                if any(sc<=0)
                    sc(sc<=0) = 0;
                    tc(sc<=0) = d(sc<=0) ./ b(sc<=0);
                end
                if any(sc>=1)
                    sc(sc>=1) = 1;
                    tc(sc>=1) = ( a(sc>=1) + d(sc>=1) )./ b(sc>=1);
                end
            otherwise
                %do nothing different
        end
    end
    
    %Q(t) is ray or seg
    if exist('option2','var')
        switch lower(option2)
            case 'line'
                %do nothing different
            case {'ray','rays'}
                if any(tc<=0)
                    tc(tc<=0) = 0;
                    chng = tc<=0 & tc==tc2;
                    sc(chng) = -d(chng) ./ a(chng);
                end
            case {'seg','segs','segment','segments'}
                if any(tc<=0)
                    tc(tc<=0) = 0;
                    chng = tc<=0 & tc==tc2;
                    sc(chng) = -d(chng) ./ a(chng);
                end
                if any(tc>=1)
                    tc(tc>=1) = 1;
                    chng = tc>=1 & tc==tc2;
                    sc(chng) = ( b(chng) - d(chng) )./ a(chng);
                end
            otherwise
                %do nothing different
        end
    end
    
    sc_M = repmat(sc,1,3);
    tc_M = repmat(tc,1,3);
    
    w = w0 + (sc_M.*u - tc_M.*v);
    near = sqrt(sum(w.^2,2));
    
    %default case is 2 lines
    dst = near;
    Pc = P0 + sc_M.*u;
    Qc = Q0 + tc_M.*v;
    
%     s_nearest = sc;
%     s_dst = zeros(size(dst));
%     t_nearest = tc;
%     t_dst = zeros(size(dst));

end