function [dst,Rc,Bc] = distRay2Border3D(R0,R1,B0,B1)
    %Rays: N rays (R0:Nx3, R1:Nx3), Borders: M segs (B0:Mx3,B1:Mx3) M=12 for cube
    
    %Input Error Trap
    if size(R0,1)~=size(R1,1) || size(B0,1)~=size(B1,1)
        error('Missmatch of ray or vector endpoints')
    end
    
    %initialize nRxmB array of possible dist of rays to each segs
    dist_all = zeros(size(R0,1),size(B0,1));
    %initialize nRx3xmB arrays of possible closest points on rays and segs
    Rc_all = zeros(size(R0,1),3,size(B0,1)); %closest points on ray
    Bc_all = zeros(size(R0,1),3,size(B0,1)); %closest points on box
    
    %loop through each edge (i in M), comparing rays
    for i = 1:size(B0,1)
        %call distRay2Seg3D
        %store in possible arrays
        [dist_all(:,i),Rc_all(:,:,i),Bc_all(:,:,i)] = distRay2Seg3D(R0,R1,repmat(B0(i,:),size(R0,1),1),repmat(B1(i,:),size(R0,1),1));
%         dist_all(:,i) = tmp_dst;
%         Rc_all(:,:,i) = tmp_Rc;
%         Bc_all(:,:,i) = tmp_Bc;
    end
    %find minimum dist of each segments
    [dst,Bi_min] = min(dist_all,[],2);
    %collapse possible dist/closest points to shortest distance/pts
    Rc = Rc_all(:,:,Bi_min);
    Bc = Bc_all(:,:,Bi_min);
end

function [dst,Pc,Qc] = distRay2Seg3D(P0,P1,Q0,Q1) %subset implementation of distLn2Ln3D using only ray and seg input
    u = P1-P0;
    v = Q1-Q0;
    w0 = P0-Q0;
    
    %find linear equation parameters
    a = dot(u,u,2);
    b = dot(u,v,2);
    c = dot(v,v,2);
    d = dot(u,w0,2);
    e = dot(v,w0,2);    
    
    %find closest point in parametric form P(s) and Q(t)
    sc = (b.*e - c.*d)./(a.*c - b.^2);
    tc = (a.*e - b.*d)./(a.*c - b.^2);
    
    %case {'ray','rays'}
    if any(sc<=0)
        sc(sc<=0) = 0;
        %tc(sc<=0) = d(sc<=0) ./ b(sc<=0); %failed version
        tc(sc<=0) = e(sc<=0) ./ c(sc<=0);
        %s_lock_low = s_lock_low | sc<=0;
        %fprintf('MODE: 1 (Sc<0)\n');
    end
    
    %case {'seg','segs','segment','segments'}
    if any(tc<=0)
        tc(tc<=0) = 0;
        chng = tc<=0; %& ~s_lock;
        sc(chng) = -d(chng) ./ a(chng);
    end
    if any(tc>=1)
        tc(tc>=1) = 1;
        chng = tc>=1; %& ~s_lock;
        sc(chng) = ( b(chng) - d(chng) )./ a(chng);
    end
    
    sc_M = repmat(sc,1,3);
    tc_M = repmat(tc,1,3);
    
    w = w0 + (sc_M.*u - tc_M.*v);
    near = sqrt(sum(w.^2,2));
    
    %default case is 2 lines
    dst = near;
    Pc = P0 + sc_M.*u;
    Qc = Q0 + tc_M.*v;

end