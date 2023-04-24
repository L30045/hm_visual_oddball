function behav_struct = cal_behav(epoch_struct)
%% calculate head rotation velocity and eye rotation velocity
be_lib = struct('headAng',[],'headRot',[],'headAngDiff',[],...
                'eyeAng',[],'eyeRot',[],'eyeAngDiff',[],...
                'gipAng',[],'gipRot',[],'gipAngDiff',[],...
                'dist',[]);
behav_struct = struct('std_epoch',[],'dev_epoch',[],...
                      'gip_std',[],'gip_dev',[],...
                      'fix_std',[],'fix_dev',[]);
                  
% find target cube location
upLoc = epoch_struct.event_time.upLoc;
downLoc = epoch_struct.event_time.downLoc;
leftLoc = epoch_struct.event_time.leftLoc;
rightLoc = epoch_struct.event_time.rightLoc;
tar_lib = [reshape(upLoc,1,3);...
           reshape(downLoc,1,3);...
           reshape(leftLoc,1,3);...
           reshape(rightLoc,1,3)];
ev_list = {'std_epoch','dev_epoch','gip_std','gip_dev','fix_std','fix_dev'};
ev_direct = {[epoch_struct.event_time.std_up;epoch_struct.event_time.std_down;...
              epoch_struct.event_time.std_left; epoch_struct.event_time.std_right];...
             [epoch_struct.event_time.dev_up;epoch_struct.event_time.dev_down;...
              epoch_struct.event_time.dev_left; epoch_struct.event_time.dev_right]};
dir_idx = [1,2,1,2,1,2];

for e_i = 1:length(ev_list)
    tmp_be = be_lib;
    tar_epoch = epoch_struct.(ev_list{e_i});
    if isempty(tar_epoch)
        break
    end
    idx_headLoc = find(ismember({tar_epoch.chanlocs.labels},{'HeadLoc_x','HeadLoc_y','HeadLoc_z'}));
    idx_headDirect = find(ismember({tar_epoch.chanlocs.labels},{'HeadDirect_x','HeadDirect_y','HeadDirect_z'}));
    idx_gip = find(ismember({tar_epoch.chanlocs.labels},{'GIP_x','GIP_y','GIP_z'}));
    idx_eyeDirect = find(ismember({tar_epoch.chanlocs.labels},{'EyeDirect_x','EyeDirect_y','EyeDirect_z'}));
    tar_direct = ev_direct{dir_idx(e_i)};
    if e_i == 3
        s_t = isnan(epoch_struct.event_time.diff_gip_std);
    elseif e_i == 4
        s_t = isnan(epoch_struct.event_time.diff_gip_dev);
    elseif e_i == 5
        s_t = isnan(epoch_struct.event_time.fixStd_time - epoch_struct.event_time.std_time);
    elseif e_i == 6
        s_t = isnan(epoch_struct.event_time.fixDev_time - epoch_struct.event_time.dev_time);
    else
        s_t = false(1,size(tar_direct,2));
    end
    tar_direct(:,s_t) = [];
    headLoc = tar_epoch.data(idx_headLoc,:,:);
    headDirect = tar_epoch.data(idx_headDirect,:,:);
    GIP = tar_epoch.data(idx_gip,:,:);
    eyeDirect = tar_epoch.data(idx_eyeDirect,:,:);
    % calculate GIP distance
    tmp_be.dist = dist2Box(GIP, headLoc, tar_lib, tar_direct);

    % calculate rotation
    [headAng, headRot, headAngDiff, ~] = cal_rot(headDirect, tar_direct, tar_epoch.srate);
    [eyeAng, eyeRot, eyeAngDiff, ~] = cal_rot(eyeDirect, tar_direct, tar_epoch.srate);
    [gipAng, gipRot, gipAngDiff, ~] = cal_rot(GIP, tar_direct, tar_epoch.srate);

    % save data to lib
    tmp_be.headAng = headAng;
    tmp_be.headRot = headRot;
    tmp_be.headAngDiff = headAngDiff;
    tmp_be.eyeAng = eyeAng;
    tmp_be.eyeRot = eyeRot;
    tmp_be.eyeAngDiff = eyeAngDiff;
    tmp_be.gipAng = gipAng;
    tmp_be.gipRot = gipRot;
    tmp_be.gipAngDiff = gipAngDiff;
    
    behav_struct.(ev_list{e_i}) = tmp_be;
end

end