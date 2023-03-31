function [gip_loc, gip_onset, fix_onset, if_early] = cal_proj_gip_loc(epoch_struct, event_time, center_align)
%% visualize GIP location
gip_loc = zeros(2,epoch_struct.pnts,epoch_struct.trials);
gip_onset = nan(1,epoch_struct.trials);
fix_onset = nan(size(gip_onset));
if_early = nan(size(gip_onset));

%% Get GIP from epoch_struct
plt_t = epoch_struct.times;
gip_xyz_ori = epoch_struct.data(ismember({epoch_struct.chanlocs.labels},{'GIP_x','GIP_y','GIP_z'}),:,:);
% get head location
head_xyz_ori = epoch_struct.data(ismember({epoch_struct.chanlocs.labels},{'HeadLoc_x','HeadLoc_y','HeadLoc_z'}),:,:);

% get target location
upLoc = event_time.upLoc(1:2);
downLoc = event_time.downLoc(1:2);
rightLoc = event_time.rightLoc(1:2);
leftLoc = event_time.leftLoc(1:2);
centerLoc = mean([upLoc,downLoc,rightLoc,leftLoc],2);
z_plane = event_time.upLoc(3);

% align center
offset = center_align-centerLoc;
% centerLoc = centerLoc+offset;
% upLoc = upLoc+offset;
% downLoc = downLoc+offset;
% rightLoc = rightLoc+offset;
% leftLoc = leftLoc+offset;
gip_xyz_ori(1:2,:,:) = gip_xyz_ori(1:2,:,:)+offset;
head_xyz_ori(1:2,:,:) = head_xyz_ori(1:2,:,:)+offset;

%% plot 2D
% recenter at zplane
plt_gip = gip_xyz_ori;
plt_gip(3,:,:) = plt_gip(3,:,:)-z_plane;
plt_head = head_xyz_ori;
plt_head(3,:,:) = plt_head(3,:,:)-z_plane;
% ignore gip that is not on the plane
% plt_idx = squeeze(abs(plt_gip(3,:,:)) <= len_cube);

% draw line between head location and GIP
plt_x = cat(1,plt_gip(1,:,:),plt_head(1,:,:));
plt_y = cat(1,plt_gip(2,:,:),plt_head(2,:,:));
plt_z = cat(1,plt_gip(3,:,:),plt_head(3,:,:));
% extract xy value on zplane
scale = -plt_z(2,:,:)./(sum(abs(plt_z),1));
x_scaled = squeeze(-diff(plt_x,1,1).*scale+plt_x(2,:,:));
y_scaled = squeeze(-diff(plt_y,1,1).*scale+plt_y(2,:,:));
gip_loc(1,:,:) = x_scaled;
gip_loc(2,:,:) = y_scaled;

%% visualize all trials
diff_fix_std = event_time.fixStd_time - event_time.std_time;
for i = 1:size(x_scaled,2)
    loc_gip_onset = event_time.diff_gip_std(i);
    loc_fix_onset = diff_fix_std(i);
    if ~isnan(loc_gip_onset)
        [~,plt_onset_gip] = min(abs(plt_t-loc_gip_onset));
        gip_onset(i) = plt_onset_gip;
        if ~isnan(loc_fix_onset)
            [~,plt_onset_fix] = min(abs(plt_t-loc_fix_onset));
            fix_onset(i) = plt_onset_fix;
            if loc_fix_onset <= loc_gip_onset
                if_early(i) = true;
            else
                if_early(i) = false;
            end
        end
    end
end

end

