function vis_gip_loc(epoch_struct, event_time)
%% visualize GIP location
% Input:
%   epoch_struct: EEG structure to plot. GIP locations are required in the
%   channels.
%   event_time: Corresponding event_time structure which contain
%   information of target locations and trial directions

%% Get GIP from epoch_struct
plt_3D = true;
test_flag = true;
len_cube = 0.110433;
plt_t = epoch_struct.times;
gip_xyz_ori = epoch_struct.data(ismember({epoch_struct.chanlocs.labels},{'GIP_x','GIP_y','GIP_z'}),:,:);
% get head location
head_xyz_ori = epoch_struct.data(ismember({epoch_struct.chanlocs.labels},{'HeadLoc_x','HeadLoc_y','HeadLoc_z'}),:,:);

% get target location
if plt_3D
    upLoc = event_time.upLoc;
    downLoc = event_time.downLoc;
    rightLoc = event_time.rightLoc;
    leftLoc = event_time.leftLoc;
    centerLoc = mean([upLoc,downLoc,rightLoc,leftLoc],2);
else
    upLoc = event_time.upLoc(1:2);
    downLoc = event_time.downLoc(1:2);
    rightLoc = event_time.rightLoc(1:2);
    leftLoc = event_time.leftLoc(1:2);
    centerLoc = mean([upLoc,downLoc,rightLoc,leftLoc],2);
end
z_plane = event_time.upLoc(3);

%% visualize all trials
plt_gip = gip_xyz_ori;
plt_idx = squeeze(abs(plt_gip(3,:,:)-z_plane) <= len_cube);
figure
% plot out target location
plt_cube(centerLoc,len_cube,'b');
hold on
grid on
plt_cube(upLoc,len_cube,'g');
plt_cube(downLoc,len_cube,'g');
plt_cube(leftLoc,len_cube,'g');
plt_cube(rightLoc,len_cube,'g');
axis('equal')

% plot3(x_scaled(plt_idx),y_scaled(plt_idx),'bo','linewidth',1);
for i = 1:size(plt_gip,3)
    gip_onset = event_time.diff_gip_std(i);
    if ~isnan(gip_onset)
        [~,plt_onset] = min(abs(plt_t-gip_onset));
        plot3(plt_gip(1,plt_onset,i),plt_gip(2,plt_onset,i),plt_gip(3,plt_onset,i),'rx','linewidth',3,'markersize',15);
        % head location
        plt_cube(head_xyz_ori(:,plt_onset,i),len_cube,'r');
        % draw line between head and GIP
        plot3([plt_gip(1,plt_onset,i),head_xyz_ori(1,plt_onset,i)],...
              [plt_gip(2,plt_onset,i),head_xyz_ori(2,plt_onset,i)],...
              [plt_gip(3,plt_onset,i),head_xyz_ori(3,plt_onset,i)],...
              'linewidth',3,'color','b');
%         pause()
    end
end

%% plot 3D first trial
if test_flag
    plt_gip = gip_xyz_ori;
    plt_idx = squeeze(abs(plt_gip(3,:,:)-z_plane) <= len_cube);
    figure
    % plot out target location
    plt_cube(centerLoc,len_cube,'b');
    hold on
    grid on
    plt_cube(upLoc,len_cube,'g');
    plt_cube(downLoc,len_cube,'g');
    plt_cube(leftLoc,len_cube,'g');
    plt_cube(rightLoc,len_cube,'g');
    axis('equal')
    gip_onset = event_time.diff_gip_std(1);
%     gip_onset = 0;
    [~,plt_onset] = min(abs(plt_t-gip_onset));
    plot3(plt_gip(1,plt_onset,1),plt_gip(2,plt_onset,1),plt_gip(3,plt_onset,1),'rx','linewidth',3,'markersize',15);
    
    pause();
    for i = 1:size(plt_gip,2)           
        % head location
        head_obj = plt_cube(head_xyz_ori(:,i,1),len_cube,'r');
        % plot gip point by point
        plot3(plt_gip(1,i,1),plt_gip(2,i,1),plt_gip(3,i,1),'bo','linewidth',3,'markersize',15);
        pause(0.01);
        delete(head_obj);
    end
    
end


%% plot 2D
% % recenter at zplane
% upLoc = event_time.upLoc(1:2);
% downLoc = event_time.downLoc(1:2);
% rightLoc = event_time.rightLoc(1:2);
% leftLoc = event_time.leftLoc(1:2);
% centerLoc = mean([upLoc,downLoc,rightLoc,leftLoc],2);
% plt_gip = gip_xyz_ori;
% plt_gip(3,:,:) = plt_gip(3,:,:)-z_plane;
% plt_head = head_xyz_ori;
% plt_head(3,:,:) = plt_head(3,:,:)-z_plane;
% % ignore gip that is not on the plane
% plt_idx = squeeze(abs(plt_gip(3,:,:)) <= len_cube);
% 
% % draw line between head location and GIP
% plt_x = cat(1,plt_gip(1,:,:),plt_head(1,:,:));
% plt_y = cat(1,plt_gip(2,:,:),plt_head(2,:,:));
% plt_z = cat(1,plt_gip(3,:,:),plt_head(3,:,:));
% % extract xy value on zplane
% scale = -plt_z(2,:,:)./(sum(abs(plt_z),1));
% x_scaled = squeeze(-diff(plt_x,1,1).*scale+plt_x(2,:,:));
% y_scaled = squeeze(-diff(plt_y,1,1).*scale+plt_y(2,:,:));
% 
% %% visualize all trials
% figure
% % plot out target location
% rectangle('Position',[centerLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
% hold on
% rectangle('Position',[upLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
% rectangle('Position',[downLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
% rectangle('Position',[rightLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
% rectangle('Position',[leftLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
% axis('equal')
% plot(centerLoc(1),centerLoc(2),'ko','linewidth',3,'markersize',15);
% plot(upLoc(1),upLoc(2),'ko','linewidth',3,'markersize',15);
% plot(downLoc(1),downLoc(2),'ko','linewidth',3,'markersize',15);
% plot(rightLoc(1),rightLoc(2),'ko','linewidth',3,'markersize',15);
% plot(leftLoc(1),leftLoc(2),'ko','linewidth',3,'markersize',15);
% 
% plot(x_scaled(plt_idx),y_scaled(plt_idx),'bo','linewidth',1);
% for i = 1:size(x_scaled,2)
%     gip_onset = event_time.diff_gip_std(i);
%     if ~isnan(gip_onset)
%         [~,plt_onset] = min(abs(plt_t-gip_onset));
%         plot(x_scaled(plt_onset,i),y_scaled(plt_onset,i),'rx','linewidth',3,'markersize',15);
%     end
% %     pause()
% end
% 
% 
% 
% 
% 
% %% visualize trials one by one
% figure
% % plot out target location
% rectangle('Position',[centerLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
% hold on
% rectangle('Position',[upLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
% rectangle('Position',[downLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
% rectangle('Position',[rightLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
% rectangle('Position',[leftLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
% axis('equal')
% plot(centerLoc(1),centerLoc(2),'ko','linewidth',3,'markersize',15);
% plot(upLoc(1),upLoc(2),'ko','linewidth',3,'markersize',15);
% plot(downLoc(1),downLoc(2),'ko','linewidth',3,'markersize',15);
% plot(rightLoc(1),rightLoc(2),'ko','linewidth',3,'markersize',15);
% plot(leftLoc(1),leftLoc(2),'ko','linewidth',3,'markersize',15);
% % plot gip on zplane
% for i = 1:size(x_scaled,2)
%     plot(x_scaled(:,i),y_scaled(:,i));
%     gip_onset = event_time.diff_gip_std(i);
%     [~,plt_onset] = min(abs(plt_t-gip_onset));
%     plot(x_scaled(plt_onset,i),y_scaled(plt_onset,i),'rx','linewidth',3,'markersize',15);
%     fprintf('Trial %d\n',i);
%     pause()
% end



end

