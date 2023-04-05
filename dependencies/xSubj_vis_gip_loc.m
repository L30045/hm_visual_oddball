%% Cross subjects GIP visualization
% load eopch_lib
filepath = 'D:/Research/oddball_epoch/';
subj_list = {dir([filepath, 'rmPreStim*']).name};

epoch_lib = cell(2,length(subj_list));
for i = 1:length(subj_list)
    fprintf('Current Subject: %d / %d\n',i, length(subj_list));
    load([filepath,subj_list{i}]);
    epoch_lib{1,i} = epoch_struct_noHm;
    epoch_lib{2,i} = epoch_struct_Hm;
end
% savepath = filepath;
disp('Done.\n')

%% parameter
cond_i = 2;
len_cube = 0.110433;
center_align = [0;0];
proj_x = zeros(1,epoch_lib{1}.std_epoch.pnts);
proj_y = zeros(size(proj_x));
if_early_lib = [];
gip_onset_lib = [];
fix_onset_lib = [];


for subj_i = 1:size(epoch_lib,2)
    [gip_loc, gip_onset, fix_onset, if_early] = cal_proj_gip_loc(epoch_lib{cond_i,subj_i}.std_epoch,...
                                                                 epoch_lib{cond_i,subj_i}.event_time, center_align);
    proj_x = [proj_x;squeeze(gip_loc(1,:,:))'];
    proj_y = [proj_y;squeeze(gip_loc(2,:,:))'];
    gip_onset_lib = [gip_onset_lib,gip_onset];
    fix_onset_lib = [fix_onset_lib,fix_onset];
    if_early_lib = [if_early_lib,if_early];    
    
end
proj_x(1,:) = [];
proj_y(1,:) = [];
proj_x = proj_x';
proj_y = proj_y';

%% visualization
upLoc = epoch_lib{cond_i,1}.event_time.upLoc(1:2);
downLoc = epoch_lib{cond_i,1}.event_time.downLoc(1:2);
rightLoc = epoch_lib{cond_i,1}.event_time.rightLoc(1:2);
leftLoc = epoch_lib{cond_i,1}.event_time.leftLoc(1:2);
centerLoc = mean([upLoc,downLoc,rightLoc,leftLoc],2);
offset = center_align-centerLoc;
upLoc = upLoc+offset;
downLoc = downLoc+offset;
rightLoc = rightLoc+offset;
leftLoc = leftLoc+offset;
centerLoc = center_align;
val_alpha = 0.5;

figure
% plot out target location
rectangle('Position',[centerLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
hold on
rectangle('Position',[upLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
rectangle('Position',[downLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
rectangle('Position',[rightLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
rectangle('Position',[leftLoc'-len_cube/2, len_cube, len_cube],'linewidth',3);
axis('equal')
plot(centerLoc(1),centerLoc(2),'ko','linewidth',3,'markersize',15);
plot(upLoc(1),upLoc(2),'ko','linewidth',3,'markersize',15);
plot(downLoc(1),downLoc(2),'ko','linewidth',3,'markersize',15);
plot(rightLoc(1),rightLoc(2),'ko','linewidth',3,'markersize',15);
plot(leftLoc(1),leftLoc(2),'ko','linewidth',3,'markersize',15);

% only plot until gip onset
plt_x = proj_x;
plt_y = proj_y;
plt_x_f = nan(size(plt_x));
plt_y_f = nan(size(plt_y));
gip_x = nan(size(gip_onset_lib));
gip_y = nan(size(gip_x));
fix_x = nan(size(gip_x));
fix_y = nan(size(gip_x));
for i = 1:size(plt_x,2)
    if ~isnan(gip_onset_lib(i))
        plt_x(gip_onset_lib(i)+1:end,i) = nan;
        plt_y(gip_onset_lib(i)+1:end,i) = nan;
        if ~isnan(fix_onset_lib(i))
            if ~if_early_lib(i)
                plt_x_f(gip_onset_lib(i):fix_onset_lib(i),i) = proj_x(gip_onset_lib(i):fix_onset_lib(i),i);
                plt_y_f(gip_onset_lib(i):fix_onset_lib(i),i) = proj_y(gip_onset_lib(i):fix_onset_lib(i),i);
            else
                plt_x_f(fix_onset_lib(i):gip_onset_lib(i),i) = proj_x(fix_onset_lib(i):gip_onset_lib(i),i);
                plt_y_f(fix_onset_lib(i):gip_onset_lib(i),i) = proj_y(fix_onset_lib(i):gip_onset_lib(i),i);
            end
            fix_x(i) = proj_x(fix_onset_lib(i),i);
            fix_y(i) = proj_y(fix_onset_lib(i),i);
        end
        gip_x(i) = proj_x(gip_onset_lib(i),i);
        gip_y(i) = proj_y(gip_onset_lib(i),i);
    else
        plt_x(:,i) = nan;
        plt_y(:,i) = nan;
    end
    
end
plot(plt_x,plt_y,'-','linewidth',1,'color',[0,0,1,val_alpha/10]);
% scatter(plt_x,plt_y,'MarkerFaceColor','b','MarkerFaceAlpha',val_alpha,...
%     'MarkerEdgeColor','b','MarkerEdgeAlpha',val_alpha);
% color map
plot(plt_x_f(:,if_early_lib==1),plt_y_f(:,if_early_lib==1),'-','linewidth',1,'color',[1,0,0,val_alpha]);
plot(plt_x_f(:,if_early_lib==0),plt_y_f(:,if_early_lib==0),'-','linewidth',1,'color',[0,1,0,val_alpha]);
% scatter(plt_x_f,plt_y_f,'MarkerFaceColor','g','MarkerFaceAlpha',val_alpha,...
%     'MarkerEdgeColor','g','MarkerEdgeAlpha',val_alpha);
% plot gip onset and fixation onset
scatter(gip_x,gip_y,'MarkerFaceColor','b','MarkerFaceAlpha',val_alpha,...
    'MarkerEdgeColor','b','MarkerEdgeAlpha',val_alpha,'marker','*');
scatter(fix_x,fix_y,'MarkerFaceColor','r','MarkerFaceAlpha',val_alpha,...
    'MarkerEdgeColor','r','MarkerEdgeAlpha',val_alpha,'marker','d');
