function [ang, v_ang, ang_diff, ang_cumsum] = cal_rot(data, tar_direct, srate)
% project 3D data onto 2D target direction and estimate the rotation
% velocity. The positive values mean toward the target box.
% WARNING: this is just an estimation
velocity_smooth_win_len = 40; % msec
ang = nan(size(data,2),size(data,3)); % time by trials
v_ang = nan(size(ang)); % time by trials
ang_diff = nan(size(ang)); % time by trials
ang_cumsum = nan(size(ang)); % time by trials
% reformat data
for d_i = 1:size(tar_direct,2)
    switch find(tar_direct(:,d_i))
        case 1 %up
            loc_d = data([2,3],:,d_i);
            pos_direct = 1;
        case 2 %down
            loc_d = data([2,3],:,d_i);
            pos_direct = -1;
        case 3 %left
            loc_d = data([1,3],:,d_i);
            pos_direct = -1;
        case 4 %right
            loc_d = data([1,3],:,d_i);
            pos_direct = 1;        
    end
    int_p = diff([0 loc_d(1,:)],1,2)~=0;
    loc_d = loc_d(:,int_p);
    loc_srate = 1/(mean(diff(find(int_p)))/srate);
    if ~isnan(loc_srate)
        [loc_angdiff, loc_v_ang, loc_ang, pos_d] = smooth_cal_ang(loc_d,loc_srate,velocity_smooth_win_len);
        % interp
        loc_angcumsum = cumsum(loc_angdiff,'omitnan');
        tmp = const_interp_data([loc_angdiff;loc_v_ang;loc_ang;loc_angcumsum],int_p,[]);
        ang_diff(:,d_i) = pos_direct*tmp(1,:);
        v_ang(:,d_i) = pos_direct*tmp(2,:);
        ang(:,d_i) = pos_direct*tmp(3,:);
        ang_cumsum(:,d_i) = pos_direct*tmp(4,:);
    end
end


end


function [angdiff, v_ang, ang, pos_d] = smooth_cal_ang(loc_d, srate, velocity_smooth_win_len)
angdiff = nan(1,size(loc_d,2));
v_ang = nan(1,size(loc_d,2));
ang = nan(1,size(loc_d,2));
pos_d = nan(1,size(loc_d,2));
vel_win_len = round(0.001*velocity_smooth_win_len*srate/2);
for v_i = vel_win_len+1:size(loc_d,2)-vel_win_len
    if ~isnan(loc_d(1,v_i+[-vel_win_len,vel_win_len]))
        nomi = double(loc_d(:,v_i-vel_win_len)'*loc_d(:,v_i+vel_win_len));
        denomi = double((norm(loc_d(:,v_i-vel_win_len))*norm(loc_d(:,v_i+vel_win_len))));
        tolerance = 1 - cos(pi/180*0.1); % adding a 0.1 degree tolerance when calculating acos
        if nomi/denomi > 1 + tolerance
            error('v_i = %d',v_i)
        end
        pos_d(v_i) = ((loc_d(1,v_i+vel_win_len)-loc_d(1,v_i-vel_win_len))>=0)*2-1;
        angdiff(v_i) = acos(nomi/denomi - tolerance)/pi*180 * pos_d(v_i);
%         angdiff(v_i) = angle(nomi/denomi)/pi*180 * pos_d(v_i);
        v_ang(v_i) = angdiff(v_i) / (0.001*velocity_smooth_win_len);
        ang(v_i) = acos(double(loc_d(2,v_i)/norm(loc_d(:,v_i))) - tolerance)/pi*180 * sign(loc_d(1,v_i));
%         ang(v_i) = angle(double(loc_d(2,v_i)/norm(loc_d(:,v_i))))/pi*180 * sign(loc_d(1,v_i));
    end
end

end