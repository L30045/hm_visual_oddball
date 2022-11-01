function output = find_ev_time(epoch_struct)
ev_time = struct('diff_gip_stim',[],'diff_grab_stim',[],'diff_fix_stim',[]);
output = struct('cir', ev_time, 'tri', ev_time);


for lock_i = 1:2
    switch lock_i
        case 1
            EEG = epoch_struct.std_epoch;
            gip_marker = 'circle_gip_start';
            fix_marker = 'circle_fix_start';
            field_name = 'cir';
        case 2
            EEG = epoch_struct.dev_epoch;
            gip_marker = 'triangle_gip_start';
            fix_marker = 'triangle_fix_start';
            field_name = 'tri';
    end
    
    gip_idx = cellfun(@(x) find(cellfun(@(y) strcmp(y,gip_marker),x),1), {EEG.epoch.eventtype},'uniformoutput',0);
    fix_idx = cellfun(@(x) find(cellfun(@(y) strcmp(y,fix_marker),x),1), {EEG.epoch.eventtype},'uniformoutput',0);   
    % remove missing event
    rm_idx = cellfun(@isempty,gip_idx)|cellfun(@isempty,fix_idx);
    if lock_i == 1
        grab_idx = cellfun(@(x) find(cellfun(@(y) strcmp(y,'grab'),x),1), {EEG.epoch.eventtype},'uniformoutput',0);
        rm_idx = rm_idx|cellfun(@isempty,grab_idx);
        grab_idx = grab_idx(~rm_idx);
        grab_t = zeros(1,length(grab_idx));
    end
    gip_idx = gip_idx(~rm_idx);
    fix_idx = fix_idx(~rm_idx);
    gip_t = zeros(1,length(gip_idx));
    fix_t = zeros(1,length(gip_idx));
    
    count = 1;
    for ev_i = find(~rm_idx)
        gip_t(count) = EEG.epoch(ev_i).eventlatency{gip_idx{count}}; % ms
        fix_t(count) = EEG.epoch(ev_i).eventlatency{fix_idx{count}}; % ms
        if lock_i == 1
            grab_t(count) = EEG.epoch(ev_i).eventlatency{grab_idx{count}}; % ms
        end
        count = count+1;
    end
    % remove event with gip - stim < 100 ms
    rm_idx = gip_t < 100;
    gip_t = gip_t(~rm_idx);
    fix_t = fix_t(~rm_idx);
    output.(field_name).diff_gip_stim = gip_t;
    output.(field_name).diff_fix_stim = fix_t;
    if lock_i == 1
        grab_t = grab_t(~rm_idx);
        output.(field_name).diff_grab_stim = grab_t;
    end
end