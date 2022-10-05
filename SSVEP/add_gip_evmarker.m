function [EEG, stimTime, gipTime] = add_gip_evmarker(EEG, ev_stim, ev_gip, marker_name)
% find out stimulus onset time
f_stim= find(ev_stim);
gip_idx = find(ev_gip);
stimTime = zeros(1,length(f_stim)); 
gipTime = zeros(1,length(f_stim));
for i = 1:length(f_stim)
    stimTime(i) = EEG.event(f_stim(i)).latency/EEG.srate*1000; % msec
    t_f = gip_idx(find(gip_idx>f_stim(i),1));
    if ~isempty(t_f)
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = [marker_name,'_start'];
        gipTime(i) = EEG.event(t_f).latency/EEG.srate*1000; % msec
        EEG.event(len_event+1).latency = EEG.event(t_f).latency; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
    else
        gipTime(i) = NaN;
        len_event = length(EEG.event);
        EEG.event(len_event+1).type = [marker_name,'_missing'];
        EEG.event(len_event+1).latency = EEG.event(f_stim(i)).latency; % change to sample point
        EEG.event(len_event+1).urevent = len_event+1;
    end
end
end