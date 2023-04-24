%% This function epoch data

%% load task data
filepath = '/home/yuan/Documents/2021 HM_visual_oddball/preproc_data/';
EEG_amp = pop_loadset([filepath,'1145_ssvep_amp.set']);
EEG_halo = pop_loadset([filepath,'1145_ssvep_halo.set']);

%% extract stimulus onset epoch
% tarCh = {'O1','O2','Oz','POz','PO4','PO3','PO7','PO8'};
tarCh = {'O1','O2','Oz','POz','PO4','PO3'};
tarFreq = 1:20;
epoch_len = [0 3];
psd_lib_amp = zeros(2,4,length(tarFreq)); % ring by direct by freq
psd_lib_halo = zeros(2,4,length(tarFreq)); % ring by direct by freq
idx_check_halo = false(2,4,length(EEG_halo.event));
idx_check_amp = false(2,4,length(EEG_amp.event));

for ring_txt = 1:2
    for direction = 0:3
        ring_ev = cellfun(@(x) ~isempty(regexp(x,sprintf('Ring %d',ring_txt),'ONCE')),{EEG_halo.event.type});
        dir_ev = cellfun(@(x) ~isempty(regexp(x,['index\(.\)\:\s',sprintf('%d',direction)],'ONCE')),{EEG_halo.event.type});
        tar_ev_halo = ring_ev & dir_ev;
        ring_ev = cellfun(@(x) ~isempty(regexp(x,sprintf('Ring %d',ring_txt),'ONCE')),{EEG_amp.event.type});
        dir_ev = cellfun(@(x) ~isempty(regexp(x,['index\(.\)\:\s',sprintf('%d',direction)],'ONCE')),{EEG_amp.event.type});
        tar_ev_amp = ring_ev & dir_ev;
        idx_check_halo(ring_txt, direction+1, :) = tar_ev_halo;
        idx_check_amp(ring_txt, direction+1, :) = tar_ev_amp;
        % epoch
        tar_halo = pop_epoch(EEG_halo,{},epoch_len,'eventindices',find(tar_ev_halo));
        tar_amp = pop_epoch(EEG_amp,{},epoch_len,'eventindices',find(tar_ev_amp));
        % calculate PSD
        tar_halo_data = reshape(tar_halo.data(ismember({tar_halo.chanlocs.labels},tarCh),:,:),length(tarCh),[]);
        tar_amp_data = reshape(tar_amp.data(ismember({tar_amp.chanlocs.labels},tarCh),:,:),length(tarCh),[]);
        [spec, freq] = spectopo(tar_halo_data,0,tar_halo.srate,'plot','off');
        psd_lib_halo(ring_txt,direction+1,:) = mean(spec(:,ismember(freq,tarFreq)));
        [spec, freq] = spectopo(tar_amp_data,0,tar_amp.srate,'plot','off');
        psd_lib_amp(ring_txt,direction+1,:) = mean(spec(:,ismember(freq,tarFreq)));
    end
end


%% visualization
ring_i = 1;
cmap = {'b','r','g','m'};
disname = {'8Hz','9Hz','10Hz','11Hz'};

figure
grid on
hold on
% Halo
for dir_i = 1:4
    plot(tarFreq, squeeze(psd_lib_halo(ring_i,dir_i,:)),'-o',...
        'color',cmap{dir_i},'linewidth',3,'DisplayName',['Halo ',disname{dir_i}]);
end
% Amp
% for dir_i = 1:4
%     plot(tarFreq, squeeze(psd_lib_amp(ring_i,dir_i,:)),'-x',...
%         'color',cmap{dir_i},'linewidth',3,'DisplayName',['Amp ',disname{dir_i}]);
% end
legend
xlabel('Frequency(Hz)')
ylabel('Power (\muV^2)')
set(gca,'fontsize',20)
set(gcf,'color','w')
if ring_i == 1
    title('Inner Ring')
else
    title('Outer Ring')
end

%% extract gip onset epoch
tarCh = {'O1','O2','Oz','POz','PO4','PO3'};
tarFreq = 1:20;
epoch_len = [-1 2];
psd_lib_amp = zeros(2,4,length(tarFreq)); % ring by direct by freq
psd_lib_halo = zeros(2,4,length(tarFreq)); % ring by direct by freq
idx_check_halo = false(2,4,length(EEG_halo.event));
idx_check_amp = false(2,4,length(EEG_amp.event));
time_signal_halo = cell(2,4); % ring by direct
time_signal_amp = cell(2,4); % ring by direct

for ring_txt = 1:2
    % halo
    ring_ev = cellfun(@(x) ~isempty(regexp(x,sprintf('Ring %d',ring_txt),'ONCE')),{EEG_halo.event.type});
    % find 1st gip after stimulus onset
    [gip_up_idx, gip_down_idx, gip_left_idx, gip_right_idx,miss_idx] = find_1st_gip(EEG_halo,ring_ev);
    % epoch
    for d_i = 1:4
        switch d_i 
            case 1
                tar_ev = gip_right_idx;
            case 2
                tar_ev = gip_up_idx;
            case 3
                tar_ev = gip_left_idx;
            case 4
                tar_ev = gip_down_idx;
        end
        tar_halo = pop_epoch(EEG_halo,{},epoch_len,'eventindices',find(tar_ev));
        time_signal_halo{ring_txt,d_i} = tar_halo.data(ismember({tar_halo.chanlocs.labels},tarCh),:,:);
        cal_spec_time = tar_halo.times>=0;
        % calculate PSD
        tar_halo_data = reshape(tar_halo.data(ismember({tar_halo.chanlocs.labels},tarCh),cal_spec_time,:),length(tarCh),[]);
        [spec, freq] = spectopo(tar_halo_data,0,tar_halo.srate,'plot','off');
        psd_lib_halo(ring_txt,d_i,:) = mean(spec(:,ismember(freq,tarFreq)));
    end
    % amp    
    ring_ev = cellfun(@(x) ~isempty(regexp(x,sprintf('Ring %d',ring_txt),'ONCE')),{EEG_amp.event.type});
    [gip_up_idx, gip_down_idx, gip_left_idx, gip_right_idx,miss_idx] = find_1st_gip(EEG_amp,ring_ev);
    % epoch
    for d_i = 1:4
        switch d_i 
            case 1
                tar_ev = gip_right_idx;
            case 2
                tar_ev = gip_up_idx;
            case 3
                tar_ev = gip_left_idx;
            case 4
                tar_ev = gip_down_idx;
        end
        tar_amp = pop_epoch(EEG_amp,{},epoch_len,'eventindices',find(tar_ev));
        time_signal_amp{ring_txt,d_i} = tar_amp.data(ismember({tar_amp.chanlocs.labels},tarCh),:,:);
        cal_spec_time = tar_amp.times>=0;
        % calculate PSD
        tar_amp_data = reshape(tar_amp.data(ismember({tar_amp.chanlocs.labels},tarCh),cal_spec_time,:),length(tarCh),[]);
        [spec, freq] = spectopo(tar_amp_data,0,tar_amp.srate,'plot','off');
        psd_lib_amp(ring_txt,d_i,:) = mean(spec(:,ismember(freq,tarFreq)));
    end
end

%% visualization
ring_i = 2;
cmap = {'b','r','g','m'};
disname = {'8Hz','9Hz','10Hz','11Hz'};

figure
grid on
hold on
% Halo
% for dir_i = 1:4
%     plot(tarFreq, squeeze(psd_lib_halo(ring_i,dir_i,:)),'-o',...
%         'color',cmap{dir_i},'linewidth',3,'DisplayName',['Halo ',disname{dir_i}]);
% end
% Amp
for dir_i = 1:4
    plot(tarFreq, squeeze(psd_lib_amp(ring_i,dir_i,:)),'-x',...
        'color',cmap{dir_i},'linewidth',3,'DisplayName',['Amp ',disname{dir_i}]);
end
legend
xlabel('Frequency(Hz)')
ylabel('Power (\muV^2)')
set(gca,'fontsize',20)
set(gcf,'color','w')
if ring_i == 1
    title('Inner Ring')
else
    title('Outer Ring')
end


%% plot time course 
ring_i = 1;
dir_i = 1;
ch_idx = 3; % Oz
plt_data = time_signal_amp{ring_i,dir_i};

figure
plot(tar_amp.times, mean(plt_data(ch_idx,:,:),3),'linewidth',1);
hold on
grid on
xline(0,'k--')
