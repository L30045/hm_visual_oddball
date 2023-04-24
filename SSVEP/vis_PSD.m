function [psd_lib, EEG_epoch, pars] = vis_PSD(EEG,timelock,epoch_len,varargin)
%% This function epoch data into SSVEP epoch by timelock and epoch length,
%  and plot the PSD of SSVEP Epoch.

%% parameter setting
p = inputParser;
p.KeepUnmatched = true;
addRequired(p,'EEG');
addRequired(p,'timelock');
addRequired(p,'epoch_len') % epoch length
addOptional(p,'tarCh',{'O1','O2','Oz','POz','PO4','PO3'}) % target channel for calculating SSVEP
addOptional(p,'tarFreq',1:20) % target frequency for calculating PSD
addOptional(p,'cal_spec_time',[]) % time point to start calculate PSD. (Default: 0 to end of the trial)
parse(p,EEG,timelock,epoch_len,varargin{:})
tarCh = p.Results.tarCh;
tarFreq = p.Results.tarFreq;
cal_spec_time = p.Results.cal_spec_time;
pars = p.Results;

psd_lib = zeros(2,4,length(tarFreq)); % ring by direct by freq
time_signal = cell(2,4); % ring by direct
dir_lib = cell(2,1); % ring

%% epoch data
switch timelock
    case 'stim'
        for ring_txt = 1:2
            for direction = 0:3
                ring_ev = cellfun(@(x) ~isempty(regexp(x,sprintf('Ring %d',ring_txt),'ONCE')),{EEG.event.type});
                dir_ev = cellfun(@(x) ~isempty(regexp(x,['index\(.\)\:\s',sprintf('%d',direction)],'ONCE')),{EEG.event.type});
                tar_ev = ring_ev & dir_ev;
                dir_lib{ring_txt} = [dir_lib{ring_txt};tar_ev];
                % epoch
                EEG_epoch = pop_epoch(EEG,{},epoch_len,'eventindices',find(tar_ev));
                time_signal{ring_txt,direction+1} = EEG_epoch.data(ismember({EEG_epoch.chanlocs.labels},tarCh),:,:);
                if isempty(cal_spec_time)
                    cal_spec_time = EEG_epoch.times>=0;
                end
                % calculate PSD
                tar_data = reshape(EEG_epoch.data(ismember({EEG_epoch.chanlocs.labels},tarCh),cal_spec_time,:),length(tarCh),[]);
                [spec, freq] = spectopo(tar_data,0,EEG_epoch.srate,'plot','off');
                psd_lib(ring_txt,direction+1,:) = mean(spec(:,ismember(freq,tarFreq)));
            end
        end

    case 'gip'
        gip_miss_idx = cell(2,1); % ring
        for ring_txt = 1:2
            % halo
            ring_ev = cellfun(@(x) ~isempty(regexp(x,sprintf('Ring %d',ring_txt),'ONCE')),{EEG.event.type});
            % find 1st gip after stimulus onset
            [gip_up_idx, gip_down_idx, gip_left_idx, gip_right_idx,miss_idx] = find_1st_gip(EEG,ring_ev);
            dir_lib{ring_txt} = [gip_right_idx;gip_up_idx;gip_left_idx;gip_down_idx];
            gip_miss_idx{ring_txt} = miss_idx;
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
                EEG_epoch = pop_epoch(EEG,{},epoch_len,'eventindices',find(tar_ev));
                time_signal{ring_txt,d_i} = EEG_epoch.data(ismember({EEG_epoch.chanlocs.labels},tarCh),:,:);
                if isempty(cal_spec_time)
                    cal_spec_time = EEG_epoch.times>=0;
                end
                % calculate PSD
                tar_data = reshape(EEG_epoch.data(ismember({EEG_epoch.chanlocs.labels},tarCh),cal_spec_time,:),length(tarCh),[]);
                [spec, freq] = spectopo(tar_data,0,EEG_epoch.srate,'plot','off');
                psd_lib(ring_txt,d_i,:) = mean(spec(:,ismember(freq,tarFreq)));
            end
        end
        pars.gip_miss_idx = gip_miss_idx;
end
pars.cal_spec_time = cal_spec_time;
pars.time_signal = time_signal;
pars.dir_lib = dir_lib;

%% visualization
cmap = {'b','r','g','m'};
disname = {'8Hz','9Hz','10Hz','11Hz'};

for ring_i = 1:2
    figure
    grid on
    hold on
    for dir_i = 1:4
        plot(tarFreq, squeeze(psd_lib(ring_i,dir_i,:)),'-o',...
            'color',cmap{dir_i},'linewidth',3,'DisplayName',disname{dir_i});
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
end
