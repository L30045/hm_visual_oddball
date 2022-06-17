function [selectCh,chweights] = selectCh_CCA(EEG,top_k)
%% this function use Canoical Correlation Analysis to weight the channels for calculating SSVEP
% Input:
%   EEG: preprocessed EEG with SSVEP event markers.
%   top_k: top k channels with highest weight in the first CCA component
%   to be return. If it is empty, return all channels.
% Output:
%   selectCh: selected channels. (Default: all)
%   chweights: CCA weights.

%% parameter setting
if isempty(top_k)
    top_k = EEG.nbchan;
end

%% find the event markers
tarFreq = 8;
if ismember({EEG.event.type},'SSVEP baseline: Trial 0')==0
    % process as experiment file
    idx_ev = cellfun(@(x) ~isempty(regexp(x,'Ring 1.*Frequency: 8','once')),{EEG.event.type});
    if sum(idx_ev)==0
        idx_ev = cellfun(@(x) ~isempty(regexp(x,'Ring 1.*Cube index\(.\)..0','ONCE')),{EEG.event.type});
    end
    % epoch EEG
    len_epoch = [0 3]; % sec
    EEG_epoch = pop_epoch(EEG,{},len_epoch,'eventindices',find(idx_ev),'epochinfo', 'yes');
else
    % process as baseline file
    EEG_epoch = pop_epoch(EEG,{'SSVEP baseline: Trial 0',...
                        'SSVEP baseline: Trial 1',...
                        'SSVEP baseline: Trial 2',...
                        'SSVEP baseline: Trial 3'}, [0 3]);
end

%% CCA
% create sin and cos wave
t = EEG_epoch.times/1000;
sin_tar = sin(tarFreq*2*pi*t)';
cos_tar = cos(tarFreq*2*pi*t)';
Y = [sin_tar,cos_tar];

% cca with epoch data for each trial
cca_weight_lib = zeros(EEG_epoch.nbchan,size(EEG_epoch.data,3));
for t_i = 1:size(EEG_epoch.data,3)
    [X_cca,Y_cca] = canoncorr(EEG_epoch.data(:,:,t_i)',Y);
    cca_weight_lib(:,t_i) = X_cca(:,1);
end

% normalize cca_weight_lib
cca_weight_lib = cca_weight_lib./sqrt(sum(cca_weight_lib.^2,1));
% mean weight across trials
chweights = mean(cca_weight_lib,2);

% select channels
[~, sort_i] = sort(chweights,'descend');
selectCh = {EEG_epoch.chanlocs(sort_i(1:top_k)).labels};

% visualization
tmp = pop_select(EEG_epoch,'channel',selectCh);
figure
topoplot([],tmp.chanlocs,'style','blank','electrodes','labelpoint','chaninfo',tmp.chaninfo);

end
