%% baseline
subj_i = 5;
filename = sprintf('s%02d_ssvep_baseline.xdf',subj_i);
filepath = 'D:\ssvep dataset new\';
eegpath = which('eeglab');
eegpath = eegpath(1:end-8);
EEG = pop_loadxdf([filepath,filename]);

%% preprocessing
% band pass
EEG = pop_eegfiltnew(EEG,1,20);
% look up channel location
% Rename channels label
chLabel_new = {'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3',...
               'FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1',...
               'P3','P5','P7','P9','PO7','PO3','O1','Iz','Oz','POz','Pz','CPz',...
               'Fpz','Fp2','AF8','AF4','AFz','Fz','F2','F4','F6','F8',...
               'FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8',...
               'TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2'};
EEG = pop_select(EEG,'nochannel',{'Trig1','EX3','EX4','EX5','EX6','EX7','EX8'});
EEG = pop_reref(EEG,{'EX1','EX2'});
[EEG.chanlocs.labels] = deal(chLabel_new{:});
% Add channel location    
EEG = pop_chanedit(EEG, 'lookup',[eegpath,'plugins/dipfit/standard_BEM/elec/standard_1005.elc']); % MNI
thresFlatChannel = 5;
highPassBand = -1;
thresPoorCorrChannel = 0.7;
thresLineNoiseChannel = 4;
thresASR = 10;
thresWindow = -1;
EEG = clean_rawdata(EEG, thresFlatChannel, highPassBand, thresPoorCorrChannel, thresLineNoiseChannel, thresASR, thresWindow);

%% extract amp epoch
% ssvep epoch baseline
ssvep8 = pop_epoch(EEG,{'SSVEP baseline: Trial 0',...
                        'SSVEP baseline: Trial 1',...
                        'SSVEP baseline: Trial 2',...
                        'SSVEP baseline: Trial 3'}, [0 3]);
ssvep9 = pop_epoch(EEG,{'SSVEP baseline: Trial 4',...
                        'SSVEP baseline: Trial 5',...
                        'SSVEP baseline: Trial 6',...
                        'SSVEP baseline: Trial 7'}, [0 3]);
ssvep10 = pop_epoch(EEG,{'SSVEP baseline: Trial 8',...
                        'SSVEP baseline: Trial 9',...
                        'SSVEP baseline: Trial 10',...
                        'SSVEP baseline: Trial 11'}, [0 3]);
ssvep11 = pop_epoch(EEG,{'SSVEP baseline: Trial 12',...
                        'SSVEP baseline: Trial 13',...
                        'SSVEP baseline: Trial 14',...
                        'SSVEP baseline: Trial 15'}, [0 3]);
% ssvep8_h = pop_epoch(EEG,{'SSVEP baseline: Trial 0',...
%                         'SSVEP baseline: Trial 1',...
%                         'SSVEP baseline: Trial 2',...
%                         'SSVEP baseline: Trial 3'}, [0 3]);
% ssvep9_h = pop_epoch(EEG,{'SSVEP baseline: Trial 4',...
%                         'SSVEP baseline: Trial 5',...
%                         'SSVEP baseline: Trial 6',...
%                         'SSVEP baseline: Trial 7'}, [0 3]);
% ssvep10_h = pop_epoch(EEG,{'SSVEP baseline: Trial 8',...
%                         'SSVEP baseline: Trial 9',...
%                         'SSVEP baseline: Trial 10',...
%                         'SSVEP baseline: Trial 11'}, [0 3]);
% ssvep11_h = pop_epoch(EEG,{'SSVEP baseline: Trial 12',...
%                         'SSVEP baseline: Trial 13',...
%                         'SSVEP baseline: Trial 14',...
%                         'SSVEP baseline: Trial 15'}, [0 3]);
                    
%% calculate PSD
tarCh = {'O1','O2','Oz','POz','PO4','PO3','PO7','PO8'};
[spec8, freqs] = spectopo(ssvep8.data,0,ssvep8.srate,'plot','off');
[spec9, ~] = spectopo(ssvep9.data,0,ssvep9.srate,'plot','off');
[spec10, ~] = spectopo(ssvep10.data,0,ssvep10.srate,'plot','off');
[spec11, ~] = spectopo(ssvep11.data,0,ssvep11.srate,'plot','off');
% [spec8_h, freqs] = spectopo(ssvep8_h.data,0,ssvep8_h.srate,'plot','off');
% [spec9_h, ~] = spectopo(ssvep9_h.data,0,ssvep9_h.srate,'plot','off');
% [spec10_h, ~] = spectopo(ssvep10_h.data,0,ssvep10_h.srate,'plot','off');
% [spec11_h, ~] = spectopo(ssvep11_h.data,0,ssvep11_h.srate,'plot','off');
plt_f = 0:20;
plt8 = mean(spec8(ismember({ssvep8.chanlocs.labels},tarCh),plt_f+1),1);
plt9 = mean(spec9(ismember({ssvep8.chanlocs.labels},tarCh),plt_f+1),1);
plt10 = mean(spec10(ismember({ssvep8.chanlocs.labels},tarCh),plt_f+1),1);
plt11 = mean(spec11(ismember({ssvep8.chanlocs.labels},tarCh),plt_f+1),1);
% plt8_h = mean(spec8_h(:,plt_f+1),1);
% plt9_h = mean(spec9_h(:,plt_f+1),1);
% plt10_h = mean(spec10_h(:,plt_f+1),1);
% plt11_h = mean(spec11_h(:,plt_f+1),1);
figure
plot(plt_f, plt8,'b-o','linewidth',3,'DisplayName','8Hz');
hold on
grid on
plot(plt_f, plt9,'r-o','linewidth',3,'DisplayName','9Hz');
plot(plt_f, plt10,'g-o','linewidth',3,'DisplayName','10Hz');
plot(plt_f, plt11,'m-o','linewidth',3,'DisplayName','11Hz');
legend
xlabel('Frequency(Hz)')
ylabel('Power (\muV^2)')
set(gca,'fontsize',20)
set(gcf,'color','w')
title('Amp')
% figure
% plot(plt_f, plt8_h,'b-o','linewidth',3,'DisplayName','8Hz halo');
% hold on
% grid on
% plot(plt_f, plt9_h,'r-o','linewidth',3,'DisplayName','9Hz halo');
% plot(plt_f, plt10_h,'g-o','linewidth',3,'DisplayName','10Hz halo');
% plot(plt_f, plt11_h,'m-o','linewidth',3,'DisplayName','11Hz halo');
% legend
% xlabel('Frequency(Hz)')
% ylabel('Power (\muV^2)')
% set(gca,'fontsize',20)
% set(gcf,'color','w')
% title('Halo')

%% histogram
figure
hold on
grid on
histogram(reshape(ssvep8.data(1,:,:),1,[]),'binwidth',1,'DisplayName','8');
histogram(reshape(ssvep9.data(1,:,:),1,[]),'binwidth',1,'DisplayName','9');
histogram(reshape(ssvep10.data(1,:,:),1,[]),'binwidth',1,'DisplayName','10');
histogram(reshape(ssvep11.data(1,:,:),1,[]),'binwidth',1,'DisplayName','11');
legend
title('Amp')
% figure
% hold on
% grid on
% histogram(reshape(ssvep8_h.data(1,:,:),1,[]),'binwidth',1,'DisplayName','8');
% histogram(reshape(ssvep9_h.data(1,:,:),1,[]),'binwidth',1,'DisplayName','9');
% histogram(reshape(ssvep10_h.data(1,:,:),1,[]),'binwidth',1,'DisplayName','10');
% histogram(reshape(ssvep11_h.data(1,:,:),1,[]),'binwidth',1,'DisplayName','11');
% legend
% title('Halo')

%% load task data
amp_EEG = pop_loadxdf('s05_ssvep_amp.xdf');
halo_EEG = pop_loadxdf('s05_ssvep_halo.xdf');

%% preprocessing
% band pass
amp_EEG = pop_eegfiltnew(amp_EEG,1,20);
halo_EEG = pop_eegfiltnew(halo_EEG,1,20);
% look up channel location
% Rename channels label
chLabel_new = {'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3',...
               'FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1',...
               'P3','P5','P7','P9','PO7','PO3','O1','Iz','Oz','POz','Pz','CPz',...
               'Fpz','Fp2','AF8','AF4','AFz','Fz','F2','F4','F6','F8',...
               'FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8',...
               'TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2'};
amp_EEG = pop_select(amp_EEG,'nochannel',{'Trig1','EX3','EX4','EX5','EX6','EX7','EX8'});
amp_EEG = pop_reref(amp_EEG,{'EX1','EX2'});
[amp_EEG.chanlocs.labels] = deal(chLabel_new{:});
% Add channel location    
amp_EEG = pop_chanedit(amp_EEG, 'lookup',[eegpath,'plugins/dipfit/standard_BEM/elec/standard_1005.elc']); % MNI
thresFlatChannel = 5;
highPassBand = -1;
thresPoorCorrChannel = 0.7;
thresLineNoiseChannel = 4;
thresASR = 10;
thresWindow = -1;
EEG_amp = clean_rawdata(amp_EEG, thresFlatChannel, highPassBand, thresPoorCorrChannel, thresLineNoiseChannel, -1, thresWindow);

%% halo
halo_EEG = pop_select(halo_EEG,'nochannel',{'Trig1','EX3','EX4','EX5','EX6','EX7','EX8'});
halo_EEG = pop_reref(halo_EEG,{'EX1','EX2'});
[halo_EEG.chanlocs.labels] = deal(chLabel_new{:});
% Add channel location    
halo_EEG = pop_chanedit(halo_EEG, 'lookup',[eegpath,'plugins/dipfit/standard_BEM/elec/standard_1005.elc']); % MNI
EEG_halo = clean_rawdata(halo_EEG, thresFlatChannel, highPassBand, thresPoorCorrChannel, thresLineNoiseChannel, -1, thresWindow);

%% extract epoch
% tarCh = {'O1','O2','Oz','POz','PO4','PO3','PO7','PO8'};
tarCh = {'O1','O2','Oz','POz','PO4','PO3'};
tarFreq = 1:20;
epoch_len = [0 3];
psd_lib_amp = zeros(2,4,length(tarFreq)); % ring by direct by freq
psd_lib_halo = zeros(2,4,length(tarFreq)); % ring by direct by freq

for ring_txt = 1:2
    for direction = 0:3
        ring_ev = cellfun(@(x) ~isempty(regexp(x,sprintf('Ring %d',ring_txt),'ONCE')),{EEG_halo.event.type});
        dir_ev = cellfun(@(x) ~isempty(regexp(x,['index\(L\)\:\s',sprintf('%d',direction)],'ONCE')),{EEG_halo.event.type});
        tar_ev_halo = ring_ev & dir_ev;
        ring_ev = cellfun(@(x) ~isempty(regexp(x,sprintf('Ring %d',ring_txt),'ONCE')),{EEG_amp.event.type});
        dir_ev = cellfun(@(x) ~isempty(regexp(x,['index\(L\)\:\s',sprintf('%d',direction)],'ONCE')),{EEG_amp.event.type});
        tar_ev_amp = ring_ev & dir_ev;
        tar_halo = pop_epoch(EEG_halo,{},epoch_len,'eventindices',find(tar_ev_halo));
        tar_amp = pop_epoch(EEG_amp,{},epoch_len,'eventindices',find(tar_ev_amp));
        tar_halo_data = reshape(tar_halo.data(ismember({tar_halo.chanlocs.labels},tarCh),:,:),length(tarCh),[]);
        tar_amp_data = reshape(tar_amp.data(ismember({tar_amp.chanlocs.labels},tarCh),:,:),length(tarCh),[]);
        [spec, freq] = spectopo(tar_halo_data,0,tar_halo.srate,'plot','off');
        psd_lib_halo(ring_txt,direction+1,:) = mean(spec(:,ismember(freq,tarFreq)));
        [spec, freq] = spectopo(tar_amp_data,0,tar_amp.srate,'plot','off');
        psd_lib_amp(ring_txt,direction+1,:) = mean(spec(:,ismember(freq,tarFreq)));
    end
end

%%
plt_ring = 1;
plt_dir = 4;
figure
plot(tarFreq,squeeze(psd_lib_amp(plt_ring,plt_dir,:)))


%% artifact removal for Ring 2
EEG = pop_loadset('s04_halo_stimLock_left.set');
% before artifact cleaning
pop_prop(EEG, 1, 1,NaN, {'freqrange',[1 20]});
% my_rmEpoch
[rm_Amp, rm_idx] = my_rmEpoch(EEG);
pop_prop(rm_Amp, 1, 1,NaN, {'freqrange',[1 20]});

%% Finding PC regressors
% perform moving average on single epoch matrix to reduce noises for the
% following PCA process

%% parameter setting
tarCh = 'O1';
mvavg_winlen = 3; % trial
t_interest = [500 2000]; % msec
nb_PC = 1;

%% 
epoch = squeeze(EEG.data(ismember({EEG.chanlocs.labels},tarCh),:,:))';
random_group = randperm(size(epoch,1));
sort_epoch = epoch(random_group,:);
mvavg_idx = discretize(1:size(sort_epoch,1),ceil(size(sort_epoch,1)/mvavg_winlen));
mvavg_epoch = zeros(max(mvavg_idx), size(sort_epoch,2));
for i = 1:max(mvavg_idx)
    mvavg_epoch(i,:) = mean(sort_epoch(mvavg_idx==i,:));
end

% segement out data within time period of interest
% reg_mat = zeros(size(epoch,2),nb_PC);
pca_epoch = mvavg_epoch;
% PCA
% trial as features, time as observation
% built covariance matrix (channel by channel)
xcov = pca_epoch * pca_epoch';
% eigenvalue decomposition
[V, ~] = eig(xcov);
% extract nb_PC
V_pick = V(:,end-nb_PC+1:end);
% reconstruct latent component (comp by time)
pc_act = V_pick'*pca_epoch;
% [pc_weight,pc_act] = pca(pca_epoch,nb_PC);
% zero paddling
pc_regressors = pc_act;
% regressors matrix
reg_mat = pc_regressors';

reg_mat = [ones(size(epoch,2),1), reg_mat];

%% Multiple Linear Regression with dispertion term
reg_epoch = zeros(size(epoch));
for i = 1:size(epoch,1)
    single_epoch = epoch(i,:)';
	b = regress(single_epoch, reg_mat);
    reg_epoch(i,:) = reg_mat*b;
end

% substract regression results
rmReg_epoch = epoch - reg_epoch;

%% visualization
[spec,freq] = spectopo(rmReg_epoch,0,EEG.srate,'plot','off');
plt_f = 0:20;
plt_spec = spec(:,plt_f+1);
figure
plot(plt_f, mean(plt_spec))
pseudoEEG = EEG;
pseudoEEG.data(1,:,:) = rmReg_epoch';
pop_prop(pseudoEEG, 1, 1,NaN, {'freqrange',[1 20]});

%% back project PC
[coeff, score, latent, tsquared, explained] = pca(epoch');
latent_norm = latent./sum(latent);
nb_PC = find(cumsum(latent_norm)>=0.9,1);

figure
plot(score(:,1:nb_PC))
grid on
legend({'PC1','PC2','PC3'})
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')

preserved_PC = score;
preserved_PC(:,1:nb_PC) = 0;
proj_data = pinv(coeff')* preserved_PC';

figure
plot(0:749, proj_data)
hold on
grid on
plot(0:749, mean(proj_data), 'k-', 'linewidth',3)
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')

[spec, freq] = spectopo(proj_data,0,EEG.srate,'plot','off');
plt_f = 0:20;
plt_spec = spec(:,plt_f+1);
figure
plot(plt_f, plt_spec)
hold on
grid on
plot(plt_f, mean(plt_spec), 'k-','linewidth',3)
xlabel('Frequency (Hz)')
ylabel('Power (\muV^2)')

[spec, freq] = spectopo(epoch,0,EEG.srate,'plot','off');
plt_f = 0:20;
plt_spec = spec(:,plt_f+1);
figure
plot(plt_f, plt_spec)
hold on
grid on
plot(plt_f, mean(plt_spec), 'k-','linewidth',3)
xlabel('Frequency (Hz)')
ylabel('Power (\muV^2)')









