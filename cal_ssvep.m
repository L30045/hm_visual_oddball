filename = 's04_amp_baseline.xdf';
amp_EEG = pop_loadxdf([filepath,filesep,filename]);
halo_EEG = pop_loadxdf([filepath,filesep,'s04_halo_baseline.xdf']);

%% preprocessing
% bandpass filter
EEG_bp_amp = pop_eegfiltnew(amp_EEG,0.5,20);
% preserve only channls O1, O2, POz
EEG_rmCh_amp = pop_select(EEG_bp_amp,'channel',{'O1','O2','POz'});
% reref
EEG_amp = pop_reref(EEG_rmCh_amp,[]);
% bandpass filter
EEG_bp_halo = pop_eegfiltnew(halo_EEG,0.5,20);
% preserve only channls O1, O2, POz
EEG_rmCh_halo= pop_select(EEG_bp_halo,'channel',{'O1','O2','POz'});
% reref
EEG_halo = pop_reref(EEG_rmCh_halo,[]);

pop_saveset(EEG_amp,'s04_amp_baseline.set');
pop_saveset(EEG_halo,'s04_halo_baseline.set');

%% extract amp epoch
% ssvep epoch baseline
ssvep8 = pop_epoch(EEG_amp,{'SSVEP baseline: Trial 0',...
                        'SSVEP baseline: Trial 1',...
                        'SSVEP baseline: Trial 2',...
                        'SSVEP baseline: Trial 3'}, [0 3]);
ssvep9 = pop_epoch(EEG_amp,{'SSVEP baseline: Trial 4',...
                        'SSVEP baseline: Trial 5',...
                        'SSVEP baseline: Trial 6',...
                        'SSVEP baseline: Trial 7'}, [0 3]);
ssvep10 = pop_epoch(EEG_amp,{'SSVEP baseline: Trial 8',...
                        'SSVEP baseline: Trial 9',...
                        'SSVEP baseline: Trial 10',...
                        'SSVEP baseline: Trial 11'}, [0 3]);
ssvep11 = pop_epoch(EEG_amp,{'SSVEP baseline: Trial 12',...
                        'SSVEP baseline: Trial 13',...
                        'SSVEP baseline: Trial 14',...
                        'SSVEP baseline: Trial 15'}, [0 3]);
ssvep8_h = pop_epoch(EEG_halo,{'SSVEP baseline: Trial 0',...
                        'SSVEP baseline: Trial 1',...
                        'SSVEP baseline: Trial 2',...
                        'SSVEP baseline: Trial 3'}, [0 3]);
ssvep9_h = pop_epoch(EEG_halo,{'SSVEP baseline: Trial 4',...
                        'SSVEP baseline: Trial 5',...
                        'SSVEP baseline: Trial 6',...
                        'SSVEP baseline: Trial 7'}, [0 3]);
ssvep10_h = pop_epoch(EEG_halo,{'SSVEP baseline: Trial 8',...
                        'SSVEP baseline: Trial 9',...
                        'SSVEP baseline: Trial 10',...
                        'SSVEP baseline: Trial 11'}, [0 3]);
ssvep11_h = pop_epoch(EEG_halo,{'SSVEP baseline: Trial 12',...
                        'SSVEP baseline: Trial 13',...
                        'SSVEP baseline: Trial 14',...
                        'SSVEP baseline: Trial 15'}, [0 3]);
                    
%% calculate PSD
[spec8, freqs] = spectopo(ssvep8.data,0,ssvep8.srate,'plot','off');
[spec9, ~] = spectopo(ssvep9.data,0,ssvep9.srate,'plot','off');
[spec10, ~] = spectopo(ssvep10.data,0,ssvep10.srate,'plot','off');
[spec11, ~] = spectopo(ssvep11.data,0,ssvep11.srate,'plot','off');
[spec8_h, freqs] = spectopo(ssvep8_h.data,0,ssvep8_h.srate,'plot','off');
[spec9_h, ~] = spectopo(ssvep9_h.data,0,ssvep9_h.srate,'plot','off');
[spec10_h, ~] = spectopo(ssvep10_h.data,0,ssvep10_h.srate,'plot','off');
[spec11_h, ~] = spectopo(ssvep11_h.data,0,ssvep11_h.srate,'plot','off');
plt_f = 0:20;
plt8 = mean(spec8(:,plt_f+1),1);
plt9 = mean(spec9(:,plt_f+1),1);
plt10 = mean(spec10(:,plt_f+1),1);
plt11 = mean(spec11(:,plt_f+1),1);
plt8_h = mean(spec8_h(:,plt_f+1),1);
plt9_h = mean(spec9_h(:,plt_f+1),1);
plt10_h = mean(spec10_h(:,plt_f+1),1);
plt11_h = mean(spec11_h(:,plt_f+1),1);
figure
plot(plt_f, plt8,'b-o','linewidth',3,'DisplayName','8Hz amp');
hold on
grid on
plot(plt_f, plt9,'r-o','linewidth',3,'DisplayName','9Hz amp');
plot(plt_f, plt10,'g-o','linewidth',3,'DisplayName','10Hz amp');
plot(plt_f, plt11,'m-o','linewidth',3,'DisplayName','11Hz amp');
legend
xlabel('Frequency(Hz)')
ylabel('Power (\muV^2)')
set(gca,'fontsize',20)
set(gcf,'color','w')
title('Amp')
figure
plot(plt_f, plt8_h,'b-o','linewidth',3,'DisplayName','8Hz halo');
hold on
grid on
plot(plt_f, plt9_h,'r-o','linewidth',3,'DisplayName','9Hz halo');
plot(plt_f, plt10_h,'g-o','linewidth',3,'DisplayName','10Hz halo');
plot(plt_f, plt11_h,'m-o','linewidth',3,'DisplayName','11Hz halo');
legend
xlabel('Frequency(Hz)')
ylabel('Power (\muV^2)')
set(gca,'fontsize',20)
set(gcf,'color','w')
title('Halo')

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
figure
hold on
grid on
histogram(reshape(ssvep8_h.data(1,:,:),1,[]),'binwidth',1,'DisplayName','8');
histogram(reshape(ssvep9_h.data(1,:,:),1,[]),'binwidth',1,'DisplayName','9');
histogram(reshape(ssvep10_h.data(1,:,:),1,[]),'binwidth',1,'DisplayName','10');
histogram(reshape(ssvep11_h.data(1,:,:),1,[]),'binwidth',1,'DisplayName','11');
legend
title('Halo')

%% load task data
amp_EEG = pop_loadxdf('s04_amp.xdf');
halo_EEG = pop_loadxdf('s04_halo.xdf');

%% preprocessing
% bandpass filter
EEG_bp_amp = pop_eegfiltnew(amp_EEG,0.5,20);
% preserve only channls O1, O2, POz
EEG_rmCh_amp = pop_select(EEG_bp_amp,'channel',{'O1','O2','POz'});
% reref
EEG_amp = pop_reref(EEG_rmCh_amp,[]);
% bandpass filter
EEG_bp_halo = pop_eegfiltnew(halo_EEG,0.5,20);
% preserve only channls O1, O2, POz
EEG_rmCh_halo= pop_select(EEG_bp_halo,'channel',{'O1','O2','POz'});
% reref
EEG_halo = pop_reref(EEG_rmCh_halo,[]);
pop_saveset(EEG_amp,'s04_amp.set');
pop_saveset(EEG_halo,'s04_halo.set');

%%
figure; histogram(EEG_halo.data(1,:),'binwidth',1);
figure; histogram(EEG_amp.data(1,:),'binwidth',1);

%% extract epoch
psd_lib_amp = zeros(2,4,22); % ring by direct by freq
psd_lib_halo = zeros(2,4,22); % ring by direct by freq


% for ring_txt = 1:2
%     for direction = 0:3
%         tar_ev_halo = cellfun(@(x) ~isempty(regexp(x,[sprintf('Ring %d\;',ring_txt),'\sTrial .; ',sprintf('Cube index(L): %d',direction)]),'ONCE'));
%         tar_ev_amp = 
%         tar_halo = pop_epoch(EEG_amp,tar_ev_halo,[0,3]);
%         tar_amp  = pop_epoch(EEG_amp,tar_ev_amp,[0,3]);

























