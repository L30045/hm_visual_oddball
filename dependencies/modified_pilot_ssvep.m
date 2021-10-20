%% SSVEP pilot analysis
EEG_baseline = pop_loadxdf('hm_oddball_ssvep_baseline_ec.xdf','streamtype','EEG');
EEG_ssvep = pop_loadxdf('hm_oddball_ssvep_1.xdf','streamtype','EEG');
% addpath('xdfimport1.13')
% addpath('related_function')
% EEG_baseline = pop_loadxdf('hm_oddball_ssvep_2_baseline_eo.xdf','streamtype','EEG');
% EEG_ssvep = pop_loadxdf('hm_oddball_ssvep_2.xdf','streamtype','EEG');


%% data preprocessing
[EEG_ssvep_prep,EEG_baseline_prep ] = preproc_EEG_hm(EEG_ssvep,EEG_baseline);
% EEG_ssvep_prep = preproc_EEG_hm(EEG_ssvep);
% EEG_baseline_prep = preproc_EEG_hm(EEG_baseline);

%% extract first cube without baseline
idx_tar_cube = cellfun(@(x) ~isempty(regexp(x,'Ring\w*','match','ONCE')),{EEG_ssvep_prep.event.type});
idx_cube_2 = cellfun(@(x) ~isempty(regexp(x,'\w*Second cube\w*','match','ONCE')),{EEG_ssvep_prep.event.type});
idx_cube_1 = xor(idx_tar_cube,idx_cube_2);

%% NEW PART USING THE FUNCTION find_trails function
%boolean array where trail starts and ends
bool_start = idx_cube_1; % You can use the second cube as start points
bool_end = cellfun(@(x) ~isempty(regexp(x,'End Trial','match','ONCE')),{EEG_ssvep_prep.event.type});

%boolean arrays where the LEFT/RIGHT or UP/DOWN event markers occured (I assumed
%they were mutually trail exclusive i.e. both do not occur in the same trail)
condition_1 = cellfun(@(x) ~isempty(regexp(x,'Up','match','ONCE')),{EEG_ssvep_prep.event.type}) + ...
    cellfun(@(x) ~isempty(regexp(x,'Bottom','match','ONCE')),{EEG_ssvep_prep.event.type});
condition_2 = cellfun(@(x) ~isempty(regexp(x,'Left','match','ONCE')),{EEG_ssvep_prep.event.type}) + ...
    cellfun(@(x) ~isempty(regexp(x,'Right','match','ONCE')),{EEG_ssvep_prep.event.type});
idx_cube_1 = find_trails(bool_start,bool_end,condition_1); %This is a boolean array denoting trails going Up or Down
idx_cube_2 = find_trails(bool_start,bool_end,condition_2); %This is a boolean array denoting trails going Left or Right

%% Select Events
cube_1 = pop_selectevent(EEG_ssvep_prep,'event',find(idx_cube_1),'deleteevents','on');
cube_2 = pop_selectevent(EEG_ssvep_prep,'event',find(idx_cube_2),'deleteevents','on');

% epoch data
len_epoch = [140 1140];
cube_1 = pop_epoch(cube_1,{},len_epoch/1000,'epochinfo','yes');
cube_2 = pop_epoch(cube_2,{},len_epoch/1000,'epochinfo','yes');

% select channels

%% tar_Ch = {'O1','O2','Cz','P3','P4','CPz','Pz','P8','P7','POz'};
tar_Ch = 'POz';
idx_tar_Ch = ismember({cube_1.chanlocs.labels},tar_Ch);
plt_data = cube_1.data;

% plot psd
figure
[cross_spectra, cross_freqs] = spectopo(reshape(squeeze(plt_data(idx_tar_Ch,:,:)),1,[]), 0, cube_1.srate, 'freqrange', [1 50]);
grid on 
title('Cross Trials PSD')
[spectra, freqs] = spectopo(squeeze(plt_data(idx_tar_Ch,:,:)), 0, cube_1.srate, 'freqrange', [1 50], 'plot', 'off');

plt_freq = freqs(1:find(freqs>50,1)-1);
plt_spectra = spectra(:,1:length(plt_freq));

figure
grid on
hold on
fs = shadedErrorBar(plt_freq,plt_spectra,{@mean @std},'lineprops','r-');
fs.mainLine.LineWidth = 3;
xlabel('Frequency (Hz)')
ylabel('Log Power Spectral Density 10*log_{10} (\muV^2/Hz)')
set(gca,'fontsize',10)
title('Single Trial PSD')

%% time course
figure
grid on
hold on
plot(cube_1.times, mean(plt_data(idx_tar_Ch,:,:),3), 'b-','linewidth',3);
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
set(gca,'fontsize',10)
