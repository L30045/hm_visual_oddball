%% This script epoch baseline data and visualize SSVEP
function ssvep_struct = vis_baseline(EEG,varargin)
if isempty(varargin)
    tarCh = {'O1','O2','Oz','POz','PO4','PO3','PO7','PO8'};
else
    tarCh = varargin{1};
end

%% load data
% filepath = '/home/yuan/Documents/2021 HM_visual_oddball/preproc_data/';
% filename = '1145_ssvep_baseline_halo.set';
% EEG = pop_loadset([filepath,filename]);

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
                    
%% calculate PSD
[spec8, freqs] = spectopo(ssvep8.data,0,ssvep8.srate,'plot','off');
[spec9, ~] = spectopo(ssvep9.data,0,ssvep9.srate,'plot','off');
[spec10, ~] = spectopo(ssvep10.data,0,ssvep10.srate,'plot','off');
[spec11, ~] = spectopo(ssvep11.data,0,ssvep11.srate,'plot','off');
plt_f = 0:20;
plt8 = mean(spec8(ismember({ssvep8.chanlocs.labels},tarCh),plt_f+1),1);
plt9 = mean(spec9(ismember({ssvep8.chanlocs.labels},tarCh),plt_f+1),1);
plt10 = mean(spec10(ismember({ssvep8.chanlocs.labels},tarCh),plt_f+1),1);
plt11 = mean(spec11(ismember({ssvep8.chanlocs.labels},tarCh),plt_f+1),1);

%% plotting
figure
plot(plt_f, plt8,'b-o','linewidth',3,'DisplayName','8Hz','markersize',10);
hold on
grid on
plot(plt_f, plt9,'r-o','linewidth',3,'DisplayName','9Hz','markersize',10);
plot(plt_f, plt10,'g-o','linewidth',3,'DisplayName','10Hz','markersize',10);
plot(plt_f, plt11,'m-o','linewidth',3,'DisplayName','11Hz','markersize',10);
legend
xlabel('Frequency(Hz)')
ylabel('Power (\muV^2)')
set(gca,'fontsize',20)
set(gcf,'color','w')
if ~isempty(regexp(EEG.filename,'amp','ONCE'))
    title('Amp')
else
    title('Halo')
end

%% output epoch structure
ssvep_struct = struct('spec8',spec8,'spec9',spec9,'spec10',spec10,'spec11',spec11,'freqs',freqs);
ssvep_struct.chanlocs = {ssvep8.chanlocs.labels};

end