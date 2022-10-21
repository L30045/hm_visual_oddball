%% New file (Biosemi)
filepath = '/home/yuan/Documents/2021 HM_visual_oddball/dataset/preproc_data/';
loadpath = '/home/yuan/Documents/2021 HM_visual_oddball/dataset/oddball/';
savepath = '/home/yuan/Documents/2021 HM_visual_oddball/dataset/new epoch/'; % Correct the synchronization issue in epoch_ez and create new epoch.

filename_list = reshape({dir([loadpath,'*Oddball*']).name},2,[])';
inner_list = filename_list(:,1);
outer_list = filename_list(:,2);
inner_data = cellfun(@(x) [x(1:end-4),'_resample_250Hz.set'],filename_list(:,1),'uniformoutput',0);
outer_data = cellfun(@(x) [x(1:end-4),'_resample_250Hz.set'],filename_list(:,2),'uniformoutput',0);

parfor i = 1:size(filename_list,1)
    EEG_noHm = pop_loadset([filepath,inner_data{i}]);
    EEG_Hm = pop_loadset([filepath,outer_data{i}]);
    [~, ~, ~, ~, ~, epoch_struct_noHm, ~]...
    = epoch_ez([loadpath,inner_list{i}], EEG_noHm);
    [~, ~, ~, ~, ~, epoch_struct_Hm, ~]...
    = epoch_ez([loadpath,outer_list{i}], EEG_Hm);
    parsave([savepath,sprintf('rmPreStim_new_s%02d_epoch_%s.mat',i,inner_list{i}(1:4))],epoch_struct_noHm,epoch_struct_Hm);
end

disp('Done')

%% Old file (Smarting)
filepath = '/home/yuan/Documents/2021 HM_visual_oddball/dataset/preproc_data/';
loadpath = '/home/yuan/Documents/2021 HM_visual_oddball/dataset/';
savepath = '/home/yuan/Documents/2021 HM_visual_oddball/dataset/new epoch/';

subj_list = [1,4,6,8:10];
parfor j = 1:length(subj_list)
    subj_i = subj_list(j);
    EEG_noHm = pop_loadset([filepath,sprintf('s%02d_Inner_resampled_250.set',subj_i)]);
    EEG_Hm = pop_loadset([filepath,sprintf('s%02d_Outer_resampled_250.set',subj_i)]);
    [~, ~, ~, ~, ~, epoch_struct_noHm, ~]...
    = epoch_ez([loadpath,sprintf('hm_visual_oddball_s%02d_cond1.xdf',subj_i)], EEG_noHm);
    [~, ~, ~, ~, ~, epoch_struct_Hm, ~]...
    = epoch_ez([loadpath,sprintf('hm_visual_oddball_s%02d_cond2.xdf',subj_i)], EEG_Hm);
    parsave([savepath,sprintf('rmPreStim_s%02d_epoch.mat',subj_i)],epoch_struct_noHm,epoch_struct_Hm);
end
disp('Done')