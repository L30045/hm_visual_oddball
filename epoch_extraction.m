filepath = '/home/yuan/Documents/2021 HM_visual_oddball/dataset/preproc_data/';
loadpath = '/home/yuan/Documents/2021 HM_visual_oddball/dataset/oddball/';
savepath = '/home/yuan/Documents/2021 HM_visual_oddball/dataset/new epoch/'; % Correct the synchronization issue in epoch_ez and create new epoch.

filename_list = reshape({dir([filepath,'*Oddball*']).name},2,[])';
inner_list = filename_list(:,1);
outer_list = filename_list(:,2);

parfor i = 1:size(filename_list,1)
    EEG_noHm = pop_loadset([filepath,inner_list{i}]);
    EEG_Hm = pop_loadset([filepath,outer_list{i}]);
    [~, ~, ~, ~, ~, epoch_struct_noHm, ~]...
    = epoch_ez([loadpath,sprintf('%s.xdf',inner_list{i}(1:end-4))], EEG_noHm);
    [~, ~, ~, ~, ~, epoch_struct_Hm, ~]...
    = epoch_ez([loadpath,sprintf('%s.xdf',outer_list{i}(1:end-4))], EEG_Hm);
    parsave([savepath,sprintf('new_s%02d_epoch_%s.mat',i,inner_list{i}(1:4))],epoch_struct_noHm,epoch_struct_Hm);
end

disp('Done')