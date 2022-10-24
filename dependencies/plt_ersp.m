function fig = plt_ersp(plt_EEG,tarCh,lock_name,ev_name,smoothing)
%% plot ERPImage
plt_EEG = pop_select(plt_EEG,'channel',{tarCh});
% thres_amp = 100;
% [plt_EEG, rm_idx] = pop_eegthresh(plt_EEG,1,1,-thres_amp,thres_amp,plt_EEG.xmin,plt_EEG.xmax,0,0);
% plt_EEG = pop_rejepoch(plt_EEG, rm_idx, 0);
idx_cir_1 = cellfun(@(x) ~isempty(regexp(x,'Standard','once')),{plt_EEG.event.type});
idx_cir_2 = cellfun(@(x) ~isempty(regexp(x,'\(L\)','once')),{plt_EEG.event.type});
idx_cir = idx_cir_1 | idx_cir_2;
idx_tri_1 = cellfun(@(x) ~isempty(regexp(x,'Deviant','once')),{plt_EEG.event.type});
idx_tri_2 = cellfun(@(x) ~isempty(regexp(x,'\(R\)','once')),{plt_EEG.event.type});
idx_tri = idx_tri_1 | idx_tri_2;
switch lock_name
    case 'stim'
        if strcmp(ev_name,'cir')
            sort_name = {'circle_gip_start'};
        else
            sort_name = {'triangle_gip_start'};
        end
    case 'gip'
        if strcmp(ev_name,'cir')
            sort_name = {plt_EEG.event(idx_cir).type};
        else
            sort_name = {plt_EEG.event(idx_tri).type};
        end
end
% fig = figure('units','normalized','outerposition',[0 0 1 1]);
fig = figure;
pop_erpimage(plt_EEG,1, [1],[[]],tarCh,smoothing,1,sort_name,[],'latency','yerplabel','\muV','erp','on','cbar','on','topo', { [1] plt_EEG.chanlocs plt_EEG.chaninfo } );
end