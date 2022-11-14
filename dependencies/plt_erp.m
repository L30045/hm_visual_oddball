function fig = plt_erp(cir_epoch, tri_epoch, tarCh, ev_name, shaded_method)
ch_idx = ismember({cir_epoch.chanlocs.labels},tarCh);
plt_t = cir_epoch.times;
cir_data = squeeze(cir_epoch.data(ch_idx,:,:))';
tri_data = squeeze(tri_epoch.data(ch_idx,:,:))';

[~,p] = ttest2(cir_data,tri_data);
[~,plt_p] = bonf_holm(p);
% plt_p = p<=0.05/sqrt(size(cir_data,2));
% plt_p = p<=0.05;

% plot ERP
fig = figure('units','normalized','outerposition',[0 0 1 1]);
% >> Triangle
ht = shadedErrorBar(plt_t, tri_data, shaded_method,'lineprops',...
    {'color','b','linewidth',3,'DisplayName','Standard'});
ht.patch.FaceAlpha = 0.1;
grid on
hold on
% >> Circle
hc = shadedErrorBar(plt_t, cir_data, shaded_method,'lineprops',...
    {'color','r','linewidth',3,'DisplayName','Deviant'});
hc.patch.FaceAlpha = 0.1;
% sig
plot(plt_t(plt_p), mean(cir_data(:,plt_p),1), 'ro','linewidth',3,'markersize',15)

xline(0,'k-','DisplayName',ev_name,'linewidth',3)
legend(findobj(gca,'-regexp','DisplayName', '[^'']'),'location','northwest')
set(gca,'fontsize',30)
set(gcf,'color',[1 1 1])
set(gca,'xtick',round(plt_t(1):100:plt_t(end)))
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
end