function fig = plt_erp_diff(cir_epoch, tri_epoch, tarCh, ev_name)
ch_idx = ismember({cir_epoch.chanlocs.labels},tarCh);
plt_t = cir_epoch.times;
cir_data = squeeze(cir_epoch.data(ch_idx,:,:));
tri_data = squeeze(tri_epoch.data(ch_idx,:,:));
% plot ERP
fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(plt_t, mean(cir_data,2) - mean(tri_data,2),'b-','linewidth',3,'displayname','Dev - Std')
hold on
grid on
xline(0,'k--','DisplayName',ev_name,'linewidth',3)
yline(0,'k-','linewidth',3)
legend(findobj(gca,'-regexp','DisplayName', '[^'']'),'location','northwest')
set(gca,'fontsize',30)
set(gcf,'color',[1 1 1])
set(gca,'xtick',round(plt_t(1):100:plt_t(end)))
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
end