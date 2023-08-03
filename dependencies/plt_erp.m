function fig = plt_erp(cir_epoch, tri_epoch, tarCh, ev_name, shaded_method, is_correct)
ch_idx = ismember({cir_epoch.chanlocs.labels},tarCh);
plt_t = cir_epoch.times;
cir_data = squeeze(cir_epoch.data(ch_idx,:,:))';
tri_data = squeeze(tri_epoch.data(ch_idx,:,:))';

% further remove trial based on high variance
var_thres = 0.95;
var_dist_cir = reshape(var(cir_data,[],2),[],1);
var_dist_tri = reshape(var(tri_data,[],2),[],1);
cir_data(var_dist_cir > quantile(var_dist_cir,var_thres),:) = [];
tri_data(var_dist_tri > quantile(var_dist_tri,var_thres),:) = [];

% [~,p] = ttest2(cir_data,tri_data);
% [~,plt_p] = bonf_holm(p);
% plt_p = p<=0.05/sqrt(size(cir_data,2));
p = rowfun(@ranksum,table(cir_data',tri_data'));
p = p{:,:}';
if is_correct
    time_test = [200,600]; % ms, time period for statistical test
    idx_test = plt_t >= time_test(1) & plt_t<= time_test(2);
    [~, h] = fdr(p(idx_test), 0.05, 'nonparametric');
    plt_p = false(size(p));
    plt_p(idx_test) = h;
else
    plt_p = p <= 0.05;
end

% plot ERP
fig = figure('units','normalized','outerposition',[0 0 1 1]);
% >> Triangle
ht = shadedErrorBar(plt_t, tri_data, shaded_method,'lineprops',...
    {'color','b','linewidth',3,'DisplayName',sprintf('Standard (%d)',size(tri_data,1))});
ht.patch.FaceAlpha = 0.1;
grid on
hold on
% >> Circle
% disp(size(cir_data))
hc = shadedErrorBar(plt_t, cir_data, shaded_method,'lineprops',...
    {'color','r','linewidth',3,'DisplayName',sprintf('Deviant (%d)',size(cir_data,1))});
hc.patch.FaceAlpha = 0.1;
% sig
plot(plt_t(plt_p), shaded_method{1}(cir_data(:,plt_p)), 'kx','linewidth',3,'markersize',15)

xline(0,'k-','DisplayName',ev_name,'linewidth',3)
legend(findobj(gca,'-regexp','DisplayName', '[^'']'),'location','northwest')
set(gca,'fontsize',30)
set(gcf,'color',[1 1 1])
set(gca,'xtick',round(plt_t(1):100:plt_t(end)))
set(gca,'XTickLabelRotation',-45)
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
end