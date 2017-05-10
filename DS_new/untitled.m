color = 'brgkcmy';

figure

for ct = 1:4
    plot_rf_summaries(datawn, intersect(id_dir{ct}, datawn.cell_ids), 'fit_color', color(ct), 'clear', false);
    hold on
end

figure
plot_rf_summaries(datawn, intersect(id_dir{2}, datawn.cell_ids), 'fit_color', color(2), 'clear', false);
plot_rf_summaries(datawn, 'ON transient', 'fit_color', color(5), 'clear', false);
plot_rf_summaries(datawn, 'ON sustained', 'fit_color', color(6), 'clear', false);
plot_rf_summaries(datawn, 'OFF transient slow', 'fit_color', color(7), 'clear', false);

figure
plot_rf_summaries(datawn, intersect(id_dir{2}, datawn.cell_ids), 'fit_color', color(2), 'clear', false);
plot_rf_summaries(datawn, 'OFF transient large', 'fit_color', color(3), 'clear', false);

figure
plot_rf_summaries(datawn, 'OFF transient large', 'clear', false, 'foa', 0, 'label', true);


color = 'brgkcmy';


for ct = 1:4
    plot_rf_summaries(datawn, intersect(id_dir{ct}, datawn.cell_ids), 'clear', false, 'foa', 0, 'label', true);
    hold on
end
plot_rf_summaries(datawn, 'ON transient', 'clear', false, 'foa', 0, 'label', true);
plot_rf_summaries(datawn, 'ON sustained', 'clear', false, 'foa', 0, 'label', true);
plot_rf_summaries(datawn, 'OFF transient slow', 'clear', false, 'foa', 0, 'label', true);
%%
id = [2628 4970 3336 4097 4623 5088 3556 4488];
figure
for i = 1:length(id)
    plot_rf_summaries(datawn, id(i), 'clear', false, 'fit_color', color(mod(i,7)+1))
end
plot(array_location_display(:,1)/10, array_location_display(:,2)/10, 'k')

%%
figure
errorbar(ctr_x, mean(pd_ds{1}), std(pd_ds{1})/sqrt(size(pd_ds{1}, 1)), 'k');
hold on
errorbar(ctr_x, mean(pd_ds{2}), std(pd_ds{2})/sqrt(size(pd_ds{2}, 1)), 'k--');
errorbar(ctr_x, mean(pd_ds{3}), std(pd_ds{3})/sqrt(size(pd_ds{3}, 1)), 'k:');

% set(gca, 'Xscale', 'log')
% xlabel('% contrast')
% ylabel('spike count')
% xlim([3 500])
% 
% figure
errorbar(ctr_x, mean(pd_ct_spikes{1}{2}), std(pd_ct_spikes{1}{2})/sqrt(size(pd_ct_spikes{1}{2}, 1)), 'r');
errorbar(ctr_x, mean(pd_ct_spikes{2}{2}), std(pd_ct_spikes{2}{2})/sqrt(size(pd_ct_spikes{2}{2}, 1)), 'r--');
errorbar(ctr_x, mean(pd_ct_spikes{3}{2}), std(pd_ct_spikes{3}{2})/sqrt(size(pd_ct_spikes{3}{2}, 1)), 'r:');

set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike count')
xlim([3 500])
legend('DS control',  'DS hex',  'DS wash','ON 2 control','ON 2 hex', 'ON 2 wash');
