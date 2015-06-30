function plot_single_cell_speed_tuning(datadg, MAG_all_norm_dg, idx, ll)

f = 4; % 1 pixel = 4 micron;
n = length(idx);
a = ceil(sqrt(n)); b = ceil(n/a);
figure
for i = 1:n
    v = datadg{idx(i)}.stimulus.params.SPATIAL_PERIOD./datadg{idx(i)}.stimulus.params.TEMPORAL_PERIOD*f;
    subplot(b, a, i)
    semilogx(v, exciseColumn(MAG_all_norm_dg{idx(i)}), 'b')
    xlabel('micron/second')
    ylabel('Response')
    title(ll{idx(i)})
    xlim([v(end) v(1)])
end
end

