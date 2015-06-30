function plot_classified_speed_tuning(datadg, MAG_all_norm_dg, idx, idx_sub, ll)

f = 4; % 1 pixel = 4 micron;
n = length(idx);
n_cluster = length(idx_sub);
a = ceil(sqrt(n)); b = ceil(n/a);
color = 'brgk';
figure
for i = 1:n
    v = datadg{idx(i)}.stimulus.params.SPATIAL_PERIOD./datadg{idx(i)}.stimulus.params.TEMPORAL_PERIOD*f;
    subplot(b, a, i)
    for j = 1:n_cluster
        semilogx(v, exciseColumn(MAG_all_norm_dg{idx(i)}(:, idx_sub{j})), color(j))
        hold on
    end
    xlabel('micron/second')
    ylabel('Response')
    title(ll{idx(i)})
    xlim([v(end) v(1)])
end
end

