function plot_mfs_psth(field_width, field_height, bin_size, histg, cell_idx, fig_n)


H = figure(fig_n);
set(H, 'Position', [1 1 1080 1080])

h = subplot(field_height, field_width, 1); p1 = get(h, 'pos');
h = subplot(field_height, field_width, 2); p2 = get(h, 'pos');
width = (p2(1) - p1(1))*0.8;

h = subplot(field_height, field_width, 1); p1 = get(h, 'pos');
h = subplot(field_height, field_width, field_width+1); p2 = get(h, 'pos');
height = p1(2) - p2(2);

xx = [bin_size/2:bin_size:2-bin_size/2];

for p = 1:field_width*field_height
    h = subplot(field_width, field_height, p);
    l = get(h, 'pos');
    l(3) = width; l(4) = height;
    set(h, 'pos', l);
    n = 1/bin_size;
    plot(xx(1:n), histg{cell_idx}{p}(1:n), 'b')
    hold on
    plot(xx(n+1:end), histg{cell_idx}{p}(n+1:end), 'r')
    ylim([0 max(max(cell2mat(histg{cell_idx}')))])
    axis off
end
