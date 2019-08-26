%% load data
rf0829 = analyzeRF160829;
rf1115 = analyzeRF161115;
[rf1122, histg1122, rf_mean1122, rf_ste1122] = analyzeRF161122;

%% fig 2A

field_width = 17; field_height = 17; bin_size = 0.02;

h = subplot(field_height, field_width, 1); p1 = get(h, 'pos');
h = subplot(field_height, field_width, 2); p2 = get(h, 'pos');
width = (p2(1) - p1(1))*0.8;

h = subplot(field_height, field_width, 1); p1 = get(h, 'pos');
h = subplot(field_height, field_width, field_width+1); p2 = get(h, 'pos');
height = p1(2) - p2(2);

cell_id = 80;
cell_idx = 2;
ll = 3;

H = figure(1);
set(H, 'Position', [1 1 1080 1080])
for p = 1:289
    h = subplot(field_width, field_height, p);
    l = get(h, 'pos');
    l(3) = width; l(4) = height;
    set(h, 'pos', l);
    n = 1/bin_size;
    plot(xx(1:n), histg1122{ll}{cell_idx}{p}(1:n), 'b')
    hold on
    plot(xx(n+1:end), histg1122{ll}{cell_idx}{p}(n+1:end), 'r')
    ylim([0 max(max(cell2mat(histg1122{ll}{cell_idx}')))])
    axis off
end

%% fig 2B-2C, same as fig 2D 20 R*/rod/s

%% fig 2D
figure(2)
set(gcf, 'Position', [1 1 400 600])
for i = 1:3
    rf = rf1122{i*2-1}{cell_idx};
    subplot(3,2,i*2-1)
    imagesc(rf(:,:,1))
    colormap gray
    axis image
    axis off

    subplot(3,2,i*2)
    imagesc(rf(:,:,2))
    colormap gray
    axis image
    axis off
end

%% fig 2E
cell_id = 618;
cell_idx = 7;
figure(3)
set(gcf, 'Position', [1 1 400 600])
for i = 1:3
    rf = rf0829{i}{cell_idx};
    subplot(3,2,i*2-1)
    imagesc(rf(:,:,1))
    colormap gray
    axis image
    axis off

    subplot(3,2,i*2)
    imagesc(rf(:,:,2))
    colormap gray
    axis image
    axis off
end

%% fig 2G
ONOFF = {'ON', 'OFF'};
dir = 1; % superior cells

figure
for onoff = 1:2
    errorbar(0:4, rf_mean1122{onoff}(:, dir), rf_ste1122{onoff}(:, dir), 'color', color(onoff))
    hold on
    ylim([0 0.14])
    xlim([-1 5])
    xlabel('log(background intensity)')
    ylabel('RF area (mm^2)')
    legend(ONOFF)
end