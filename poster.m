[raster_temp, ~, ~] = deal(cell(4, 1));
for i = 1:4    
    [NumSpikesCell, StimComb] = get_spikescellstim(datarun{i},19,0);
    nDS{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster_temp{i} = get_ds_raster(datarun{i}, 19);
end

for cc = 1:1 %length(id)
    plot_ds_raster(nDS, raster_temp, cc, 47, 'title', 2, 2, 0)
end

figure
x = 0:45:315;
plot(x, dstuning)
hold on
plot(x, ndstuning, 'r')
xlim([0 315])

%% ds frequency doubling
% on-off example: 2014-07-07-0 id:5731 idx:36 file:5731_5 
% on example: 2014-07-07-0 id:3436 idx:14 file:3436_5

XX = 0.2:0.4:7.8;
FigHandle = figure;
set(FigHandle, 'Position', [1 1 500 400])

subplot(2, 2, 1)
plot_raster(squeeze(raster{4}{36}(1, 5, p_idx{4}(36), :)), 0, 8)
hold on
h = hist(raster_p_sum{4}{36}{5}, XX);
h = h/max(h)*3;
plot(XX, h+4)
ylim([0 7.5])
xlabel('seconds')

subplot(2, 2, 2)
plot_raster(squeeze(raster{4}{36}(1, 5, n_idx{4}(36), :)), 0, 8)
hold on
h = hist(raster_n_sum{4}{36}{5}, XX);
h = h/max(h)*3;
plot(XX, h+4)
ylim([0 7.5])

subplot(2, 2, 3)
plot_raster(squeeze(raster{2}{36}(1, 5, p_idx{2}(36), :)), 0, 8)
hold on
h = hist(raster_p_sum{2}{36}{5}, XX);
h = h/max(h)*3;
plot(XX, h+4)
ylim([0 7.5])

subplot(2, 2, 4)
plot_raster(squeeze(raster{2}{36}(1, 5, n_idx{2}(36), :)), 0, 8)
hold on
h = hist(raster_n_sum{2}{36}{5}, XX);
h = h/max(h)*3;
plot(XX, h+4)
ylim([0 7.5])

% on
XX = 0.2:0.4:7.8;
FigHandle = figure;
set(FigHandle, 'Position', [1 1 500 400])

subplot(2, 2, 1)
plot_raster(squeeze(raster{4}{14}(1, 5, p_idx{4}(14), :)), 0, 8)
hold on
h = hist(raster_p_sum{4}{14}{5}, XX);
h = h/max(h)*3;
plot(XX, h+4)
ylim([0 7.5])
xlabel('seconds')

subplot(2, 2, 2)
plot_raster(squeeze(raster{4}{14}(1, 5, n_idx{4}(14), :)), 0, 8)
hold on
h = hist(raster_n_sum{4}{14}{5}, XX);
h = h/max(h)*3;
plot(XX, h+4)
ylim([0 7.5])

subplot(2, 2, 3)
plot_raster(squeeze(raster{2}{14}(1, 5, p_idx{2}(14), :)), 0, 8)
hold on
h = hist(raster_p_sum{2}{14}{5}, XX);
h = h/max(h)*3;
plot(XX, h+4)
ylim([0 7.5])

subplot(2, 2, 4)
plot_raster(squeeze(raster{2}{14}(1, 5, n_idx{2}(14), :)), 0, 8)
hold on
h = hist(raster_n_sum{2}{14}{5}, XX);
h = h/max(h)*3;
plot(XX, h+4)
ylim([0 7.5])

%% response of example DSGC at NDF 4
% 2014-07-07-0 id:4832 idx:28 file:4832_3
for cc = 28:28
    plot_ds_raster(DS, raster, cc, id(cc), ll, 2, 2, 0)
end

x = 0:45:360;
tuning = repmat(DS{1}.rho{6}(28, :), 1, 2);
tuning = tuning(4:12);
figure
plot(x, tuning, 'k')
ylim([0 1])
xlim([0 360])
xlabel('degree')
ylabel('normalized spike number')
