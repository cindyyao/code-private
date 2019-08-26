cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

datadg = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-03-16-0/data004-sorted/data004-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2018-03-16-0/stimuli/s04.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

params_idx = [3 2]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

load('DS180316.mat')

% datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-03-16-0/datawn-map/datawn-map', opt);
% datarun = load_ei(datarun, id_dir{1}, 'array_id',1501);
datawn{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-03-16-0/data001-map/data001-map', opt);
datawn{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-03-16-0/data003-map/data003-map', opt);
datawn{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-03-16-0/data006-map/data006-map', opt);
datawn{3} = load_ei(datawn{3}, id_dir{1}, 'array_id',1501);

datawn_alpha{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-03-16-0/data001-map-from-006/data001-map-from-006', opt);
datawn_alpha{1} = load_sta(datawn_alpha{1});

datawn_alpha{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-03-16-0/data003-map-from-006/data003-map-from-006', opt);
datawn_alpha{2} = load_sta(datawn_alpha{2});

datawn_alpha{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-03-16-0/data006/data006', opt);
datawn_alpha{3} = load_sta(datawn_alpha{3});
%% dg
% ds_id = datadg.cell_ids;
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);

delta_p = 1; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(DG{i});
rep = datadg.stimulus.repetitions;
%% classification based on speed tunning
L = 1;
mag_pca = MAG_all_norm_dg{L};
mag_pca = mag_pca';
[id_sub, idx_sub] = deal(cell(2, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(mag_pca);
pc1 = 1; pc2 = 3;
plot(scores(:, pc1), scores(:, pc2), 'o')
hold on
for i = 1:2
    [x, y] = ginput;
    plot(x, y)
    IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
    [~, idx_sub{i}] = find(IN' == 1);
    id_sub{i} = ds_id(idx_sub{i});
end
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
title('NDF 0')

figure
plot(scores(idx_sub{1}, pc1), scores(idx_sub{1}, pc2), 'ro', scores(idx_sub{2}, pc1), scores(idx_sub{2}, pc2), 'bo')
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
title('NDF 0')

v = 240*datadg.stimulus.params.SPATIAL_PERIOD./datadg.stimulus.params.TEMPORAL_PERIOD;
figure
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
xlim([v(end) v(1)])

t = 2;
figure
compass(DG{1}.U{t}(idx_sub{1}), DG{1}.V{t}(idx_sub{1}), 'r')
hold on
compass(DG{1}.U{t}(idx_sub{2}), DG{1}.V{t}(idx_sub{2}), 'b')

%%
d = 1;
t = 1;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(DG{d}.U{t}(idx_sub{2}), DG{d}.V{t}(idx_sub{2}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}(idx_sub{2}), DG{d}.V{t}(idx_sub{2}), x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = idx_sub{2}(I);
    id_dir{i} = ds_id(idx_dir{i});
end

d = 1;
t = 1;
h = figure;
dirn = 3;
set(h, 'Position', [1 1 1080 500])
compass(DG{d}.U{t}(idx_sub{1}), DG{d}.V{t}(idx_sub{1}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}(idx_sub{1}), DG{d}.V{t}(idx_sub{1}), x, y);
    [~, I] = find(IN == 1);
    idx_dir_on{i} = idx_sub{1}(I);
    id_dir_on{i} = ds_id(idx_dir_on{i});
end


%% plot cell summary
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [4 5]);

for cc = 1:length(ds_id)
    plot_ds_raster(DG_cut, raster_dg_cut, cc, ds_id(cc), 1, 1, 1, 1)
end


%% neighboring pairs
pos = datawn{3}.ei.position;
mode = 'neg';
neighbors = [];
ct = 1;
id_wn = id_dir{ct};
for cc1 = 1:length(id_wn)
    for cc2 = cc1+1:length(id_wn)
        id1 = id_wn(cc1);
        idx1 = get_cell_indices(datawn{3}, id1);
        ei1 = datawn{3}.ei.eis{idx1};
        com1 = ei_com_xy(ei1, pos, 30*3, mode);
        id2 = id_wn(cc2);
        idx2 = get_cell_indices(datawn{3}, id2);
        ei2 = datawn{3}.ei.eis{idx2};
        com2 = ei_com_xy(ei2, pos, 30*3, mode);
        if pdist([com1;com2]) < 150
            neighbors = [neighbors; id1 id2];
        end
    end
end

cn = 0;
ct = 1;
celltype = {'superior', 'anterior', 'inferior', 'posterior'};
coms = [];
for cc = 1:length(id_wn)
    id = id_wn(cc);
    idx = get_cell_indices(datawn{3}, id);
    ei = datawn{3}.ei.eis{idx};
    com = ei_com_xy(ei, pos, 30*3, mode);
    coms = [coms; com];
end

corner_i = [4 126 195 264 386 455 4];
corner_position = datawn{3}.ei.position(corner_i, :);
figure
for cc = 1:length(id_wn)
    plot(coms(cc, 1), coms(cc, 2),'ko')
    hold on
    text(coms(cc, 1)+5, coms(cc, 2)+5, num2str(id_wn(cc)), 'FontSize', 10)
    
end

plot(corner_position(:, 1), corner_position(:, 2), 'color', [.5 .5 .5])
axis off
title(celltype{ct})

% ct = 1;
cp_i = [];
for c1 = 1:length(id_wn)-1
    for c2 = c1+1:length(id_wn)
        if norm([coms(c1, :) - coms(c2, :)]) < 150
            cn = cn + 1;
            cp_i = [cp_i; c1 c2];
        end
    end
end

%%
duration = 2400;
bin_size = 0.00025;
max_lag =40;
ct = 1;
N = 10000;
LL = {'NDF 4', 'NDF 2', 'NDF 0'};
firing_rate = zeros(length(id_dir{ct}), length(LL));
for cc = 1:length(id_dir{ct})
    for ll = 1:3
        idx = get_cell_indices(datawn{ll}, id_dir{ct}(cc));
        spikes = datawn{ll}.spikes{idx};
        firing_rate(cc, ll) = length(spikes)/duration;
    end
end

figure(1)
subplot(2,3,1)
plot(firing_rate')
xlabel('light level')
ylabel('Hz')
title('mean firing rate')


firing_rate_mean = mean(firing_rate);
firing_rate_ste = std(firing_rate)/sqrt(size(firing_rate, 1));
subplot(2,3,4)
errorbar(1:3,firing_rate_mean, firing_rate_ste)
ylim([0 5])
xlabel('light level')
ylabel('Hz')
title('mean firing rate')

A_all = []; A_all_coeff = [];

for cp =21:21 %size(neighbors, 1) 
    xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
    FigHandle = figure;
    set(FigHandle, 'Position', [1 1 300 800])
    for ll = 1:3
    %     id1 = corr_cells(cp, 1);
    %     id2 = corr_cells(cp, 2);

        id1 = neighbors(cp, 1);
        id2 = neighbors(cp, 2);

        idx1 = get_cell_indices(datawn{ll}, id1);
        idx2 = get_cell_indices(datawn{ll}, id2);
        spikes1 = datawn{ll}.spikes{idx1};
        spikes1_TF= ceil(spikes1/bin_size);
        spikes1 = zeros(duration/bin_size, 1);
        spikes1(spikes1_TF) = 1;

        spikes2 = datawn{ll}.spikes{idx2};
        spikes2_TF= ceil(spikes2/bin_size);
        spikes2 = zeros(duration/bin_size, 1);
        spikes2(spikes2_TF) = 1;

        A = xcorr(spikes1, spikes2, max_lag)/duration;
        A_all(cp, :, ll) = A;
        
        A = xcorr(spikes1, spikes2, max_lag, 'coeff');
        A_all_coeff(cp, :, ll) = A;
        
        subplot(3,2,ll*2-1)
        bar(xx, A, 'k')
        title([num2str(id1) '  ' num2str(id2) ' ' LL{ll} ])
        xlim([-0.01 0.01])
    end
%     print_close(1, [12,4], [num2str(id1) '  ' num2str(id2)])
end

A_all_coeff_clean = A_all_coeff;
A_all_coeff_clean([5 20], :, :) = [];
t2p_coeff = squeeze(max(A_all_coeff_clean, [], 2)) - squeeze(A_all_coeff_clean(:, 41, :));

subplot(2,3,2)
plot(t2p_coeff')
xlabel('light level')
ylabel('Hz')
title('correlation coefficient')

t2p_coeff_mean = mean(t2p_coeff);
t2p_coeff_ste = std(t2p_coeff)/sqrt(size(t2p_coeff,1));
subplot(2,3,5)
errorbar(1:3, t2p_coeff_mean, t2p_coeff_ste);
xlabel('light level')
ylabel('Hz')
title('correlation coefficient')

A_all_clean = A_all;
% A_all_clean([15 17], :, :) = [];
t2p = squeeze(max(A_all_clean, [], 2)) - squeeze(A_all_clean(:, 41, :));

subplot(2,3,3)
plot(t2p')
xlabel('light level')
ylabel('Hz')
title('correlated firing rate')

t2p_mean = mean(t2p);
t2p_ste = std(t2p)/sqrt(size(t2p,1));
subplot(2,3,6)
errorbar(1:3, t2p_mean, t2p_ste);
xlabel('light level')
ylabel('Hz')
title('correlated firing rate')

% paper figure
figure
for cc = 1:size(t2p_coeff, 1)
    plot(1:3, t2p_coeff(cc, :), 'ko')
    hold on
end
t2p_coeff_mean = mean(t2p_coeff);
t2p_coeff_ste = std(t2p_coeff)/sqrt(size(t2p_coeff,1));
errorbar(1:3, t2p_coeff_mean, t2p_coeff_ste, 'rd-', 'MarkerSize', 10);
xlabel('light level')
ylabel('Hz')
title('correlation coefficient')


% -5 to 5 ms above baseline
A_all_coeff_clean = A_all_coeff;
base = repmat(A_all_coeff_clean(:,41,:),1,41,1);

double_peaks = A_all_coeff_clean(:,21:61,:) - base;
double_peaks = squeeze(sum(double_peaks, 2));
% double_peaks = double_peaks(~any(double_peaks'<0), :);
% double_peaks = double_peaks(~any(double_peaks'>0.4), :);

figure
for cc = 1:size(double_peaks, 1)
    plot(1:3, double_peaks(cc, :), 'color', [1 1 1]*.5, 'Marker', 'o')
    hold on
end
double_peaks_mean = mean(double_peaks);
double_peaks_ste = std(double_peaks)/sqrt(size(double_peaks,1));
errorbar(1:3, double_peaks_mean, double_peaks_ste, 'rd-', 'MarkerSize', 10);
xlabel('light level')
ylabel('area under CCF')


%%
duration = 2400;
bin_size = 0.00025;
max_lag = 40;
ids = get_cell_ids(datawn_alpha{3}, 'OFF type3');
for cc1 = 2:2%length(ids)-1
    xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
    FigHandle = figure;
    set(FigHandle, 'Position', [1 1 2000 1500])
    for cc2 = 1:length(ids)
        id1 = ids(cc1);
        id2 = ids(cc2);

        idx1 = get_cell_indices(datawn_alpha{3}, id1);
        idx2 = get_cell_indices(datawn_alpha{3}, id2);
        spikes1 = datawn_alpha{3}.spikes{idx1};
        spikes1_TF= ceil(spikes1/bin_size);
        spikes1 = zeros(duration/bin_size, 1);
        spikes1(spikes1_TF) = 1;

        spikes2 = datawn_alpha{3}.spikes{idx2};
        spikes2_TF= ceil(spikes2/bin_size);
        spikes2 = zeros(duration/bin_size, 1);
        spikes2(spikes2_TF) = 1;
        A = xcorr(spikes1, spikes2, max_lag)/duration;
        subplot(6,6,cc2)
        bar(xx, A, 'k')
        title([num2str(id1) '  ' num2str(id2)])
        xlim([-0.01 0.01])
    end
    print_close(1, [24,12], num2str(id1))
end

%% OFF alpha cells
duration = 2400;
bin_size = 0.00025;
max_lag = 40;
    
for cp = 1:length(neighbors)
    xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
%     FigHandle = fiagure;
%     set(FigHandle, 'Position', [1 1 2000 1500])
    id1 = neighbors(cp, 1);
    id2 = neighbors(cp, 2);
    for ll = 1:length(datawn_alpha)

        idx1 = get_cell_indices(datawn_alpha{ll}, id1);
        idx2 = get_cell_indices(datawn_alpha{ll}, id2);
        spikes1 = datawn_alpha{ll}.spikes{idx1};
        spikes1_TF= ceil(spikes1/bin_size);
        spikes1 = zeros(duration/bin_size, 1);
        spikes1(spikes1_TF) = 1;

        spikes2 = datawn_alpha{ll}.spikes{idx2};
        spikes2_TF= ceil(spikes2/bin_size);
        spikes2 = zeros(duration/bin_size, 1);
        spikes2(spikes2_TF) = 1;

        A = xcorr(spikes1, spikes2, max_lag)/duration;
        A_all(cp, :, ll) = A;

%             subplot(2,3,ll)
%             bar(xx, A, 'k')
%             title('correlated firing rate')
%             xlim([-0.01 0.01])

        A = xcorr(spikes1, spikes2, max_lag, 'coeff');
        A_all_coeff(cp, :, ll) = A;

%             subplot(2,3,ll+3)
%             bar(xx, A, 'k')
%             title('correlation coefficient')
%             xlim([-0.01 0.01])
    end
%     print_close(1, [24,12], [num2str(id1) ' ' num2str(id2)])
end

firing_rate = zeros(length(off_alpha_id), length(datawn_alpha));
for cc = 1:length(off_alpha_id)
    for ll = 1:3
        idx = get_cell_indices(datawn_alpha{ll}, off_alpha_id(cc));
        spikes = datawn_alpha{ll}.spikes{idx};
        firing_rate(cc, ll) = length(spikes)/duration;
    end
end

figure(1)
subplot(2,3,1)
plot(firing_rate')
xlabel('light level')
ylabel('Hz')
title('mean firing rate')


firing_rate_mean = mean(firing_rate);
firing_rate_ste = std(firing_rate)/sqrt(size(firing_rate, 1));
subplot(2,3,4)
errorbar(1:3,firing_rate_mean, firing_rate_ste)
xlabel('light level')
ylabel('Hz')
title('mean firing rate')


A_all_coeff_clean = A_all_coeff;
t2p_coeff = squeeze(max(A_all_coeff_clean, [], 2)) - squeeze(A_all_coeff_clean(:, 41, :));

subplot(2,3,2)
plot(t2p_coeff')
xlabel('light level')
ylabel('Hz')
title('correlation coefficient')

t2p_coeff_mean = mean(t2p_coeff);
t2p_coeff_ste = std(t2p_coeff)/sqrt(size(t2p_coeff,1));
subplot(2,3,5)
errorbar(1:3, t2p_coeff_mean, t2p_coeff_ste);
xlabel('light level')
ylabel('Hz')
title('correlation coefficient')


A_all_clean = A_all;
t2p = squeeze(max(A_all_clean, [], 2)) - squeeze(A_all_clean(:, 41, :));

subplot(2,3,3)
plot(t2p')
xlabel('light level')
ylabel('Hz')
title('correlated firing rate')

t2p_mean = mean(t2p);
t2p_ste = std(t2p)/sqrt(size(t2p,1));
subplot(2,3,6)
errorbar(1:3, t2p_mean, t2p_ste);
xlabel('light level')
ylabel('Hz')
title('correlated firing rate')

%%
figure
for ll = 1:3
    subplot(3,2,ll*2-1)
    bar(xx, A_all_coeff(21, :, ll), 'k')
    xlim([-0.01 0.01])
end
subplot(3,2,[2 4 6])
for i = 1:size(t2p_coeff, 1)
    hold on
    if i == 19
        plot(t2p_coeff(i, :), 'ro-')
    else
        plot(t2p_coeff(i, :), 'ko-')
    end
end
t2p_coeff_mean = mean(t2p_coeff);
t2p_coeff_ste = std(t2p_coeff)/sqrt(size(t2p_coeff,1));
errorbar(1:3, t2p_coeff_mean, t2p_coeff_ste, 'bd-', 'markersize', 10);
xlabel('light level')
ylabel('Hz')
title('correlation coefficient')
xlim([0.5 3.5])

t2p_coeff_0316 = t2p_coeff;
load('DS180601.mat', 't2p_coeff')
t2p_coeff_0601 = t2p_coeff;
t2p_coeff_all = [t2p_coeff_0316; t2p_coeff_0601];
figure
for i = 1:size(t2p_coeff_all, 1)
    hold on
    if i == 19
        plot(t2p_coeff_all(i, :), 'ro-')
    else
        plot(t2p_coeff_all(i, :), 'ko-')
    end
end
t2p_coeff_all_mean = mean(t2p_coeff_all);
t2p_coeff_all_ste = std(t2p_coeff_all)/sqrt(size(t2p_coeff_all,1));
errorbar(1:3, t2p_coeff_all_mean, t2p_coeff_all_ste, 'bd-', 'markersize', 10);
xlabel('light level')
ylabel('Hz')
title('correlation coefficient')
xlim([0.5 3.5])

