cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);


datadg{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data000-map/data000-map', opt);
datadg{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s00.mat';
datadg{1} = load_stim_matlab(datadg{1}, 'user_defined_trigger_interval', 10);
datadg{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data003-map/data003-map', opt);
% datadg{2} = load_data('/Analysis/xyao/2015-07-03-0/data003/data003', opt);
datadg{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s03.mat';
datadg{2} = load_stim_matlab(datadg{2}, 'user_defined_trigger_interval', 10);
datadg{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data006-map/data006-map', opt);
datadg{3}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s06.mat';
datadg{3} = load_stim_matlab(datadg{3}, 'user_defined_trigger_interval', 10);
datadg{4} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data009-map/data009-map', opt);
datadg{4}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s09.mat';
datadg{4} = load_stim_matlab(datadg{4}, 'user_defined_trigger_interval', 10);
datadg{5} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data012-map/data012-map', opt);
datadg{5}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s12.mat';
datadg{5} = load_stim_matlab(datadg{5}, 'user_defined_trigger_interval', 10);

%% fig 1C 
n = 5;
i = 5;
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg{i},datadg{i}.cell_ids,0, 1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg{i});

params_idx = [4 5]; % which parameters to use for classification
[ds_id, nonds_id] = classify_ds(datadg{i}, ds_struct, params_idx);
title('Normalized Vector Sum Amplitude')
xlabel('High speed - 480 micron/s')
ylabel('Low speed - 240 micron/s')

%% drifting grating

load('DS150703-1.mat')
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg{i},ds_id,0,1);
    DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg{i}));
    raster_dg{i} = get_ds_raster(datadg{i}, ds_id);
    for j = 1:length(raster_dg{i})
        if(dg_idx(j, i))
            raster_dg{i}{j} = [];
        end
    end
end

[raster_dg_nds, DG_nds] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg{i},nonds_id,0,1);
    DG_nds{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg{i}));
    raster_dg_nds{i} = get_ds_raster(datadg{i}, nonds_id);
end

delta_p = 4; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

for i = 1:n
    if ismember(i, [2 5])
        [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
    else
        [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p-3});
    end
    MAG_all_norm_dg{i} = normalize_MAG(DG{i});
    rep = datadg{i}.stimulus.repetitions;
end

%% fig 1A
ll = 5; t = 4; cc = 23;
[idx, xx, yy] = subplot_idx(1, 1);
tt = DG{1}.theta{1}(1, :);
FigHandle = figure;
set(FigHandle, 'Position', [1 1 400 400])
h = subplot(xx, yy, idx(1)); 
u_temp = DG{ll}.U{t}(cc);
v_temp = DG{ll}.V{t}(cc);
P = polar(0, 3);
set(P, 'Visible', 'off')
hold on
compass(DG{ll}.U{t}(cc), DG{ll}.V{t}(cc), 'r');
polar(tt, DG{ll}.rho{t}(cc, :), 'b');
polar_theta_off(h)
for i = 2:9
    subplot(xx, yy, idx(i)); plot_raster(squeeze(raster_dg{ll}{cc}(1, t, i-1, :)), 0, 8);
    if idx(i) == 7
        ylabel('trials')
        xlabel('time (s)')
    end
end
subplot(3, 3, 5)
tt = linspace(0, 2*pi, 100);
polar(tt, ones(1, 100), 'k');
polar(tt, 2*ones(1, 100), 'k');
polar(tt, 3*ones(1, 100), 'k');

%% fig 1B
ll = 5; t = 4; cc = 60;
[idx, xx, yy] = subplot_idx(1, 1);
tt = DG_nds{1}.theta{1}(1, :);
FigHandle = figure;
set(FigHandle, 'Position', [1 1 400 400])
h = subplot(xx, yy, idx(1)); 
u_temp = DG_nds{ll}.U{t}(cc);
v_temp = DG_nds{ll}.V{t}(cc);
P = polar(0, 3);
set(P, 'Visible', 'off')
hold on
compass(DG_nds{ll}.U{t}(cc), DG_nds{ll}.V{t}(cc), 'r');
polar(tt, DG_nds{ll}.rho{t}(cc, :), 'b');
polar_theta_off(h)
for i = 2:9
    subplot(xx, yy, idx(i)); plot_raster(squeeze(raster_dg_nds{ll}{cc}(1, t, i-1, :)), 0, 8);
    if idx(i) == 7
        ylabel('trials')
        xlabel('time (s)')
    end
end
subplot(3, 3, 5)
tt = linspace(0, 2*pi, 100);
polar(tt, ones(1, 100), 'k');
polar(tt, 2*ones(1, 100), 'k');
polar(tt, 3*ones(1, 100), 'k');

%% fig 1C inset
ll = 5; t = 5;
ds_vs = DG{ll}.mag{t};
nds_vs = DG_nds{ll}.mag{t};
max_vs = max([ds_vs nds_vs]);
min_vs = min([ds_vs nds_vs]);

X = 0.05:0.1:2.15;
ds_vs_hist = hist(ds_vs, X);
nds_vs_hist = hist(nds_vs, X);
bar(X, [ds_vs_hist' nds_vs_hist'], 1, 'stacked')
xlim([0 2.2])

%% fig 1D
clear X
mag_pca = MAG_all_norm_dg{5};
mag_pca = mag_pca';
[id_sub, idx_sub] = deal(cell(2, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(mag_pca);
pc1 = 1; pc2 = 3;
plot(scores(:, pc1), scores(:, pc2), 'o')
hold on
[x, y] = ginput;
plot(x, y)
IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
[~, idx_init] = find(IN' == 1);
id_init = ds_id(idx_init);

[~, ~, ib] = intersect(id_init, ds_id);
vc = ones(length(ds_id),1);
vc(ib) = 2; %initializing ds cells to cluster 2, everything else cluster 1

X(:,1) = scores(:, pc1);
X(:,2) = scores(:, pc2);
[idx, ~, ~] = clustering_analysis_plots(X, 0,1, 2, 0, 1, 0, 0, 0,0, vc);

oods_id = datarun.cell_ids(idx==2);
onds_id = datarun.cell_ids(idx==1);
% ds_id = id_init;

xlabel('1st Principal Component')
ylabel('3rd Principal Component')
title('NDF 0')
% xlim([-1.5 0.7])
% ylim([-0.7 0.5])

figure
v = 4*datadg{2}.stimulus.params.SPATIAL_PERIOD./datadg{2}.stimulus.params.TEMPORAL_PERIOD;
subplot(1, 2, 1)
semilogx(v, exciseColumn(MAG_all_norm_dg{5}(:, idx_sub{1})), 'r')
xlabel('micron/second')
ylabel('Normalized Response')
xlim([v(end) v(1)])
subplot(1, 2, 2)
semilogx(v, exciseColumn(MAG_all_norm_dg{5}(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Normalized Response')
xlim([v(end) v(1)])

t = 4;
figure
subplot(1, 2, 1)
compass(DG{5}.U{t}(idx_sub{1}), DG{5}.V{t}(idx_sub{1}), 'b')
subplot(1, 2, 2)
compass(DG{5}.U{t}(idx_sub{2}), DG{5}.V{t}(idx_sub{2}), 'r')

%% fig 1G
stixel_size_ = 10;
corner_i = [4 126 195 264 386 455 4];
corner_position = datadg{5}.ei.position(corner_i, :);
center_ds = tforminv(Tform, center_corrected * stixel_size_);
radius = 75;

% on-off DSGC
figure
for ct = 1:4
    subplot(2, 2, ct)
    plot(corner_position(:, 1), corner_position(:, 2))
    hold on
    for i = 1:length(idx_dir{ct})
        [x, y] = circle(center_ds(idx_dir{ct}(i), 1), center_ds(idx_dir{ct}(i), 2), radius);
        plot(x, y, 'k')
        hold on
    end
    xlim([-700 700])
    ylim([-500 500])
    axis off
%     title(oo_ds{ct})
end

% on DSGC
figure
for ct = 1:3
    subplot(2, 2, ct)
    plot(corner_position(:, 1), corner_position(:, 2))
    hold on
    for i = 1:length(idx_dir_on{ct})
        [x, y] = circle(center_ds(idx_dir_on{ct}(i), 1), center_ds(idx_dir_on{ct}(i), 2), radius);
        plot(x, y, 'k')
        hold on
    end
    xlim([-700 700])
    ylim([-500 500])
    axis off
%     title(on_ds{ct})
end


%% fig 4A
load('DS140707.mat')
bin_size = 0.05;
XX = bin_size/2:bin_size:8-bin_size/2;
FigHandle = figure;
set(FigHandle, 'Position', [1 1 500 400])
h_all = [];

subplot(4, 2, 3)
plot_raster(squeeze(raster{4}{36}(1, 5, p_idx{4}(36), :)), 0, 8)
subplot(4, 2, 1)
h = hist(raster_p_sum{4}{36}{5}, XX)/(bin_size*3);
h_all = [h_all h];
plot(XX, h)
xlim([0 8])

subplot(4, 2, 4)
plot_raster(squeeze(raster{4}{36}(1, 5, n_idx{4}(36), :)), 0, 8)
subplot(4, 2, 2)
h = hist(raster_n_sum{4}{36}{5}, XX)/(bin_size*3);
h_all = [h_all h];
plot(XX, h)
xlim([0 8])

subplot(4, 2, 7)
plot_raster(squeeze(raster{2}{36}(1, 5, p_idx{2}(36), :)), 0, 8)
subplot(4, 2, 5)
h = hist(raster_p_sum{2}{36}{5}, XX)/(bin_size*3);
h_all = [h_all h];
plot(XX, h)
xlim([0 8])

subplot(4, 2, 8)
plot_raster(squeeze(raster{2}{36}(1, 5, n_idx{2}(36), :)), 0, 8)
subplot(4, 2, 6)
h = hist(raster_n_sum{2}{36}{5}, XX)/(bin_size*3);
h_all = [h_all h];
plot(XX, h)
xlim([0 8])
h_max = max(h_all);
subplot(4, 2, 1); ylim([0 h_max])
subplot(4, 2, 2); ylim([0 h_max])
subplot(4, 2, 5); ylim([0 h_max])
subplot(4, 2, 6); ylim([0 h_max])

%
FigHandle = figure;
set(FigHandle, 'Position', [1 1 500 400])
h_all = [];
idx = 14;
subplot(4, 2, 3)
plot_raster(squeeze(raster{4}{idx}(1, 5, p_idx{4}(idx), :)), 0, 8)
subplot(4, 2, 1)
h = hist(raster_p_sum{4}{idx}{5}, XX)/(bin_size*3);
h_all = [h_all h];
plot(XX, h)
xlim([0 8])

subplot(4, 2, 4)
plot_raster(squeeze(raster{4}{idx}(1, 5, n_idx{4}(idx), :)), 0, 8)
subplot(4, 2, 2)
h = hist(raster_n_sum{4}{idx}{5}, XX)/(bin_size*3);
h_all = [h_all h];
plot(XX, h)
xlim([0 8])

subplot(4, 2, 7)
plot_raster(squeeze(raster{2}{idx}(1, 5, p_idx{2}(idx), :)), 0, 8)
subplot(4, 2, 5)
h = hist(raster_p_sum{2}{idx}{5}, XX)/(bin_size*3);
h_all = [h_all h];
plot(XX, h)
xlim([0 8])

subplot(4, 2, 8)
plot_raster(squeeze(raster{2}{idx}(1, 5, n_idx{2}(idx), :)), 0, 8)
subplot(4, 2, 6)
h = hist(raster_n_sum{2}{idx}{5}, XX)/(bin_size*3);
h_all = [h_all h];
plot(XX, h)
xlim([0 8])
h_max = max(h_all);
subplot(4, 2, 1); ylim([0 h_max])
subplot(4, 2, 2); ylim([0 h_max])
subplot(4, 2, 5); ylim([0 h_max])
subplot(4, 2, 6); ylim([0 h_max])


%% fig 4d
% id = 5731
idx = 36; 
bin_size = 0.05; d = 3;
XX = bin_size/2:bin_size:4*d-bin_size/2;
for i = 1:2
    figure
    subplot(2, 1, 2)
    a = [0 d/2 d/2 3*d/2 3*d/2 5*d/2 5*d/2 7*d/2 7*d/2 4*d]';
    b = [1 1 2 2 1 1 0 0 1 1]';
    trialnn = length(raster_ff{i}{idx});
    trialn = 5;
    for j = 1:trialn
        SpikeTime = raster_ff{i}{idx}{j};
        SpikeTime = SpikeTime';
        X = [SpikeTime; SpikeTime];
        Y = [ones(1, length(SpikeTime))*(j-0.5); ones(1, length(SpikeTime))*j];
        line(X, Y, 'color', 'b');
        xlim([0 4*d])
        hold on
    end   
    line(a, b+trialn+2)
    xlabel('time/sec')
    subplot(2, 1, 1)
    h = hist(raster_ff_all{i}{idx}, XX)/bin_size/trialnn;
    plot(XX, h)
    xlim([0 4*d])
end

%% fig 4b
load('DS150603.mat')
t = 1;
figure
for ct = 1:4
    errorbar(ratio_dir_mean(:, ct, t), ratio_dir_ste(:, ct, t), '--')
    hold on 
end
errorbar(ratio_oo_mean(:, t), ratio_oo_ste(:, t))

legend('all', 'S', 'P', 'I', 'A')
% ylim([0 3])
% xlim([0.5 5.5])

%% fig 4e
load('DS141212.mat')
figure
errorbar(ratio_ffp_avg, ratio_ffp_ste)
title('ffp off on ratio')
xlabel('light level')
ylabel('off on ratio')


%% plot example ei
pos = datadg{5}.ei.position;
figure
plot_electrodes(pos, 'scale', 0.05, 'alpha', 0);
axis off

cell_id = [1205];
for i = 1:1
    ds_idx = find(id_dir{1} == cell_id(i));
    cell_idx = get_cell_indices(datadg{5}, cell_id(i));
    ei = datadg{5}.ei.eis{cell_idx};
    figure
    plot_ei_(ei, pos, 0, 'scale', 3.5);%, 'Alpha', 0);
    hold on
%     plot(com_oo{1}(ds_idx, 1), com_oo{1}(ds_idx, 2), 'ro')
    plot(pos(455, 1), pos(455, 2), 'bo')
    axis off
end

%% NNND
for dir = 1:4
    com = center_ds(idx_dir{dir}, :);
    nnnd_oo{dir} = get_nnnds_xy(com, repmat(radius,1,2));
end
for dir = 1:3
    com = center_ds(idx_dir_on{dir}, :);
    nnnd_on{dir} = get_nnnds_xy(com, repmat(radius,1,2));
end

figure
for dir = 1:4
    subplot(2,4,dir)
    hist(nnnd_oo{dir}, [0.05:0.1:2.95])
    xlim([0 inf])
end
for dir = 1:3
    subplot(2,4,dir+4)
    hist(nnnd_on{dir}, [0.05:0.1:2.95])
    xlim([0 inf])
end

nnnd_all = [cell2mat(nnnd_oo'); cell2mat(nnnd_on')];
figure
a = hist(nnnd_all, [0.05:0.1:3.45]);
bar([0.05:0.1:3.45], a, 1)
xlabel('nnnd')
ylabel('cell number')