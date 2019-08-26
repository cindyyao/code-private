cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

datadg = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-06-01-0/data005-xy/data005-xy', opt);
datadg.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2018-06-01-0/stimuli/s05.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);
datadg = load_ei(datadg, ds_id, 'array_id',1501);

% datadg = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-06-01-0/data003-sorted/data003-sorted', opt);
% datadg.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2018-06-01-0/stimuli/s03.txt';
% datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);
% datadg = load_ei(datadg, ds_id);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

params_idx = [1 2]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

datawn{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-06-01-0/data000-map/data000-map', opt);
datawn{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-06-01-0/data001-map/data001-map', opt);
datawn{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-06-01-0/data002-map/data002-map', opt);
datawn{3} = load_ei(datawn{3}, id_dir{1}, 'array_id',1501);


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
pc1 = 1; pc2 = 2;
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
t = 2;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(DG{d}.U{t}, DG{d}.V{t});
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}, DG{d}.V{t}, x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = I;
    id_dir{i} = ds_id(idx_dir{i});
end

%% plot cell summary
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [5 6]);

for cc = 1:length(ds_id)
    plot_ds_raster(DG_cut, raster_dg_cut, cc, ds_id(cc), 1, 1, 1, 1)
end

%% neighboring pairs
pos = datadg.ei.position;
mode = 'neg';
neighbors = [];
ct = 4;
id_wn = id_dir{ct};
for cc1 = 1:length(id_wn)
    for cc2 = cc1+1:length(id_wn)
        id1 = id_wn(cc1);
        idx1 = get_cell_indices(datadg, id1);
        ei1 = datadg.ei.eis{idx1};
        com1 = ei_com_xy(ei1, pos, 30*3, mode);
        id2 = id_wn(cc2);
        idx2 = get_cell_indices(datadg, id2);
        ei2 = datadg.ei.eis{idx2};
        com2 = ei_com_xy(ei2, pos, 30*3, mode);
        if pdist([com1;com2]) < 150
            neighbors = [neighbors; id1 id2];
        end
    end
end

cn = 0;
celltype = {'superior', 'anterior', 'inferior', 'posterior'};
coms = [];
for cc = 1:length(id_wn)
    id = id_wn(cc);
    idx = get_cell_indices(datadg, id);
    ei = datadg.ei.eis{idx};
    com = ei_com_xy(ei, pos, 30*3, mode);
    coms = [coms; com];
end

corner_i = [4 126 195 264 386 455 4];
corner_position = datadg.ei.position(corner_i, :);
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
duration = 2700;
bin_size = 0.00025;
max_lag = 40;
ids = id_wn;
for cc1 = 1:length(ids)-1
    xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
    FigHandle = figure;
    set(FigHandle, 'Position', [1 1 2000 1500])
    for cc2 = cc1+1:length(ids)
        id1 = ids(cc1);
        id2 = ids(cc2);

        idx1 = get_cell_indices(datawn{3}, id1);
        idx2 = get_cell_indices(datawn{3}, id2);
        spikes1 = datawn{3}.spikes{idx1};
        spikes1_TF= ceil(spikes1/bin_size);
        spikes1 = zeros(duration/bin_size, 1);
        spikes1(spikes1_TF) = 1;

        spikes2 = datawn{3}.spikes{idx2};
        spikes2_TF= ceil(spikes2/bin_size);
        spikes2 = zeros(duration/bin_size, 1);
        spikes2(spikes2_TF) = 1;
        A = xcorr(spikes1, spikes2, max_lag, 'coeff');
        subplot(3,5,cc2)
        bar(xx, A, 'k')
        title([num2str(id1) '  ' num2str(id2)])
        xlim([-0.01 0.01])
    end
%     pause
    print_close(1, [24,12], num2str(id1))
end

%%
LL = {'NDF 4', 'NDF 2', 'NDF 0'};
duration = 2700;
bin_size = 0.0005;
max_lag = 20;

neighbors = [1789 2134;2134 2629; 2134 2673];
for cp =1:size(neighbors, 1) 
    xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
    FigHandle = figure;
    set(FigHandle, 'Position', [1 1 800 300])
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
%         A_all(cp, :, ll) = A;
        
        A = xcorr(spikes1, spikes2, max_lag, 'coeff');
%         A_all_coeff(cp, :, ll) = A;
        
        subplot(1,3,ll)
        bar(xx, A, 'k')
        title([num2str(id1) '  ' num2str(id2) ' ' LL{ll} ])
        xlim([-0.01 0.01])
    end
%     print_close(1, [12,4], [num2str(id1) '  ' num2str(id2)])
end

t2p_coeff = squeeze(max(A_all_coeff, [], 2)) - squeeze(A_all_coeff(:, 41, :));

