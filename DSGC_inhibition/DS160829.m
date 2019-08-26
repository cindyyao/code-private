%% load drifting grating data
cd /Volumes/dusom_fieldlab/All_Staff/lab/Development/matlab/private/xyao/matlab/code-private/DS_new
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);
path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2016-08-29-0/';

datadg = load_data(strcat(path, 'data003-sorted/data003-sorted'), opt);
datadg.names.stimulus_path = strcat(path, 'stimuli/s03.mat');
datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);

%% identify DSGCs
[NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg, datadg.cell_ids, 0, 1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb, datadg);
params_idx = [2 5]; % which parameters to use for classification

[ds_id, nonds_id, id_init] = classify_ds(datadg, ds_struct, params_idx, 'manual', true);

%% classify ON-OFF vs ON DSGCs based on speed tunning
[NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg,ds_id,0,1);
DG = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg = get_ds_raster(datadg, ds_id);

delta_p = 4; % choose which params to use to calculate prefer direction indices 
[raster_p_sum, p_idx] = get_pdirection_raster(raster_dg, DG.angle{delta_p});
MAG_all_norm_dg = normalize_MAG(DG);
rep = datadg.stimulus.repetitions;

mag_pca = MAG_all_norm_dg;
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
ylabel('2nd Principal Component')
title('NDF 0')

v = 4*datadg.stimulus.params.SPATIAL_PERIOD./datadg.stimulus.params.TEMPORAL_PERIOD;
figure
semilogx(v, exciseColumn(MAG_all_norm_dg(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
xlim([v(end) v(1)])

t = 6;
figure
compass(DG.U{t}(idx_sub{1}), DG.V{t}(idx_sub{1}), 'r')
hold on
compass(DG.U{t}(idx_sub{2}), DG.V{t}(idx_sub{2}), 'b')

%% Classify DSGCs by direction

% ON-OFF
t = 5;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(DG.U{t}(idx_sub{2}), DG.V{t}(idx_sub{2}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG.U{t}(idx_sub{2}), DG.V{t}(idx_sub{2}), x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = idx_sub{2}(I);
    id_dir{i} = ds_id(idx_dir{i});
end

% ON
t = 3;
h = figure;
dirn = 3;
set(h, 'Position', [1 1 1080 500])
compass(DG.U{t}(idx_sub{1}), DG.V{t}(idx_sub{1}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG.U{t}(idx_sub{1}), DG.V{t}(idx_sub{1}), x, y);
    [~, I] = find(IN == 1);
    idx_dir_on{i} = idx_sub{1}(I);
    id_dir_on{i} = ds_id(idx_dir_on{i});
end

% save('DS160829.mat', 'ds_id', 'id_dir', 'id_dir_on', 'idx_dir', 'idx_dir_on'. '-append')