%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);
path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-06-27-0/';
% load drifting grating data
datadg = load_data(strcat(path, 'data005-sorted/data005-sorted'), opt);
datadg.names.stimulus_path = strcat(path, 'stimuli/s05.txt');
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

% identify DS cells
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [4 5]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

datarun = load_data(strcat(path, 'datamfs-map/datamfs-map'), opt);
time_points = [2400 5100];
datafs(1:3) = split_datarun(datarun, time_points);
datafs{1}.names.stimulus_path = strcat(path, 'stimuli/s00.mat');
datafs{1} = load_stim_mfs(datafs{1});

datafs{2}.names.stimulus_path = strcat(path, 'stimuli/s01.mat');
datafs{2} = load_stim_mfs(datafs{2});

datafs{3}.names.stimulus_path = strcat(path, 'stimuli/s03.mat');
datafs{3} = load_stim_mfs(datafs{3});

%% dg
n = 1;
i = 1;
[raster_dg, dg, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
dg{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);


delta_p = 1; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, dg{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(dg{i});
rep = datadg.stimulus.repetitions;
%%
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

v = 4*datadg.stimulus.params.SPATIAL_PERIOD./datadg.stimulus.params.TEMPORAL_PERIOD;
figure
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
xlim([v(end) v(1)])

t = 3;
figure
compass(dg{1}.U{t}(idx_sub{1}), dg{1}.V{t}(idx_sub{1}), 'r')
hold on
compass(dg{1}.U{t}(idx_sub{2}), dg{1}.V{t}(idx_sub{2}), 'b')

%%
d = 1;
t = 3;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(dg{d}.U{t}(idx_sub{2}), dg{d}.V{t}(idx_sub{2}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(dg{d}.U{t}(idx_sub{2}), dg{d}.V{t}(idx_sub{2}), x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = idx_sub{2}(I);
    id_dir{i} = ds_id(idx_dir{i});
end
%% plot rfs
cell_type = {'superior', 'anterior', 'inferior', 'posterior'};
field_width = 17; field_height = 17;
field_width_sta = 40; field_height_sta = 40;
subregion = 0;
stop = 0.5; %second
for i = 1:3
%     fs_raster{i} = get_fs_raster(datafs{i}, ds_id, 'stop', 0.5);
    fs_raster{i} = get_fs_raster(datafs{i}, ds_id);
    for cc = 1:length(ds_id)
        if fs_idx(cc)
            fs_raster{i}{cc} = [];
        end
    end
    fs_spike{i} = fs_get_spike(fs_raster{i});
    [rf_all{i}, rf_std{i}] = get_fs_rf(fs_spike{i}, field_width, field_height,subregion);
end

%%
for cc = 1:length(ds_id)
    figure(1)
    set(gcf, 'Position', [1 1 1000 600])
    id = ds_id(cc);
    for i = 1:3
        if ~isempty(rf_all{i}{cc})
%             rf = padarray(rf_all{i}{cc},[7,7]);
            rf = rf_all{i}{cc};
    
            subplot(3,3,3*(i-1)+1)
            imagesc(sum(rf,3))
            colormap gray
            axis image
            axis off

            subplot(3,3,3*(i-1)+2)
            imagesc(rf(:,:,1))
            colormap gray
            axis image
            axis off

            subplot(3,3,3*i)
            imagesc(rf(:,:,2))
            colormap gray
            axis image
            axis off

        end
    end
    print_close(1,[14 10],num2str(id))
%     pause
%     close(1)
end
