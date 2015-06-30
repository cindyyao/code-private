%% load data
addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/Classification/
addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/DS' cell analysis'/

opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load data
datadg{1} = load_data('/Analysis/xyao/2014-10-24-0/data008/data008', opt);
datadg{1}.names.stimulus_path = '/Analysis/xyao/2014-10-24-0/stimuli/s08';
datadg{1} = load_stim(datadg{1}, 'user_defined_trigger_interval', 10);


%%
dt = 1;

[NumSpikesCell, StimComb] = get_spikescellstim(datadg{dt},datadg{dt}.cell_ids,0);
dg_struct = dscellanalysis(NumSpikesCell, StimComb);

figure
plot(dg_struct.mag{1, 1}, dg_struct.mag{2, 1}, 'o')
hold on
[x, y] = ginput;
plot(x, y);
IN = inpolygon(dg_struct.mag{1, 1}, dg_struct.mag{2, 1}, x, y);
[~, I] = find(IN == 1);
id_init = datadg{dt}.cell_ids(I);

[C ia ib] = intersect(id_init, datadg{dt}.cell_ids);
vc = ones(length(datadg{dt}.cell_ids),1);
vc(ib) = 2; %initializing ds cells to cluster 2, everything else cluster 1

close all;
X = [];
N = [];
p = [];
X(:,1) = log(dg_struct.mag{1,1})';
X(:,2) = log(dg_struct.mag{2,1})';
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 0, 1, 0, 0, 0,0, vc);

ds_id = [];
ds_id = datadg{dt}.cell_ids(idx==2);
nonds_id = datadg{dt}.cell_ids(idx==1);
ds_idx = get_cell_indices(datadg{dt}, ds_id);
[NumSpikesCell, StimComb] = get_spikescellstim(datadg{dt},ds_id,0);
dg_struct = dscellanalysis(NumSpikesCell, StimComb);

%% raster

n_dg = 1;
[raster, raster_p_sum, p_idx] = deal(cell(n_dg, 1));
for i = 1:n_dg   
    [NumSpikesCell, StimComb] = get_spikescellstim(datadg{i},ds_id,0);
    DS{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster{i} = get_ds_raster(datadg{i}, ds_id);
%     for j = 1:length(raster{i})
%         if(ds_idx(j, i))
%             raster{i}{j} = [];
%         end
%     end
end

param_p_t = 2; % choose which params to use to calculate prefer direction indices 
param_p_s = 2;
MAG_all_norm = cell(n_dg, 1);
max_r = cell(n_dg, 1);
norm_max_r = cell(n_dg, 1);

for i = 1:n_dg
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster{i}, DS{i}.angle{param_p_t, param_p_s});
    MAG_all_norm{i} = normalize_MAG(DS{i});
    rep = datadg{i}.stimulus.repetitions;
end

%% plot direction tuning curves

spi = 2;
tpi = 2;
dt = 1;
angle = 0:20:340;

rho = DS{dt}.rho{tpi, spi};
figure
for i = 15:15 %length(ds_id)
    plot(angle, rho(i, :), 'b')
    pause
end

p_direction = DS{dt}.angle{tpi, spi}';
xx = 0:pi/9:17*pi/9;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 18);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

figure
for i = 1:length(ds_id)
    subplot(5, 6, i)
    [x, idx] = sort(xx(i, :));
    rho_temp = rho(i, :);
    plot(x, rho_temp(idx), 'b');
end
