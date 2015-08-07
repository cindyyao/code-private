%% load data
addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/Classification/
addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/DS' cell analysis'/

opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);
cd /Users/xyao/matlab/code-private/DS_new/

% load data
datadg{1} = load_data('/Volumes/lab/analysis/2015-06-09-0/data005-map/data005-map', opt);
datadg{1}.names.stimulus_path = '/Volumes/lab/analysis/2015-06-09-0/stimuli/s05.mat';
datadg{1} = load_stim_matlab(datadg{1}, 'user_defined_trigger_interval', 10);

%%
n = 1; % number of datarun
i = 1; % which datarun to use for classification
[NumSpikesCell, StimComb] = get_spikescellstim(datadg{i},datadg{i}.cell_ids,0);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);

params_idx = [1 2]; % which parameters to use for classification
[ds_id, nonds_id] = classify_ds(datadg{i}, ds_struct, params_idx);

%
figure(1)
xlabel('TP 1')
ylabel('TP 2')
title('NDF 0')
figure(2)
legend('non DS', 'DS')
% 

%%
[raster_dg, DG, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, StimComb] = get_spikescellstim(datadg{i},ds_id,0);
    DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster_dg{i} = get_ds_raster(datadg{i}, ds_id);
end

delta_p = 1; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
    MAG_all_norm_dg{i} = normalize_MAG(DG{i});
    rep = datadg{i}.stimulus.repetitions;
end

ll = {'NDF0'};

%% plot cell summary

for cc = 1:5 %length(ds_id)
    plot_ds_raster(DG, raster_dg, cc, ds_id(cc), ll, 1, 1, 0)
end



