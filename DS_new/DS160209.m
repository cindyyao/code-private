cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);


datadg{1} = load_data('/Volumes/lab/analysis/2016-02-09-0/data002/data002', opt);
datadg{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-02-09-0/stimuli/s02';
datadg{1} = load_stim(datadg{1}, 'user_defined_trigger_interval', 10);


i = 1;
[NumSpikesCell, StimComb] = get_spikescellstim(datadg{i},datadg{i}.cell_ids,0);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);
params_idx = [5 4]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg{i}, ds_struct, params_idx);

[NumSpikesCell, StimComb] = get_spikescellstim(datadg{i},datadg{i}.cell_ids,0);

DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg{5}));
raster_dg{i} = get_ds_raster(datadg{i}, datadg{i}.cell_ids);
delta_p = 4;
[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});


[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [5 6]);
%%
% id = 962;
% cc = get_cell_indices(datadg{1}, id);
plot_ds_raster(DG_cut, raster_dg_cut, 1, 7563, 'test', 1, 1, 0)
