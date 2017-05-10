cd /Users/xyao/matlab/code-private/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);


datadg = load_data('/Volumes/lab/analysis/2016-02-09-0/data002/data002', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-02-09-0/stimuli/s02';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

% dataflash = load_data('/Volumes/lab/analysis/2016-02-19-0/data001/data001', opt);
dataflash = load_data('/Volumes/lab/analysis/2016-02-09-0/data000-001-map/data000-001-map', opt);
dataflash.DfParams.NDF =   [5,5,5,4,4,4,4,4,3,3,3,2,2,1,0] ; % on filter turret 
dataflash.DfParams.Ftime = [2,4,8,2,3,4,6,8,2,4,8,2,4,2,2] ; % ms
dataflash.DfParams.interFlashInt = [3] ; % sec
dataflash.triggers(160:161) = [];

[NumSpikesCell, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);
params_idx = [4 5]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);


ds_id_flash = ds_id(flash_idx);
ds_id_flash(24) = 6522;
ds_idx_flash = get_cell_indices(dataflash, ds_id_flash);

%% 
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell, StimComb] = get_spikescellstim(datadg,ds_id,0);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
raster_dg{i} = get_ds_raster(datadg, ds_id);
%     for j = 1:length(raster_dg{i})
%         if(dg_idx(j, i))
%             raster_dg{i}{j} = [];
%         end
%     end

delta_p = 4; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(DG{i});
rep = datadg.stimulus.repetitions;

%%
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [4 5]);

for cc = 1:length(ds_id)
    plot_ds_raster(DG_cut, raster_dg_cut, cc, ds_id(cc), '', 1, 1, 1)
end

% for cc = 2:2 %length(ds_id)
%     plot_ds_raster(DG, raster_dg, cc, ds_id(cc), '', 1, 1, 0)
% end

