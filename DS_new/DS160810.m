cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);


datadg = load_data('/Volumes/lab/analysis/2016-08-10-0/data002/data002', opt);
trigger = [0:10:2870]';
datadg.triggers = sort([datadg.triggers; trigger]);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-08-10-0/stimuli/s02.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);
% datadg = convert_stim_160810(datadg);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis_photon(NumSpikesCell, StimComb, datadg);
params_idx = [3 4]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);


datamb = load_data('/Volumes/lab/analysis/2016-08-10-0/data000-map/data000-map', opt);
datamb.names.stimulus_path = '/Volumes/lab/analysis/2016-08-10-0/stimuli/s00.txt';
datamb = load_stim(datamb, 'user_defined_trigger_set', [1:2:length(datamb.triggers)]);

