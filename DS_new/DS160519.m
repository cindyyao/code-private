cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);


datadg = load_data('/Volumes/lab/analysis/2016-05-19-0/data007-sorted/data007-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-05-19-0/stimuli/s07.mat';
datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);

[NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [1 2]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);
title('Normalized Vector Sum Amplitude')
xlabel('High speed - 480 micron/s')
ylabel('Low speed - 240 micron/s')

[NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
