cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);


datadg = load_data('/Volumes/lab/analysis/2016-06-10-0/data001/data001', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-06-10-0/stimuli/s01.mat';
datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);

datamb{1} = load_data('/Volumes/lab/analysis/2016-06-10-0/data002-003-map/data002-003-map', opt);
datamb{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-06-10-0/stimuli/s02.mat';
datamb{1} = load_stim_matlab(datamb{1});

datamb{2} = load_data('/Volumes/lab/analysis/2016-06-10-0/data004-map/data004-map', opt);
datamb{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-06-10-0/stimuli/s04.mat';
datamb{2} = load_stim_matlab(datamb{2});

datamb{3} = load_data('/Volumes/lab/analysis/2016-06-10-0/data005-map/data005-map', opt);
datamb{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-06-10-0/stimuli/s05.mat';
datamb{3} = load_stim_matlab(datamb{3});

%%
[NumSpikesCell, MaxRate,StimComb] = get_spikescellstim_mb(datamb{1},datamb{1}.cell_ids,3008.4,0.1);
ds_struct = mbcellanalysis(MaxRate, StimComb,datamb{1});
params_idx = [2 3]; % which parameters to use for classification

[id, ~] = classify_ds(datamb{1}, ds_struct, params_idx);


%
[NumSpikesCell, maxrate,StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
params_idx = [2 3]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, maxrate,StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));

n = 3;
duration = 3008.4; %sec
bin_size = 0.025; %sec
[raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},datamb{i}.cell_ids,duration,bin_size);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb{i} = get_mb_raster(datamb{i}, datamb{i}.cell_ids, duration);
%     for j = 1:length(raster_mb{i})
%         if(mb_idx(j))
%             raster_mb{i}{j} = [];
%         end
%     end
    trial_dur{i} = get_mb_trial_dur(datamb{i});
end
raster_dg = get_ds_raster(datadg, ds_id);
%%
for cc = 1:length(ds_id);
plot_mb_raster_ctr_12(MB, raster_mb, trial_dur, cc, ds_id(cc), 'NDF0', 3, 8, 1)
end

%% 
cc = 61;
plot_ds_raster({ds_struct}, {raster_dg}, cc, ds_id(cc), '',1,1, 0)
%%
figure
compass(ds_struct.U{1}, ds_struct.V{1});

[x, y] = ginput;
hold on
plot(x, y, 'r');

IN = inpolygon(ds_struct.U{1}, ds_struct.V{1}, x, y);
[~, I] = find(IN == 1);

