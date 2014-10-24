%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
datarun{1} = load_data('/Analysis/xyao/2014-09-24-0/data000-map/data000-map', opt);
datarun{1}.names.stimulus_path = '/Analysis/xyao/2014-09-24-0/stimuli/s00';
datarun{1} = load_stim(datarun{1}, 'user_defined_trigger_interval', 10);

datarun{2} = load_data('/Analysis/xyao/2014-09-24-0/data002-map/data002-map', opt);
datarun{2}.names.stimulus_path = '/Analysis/xyao/2014-09-24-0/stimuli/s02';
datarun{2} = load_stim(datarun{2}, 'user_defined_trigger_interval', 10);

datarun{3} = load_data('/Analysis/xyao/2014-09-24-0/data004-map/data004-map', opt);
datarun{3}.names.stimulus_path = '/Analysis/xyao/2014-09-24-0/stimuli/s04';
datarun{3} = load_stim(datarun{3}, 'user_defined_trigger_interval', 10);

datarun{4} = load_data('/Analysis/xyao/2014-09-24-0/data006-map/data006-map', opt);
datarun{4}.names.stimulus_path = '/Analysis/xyao/2014-09-24-0/stimuli/s06';
datarun{4} = load_stim(datarun{4}, 'user_defined_trigger_interval', 10);

[NumSpikesCell, StimComb] = get_spikescellstim(datarun{3},datarun{3}.cell_ids,0);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);

% pull out DS cells

figure
plot(ds_struct.mag{2, 1}, ds_struct.mag{3, 1}, 'o')
title('data002 after mapping')
xlabel('TP 24')
ylabel('TP 60')
hold on
[x, y] = ginput;
plot(x, y);

IN = inpolygon(ds_struct.mag{2, 1}, ds_struct.mag{3, 1}, x, y);
[~, I] = find(IN == 1);
id = datarun{3}.cell_ids(I);

I = zeros(length(id), 4);
for i = 1:4
    I(:, i) = ismember(id, datarun{i}.cell_ids);
end
I = sum(I');
id4 = id(I == 4); % DS cells that can be found in all light levels
idx4 = arrayfun(@(x) find(id == x), id4);

[raster, raster_p_sum, p_idx] = deal(cell(4, 1));
for i = 1:4    
    [NumSpikesCell, StimComb] = get_spikescellstim(datarun{i},id,0);
    DS{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster{i} = get_ds_raster(datarun{i}, id);
end

param_p = 5; % choose which params to use to calculate prefer direction indices 
MAG_all_norm = cell(4, 1);
max_r = cell(4, 1);
norm_max_r = cell(4, 1);

for i = 1:4
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster{i}, DS{i}.angle{param_p});
    MAG_all_norm{i} = normalize_MAG(DS{i});
    rep = datarun{i}.stimulus.repetitions;
    max_r{i} = max_firing_rate(raster_p_sum{i}, 0.2, 8)/rep;
    norm_max_r{i} = max_r{i}./repmat(max(max_r{i}, [], 2), 1, size(max_r{i}, 2));
end

ll = {'NDF3 SP240', 'NDF2 SP240', 'NDF1 SP240', 'NDF0 SP240'};

%% plot cell summary

for cc = 58:58 %length(id)
    plot_ds_raster(DS, raster, cc, id(cc), ll, 2, 2, 0)
end



