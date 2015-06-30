%% load data
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load data
datadg{1} = load_data('/Analysis/xyao/2014-10-28-0/data000-map/data000-map', opt);
datadg{1}.names.stimulus_path = '/Analysis/xyao/2014-10-28-0/stimuli/s00.mat';
datadg{1} = load_stim_matlab(datadg{1}, 'user_defined_trigger_interval', 10);
datadg{2} = load_data('/Analysis/xyao/2014-10-28-0/data003-map/data003-map', opt);
datadg{2}.names.stimulus_path = '/Analysis/xyao/2014-10-28-0/stimuli/s03.mat';
datadg{2} = load_stim_matlab(datadg{2}, 'user_defined_trigger_interval', 10);

datamb{1} = load_data('/Analysis/xyao/2014-10-28-0/data001-map004/data001-map004', opt);
datamb{1}.names.stimulus_path = '/Analysis/xyao/2014-10-28-0/stimuli/s01.mat';
datamb{1} = load_stim_matlab(datamb{1});
datamb{2} = load_data('/Analysis/xyao/2014-10-28-0/data004/data004', opt);
datamb{2}.names.stimulus_path = '/Analysis/xyao/2014-10-28-0/stimuli/s04.mat';
datamb{2} = load_stim_matlab(datamb{2});

dataffp{1} = load_data('/Analysis/xyao/2014-10-28-0/data002/data002', opt);
dataffp{1}.triggers = dataffp{1}.triggers(2:end);
dataffp{2} = load_data('/Analysis/xyao/2014-10-28-0/data005/data005', opt);
dataffp{2}.triggers = dataffp{2}.triggers(2:end);

%% Moving Bar

% bin_size = [1];
% for i = 1:length(bin_size)
%     [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datarun{2},datarun{2}.cell_ids,0, bin_size(i));
%     mb_struct = mbcellanalysis(MaxRate, StimComb);
%     title(['bin size ' num2str(bin_size(i))])
% end

[NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{2},datamb{2}.cell_ids,0, 1);
mb_struct = mbcellanalysis(NumSpikesCell, StimComb);

n = 2;
duration = 3653;
% pull out DS cells

figure
plot(mb_struct.mag{2, 1}, mb_struct.mag{3, 1}, 'o')
title('data002 after mapping')
xlabel('TP 24')
ylabel('TP 60')
hold on
[x, y] = ginput;
plot(x, y);

IN = inpolygon(mb_struct.mag{2, 1}, mb_struct.mag{3, 1}, x, y);
[~, I] = find(IN == 1);
id = datamb{2}.cell_ids(I);
idx = get_cell_indices(datamb{2}, id);

I = zeros(length(id), n);
for i = 1:n
    I(:, i) = ismember(id, datamb{i}.cell_ids);
end
I = sum(I');
id_all = id(I == n); % DS cells that can be found in all light levels
idx_all = arrayfun(@(x) find(id == x), id_all);

[raster, MB, trial_dur, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},id,0,1);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb));
    raster{i} = get_mb_raster(datamb{i}, id, duration);
    trial_dur{i} = get_mb_trial_dur(datamb{i});
end

param_p = 5; % choose which params to use to calculate prefer direction indices 
MAG_all_norm = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster{i}, MB{i}.angle{param_p});
    MAG_all_norm{i} = normalize_MAG(MB{i});
    rep = datamb{i}.stimulus.repetitions;
end

ll = {'NDF3', 'NDF0'};


%% plot cell summary

for cc = 28:28
    plot_mb_raster(MB, raster, trial_dur, cc, id(cc), ll, 1, 2, 1)
end



%%

for i = 1:2
    triggers_ffp{i} = dataffp{i}.triggers(1:4:length(dataffp{i}.triggers));
end


d = 1;
for cell_idx = 1:50
get_raster(dataffp{d}.spikes{cell_idx}, triggers_ffp{d});
end

%%
figure
    v = datarun{4}.stimulus.params.SPATIAL_PERIOD./datarun{4}.stimulus.params.TEMPORAL_PERIOD;
    semilogx(v, exciseColumn(MAG_all_norm{}), 'b')
    xlabel('speed')
    ylabel('Response')
    xlim([v(end) v(1)])

