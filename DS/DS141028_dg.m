%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
datadg{1} = load_data('/Analysis/xyao/2014-10-28-0/data000-map/data000-map', opt);
datadg{1}.names.stimulus_path = '/Analysis/xyao/2014-10-28-0/stimuli/s00.mat';
datadg{1} = load_stim_matlab(datadg{1}, 'user_defined_trigger_interval', 10);
datadg{2} = load_data('/Analysis/xyao/2014-10-28-0/data003-map/data003-map', opt);
datadg{2}.names.stimulus_path = '/Analysis/xyao/2014-10-28-0/stimuli/s03.mat';
datadg{2} = load_stim_matlab(datadg{2}, 'user_defined_trigger_interval', 10);

dataffp{1} = load_data('/Analysis/xyao/2014-10-28-0/data002-map/data002-map', opt);
dataffp{1}.triggers = dataffp{1}.triggers(2:end);
dataffp{2} = load_data('/Analysis/xyao/2014-10-28-0/data005-map/data005-map', opt);
dataffp{2}.triggers = dataffp{2}.triggers(2:end);

%% drifting grating

[NumSpikesCell, StimComb] = get_spikescellstim(datadg{1},datadg{1}.cell_ids,0);
dg_struct = dscellanalysis(NumSpikesCell, StimComb);

n = 2;
% pull out DS cells

figure
plot(dg_struct.mag{4, 1}, dg_struct.mag{5, 1}, 'o')
title('data002 after mapping')
xlabel('TP 24')
ylabel('TP 60')
hold on
[x, y] = ginput;
plot(x, y);

IN = inpolygon(dg_struct.mag{4, 1}, dg_struct.mag{5, 1}, x, y);
[~, I] = find(IN == 1);
id = datadg{1}.cell_ids(I);
idx = get_cell_indices(datadg{1}, id);

I = zeros(length(id), n);
for i = 1:n
    I(:, i) = ismember(id, datadg{i}.cell_ids);
end
I = sum(I');
id_all = id(I == n); % DS cells that can be found in all light levels
idx_all = arrayfun(@(x) find(id == x), id_all);

[raster, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, StimComb] = get_spikescellstim(datadg{i},id,0);
    DS{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb));
    raster{i} = get_ds_raster(datadg{i}, id);
end

param_p = 5; % choose which params to use to calculate prefer direction indices 
MAG_all_norm = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster{i}, DS{i}.angle{param_p});
    MAG_all_norm{i} = normalize_MAG(DS{i});
    rep = datadg{i}.stimulus.repetitions;
end

ll = {'NDF3', 'NDF0'};


%% plot cell summary

for cc = 1:length(id)
    plot_ds_raster(DS, raster, cc, id(cc), ll, 1, 2, 1)
end

%% plot single cell tuning curve

% use unnormalized vector sum as response
figure
for i = 1:n
    v = datadg{i}.stimulus.params.SPATIAL_PERIOD./datadg{i}.stimulus.params.TEMPORAL_PERIOD;
    subplot(1, 2, i)
    semilogx(v, exciseColumn(MAG_all_norm{i}), 'b')
    xlabel('speed')
    ylabel('Response')
    title(ll{i})
    xlim([v(end) v(1)])
end

% use maximum firing rate of preferred direction as response

figure
for i = 1:n
    v = datadg{i}.stimulus.params.SPATIAL_PERIOD./datadg{i}.stimulus.params.TEMPORAL_PERIOD;
    subplot(1, 2, i)
    semilogx(v, exciseColumn(norm_max_r{i}'), 'b')
    xlabel('speed')
    ylabel('Response')
    title(ll{i})
    xlim([v(end) v(1)])
end

%% classification based on speed tunning
%% pca
mag_pca = MAG_all_norm{1};
% mag_pca = MAG_all_norm{1}(:, idx4);
mag_pca = mag_pca';
[id_sub, idx_sub] = deal(cell(2, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(mag_pca);
pc1 = 3; pc2 = 1;
plot(scores(:, pc1), scores(:, pc2), 'o')
hold on
for i = 1:2
    [x, y] = ginput;
    plot(x, y)
    IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
    [~, idx_sub{i}] = find(IN' == 1);
    id_sub{i} = id(idx_sub{i});
end
xlabel('3rd Principal Component')
ylabel('1st Principal Component')

for j = 1:2
    figure
    for i = 1:n
        v = datadg{i}.stimulus.params.SPATIAL_PERIOD./datadg{i}.stimulus.params.TEMPORAL_PERIOD;
        subplot(1, 2, i)
        semilogx(v, exciseColumn(MAG_all_norm{i}(:, idx_sub{j})), 'b')
        xlabel('speed')
        ylabel('Response')
        title(ll{i})
        xlim([v(end) v(1)])
    end
end

%% plot average tunning curve
color = 'br';
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:n
    v = datadg{i}.stimulus.params.SPATIAL_PERIOD./datadg{i}.stimulus.params.TEMPORAL_PERIOD;
    subplot(1, 2, i)
    for ct = 1:2
        mag_temp = exciseColumn(MAG_all_norm{i}(:, idx_sub{ct}));
        tuning_avg{i}(:, ct) = mean(mag_temp, 2);
        tuning_ste{i}(:, ct) = std(mag_temp, [], 2)/sqrt(size(mag_temp, 2));
        errorbar(v, tuning_avg{i}(:, ct), tuning_ste{i}(:, ct), color(ct))
        hold on
    end
    set(gca, 'XScale', 'log')
    title(ll{i})
    xlim([min(v) max(v)])
    xlabel('speed')
    ylabel('response')
end
legend('on-off DSGC', 'on DSGC', 'location', 'southeast')

%% vector sum plot of individual type
for ct = 1:2
    [NumSpikesCell, StimComb] = get_spikescellstim(datadg{1},id_sub{ct},0);
    dscellanalysis(NumSpikesCell, StimComb);
end

%% compare across light level
SP = {'SP 240'};
CT = {'on-off DSGC', 'on-DSGC'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
v = datadg{1}.stimulus.params.SPATIAL_PERIOD./datadg{i}.stimulus.params.TEMPORAL_PERIOD;
for j = 1:2
    subplot(1, 2, j)
    errorbar(v, tuning_avg{1}(:, j), tuning_ste{1}(:, j), 'b')
    hold on
    errorbar(v, tuning_avg{2}(:, j), tuning_ste{2}(:, j), 'r')
    title([SP{1} '  ' CT{j}])
    set(gca, 'XScale', 'log')
    xlim([min(v) max(v)])
    xlabel('speed')
    ylabel('response')

end
legend('NDF 3', 'NDF 0', 'location', 'northeast')


%% full field pulses

[raster_ff, raster_ff_all] = deal(cell(2, 1));
for d = 1:2
    [raster_ff{d}, raster_ff_all{d}] = get_ffp_raster(dataffp{d}, id, 3);
end

i_empty{1} = [3 11 14 20 23 24 25 27 32 36 40];
i_empty{2} = [1 2 3 4 6 10 11 14 18 20 21 22 23 24 25 28 31 36 38 40 41 42];

for i = 1:2
    for j = 1:length(i_empty{i})
        raster_ff{i}{i_empty{i}(j)} = [];
        raster_ff_all{i}{i_empty{i}(j)} = [];
    end
end

for i = 6:length(id) 
    if ~isempty(raster_ff{1}{i}) || ~isempty(raster_ff{2}{i})
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 800 800])
        for d = 1:2
        subplot(1, 2, d)
        plot_ffp(raster_ff{d}, raster_ff_all{d}, i, 3)
        title([num2str(id(i)) ' ' ll{d}(1:4)])
        end
        
        print_close(1, [12, 12], num2str(id(i)))
    end
end

