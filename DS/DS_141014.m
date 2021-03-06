%% load data
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load data
datarun{1} = load_data('/Analysis/xyao/2014-10-14-0/data000-mapn/data000-map', opt);
datarun{1}.names.stimulus_path = '/Analysis/xyao/2014-10-14-0/stimuli/s00.mat';
datarun{1} = load_stim_matlab(datarun{1});

datarun{2} = load_data('/Analysis/xyao/2014-10-14-0/data003/data003', opt);
datarun{2}.names.stimulus_path = '/Analysis/xyao/2014-10-14-0/stimuli/s03.mat';
datarun{2} = load_stim_matlab(datarun{2});

datarun{3} = load_data('/Analysis/xyao/2014-10-14-0/data002-mapn/data002-map', opt);
datarun{3}.triggers = datarun{3}.triggers(2:end);
datarun{4} = load_data('/Analysis/xyao/2014-10-14-0/data006-mapn/data006-map', opt);
datarun{4}.triggers = datarun{4}.triggers(2:end);

datarun{5} = load_data('/Analysis/xyao/2014-10-14-0/data001-map/data001-map', opt);
datarun{5}.names.stimulus_path = '/Analysis/xyao/2014-10-14-0/stimuli/s01.mat';
datarun{5} = load_stim_matlab(datarun{5});

datarun{6} = load_data('/Analysis/xyao/2014-10-14-0/data004-map/data004-map', opt);
datarun{6}.names.stimulus_path = '/Analysis/xyao/2014-10-14-0/stimuli/s04.mat';
datarun{6} = load_stim_matlab(datarun{6});


[NumSpikesCell, StimComb] = get_spikescellstim(datarun{2},datarun{2}.cell_ids,0);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);
n = 2;
% pull out DS cells

figure
plot(ds_struct.mag{3, 1}, ds_struct.mag{4, 1}, 'o')
title('data002 after mapping')
xlabel('TP 24')
ylabel('TP 60')
hold on
[x, y] = ginput;
plot(x, y);

IN = inpolygon(ds_struct.mag{3, 1}, ds_struct.mag{4, 1}, x, y);
[~, I] = find(IN == 1);
id = datarun{2}.cell_ids(I);
idx = get_cell_indices(datarun{2}, id);
I = zeros(length(id), n);
for i = 1:n
    I(:, i) = ismember(id, datarun{i}.cell_ids);
end
I = sum(I');
id_all = id(I == n); % DS cells that can be found in all light levels
idx_all = arrayfun(@(x) find(id == x), id_all);

[raster, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, StimComb] = get_spikescellstim(datarun{i},id,0);
    DS{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster{i} = get_ds_raster(datarun{i}, id);
end

param_p = 5; % choose which params to use to calculate prefer direction indices 
MAG_all_norm = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster{i}, DS{i}.angle{param_p});
    MAG_all_norm{i} = normalize_MAG(DS{i});
    rep = datarun{i}.stimulus.repetitions;
    max_r{i} = max_firing_rate(raster_p_sum{i}, 0.2, 8)/rep;
    norm_max_r{i} = max_r{i}./repmat(max(max_r{i}, [], 2), 1, size(max_r{i}, 2));
end

ll = {'NDF3 SP240', 'NDF0 SP240'};


%% plot cell summary

for cc = 45:50 %length(id)
    plot_ds_raster(DS, raster, cc, id(cc), ll, 1, 2, 1)
end

%% plot single cell tuning curve

% use unnormalized vector sum as response
figure
for i = 1:n
    v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
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
    v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    subplot(2, 2, i)
    semilogx(v, exciseColumn(norm_max_r{i}'), 'b')
    xlabel('speed')
    ylabel('Response')
    title(ll{i})
    xlim([v(end) v(1)])
end

%% classification based on speed tunning
%% pca
mag_pca = MAG_all_norm{2};
% mag_pca = MAG_all_norm{1}(:, idx4);
mag_pca = mag_pca';
[id_sub, idx_sub] = deal(cell(2, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(mag_pca);
pc1 = 2; pc2 = 3;
plot(scores(:, pc1), scores(:, pc2), 'o')
hold on
for i = 1:3
    [x, y] = ginput;
    plot(x, y)
    IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
    [~, idx_sub{i}] = find(IN' == 1);
    id_sub{i} = id(idx_sub{i});
end
xlabel('1st Principal Component')
ylabel('2nd Principal Component')

for j = 1:3
    figure
    for i = 1:4
        v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
        subplot(2, 2, i)
        semilogx(v, exciseColumn(MAG_all_norm{i}(:, idx_sub{j})), 'b')
        xlabel('speed')
        ylabel('Response')
        title(ll{i})
        xlim([v(end) v(1)])
    end
end

%% peak_area

% % peak location
% v = datarun{1}.stimulus.params.SPATIAL_PERIOD./datarun{1}.stimulus.params.TEMPORAL_PERIOD;
% mag = MAG_all_norm{1};
% mag(:, any(isnan(mag))) = [];
% tpn = size(mag, 1);
% peak = v(mod(find(mag == 1), tpn));
% 
% % area under curve
% area = sum(mag) - mag(1, :)/2 - mag(end, :)/2;
% 
% figure
% scatter(peak, area)
% set(gca, 'XScale', 'log')
% xlabel('peak location')
% ylabel('area under curve')
% title('NDF3 SP60')
% 
%% concatenated tuning curves
% mag_all = cell2mat(MAG_all_norm);
% mag_all = mag_all(:, idx4);
% mag_all = mag_all';
% 
% [id_sub, idx_sub] = deal(cell(2, 1));
% 
% FigHandle = figure;
% set(FigHandle, 'Position', [1 1 380 400])
% 
% [~,scores,~,~] = princomp(mag_all);
% pc1 = 1; pc2 = 2;
% plot(scores(:, pc1), scores(:, pc2), 'o')
% hold on
% for i = 1:2
%     [x, y] = ginput;
%     plot(x, y)
%     IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
%     [~, idx] = find(IN' == 1);
%     idx_sub{i} = idx4(idx);
%     id_sub{i} = id4(idx);
% end
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
% 
% for j = 1:2
%     figure
%     for i = 1:4
%         v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
%         subplot(2, 2, i)
%         semilogx(v, exciseColumn(MAG_all_norm{i}(:, idx_sub{j})), 'b')
%         xlabel('speed')
%         ylabel('Response')
%         title(ll{i})
%         xlim([v(end) v(1)])
%     end
% end
% 
% 

%% plot average tunning curve
color = 'brk';
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:4
    v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    subplot(2, 2, i)
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
    [NumSpikesCell, StimComb] = get_spikescellstim(datarun{2},id_sub{ct},0);
    dscellanalysis(NumSpikesCell, StimComb);
end

%% compare across light level
SP = {'SP 60', 'SP 240'};
CT = {'on-off DSGC', 'on-DSGC'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:2
    v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    for j = 1:2
        subplot(2, 2, i+2*(j-1))
        errorbar(v, tuning_avg{i}(:, j), tuning_ste{i}(:, j), 'b')
        hold on
        errorbar(v, tuning_avg{i+2}(:, j), tuning_ste{i+2}(:, j), 'r')
        title([SP{i} '  ' CT{j}])
        set(gca, 'XScale', 'log')
        xlim([min(v) max(v)])
        xlabel('speed')
        ylabel('response')

    end
end
legend('NDF 3', 'NDF 0', 'location', 'northeast')


%% full field pulses

[raster_ff, raster_ff_all] = deal(cell(2, 1));
for d = 3:4
    [raster_ff{d-2}, raster_ff_all{d-2}] = get_ffp_raster(datarun{d}, id, 3);
end

for i = 1:length(id) 
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

%% Moving bar
[NumSpikesCell, StimComb] = get_spikescellstim_mb(datarun{6},datarun{6}.cell_ids,3653);
ds_struct = mbcellanalysis(NumSpikesCell, StimComb);
figure
plot(ds_struct.mag{2, 1}, ds_struct.mag{3, 1}, 'o')
title('data004 after mapping')
xlabel('speed 0.5')
ylabel('speed 1')
hold on
[x, y] = ginput;
plot(x, y);

IN = inpolygon(ds_struct.mag{2, 1}, ds_struct.mag{3, 1}, x, y);
[~, I] = find(IN == 1);
id = datarun{6}.cell_ids(I);
idx = get_cell_indices(datarun{6}, id);
I = zeros(length(id), n);
for i = 5:6
    I(:, i-4) = ismember(id, datarun{i}.cell_ids);
end
I = sum(I');
id_all = id(I == n); % DS cells that can be found in all light levels
idx_all = arrayfun(@(x) find(id == x), id_all);

raster_mb = cell(n, 1);
for i = 1:n    
    [NumSpikesCell, StimComb] = get_spikescellstim_mb(datarun{i+4},id,0);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb));
    raster_mb{i} = get_mb_raster(datarun{i+4}, id, 3653);
end



[idx, xx, yy] = subplot_idx(1, 1);

for cell_idx = 1:13;

for j = 1:7
    FigHandle = figure;
    set(FigHandle, 'Position', [1 1 700 700])
    for i = 2:9
        subplot(xx, yy, idx(i)); plot_raster(squeeze(raster_mb{2}{cell_idx}(1, j, i-1, :)), 0, 8)
    end
    title([num2str(id(cell_idx)) '_' num2str(j)])
    print_close(1, [12, 12], [num2str(id(cell_idx)) '_' num2str(j)])

end

end



dur = get_mb_trial_dur(datarun{5});



for cc = 45:50 %length(id)
    plot_ds_raster(DS, raster, cc, id(cc), ll, 1, 2, 1)
end

