%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
datarun{1} = load_data('/Analysis/xyao/2014-01-08-0/data000-map/data000-map', opt);
datarun{1}.names.stimulus_path = '/Analysis/xyao/2014-01-08-0/stimuli/s00';
datarun{1} = load_stim(datarun{1}, 'user_defined_trigger_interval', 10);

datarun{2} = load_data('/Analysis/xyao/2014-01-08-0/data001-map/data001-map', opt);
datarun{2}.names.stimulus_path = '/Analysis/xyao/2014-01-08-0/stimuli/s01';
datarun{2} = load_stim(datarun{2}, 'user_defined_trigger_interval', 10);

datarun{3} = load_data('/Analysis/xyao/2014-01-08-0/data002-map/data002-map', opt);
datarun{3}.names.stimulus_path = '/Analysis/xyao/2014-01-08-0/stimuli/s02';
datarun{3} = load_stim(datarun{3}, 'user_defined_trigger_interval', 10);

datarun{4} = load_data('/Analysis/xyao/2014-01-08-0/data004-map/data004-map', opt);
datarun{4}.names.stimulus_path = '/Analysis/xyao/2014-01-08-0/stimuli/s04';
datarun{4} = load_stim(datarun{4}, 'user_defined_trigger_interval', 10);

datarun{5} = load_data('/Analysis/xyao/2014-01-08-0/data005-map/data005-map', opt);
datarun{5}.names.stimulus_path = '/Analysis/xyao/2014-01-08-0/stimuli/s05';
datarun{5} = load_stim(datarun{5}, 'user_defined_trigger_interval', 10);

datarun{6} = load_data('/Analysis/xyao/2014-01-08-0/data006-map/data006-map', opt);
datarun{6}.names.stimulus_path = '/Analysis/xyao/2014-01-08-0/stimuli/s06';
datarun{6} = load_stim(datarun{6}, 'user_defined_trigger_interval', 10);

datarun{7} = load_data('/Analysis/xyao/2014-01-08-0/data003-map/data003-map', opt);
datarun{8} = load_data('/Analysis/xyao/2014-01-08-0/data007-map/data007-map', opt);


%% get DSGC ids

[NumSpikesCell, StimComb] = get_spikescellstim(datarun{3},datarun{3}.cell_ids,0);
[ds_struct] = dscellanalysis(NumSpikesCell, StimComb);

% pull out DS cells

figure
plot(ds_struct.mag{5, 1}, ds_struct.mag{6, 1}, 'o')
hold on
[x, y] = ginput;
plot(x, y);
xlabel('TP 60')
ylabel('TP 80')
title('Vector sum plot')

IN = inpolygon(ds_struct.mag{5, 1}, ds_struct.mag{6, 1}, x, y);
[~, I] = find(IN == 1);
id = datarun{3}.cell_ids(I);

I = zeros(length(id), 6);
for i = 1:6
    I(:, i) = ismember(id, datarun{i}.cell_ids);
end

I = sum(I');
id6 = id(I == 6); % DS cells that can be found in all light levels
idx6 = arrayfun(@(x) find(id == x), id6);

% I = sum(I(:, [3 6])');
% idx36 = find(I == 2);

[raster, raster_p_sum, p_idx] = deal(cell(6, 1));
for i = 1:6    
    [NumSpikesCell, StimComb] = get_spikescellstim(datarun{i},id,0);
    DS{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster{i} = get_ds_raster(datarun{i}, id);
end


param_p = 5; % choose which params to use to calculate prefer direction indices 
MAG_all_norm = cell(6, 1);
max_r = cell(6, 1);
norm_max_r = cell(6, 1);

for i = 1:6
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster{i}, DS{i}.angle{param_p});
    MAG_all_norm{i} = normalize_MAG(DS{i});
    rep = datarun{i}.stimulus.repetitions;
    max_r{i} = max_firing_rate(raster_p_sum{i}, 0.2, 8)/rep;
    norm_max_r{i} = max_r{i}./repmat(max(max_r{i}, [], 2), 1, size(max_r{i}, 2));
end

%% plot cell summary
for cc = 4:4 %length(id)
    plot_ds_raster(DS, raster, cc, id(cc), ll, 2, 3, 0)
end

%% plot single cell tuning curve
% choose the preferred direction under NDF 2, TP 30, as the preferred
% direction for plotting

% use unnormalized vector sum as response
figure
for i = 1:6
    v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    subplot(2, 3, i)
    semilogx(v, exciseColumn(MAG_all_norm{i}), 'b')
    xlabel('speed')
    ylabel('Response')
    title(ll{i})
    xlim([v(end) v(1)])
end


%% classify cells into two groups according to speed tunning
% pca
mag_pca = MAG_all_norm{3}';
% mag_pca = norm_max_r{1}(idx4, :);

[id_sub, idx_sub] = deal(cell(2, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(mag_pca);
pc1 = 1; pc2 = 2;
plot(scores(:, pc1), scores(:, pc2), 'o')
hold on
for i = 1:2
    [x, y] = ginput;
    plot(x, y)
    IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
    [~, idx_sub{i}] = find(IN' == 1);
    id_sub{i} = id(idx_sub{i});
end
xlabel('1st Principal Component')
ylabel('2nd Principal Component')

for j = 1:2
    figure
    for i = 1:6
        v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
        subplot(2, 3, i)
        semilogx(v, exciseColumn(MAG_all_norm{i}(:, idx_sub{j})), 'b')
        xlabel('speed')
        ylabel('Response')
        title(ll{i})
        xlim([v(end) v(1)])
    end
end

% vector sum plot of individual type
[NumSpikesCell, StimComb] = get_spikescellstim(datarun{3},id,0);
ds_temp = dscellanalysis(NumSpikesCell, StimComb);

%% full field pulses

dd = 8;
n = length(datarun{dd}.triggers)/4;
trigger = zeros(n, 1);
for i = 1:n
    trigger(i) = datarun{dd}.triggers(4*(i-1)+1);
end

% trigger = trigger - 1;

raster_ff = cell(length(id), 1);
for i = 1:length(id)
    if ~isempty(intersect(datarun{dd}.cell_ids, id(i)))
        idx = get_cell_indices(datarun{dd}, id(i));
        raster_ff{i, 1} = get_raster(datarun{dd}.spikes{idx}, trigger, 'plot', false);
    end
end

%% plot wn & full field flashes

a = [0 1 1 3 3 5 5 7 7 8]';
b = [1 1 2 2 1 1 0 0 1 1]';
id1ff = intersect(id1, datarun{dd}.cell_ids);
id2ff = intersect(id2, datarun{dd}.cell_ids);
idx1ff = arrayfun(@(x) find(id == x), id1ff);
idx2ff = arrayfun(@(x) find(id == x), id2ff);

for i = 7:7 %length(id1ff) 
    FigHandle = figure;
    set(FigHandle, 'Position', [1 1 400 800])
%     set(FigHandle,'yticklabel',[]);

    for j = 1:length(raster_ff{idx1ff(i)})
        SpikeTime = raster_ff{idx1ff(i)}{j};
        SpikeTime = SpikeTime';
        X = [SpikeTime; SpikeTime];
        Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
        line(X, Y, 'color', 'b');
        xlim([0 8])
        hold on
    end   
    line(a, b+57)
    xlabel('time/sec')
    title([num2str(id(idx1ff(i))) '  NDF 0'])
end

for i = 1:length(id2ff) 
    FigHandle = figure;
    set(FigHandle, 'Position', [1 1 400 800])
%     set(FigHandle,'yticklabel',[]);

    for j = 1:length(raster_ff{idx2ff(i)})
        SpikeTime = raster_ff{idx2ff(i)}{j};
        SpikeTime = SpikeTime';
        X = [SpikeTime; SpikeTime];
        Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
        line(X, Y, 'color', 'b');
        xlim([0 8])
        hold on
    end   
    line(a, b+57)
    xlabel('time/sec')
    title(num2str(id(idx2ff(i))))

end


id1wn = intersect(id1, datarun{7}.cell_ids);
id2wn = intersect(id2, datarun{7}.cell_ids);

figure;

for i = 1:length(id1wn)
    plot_autocorrelograms(datarun{7}, id1wn(i), 'foa', -1, 'clear_fig', false, 'normalize', true, 'line_color', 'b')
    hold on
end

% figure

for i = 1:length(id2wn)
    plot_autocorrelograms(datarun{7}, id2wn(i), 'foa', -1, 'clear_fig', false, 'normalize', true, 'line_color', 'r')
    hold on
end


%% ffp
[raster_ff, raster_ff_all] = deal(cell(2, 1));
for d = 7:8
    [raster_ff{d-6}, raster_ff_all{d-6}] = get_ffp_raster(datarun{d}, id, 2);
end

for i = 1:length(id) 
    if ~isempty(raster_ff{1}{i}) || ~isempty(raster_ff{2}{i})
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 800 800])
        for d = 1:2
        subplot(1, 2, d)
        plot_ffp(raster_ff{d}, raster_ff_all{d}, i, 2)
        title([num2str(id(i)) ' ' ll{3*d}(1:4)])
        end
        
        print_close(1, [12, 12], num2str(id(i)))
    end
end

%% 
d = 3;
t = 5;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
subplot(1, 2, 1)
idx_temp = idx_sub{1};
compass(DS{d}.U{t}(idx_temp), DS{d}.V{t}(idx_temp))
color = 'brgk';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DS{d}.U{t}(idx_temp), DS{d}.V{t}(idx_temp), x, y);
    [~, I] = find(IN == 1);
    idx_dir_oo{i} = idx_temp(I);
end

%%
trial_n = 50;
for cd = 1:2
[FFP.burst_spike_n{cd}, FFP.burst_spike_n_avg{cd}, FFP.burst_spike_n_ste{cd}, FFP.burst_start_time{cd}, FFP.burst_start_time_avg{cd}, FFP.burst_start_time_ste{cd}] = deal(cell(4, 1));
    for ct = 1:4
        cell_idx = idx_dir_oo{ct};
        for cc = 1:length(cell_idx)
            if ~isempty(raster_ff{cd}{cell_idx(cc)})
                for trial = 1:trial_n
                    raster_temp = raster_ff{cd}{cell_idx(cc)}{trial};
                    raster_temp = raster_temp(raster_temp > 9);
                    if ~isempty(raster_temp)
                        spike_interval = find(diff(raster_temp) > 0.2);
                        start_time_i = [1; spike_interval+1];
                        start_time = raster_temp(start_time_i);
                        spike_n = diff([0; spike_interval; length(raster_temp)]);
                        FFP.time{cd}{ct}{cc}{trial} = start_time;
                        FFP.spike_n{cd}{ct}{cc}{trial} = spike_n;
                    end
                end

                burst_n = cellfun(@(x)(length(x)), FFP.time{cd}{ct}{cc});
                max_burst_n = max(burst_n);
                for burst = 1:max_burst_n
                    FFP.burst_start_time{cd}{ct}{cc}{burst} = cellfun(@(x)(x(burst)), FFP.time{cd}{ct}{cc}(burst_n >= burst));
                    FFP.burst_spike_n{cd}{ct}{cc}{burst} = cellfun(@(x)(x(burst)), FFP.spike_n{cd}{ct}{cc}(burst_n >= burst));
                    FFP.burst_start_time_avg{cd}{ct}{cc}(burst) = mean(FFP.burst_start_time{cd}{ct}{cc}{burst});
                    FFP.burst_spike_n_avg{cd}{ct}{cc}(burst) = mean(FFP.burst_spike_n{cd}{ct}{cc}{burst});
                    FFP.burst_start_time_ste{cd}{ct}{cc}(burst) = std(FFP.burst_start_time{cd}{ct}{cc}{burst})/sqrt(length(FFP.burst_start_time{cd}{ct}{cc}{burst}));
                    FFP.burst_spike_n_ste{cd}{ct}{cc}(burst) = std(FFP.burst_spike_n{cd}{ct}{cc}{burst})/sqrt(length(FFP.burst_spike_n{cd}{ct}{cc}{burst}));
                end


            end
        end
    end
end


dir = {'Inferior', 'A(P)', 'Superior', 'P(A)'};
condition = {'NDF 3', 'NDF 0'};
color = 'brgmckybrgmcky';
figure
for cd = 1:2
    h = subplot(1, 2, cd);
    for ct = 1:4
        if ~isempty(FFP.burst_start_time_avg{cd}{ct})
            
            for cc = 1:length(FFP.burst_start_time_avg{cd}{ct});
            %     x = FFP.burst_start_time_avg{cd}{cc};
            %     x_bar = FFP.burst_start_time_ste{cd}{cc};
            %     y = FFP.burst_spike_n_avg{cd}{cc};
            %     y_bar = FFP.burst_spike_n_ste{cd}{cc};


            %     errorbar(x, y, y_bar, y_bar, 's-')
            %     hold on
                if ~isempty(FFP.burst_start_time_avg{cd}{ct}{cc})
                    x = ct;
                    y = FFP.burst_start_time_avg{cd}{ct}{cc};
                    y_bar = FFP.burst_start_time_ste{cd}{ct}{cc};

                   for i = 1:length(y)
                       errorbar(x, y(i), y_bar(i), y_bar(i), 'Marker', 's', 'color', color(i))
                       hold on
                   end
                end
            end
        end
    end
    xlim([0.5 4.5])
    set(h,'XTicklabel',dir)                
    ylabel('time (second)')  
    title(condition{cd})
end
