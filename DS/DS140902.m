%% load data
opt = struct('load_params', 1,'load_neurons', 1); %, 'load_ei', 1);

% load data
datarun{1} = load_data('/Analysis/xyao/2014-09-02-0/data002-map/data002-map', opt);
datarun{1}.names.stimulus_path = '/Analysis/xyao/2014-09-02-0/stimuli/s02';
datarun{1} = load_stim(datarun{1}, 'user_defined_trigger_interval', 10);

datarun{2} = load_data('/Analysis/xyao/2014-09-02-0/data004-map-self/data004-map-self', opt);
datarun{2}.names.stimulus_path = '/Analysis/xyao/2014-09-02-0/stimuli/s04';
datarun{2} = load_stim(datarun{2}, 'user_defined_trigger_interval', 10);

datarun{3} = load_data('/Analysis/xyao/2014-09-02-0/data003-map/data003-map', opt);
datarun{4} = load_data('/Analysis/xyao/2014-09-02-0/data005-map/data005-map', opt);


[NumSpikesCell, StimComb] = get_spikescellstim(datarun{2},datarun{2}.cell_ids,0);
[ds_struct] = dscellanalysis(NumSpikesCell, StimComb);

% pull out DS cells

figure
plot(ds_struct.mag{5}, ds_struct.mag{6}, 'o')
title('data004 after mapping')
xlabel('TP 24')
ylabel('TP 60')
hold on
[x, y] = ginput;
plot(x, y);

IN = inpolygon(ds_struct.mag{5}, ds_struct.mag{6}, x, y);
[~, I] = find(IN == 1);
id = datarun{2}.cell_ids(I);

figure
scatter(ds_struct.mag{5}, ds_struct.mag{6})
hold on
scatter(ds_struct.mag{5}(I), ds_struct.mag{6}(I), 'r')

title('data004 after mapping')
xlabel('TP 60')
ylabel('TP 120')


I = zeros(length(id), 2);
for i = 1:2
    I(:, i) = ismember(id, datarun{i}.cell_ids);
end
I = sum(I');
id2 = id(I == 2); % DS cells that can be found in all light levels
idx2 = arrayfun(@(x) find(id == x), id2);

[raster, raster_p_sum, p_idx] = deal(cell(1, 1));
for i = 1:2    
    [NumSpikesCell, StimComb] = get_spikescellstim(datarun{i},id,0);
    DS{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster{i} = get_ds_raster(datarun{i}, id);
end

param_p = 5; % choose which params to use to calculate prefer direction indices 
MAG_all_norm = cell(2, 1);
max_r = cell(2, 1);
norm_max_r = cell(2, 1);

for i = 1:2
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster{i}, DS{i}.angle{param_p});
    MAG_all_norm{i} = normalize_MAG(DS{i});
    rep = datarun{i}.stimulus.repetitions;
    max_r{i} = max_firing_rate(raster_p_sum{i}, 0.2, 8)/rep;
    norm_max_r{i} = max_r{i}./repmat(max(max_r{i}, [], 2), 1, size(max_r{i}, 2));
end

ll = {'NDF3 SP60', 'NDF0 SP60'};


%% plot cell summary

for cc = 2:2 %length(id)
    plot_ds_raster(DS, raster, cc, id(cc), ll, 1, 2, 0)
end

%% plot single cell tuning curve

% use unnormalized vector sum as response
figure
for i = 1:2
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
for i = 1:4
    v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    subplot(2, 2, i)
    semilogx(v, exciseColumn(norm_max_r{i}'), 'b')
    xlabel('speed')
    ylabel('Response')
    title(ll{i})
    xlim([v(end) v(1)])
end

%% classify cells into two groups according to speed tunning
% pca
mag_pca = MAG_all_norm{2};
% mag_pca = MAG_all_norm{1}(:, idx4);
mag_pca = mag_pca';
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

figure
scatter(scores(:, pc1), scores(:, pc2))
hold on
scatter(scores(idx_sub{2}, pc1), scores(idx_sub{2}, pc2), 'r')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')



for j = 1:2
    figure
    for i = 1:2
        v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
        subplot(1, 2, i)
        semilogx(v, exciseColumn(MAG_all_norm{i}(:, idx_sub{j})), 'b')
        xlabel('speed')
        ylabel('Response')
        title(ll{i})
        xlim([v(end) v(1)])
    end
end



%% plot average tunning curve
color = 'brk';
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:2
    v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
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
legend('on-off DSGC', 'on DSGC')

%% vector sum plot of individual type
for ct = 1:2
    [NumSpikesCell, StimComb] = get_spikescellstim(datarun{2},id_sub{ct},0);
    dscellanalysis(NumSpikesCell, StimComb);
end

%% compare across light level
SP = {'SP 60'};
CT = {'on-off DSGC', 'on-DSGC'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:1
    v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    for j = 1:2
        subplot(1, 2, j)
        errorbar(v, tuning_avg{i}(:, j), tuning_ste{i}(:, j), 'b')
        hold on
        errorbar(v, tuning_avg{i+1}(:, j), tuning_ste{i+1}(:, j), 'r')
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

for i = 1:5 %length(id) 
    if ~isempty(raster_ff{1}{i}) || ~isempty(raster_ff{2}{i})
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 800 800])
        for d = 1:2
        subplot(1, 2, d)
        plot_ffp(raster_ff{d}, raster_ff_all{d}, i, 3)
        title([num2str(id(i)) ' ' ll{d}(1:4)])
        end
        
%         print_close(1, [12, 12], num2str(id(i)))
    end
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


%% classify DSGC into subtypes (directions)
d = 2;
t = 4;
h = figure;
dirn = 4;
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
title('on DSGC')
