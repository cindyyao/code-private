%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);


datadg = load_data('/Volumes/lab/analysis/2016-09-04-0/data005-sorted/data005-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-09-04-0/stimuli/s05.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

% datarun = load_data('/Volumes/lab/analysis/2016-09-04-0/data000-003-clean/data000-003-clean', opt);
datarun = load_data('/Volumes/lab/analysis/2016-09-04-0/data000-003-map/data000-003-map', opt);
time_points = [1842 3742 5642];
datamb(1:4) = split_datarun(datarun, time_points);
datamb{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-09-04-0/stimuli/s00.txt';
datamb{1} = load_stim(datamb{1}, 'user_defined_trigger_set', [1:2:1120]);
datamb{1}.stimulus.triggers = datamb{1}.stimulus.triggers';
datamb{1} = convert_stim_160904(datamb{1});
datamb{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-09-04-0/stimuli/s01.txt';
datamb{2} = load_stim(datamb{2}, 'user_defined_trigger_set', [1:2:1120]);
datamb{2}.stimulus.triggers = datamb{2}.stimulus.triggers';
datamb{2} = convert_stim_160904(datamb{2});
datamb{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-09-04-0/stimuli/s02.txt';
datamb{3} = load_stim(datamb{3}, 'user_defined_trigger_set', [1:2:1120]);
datamb{3}.stimulus.triggers = datamb{3}.stimulus.triggers';
datamb{4}.names.stimulus_path = '/Volumes/lab/analysis/2016-09-04-0/stimuli/s03.txt';
datamb{4} = load_stim(datamb{4}, 'user_defined_trigger_set', [1:2:1120]);
datamb{4}.stimulus.triggers = datamb{4}.stimulus.triggers';


[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [3 5]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

load('DS160904.mat')
%% 
datarun = load_data('/Volumes/lab/analysis/2016-09-04-0/data004/data004', opt);
datarun = load_sta(datarun);

%% mb classification
[NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{1},datamb{1}.cell_ids,1775.2,0.025);
mbstruct = mbcellanalysis(NumSpikesCell, StimComb, datamb{1});
%% dg
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);


delta_p = 4; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(DG{i});
rep = datadg.stimulus.repetitions;

%% mb
load('DS160904.mat')
n = 4;
duration = [1775.2 1776.1 1776.2 1776.2]; %sec
bin_size = 0.025; %sec
[raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},ds_id,duration(i),bin_size);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb{i} = get_mb_raster(datamb{i}, ds_id, duration(i));
    for j = 1:length(raster_mb{i})
        if(mb_idx(j))
            raster_mb{i}{j} = [];
        end
    end
    trial_dur{i} = get_mb_trial_dur(datamb{i}, 400, 400, 0.5);
end

ctr_p = 7; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_mb = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}, raster_p_sum_all{i}] = get_pdirection_raster(raster_mb{i}, MB{1}.angle{ctr_p});
    [raster_n_sum{i}, n_idx{i}, raster_n_sum_all{i}] = get_ndirection_raster(raster_mb{i}, MB{1}.angle{ctr_p});
    MAG_all_norm_mb{i} = normalize_MAG(MB{i});
    rep = datamb{i}.stimulus.repetitions;
end
close all
%% plot cell summary
for cc = 30:30%length(ds_id)
    plot_mb_raster_ctr(MB, raster_mb, trial_dur, cc, ds_id(cc), 'NDF0', 4, 7, 1)
end

%% classification based on speed tunning
L = 1;
mag_pca = MAG_all_norm_dg{L};
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
    id_sub{i} = ds_id(idx_sub{i});
end
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
title('NDF 0')

figure
plot(scores(idx_sub{1}, pc1), scores(idx_sub{1}, pc2), 'ro', scores(idx_sub{2}, pc1), scores(idx_sub{2}, pc2), 'bo')
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
title('NDF 0')

v = 4*datadg.stimulus.params.SPATIAL_PERIOD./datadg.stimulus.params.TEMPORAL_PERIOD;
figure
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
xlim([v(end) v(1)])

t = 3;
figure
compass(DG{1}.U{t}(idx_sub{1}), DG{1}.V{t}(idx_sub{1}), 'r')
hold on
compass(DG{1}.U{t}(idx_sub{2}), DG{1}.V{t}(idx_sub{2}), 'b')

%%
d = 1;
t = 5;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(DG{d}.U{t}(idx_sub{2}), DG{d}.V{t}(idx_sub{2}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}(idx_sub{2}), DG{d}.V{t}(idx_sub{2}), x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = idx_sub{2}(I);
    id_dir{i} = ds_id(idx_dir{i});
    
%     idx_temp = idx_sub{2}(I);
%     id_temp = ds_id(idx_temp);
end

d = 1;
t = 3;
h = figure;
dirn = 3;
set(h, 'Position', [1 1 1080 500])
compass(DG{d}.U{t}(idx_sub{1}), DG{d}.V{t}(idx_sub{1}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}(idx_sub{1}), DG{d}.V{t}(idx_sub{1}), x, y);
    [~, I] = find(IN == 1);
    idx_dir_on{i} = idx_sub{1}(I);
    id_dir_on{i} = ds_id(idx_dir_on{i});
end

%% 
for dir = 1:4
    id_dir_mb{dir} = id_dir{dir}(~mb_idx(idx_dir{dir}));
    idx_dir_mb{dir} = idx_dir{dir}(~mb_idx(idx_dir{dir}));
end
%% contrast response function (spike count)

% DS cell
ctr_x = [5 10 20 40 80 150 300];
color = 'brgkc';
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'AP5+Hex', 'wash'};
for drug = 1:4
    for cc = 1:length(raster_p_sum{1})
        if ~isempty(raster_p_sum{drug}{cc})
            pd_spikes{drug}(cc,:) = cellfun('length', raster_p_sum{drug}{cc})/datamb{drug}.stimulus.repetitions;
        end
    end
    
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                spike_temp = cellfun('length', raster_p_sum{drug}{idx_dir{dir}(cc)})/datamb{drug}.stimulus.repetitions;
                if sum(spike_temp) ~= 0
                    pd_dir_spikes{drug}{dir}(CC,:) = cellfun('length', raster_p_sum{drug}{idx_dir{dir}(cc)})/datamb{drug}.stimulus.repetitions;
%                     pd_dir_spikes_nor{drug}{dir}(CC,:) = pd_dir_spikes{drug}{dir}(CC,:)/max(pd_dir_spikes{1}{dir}(CC,:));
                    pd_dir_spikes_nor{drug}{dir}(CC,:) = pd_dir_spikes{drug}{dir}(CC,:)/max(pd_dir_spikes{drug}{dir}(CC,:));
                    CC = CC + 1;
                end
            end
        end
        pd_dir_spikes_nor{drug}{dir} = nan2empty(pd_dir_spikes_nor{drug}{dir});
    end
    for dir = 1:2
        CC = 1;
        for cc = 1:length(idx_dir_on{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir_on{dir}(cc)})
                spike_temp = cellfun('length', raster_p_sum{drug}{idx_dir_on{dir}(cc)})/datamb{drug}.stimulus.repetitions;
                if sum(spike_temp) ~= 0
                    pd_dir_on_spikes{drug}{dir}(CC,:) = cellfun('length', raster_p_sum{drug}{idx_dir_on{dir}(cc)})/datamb{drug}.stimulus.repetitions;
%                     pd_dir_on_spikes_nor{drug}{dir}(CC,:) = pd_dir_on_spikes{drug}{dir}(CC,:)/max(pd_dir_on_spikes{1}{dir}(CC,:));
                    pd_dir_on_spikes_nor{drug}{dir}(CC,:) = pd_dir_on_spikes{drug}{dir}(CC,:)/max(pd_dir_on_spikes{drug}{dir}(CC,:));
                    CC = CC + 1;
                end
            end
        end
        pd_dir_on_spikes_nor{drug}{dir} = nan2empty(pd_dir_on_spikes_nor{drug}{dir});
    end

end

for drug = 1:4
    figure
    for dir = 1:4
        errorbar(ctr_x, mean(pd_dir_spikes_nor{drug}{dir}), std(pd_dir_spikes_nor{drug}{dir})/sqrt(size(pd_dir_spikes_nor{drug}{dir}, 1)), 'color', color(dir));
        hold on
    end
    legend(dscell_type)
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('normalized response')
    title(condition{drug})
end

% for drug = 1:4
%     figure
%     for dir = 1:2
%         errorbar(ctr_x, mean(pd_dir_on_spikes_nor{drug}{dir}), std(pd_dir_on_spikes_nor{drug}{dir})/sqrt(size(pd_dir_on_spikes_nor{drug}{dir}, 1)), 'color', color(dir));
%         hold on
%     end
%     legend(dscell_type)
%     set(gca, 'Xscale', 'log')
%     xlabel('% contrast')
%     ylabel('normalized response')
%     title(condition{drug})
% end

figure
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4
        errorbar(ctr_x, mean(pd_dir_spikes{drug}{dir}), std(pd_dir_spikes{drug}{dir})/sqrt(size(pd_dir_spikes{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition)
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike count')
    title(dscell_type{dir})
    xlim([3 400])
end

% for dir = 1:4
%     for cc = 1:size(pd_dir_spikes{}
%     figure(1)
%     for drug = 1:4
%         plot(ctr_x, pd_dir_on_spikes{drug}{dir}
%% CRF (Peak firing rate)
bin_size = 0.1;
xx = bin_size/2:bin_size:2.6-bin_size/2;
ctr_x = [5 10 20 40 80 150 300];
color = 'brgkc';
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'AP5+Hex', 'wash'};
for drug = 1:4    
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                a = raster_p_sum{drug}{idx_dir{dir}(cc)};
                hist_temp = cellfun(@(a) hist(a, xx), a, 'UniformOutput', false);
                hist_temp = cell2mat(squeeze(cellfun(@(hist_temp) max(hist_temp), hist_temp, 'UniformOutput', false)))/datamb{drug}.stimulus.repetitions/bin_size;
                if sum(hist_temp) ~= 0
                    pd_dir_spikes{drug}{dir}(CC,:) = hist_temp;
%                     pd_dir_spikes_nor{drug}{dir}(CC,:) = pd_dir_spikes{drug}{dir}(CC,:)/max(pd_dir_spikes{1}{dir}(CC,:));
                    pd_dir_spikes_nor{drug}{dir}(CC,:) = pd_dir_spikes{drug}{dir}(CC,:)/max(pd_dir_spikes{drug}{dir}(CC,:));
                    CC = CC + 1;
                end
            end
        end
        pd_dir_spikes_nor{drug}{dir} = nan2empty(pd_dir_spikes_nor{drug}{dir});
    end
    for dir = 1:2
        CC = 1;
        for cc = 1:length(idx_dir_on{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir_on{dir}(cc)})
                a = raster_p_sum{drug}{idx_dir_on{dir}(cc)};
                hist_temp = cellfun(@(a) hist(a, xx), a, 'UniformOutput', false);
                hist_temp = cell2mat(squeeze(cellfun(@(hist_temp) max(hist_temp), hist_temp, 'UniformOutput', false)))/datamb{drug}.stimulus.repetitions/bin_size;
                if sum(hist_temp) ~= 0
                    pd_dir_on_spikes{drug}{dir}(CC,:) = hist_temp;
                    pd_dir_on_spikes_nor{drug}{dir}(CC,:) = pd_dir_on_spikes{drug}{dir}(CC,:)/max(pd_dir_on_spikes{1}{dir}(CC,:));
                    CC = CC + 1;
                end
            end
        end
        pd_dir_on_spikes_nor{drug}{dir} = nan2empty(pd_dir_on_spikes_nor{drug}{dir});
    end

end

for drug = 1:4
    figure
    for dir = 1:4
        errorbar(ctr_x, mean(pd_dir_spikes_nor{drug}{dir}), std(pd_dir_spikes_nor{drug}{dir})/sqrt(size(pd_dir_spikes_nor{drug}{dir}, 1)), 'color', color(dir));
        hold on
    end
    legend(dscell_type)
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('normalized response')
    title(condition{drug})
end

% for drug = 1:4
%     figure
%     for dir = 1:2
%         errorbar(ctr_x, mean(pd_dir_on_spikes_nor{drug}{dir}), std(pd_dir_on_spikes_nor{drug}{dir})/sqrt(size(pd_dir_on_spikes_nor{drug}{dir}, 1)), 'color', color(dir));
%         hold on
%     end
%     legend(dscell_type)
%     set(gca, 'Xscale', 'log')
%     xlabel('% contrast')
%     ylabel('normalized response')
%     title(condition{drug})
% end

figure
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4
        errorbar(ctr_x, mean(pd_dir_spikes{drug}{dir}), std(pd_dir_spikes{drug}{dir})/sqrt(size(pd_dir_spikes{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition)
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike count')
    title(dscell_type{dir})
    xlim([3 400])
end


%% CRF (Peak firing rate, seperate ON OFF)
bin_size = 0.05;
xx = bin_size/2:bin_size:3.2-bin_size/2;
ctr_x = [5 10 20 40 80 150 300];
bar_t = datamb{1}.stimulus.params.BAR_WIDTH/datamb{1}.stimulus.params.DELTA/60;
bar_bin = ceil(bar_t/bin_size);
color = 'brgkc';
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'AP5+Hex', 'wash'};
for drug = 1:4    
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                a = raster_p_sum{drug}{idx_dir{dir}(cc)};
                hist_temp = cellfun(@(a) hist(a, xx), a, 'UniformOutput', false);
%                 if drug == 1
%                     plot(hist_temp{7})
%                     [x,~] = ginput;
%                     x = round(x);
%                     X{dir}{CC} = x;
%                 end
                hist_temp = cell2mat(squeeze(hist_temp));
%                 hist_temp_on = max(hist_temp(:, max(x-bar_bin,1):x), [], 2)/datamb{drug}.stimulus.repetitions/bin_size;
%                 hist_temp_off = max(hist_temp(:, x+1:x+bar_bin+1), [], 2)/datamb{drug}.stimulus.repetitions/bin_size;
                hist_temp_on = max(hist_temp(:, 1:X{dir}{CC}), [], 2)/datamb{drug}.stimulus.repetitions/bin_size;
                hist_temp_off = max(hist_temp(:, X{dir}{CC}+1:end), [], 2)/datamb{drug}.stimulus.repetitions/bin_size;
                
%                 if sum(hist_temp) ~= 0
                pd_dir_spikes_on{drug}{dir}(CC,:) = hist_temp_on;
                pd_dir_spikes_off{drug}{dir}(CC,:) = hist_temp_off;
%                     pd_dir_spikes_nor{drug}{dir}(CC,:) = pd_dir_spikes{drug}{dir}(CC,:)/max(pd_dir_spikes{1}{dir}(CC,:));
                pd_dir_spikes_nor_on{drug}{dir}(CC,:) = pd_dir_spikes_on{drug}{dir}(CC,:)/max(pd_dir_spikes_on{drug}{dir}(CC,:));
                pd_dir_spikes_nor_off{drug}{dir}(CC,:) = pd_dir_spikes_off{drug}{dir}(CC,:)/max(pd_dir_spikes_off{drug}{dir}(CC,:));
                CC = CC + 1;
%                 end
            end
        end
%         pd_dir_spikes_nor_on{drug}{dir} = nan2empty(pd_dir_spikes_nor_on{drug}{dir});
%         pd_dir_spikes_nor_off{drug}{dir} = nan2empty(pd_dir_spikes_nor_off{drug}{dir});
    end

end

% for drug = 1:4
%     figure
%     for dir = 1:4
%         errorbar(ctr_x, mean(pd_dir_spikes_nor{drug}{dir}), std(pd_dir_spikes_nor{drug}{dir})/sqrt(size(pd_dir_spikes_nor{drug}{dir}, 1)), 'color', color(dir));
%         hold on
%     end
%     legend(dscell_type)
%     set(gca, 'Xscale', 'log')
%     xlabel('% contrast')
%     ylabel('normalized response')
%     title(condition{drug})
% end

% for drug = 1:4
%     figure
%     for dir = 1:2
%         errorbar(ctr_x, mean(pd_dir_on_spikes_nor{drug}{dir}), std(pd_dir_on_spikes_nor{drug}{dir})/sqrt(size(pd_dir_on_spikes_nor{drug}{dir}, 1)), 'color', color(dir));
%         hold on
%     end
%     legend(dscell_type)
%     set(gca, 'Xscale', 'log')
%     xlabel('% contrast')
%     ylabel('normalized response')
%     title(condition{drug})
% end

figure
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4
        errorbar(ctr_x, mean(pd_dir_spikes_on{drug}{dir}), std(pd_dir_spikes_on{drug}{dir})/sqrt(size(pd_dir_spikes_on{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
    xlim([3 400])
end

figure
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4
        errorbar(ctr_x, mean(pd_dir_spikes_off{drug}{dir}), std(pd_dir_spikes_off{drug}{dir})/sqrt(size(pd_dir_spikes_off{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
    xlim([3 400])
end

%% CRF (1 second before minus 1 second after)
bin_size = 0.025;
xx = bin_size/2:bin_size:3.2-bin_size/2;
ctr_x = [5 10 20 40 80 150 300];
bar_t = datamb{1}.stimulus.params.BAR_WIDTH/datamb{1}.stimulus.params.DELTA/60;
bar_bin = ceil(bar_t/bin_size);
color = 'brgkc';
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'AP5+Hex', 'wash'};

for drug = 1:4    
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                t = T{dir}{CC}*bin_size;
                for ctr = 1:length(datamb{1}.stimulus.params.RGB)
                    spikes = raster_p_sum{drug}{idx_dir{dir}(cc)}{1,1,ctr};
                    spikes1 = spikes(spikes < t & spikes > t-1);
                    spikes2 = spikes(spikes > t & spikes < t+1);
%                     spikes_delta{drug}{dir}(CC, ctr) = abs(length(spikes1) - length(spikes2));
                    spikes_delta{drug}{dir}(CC, ctr) = length(spikes1);
                end
                CC = CC + 1;
            end
        end
        spikes_delta_norm{drug}{dir} = spikes_delta{drug}{dir}./repmat(max(spikes_delta{1}{dir}, [], 2), 1, length(datamb{1}.stimulus.params.RGB));
    end
end
                    
figure
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4
        errorbar(ctr_x, mean(spikes_delta_norm{drug}{dir}), std(spikes_delta_norm{drug}{dir})/sqrt(size(spikes_delta_norm{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
    xlim([3 400])
end

%% fit ds
for dir = 1:4
    for drug = 1:4
        figure
        set(gcf, 'Position', [1 1 900 800])
        for cc = 1:size(spikes_delta_norm{drug}{dir}, 1)

            ydata = spikes_delta_norm{drug}{dir}(cc, :);
            xdata = log10(ctr_x);
            [f, G] = fit_nr(xdata, ydata, 'upper', [2, 100, log10(300), min(ydata)]);
            fit_all{drug}{dir}{cc} = f;
            G_all{drug}{dir}{cc} = G;

            x = linspace(min(xdata), max(xdata), 100);
            y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;

            subplot(4,4,cc)
            plot(xdata, ydata)
            hold on
            plot(x, y)

            sigma{dir}(drug, cc) = f.sigma;
            ymax{dir}(drug, cc) = f.ymax;
            aa{dir}(drug, cc) = f.a;
            bb{dir}(drug, cc) = f.b;
        end
    end
end

%% get spontaneous activity
for drug = 1:4
    for dir = 1:4
        for cc = 1:length(id_dir_mb{dir})
            idx = get_cell_indices(datamb{drug},id_dir_mb{dir}(cc));
            spikes_temp = datamb{drug}.spikes{idx};
            bgnd_firing{dir}(drug, cc) = length(spikes_temp(spikes_temp < 1830 & spikes_temp > 1780))/50;
        end
    end
end

%%

for dir = 1:4
    for drug = 1:4
        for cc = 1:size(spikes_delta_norm{drug}{dir}, 1)
            x = linspace(sigma{1}{dir}(cc)-0.9, sigma{1}{dir}(cc)+1.1, 10);
            y = ymax{drug}{dir}(cc)*x.^aa{drug}{dir}(cc)./(x.^aa{drug}{dir}(cc) + sigma{drug}{dir}(cc)^aa{drug}{dir}(cc))+bb{drug}{dir}(cc);
            delta_y{drug}{dir}(cc, :) = y;
%             pause
        end
    end
end
        
figure
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4
        errorbar(linspace(0.1,2.1,10), mean(delta_y{drug}{dir}), std(delta_y{drug}{dir})/sqrt(size(delta_y{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
%     set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
%     xlim([3 400])
end

%% fit ds (fit before normalization)
for dir = 1:4
    for drug = 1:4
        figure
        set(gcf, 'Position', [1 1 900 800])
        for cc = 1:size(spikes_delta{drug}{dir}, 1)

            ydata = spikes_delta{drug}{dir}(cc, :);
            xdata = log10(ctr_x);
            [f, G] = fit_nr(xdata, ydata, 'upper', [500, 100, log10(300), min(ydata)]);
            fit_all{drug}{dir}{cc} = f;
            G_all{drug}{dir}{cc} = G;

            x = linspace(min(xdata), max(xdata), 100);
            y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;

            subplot(4,4,cc)
            plot(xdata, ydata)
            hold on
            plot(x, y)

            sigma{drug}{dir}(cc) = f.sigma;
            ymax{drug}{dir}(cc) = f.ymax;
            aa{drug}{dir}(cc) = f.a;
            bb{drug}{dir}(cc) = f.b;
        end
    end
end


for dir = 1:4
    for drug = 1:4
        for cc = 1:size(spikes_delta_norm{drug}{dir}, 1)
            x = linspace(sigma{1}{dir}(cc)-0.9, sigma{1}{dir}(cc)+1.1, 10);
            y = ymax{drug}{dir}(cc)*x.^aa{drug}{dir}(cc)./(x.^aa{drug}{dir}(cc) + sigma{drug}{dir}(cc)^aa{drug}{dir}(cc))+bb{drug}{dir}(cc);
            delta_y{drug}{dir}(cc, :) = y;
            delta_y_norm{drug}{dir}(cc, :) = y/max(delta_y{1}{dir}(cc, :));
%             pause
        end
    end
end
        
figure
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:3
        errorbar(linspace(0.1,2.1,10), mean(delta_y_norm{drug}{dir}), std(delta_y_norm{drug}{dir})/sqrt(size(delta_y_norm{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
%     set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
%     xlim([3 400])
end
%%
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'AP5+Hex', 'wash'};
step_size = 60;
for drug = 1:4
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                for ctr = 7:-1:1
                    a = raster_p_sum{drug}{idx_dir{dir}(cc)}{ctr};
                    hist_temp = hist(a, xx);
                    if drug == 1 && ctr == 7
                        [max_p, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
                        Max_i{dir}(CC) = max_i;
                    else
                        max_p = conv(hist_temp, ones(1,step_size), 'valid');
                        max_p = max_p(Max_i{dir}(CC));
                    end
                    a = raster_n_sum{drug}{idx_dir{dir}(cc)}{ctr};
                    mean_n = length(a)/trial_dur*bin_size*step_size;
%                     response_pn{drug}{dir}(CC, ctr) = max(max_p - mean_n, 0);
%                     response_pn{drug}{dir}(CC, ctr) = max_p - mean_n;
%                     response_pn{drug}{dir}(CC, ctr) = abs(max_p - mean_n);
                    response_pn{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
                end
                CC = CC + 1;
            end
        end
        response_pn_norm{drug}{dir} = response_pn{drug}{dir}./repmat(max(response_pn{1}{dir}, [], 2), 1, size(response_pn{drug}{dir},2));
    end
end


ctr_x = [5 10 20 40 80 150 300];
color = 'brgkc';
figure
set(gcf, 'Position', [1 1 900 800])

for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4
        errorbar(ctr_x, mean(response_pn_norm{drug}{dir}), std(response_pn_norm{drug}{dir})/sqrt(size(response_pn_norm{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
%     xlim([3 400])
end

% fit ds
for dir = 1:4
    for drug = 1:4
        figure
        set(gcf, 'Position', [1 1 900 800])
        for cc = 1:size(response_pn_norm{drug}{dir}, 1)

            ydata = response_pn_norm{drug}{dir}(cc, :);
            xdata = log10(ctr_x);
            [f, G] = fit_nr(xdata, ydata, 'upper', [1.2, 100, log10(300), 0]);
            fit_all{drug}{dir}{cc} = f;
            G_all{drug}{dir}{cc} = G;

            x = linspace(min(xdata), max(xdata), 100);
            y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;

            subplot(4,4,cc)
            plot(xdata, ydata)
            hold on
            plot(x, y)

            sigma{dir}(drug, cc) = f.sigma;
            ymax{dir}(drug, cc) = f.ymax;
            aa{dir}(drug, cc) = f.a;
            bb{dir}(drug, cc) = f.b;
        end
    end
end


for dir = 1:4
    for drug = 1:4
        for cc = 1:size(response_pn_norm{drug}{dir}, 1)
            x = linspace(sigma{dir}(1,cc)-0.9, sigma{dir}(1, cc)+1.1, 10);
            y = ymax{dir}(drug, cc)*x.^aa{dir}(drug, cc)./(x.^aa{dir}(drug, cc) + sigma{dir}(drug, cc)^aa{dir}(drug, cc))+bb{dir}(drug, cc);
            pn_y{drug}{dir}(cc, :) = y;
%             pause
        end
    end
end
        
figure
set(gcf, 'Position', [1 1 900 800])
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4 
        errorbar(linspace(0.1,2.1,10), mean(pn_y{drug}{dir}), std(pn_y{drug}{dir})/sqrt(size(pn_y{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
%     set(gca, 'Xscale', 'log')
    xlabel('normalized contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
%     xlim([3 400])
end


            

% end  


figure
set(gcf, 'Position', [1 1 900 800])
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4
        for cc = 1:size(response_pn_norm{drug}{dir}, 1)
%             plot(log10(ctr_x), response_pn_norm{drug}{dir}(cc, :), color(drug))
            plot(log10(ctr_x)-(sigma{dir}(1,cc)-1), response_pn_norm{drug}{dir}(cc, :), [color(drug) 'o'])
            hold on
        end
    end
end

for dir = 1:4
    sigma_mean{dir} = mean(10.^sigma{dir}');
    sigma_ste{dir} = std(10.^sigma{dir}')/sqrt(size(sigma{dir}, 2));
end
%%
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = dscell_type;
model_series = cell2mat(sigma_mean');
model_error = cell2mat(sigma_ste');
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('C50')
title('NDF 3')
legend('control','AP5', 'AP5+HEX', 'wash', 'location', 'northeast');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end



%%
color = {[0,0,0]/255,[0,191,255]/255, [255,20,60]/255, [0.5 0.5 0.5]};
figure
set(gcf, 'Position', [1 1 900 800])
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:3
        ydata = mean(response_pn_norm{drug}{dir});
        yste = std(response_pn_norm{drug}{dir})/sqrt(size(response_pn_norm{drug}{dir}, 1));
        xdata = log10(ctr_x);
        [f, G] = fit_nr(xdata, ydata, 'upper', [2, 100, log10(300), min(ydata)]);

        x = linspace(min(xdata), max(xdata)+0.5, 100);
        y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;
        h(drug) = plot(10.^x,y,'color', color{drug});
        hold on
        errorbar(10.^xdata, ydata, yste, 'color', color{drug}, 'marker', 'o', 'linestyle', 'none')
    end
    legend([h(1), h(2), h(3)], condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
end

%% figure
% example cells:
% Superior: 618 Anterior: 1066 Inferior: 347 Posterior: 6183
cell_id = [618 7338 347 6574];
ctr_i = [4 3 4 4];
direction = [7 8 3 5];
% bin_size = 0.1; % second
% XX = bin_size/2:bin_size:3-bin_size/2; 

h = figure;
set(h, 'position', pos)
for cc = 1:length(cell_id)
    cell_idx = find(ds_id == cell_id(cc));
    for drug = 1:4
        h = subplot(8,4,(drug-1)*4+cc);
        plot_raster(raster_mb{drug}{cell_idx}(:,:,direction(cc), ctr_i(cc),:), 0,3, 'color', color{drug})
        set(h, 'xtick', [])
        set(h, 'ytick', [])
        if drug == 1
            title(dscell_type{cc})
        end
    end
    subplot(8,4,4*[4:7]+cc)
    for drug = 1:3
        ydata = mean(response_pn_norm{drug}{cc});
        yste = std(response_pn_norm{drug}{cc})/sqrt(size(response_pn_norm{drug}{cc}, 1));
        xdata = log10(ctr_x);
        [f, G] = fit_nr(xdata, ydata, 'upper', [2, 100, log10(300), min(ydata)]);

        x = linspace(min(xdata), max(xdata)+0.5, 100);
        y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;
        h(drug) = plot(10.^x,y,'color', color{drug});
        hold on
        errorbar(10.^xdata, ydata, yste, 'color', color{drug}, 'marker', 'o', 'linestyle', 'none')
    end
    legend([h(1), h(2), h(3)], condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('contrast')
    ylabel('spike rate')
    xlim([3 1000])
    ylim([-0.1 1.2])
%     title(dscell_type{dir})
end
    
    

%%

for dir = 1:4
    for drug = 1:4
        for cc = 1:size(sigma{dir},2)
            c20 = sigma{dir}(1,cc)/(1)^(1/aa{dir}(drug,cc));
            c20r{dir}(drug, cc) = ymax{dir}(drug,cc)*c20^aa{dir}(drug,cc)/(c20^aa{dir}(drug,cc)+sigma{dir}(drug,cc)^aa{dir}(drug,cc));
        end
    end
    c20r_mean{dir} = mean(c20r{dir}');
    c20r_ste{dir} = std(c20r{dir}')/sqrt(size(sigma{dir},2));
end

cell_type = {'superior', 'anterior', 'inferior', 'posterior'};
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = cell_type;
model_series = cell2mat(c20r_mean');
model_error = cell2mat(c20r_ste');
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('c20r')
title('NDF 3')
legend('control','AP5', 'AP5+HEX', 'wash', 'location', 'northeast');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end
  
%% c20 combined DSGC

c20r_combined = cell2mat(c20r);
c20r_combined_ndf3_mean = mean(c20r_combined, 2);
c20r_combined_ndf3_ste = std(c20r_combined, [], 2)/sqrt(size(c20r_combined, 2));

errorbar(c20r_combined_ndf3_mean, c20r_combined_ndf3_ste)

%%

%% response variance
% spike count variance
for i = 1:4
    for dir = 1:4
        spike_count_var_p{i}{dir} = [];
        spike_count_mean_p{i}{dir} = [];
        for cc = 1:length(id_dir{dir})
            if ~isempty(raster_p_sum_all{i}{idx_dir{dir}(cc)})
                for ctr = 1:7
                    spike_count = cellfun(@length, squeeze(raster_p_sum_all{i}{idx_dir{dir}(cc)}(:, :, ctr, :)));
                    if sum(spike_count == 0) <= 5
                        spike_count_var_p{i}{dir} = [spike_count_var_p{i}{dir} var(spike_count)];
                        spike_count_mean_p{i}{dir} = [spike_count_mean_p{i}{dir} mean(spike_count)];
                    end
                end
            end
        end
    end
end


figure
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:3
        plot(spike_count_mean_p{drug}{dir}, spike_count_var_p{drug}{dir}, [color(drug) 'o'])
        hold on
    end
    xlabel('spike count mean')
    ylabel('spike count variance')
    legend('control', 'AP5', 'AP5+HEX')
end

% 1st spike timing variance
for i = 1:4
    for dir = 1:4
        spike_count_mean_p{i}{dir} = [];
        spike_timing_var_p{i}{dir} = [];
        for cc = 1:length(id_dir{dir})
            if ~isempty(raster_p_sum_all{i}{idx_dir{dir}(cc)})
                for ctr = 1:7
                    spike_count = cellfun(@length, squeeze(raster_p_sum_all{i}{idx_dir{dir}(cc)}(:, :, ctr, :)));
                    if sum(spike_count == 0) < 5
                        spike_count_mean_p{i}{dir} = [spike_count_mean_p{i}{dir} mean(spike_count)];
                        spike_timing = [];
                        for trial = 1:10
                            if ~isempty(raster_p_sum_all{i}{idx_dir{dir}(cc)}{:, :, ctr, trial})
                                spike_timing = [spike_timing raster_p_sum_all{i}{idx_dir{dir}(cc)}{:, :, ctr, trial}(1)];
                            end
                        end
                        spike_timing_var_p{i}{dir} = [spike_timing_var_p{i}{dir} var(spike_timing)];
                    end
                end
            end
        end
    end
end


figure
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:3
        plot(spike_count_mean_p{drug}{dir}, spike_timing_var_p{drug}{dir}, [color(drug) 'o'])
        hold on
    end
    xlabel('spike count mean')
    ylabel('first spike timing variance')

%     set(gca, 'YScale', 'log')
    legend('control', 'AP5', 'AP5+HEX')
end

%% direction tuning curves
% all ds cells

color = 'brgkcmy';
dirn = 4;
D = 1;
T = 1;
BW = 1;
CL = 7;

p_direction = MB{D}.angle{T,BW,CL}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;


%subtypes
clear rho_mb_mean rho_mb_ste dsi_mb_mean dsi_mb_ste rho_mb dsi_mb dsi_mb_mean_all dsi_mb_ste_all dsi_mb_mean_all_on dsi_mb_ste_all_on

for drug = 1:4
    for cl = 1:7
%         subplot(3, 7, (drug-1)*7+cl)
        for i = 1:4
            rho_mb{drug}{i}{cl} = [];
            RHO_mb{drug}{i}{cl} = [];
            dsi_mb{drug}{i}{cl} = [];
            for cc = 1:length(idx_dir{i})
                if ~mb_idx(idx_dir{i}(cc)) && sum(MB{drug}.RHO{T, BW,cl}(idx_dir{i}(cc), :))>0
                [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
                Y_temp = MB{drug}.RHO{T,BW,cl}(idx_dir{i}(cc), :);
                y_temp = MB{drug}.rho{T,BW,cl}(idx_dir{i}(cc), :);
%                 y_temp = MB{drug}.RHO{T,BW,cl}(idx_dir{i}(cc), :)/bin_size;
%                 plot(xsort, y_temp(seq), color(i))
%                     ylim([0 1])
    %             pause
%                 hold on
                rho_mb{drug}{i}{cl} = [rho_mb{drug}{i}{cl}; y_temp(seq)];
                RHO_mb{drug}{i}{cl} = [RHO_mb{drug}{i}{cl}; Y_temp(seq)];
                dsi_mb{drug}{i}{cl} = [dsi_mb{drug}{i}{cl}; MB{drug}.dsindex{T,BW,cl}(idx_dir{i}(cc))];
                end
            end
            if ~isempty(rho_mb{drug}{i}{cl})
                rho_mb_mean{cl}{drug}(i, :) = mean(rho_mb{drug}{i}{cl}, 1);
                rho_mb_ste{cl}{drug}(i, :) = std(rho_mb{drug}{i}{cl}, [], 1)/sqrt(size(rho_mb{drug}{i}{cl}, 1));
                RHO_mb_mean{cl}{drug}(i, :) = mean(RHO_mb{drug}{i}{cl}, 1);
                RHO_mb_ste{cl}{drug}(i, :) = std(RHO_mb{drug}{i}{cl}, [], 1)/sqrt(size(RHO_mb{drug}{i}{cl}, 1));
                dsi_mb_mean{cl}{drug}(i) = mean(dsi_mb{drug}{i}{cl});
                dsi_mb_ste{cl}{drug}(i) = std(dsi_mb{drug}{i}{cl})/sqrt(length(dsi_mb{drug}{i}{cl}));
            else
                rho_mb_mean{cl}{drug}(i, :) = nan;
                rho_mb_ste{cl}{drug}(i, :) = nan;
                RHO_mb_mean{cl}{drug}(i, :) = nan;
                RHO_mb_ste{cl}{drug}(i, :) = nan;
                dsi_mb_mean{cl}{drug}(i) = nan;
                dsi_mb_ste{cl}{drug}(i) = nan;
            end
        end
%         xlabel('direction (rad)')
%         ylabel('spike number')
%             title(ll{d})
%         xlim([-pi pi])
        dsi_mb_mean_all{cl} = cell2mat(dsi_mb_mean{cl}');
        dsi_mb_ste_all{cl} = cell2mat(dsi_mb_ste{cl}'); 

        for i = 1:2
            rho_mb_on{drug}{i}{cl} = [];
            RHO_mb_on{drug}{i}{cl} = [];
            dsi_mb_on{drug}{i}{cl} = [];
            for cc = 1:length(idx_dir_on{i})
                if ~mb_idx(idx_dir_on{i}(cc))% && sum(MB_NDF{drug, d}.RHO{T, BW,cl}(idx_dir{i}(cc), :))>0
                [xsort, seq] = sort(xx(idx_dir_on{i}(cc), :));
                Y_temp = MB{drug}.RHO{T,BW,cl}(idx_dir_on{i}(cc), :);
                y_temp = MB{drug}.rho{T,BW,cl}(idx_dir{i}(cc), :);
%                 y_temp = MB{drug}.RHO{T,BW,cl}(idx_dir{i}(cc), :)/bin_size;
%                 plot(xsort, y_temp(seq), color(i))
%                     ylim([0 1])
    %             pause
%                 hold on
                rho_mb_on{drug}{i}{cl} = [rho_mb_on{drug}{i}{cl}; y_temp(seq)];
                RHO_mb_on{drug}{i}{cl} = [RHO_mb_on{drug}{i}{cl}; Y_temp(seq)];
                dsi_mb_on{drug}{i}{cl} = [dsi_mb_on{drug}{i}{cl}; MB{drug}.dsindex{T,BW,cl}(idx_dir_on{i}(cc))];
                end
            end
            if ~isempty(rho_mb_on{drug}{i}{cl})
                rho_mb_mean_on{cl}{drug}(i, :) = mean(rho_mb_on{drug}{i}{cl},1);
                rho_mb_ste_on{cl}{drug}(i, :) = std(rho_mb_on{drug}{i}{cl},[],1)/sqrt(size(rho_mb_on{drug}{i}{cl}, 1));
                RHO_mb_mean_on{cl}{drug}(i, :) = mean(RHO_mb_on{drug}{i}{cl},1);
                RHO_mb_ste_on{cl}{drug}(i, :) = std(RHO_mb_on{drug}{i}{cl},[],1)/sqrt(size(RHO_mb_on{drug}{i}{cl}, 1));
                dsi_mb_mean_on{cl}{drug}(i) = mean(dsi_mb_on{drug}{i}{cl});
                dsi_mb_ste_on{cl}{drug}(i) = std(dsi_mb_on{drug}{i}{cl})/sqrt(length(dsi_mb_on{drug}{i}{cl}));
            end
        end
%         xlabel('direction (rad)')
%         ylabel('spike number')
%             title(ll{d})
%         xlim([-pi pi])
        dsi_mb_mean_all_on{cl} = cell2mat(dsi_mb_mean_on{cl}');
        dsi_mb_ste_all_on{cl} = cell2mat(dsi_mb_ste_on{cl}');            


    end
end
dsi_mb_mean_all = reshape(cell2mat(dsi_mb_mean_all), size(dsi_mb_mean_all{1}, 1), size(dsi_mb_mean_all{1}, 2), length(dsi_mb_mean_all));
dsi_mb_ste_all = reshape(cell2mat(dsi_mb_ste_all), size(dsi_mb_ste_all{1}, 1), size(dsi_mb_ste_all{1}, 2), length(dsi_mb_ste_all));
dsi_mb_mean_all_on = reshape(cell2mat(dsi_mb_mean_all_on), size(dsi_mb_mean_all_on{1}, 1), size(dsi_mb_mean_all_on{1}, 2), length(dsi_mb_mean_all_on));
dsi_mb_ste_all_on = reshape(cell2mat(dsi_mb_ste_all_on), size(dsi_mb_ste_all_on{1}, 1), size(dsi_mb_ste_all_on{1}, 2), length(dsi_mb_ste_all_on));

% plot average (cell type)
ct = {'superior', 'anterior', 'inferior', 'posterior'};
h = figure;
set(h, 'Position', [1 1 1520,1080])
for drug = 1:4
    for cl = 1:7
        subplot(4, 7, (drug-1)*7+cl)
        for i = 1:dirn
            if i<=size(rho_mb_mean{cl}{drug},1)
                errorbar(xsort/pi*180, rho_mb_mean{cl}{drug}(i, :), rho_mb_ste{cl}{drug}(i, :), color(i));
                hold on
            end
        end
        xlabel('degrees')
        ylabel('spike number')
%                 title(ll{d});

    end
end
legend(ct)

h = figure;
set(h, 'Position', [1 1 1520,1080])
for drug = 1:4
    for cl = 1:7
        subplot(4, 7, (drug-1)*7+cl)
        for i = 1:2
            if i<=size(rho_mb_mean_on{cl}{drug},1)
                errorbar(xsort/pi*180, rho_mb_mean_on{cl}{drug}(i, :), rho_mb_ste_on{cl}{drug}(i, :), color(i));
                hold on
            end
        end
        xlabel('degrees')
        ylabel('spike number')
%                 title(ll{d});

    end
end
legend(ct)

% plot average (cell type)
% Drug = {'control', 'AP5', 'AP5+HEX', 'wash'};
Drug = {'control', 'AP5', 'wash'};

h = figure;
set(h, 'Position', [1 1 1520,1080])
for i = 1:4
    for cl = 1:7
        subplot(4, 7, (i-1)*7+cl)
        for drug = [1 2 4]
            if i<=size(rho_mb_mean{cl}{drug},1)
                errorbar(xsort/pi*180, rho_mb_mean{cl}{drug}(i, :), rho_mb_ste{cl}{drug}(i, :), color(drug));
                hold on
            end
        end
        xlabel('degrees')

        if cl == 1
            ylabel(ct{i})
        end
        if i == 1
            title(ctr{cl})
        end
    end
    if i == 1
        legend(Drug)
    end
end

h = figure;
set(h, 'Position', [1 1 1520,1080])
for i = 1:2
    for cl = 1:7
        subplot(4, 7, (i-1)*7+cl)
        for drug = 1:3
            if i<=size(rho_mb_mean_on{cl}{drug},1)
                errorbar(xsort/pi*180, rho_mb_mean_on{cl}{drug}(i, :), rho_mb_ste_on{cl}{drug}(i, :), color(drug));
                hold on
            end
        end
        xlabel('degrees')

        if cl == 1
            title(ct{i})
        end
    end
    if i == 1
        legend(Drug)
    end
end

% control
ctr = {'5%', '10%', '20%', '40%', '80%', '150%', '300%'};
ct = {'superior', 'anterior', 'inferior', 'posterior'};
h = figure;
set(h, 'Position', [1 1 1520,1080])
drug = 1;
for i = 1:4
    subplot(2, 2, i)
    for cl = 1:7
        if i<=size(rho_mb_mean{cl}{drug},1)
            errorbar(xsort/pi*180, rho_mb_mean{cl}{drug}(i, :)/max(rho_mb_mean{cl}{drug}(i, :)), rho_mb_ste{cl}{drug}(i, :), color(cl));
            hold on
        end
    end
    xlabel('degrees')
    ylabel('spike number')
    title(ct{i})
%                 title(ll{d});

end
legend(ctr)
% subplot(3,2,5)

h = figure;
set(h, 'Position', [1 1 1520,1080])
drug = 1;
for i = 1:2
    subplot(2, 2, i)
    for cl = 1:7
        if i<=size(rho_mb_mean_on{cl}{drug},1)
            errorbar(xsort/pi*180, rho_mb_mean_on{cl}{drug}(i, :), rho_mb_ste_on{cl}{drug}(i, :), color(cl));
            hold on
        end
    end
    xlabel('degrees')
    ylabel('spike number')
    title(ct{i})
%                 title(ll{d});

end
legend(ctr)
subplot(2,2,3)

% individual cell
figure
drug = 1;
for i = 7:7%size(RHO_mb{drug}{2}{1}, 1)
%     subplot(5,6,i)
    for cl = 1:7
        plot(xsort/pi*180, RHO_mb{drug}{2}{cl}(i,:))
        hold on
    end
end
xlim([-200 200])
xlabel('degrees')
ylabel('spike number')

% DSI against contrast
drug = 1;
ctr_x = [300 150 80 40 20 10 5];
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
figure;
for i = 1:4
    errorbar(ctr_x, dsi_mb_mean_all(drug, i, :), dsi_mb_ste_all(drug, i, :))
    hold on
end
set(gca, 'Xscale', 'log')
legend(dscell_type, 'location', 'NorthEastOutside')
xlabel('contrast %')
ylabel('DSI')

%% separate ON and OFF responses for CRF

ctr = 7;
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur{1}-bin_size/2;
drug = 1;
for cc = 1:length(ds_id)
    if ~mb_idx(cc)
        h = figure(1);
        set(h, 'position', [1 1 1000 1000])
        for i = 7:-1:1
            subplot(7,1,i)
            temp = hist(raster_p_sum{drug}{cc}{ctr-7+i}, xx);
            plot(xx, temp)
            hold on
        end
        [x,~] = ginput;
        BreakIndices(cc, :) = round(x'/bin_size);
        close(1)
    else
        BreakIndices(cc, :) = nan(1, ctr);


    end
end

%% max window
% on-off DSGC
clear Max_i
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'AP5+Hex', 'wash'};
step_size = 40;
for drug = 1:4
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                for ctr = 7:-1:1
                    a = raster_p_sum{drug}{idx_dir{dir}(cc)}{ctr};
                    hist_temp = hist(a, xx);
                    OnRes = hist_temp(BreakIndices(idx_dir{dir}(cc), ctr)-40:BreakIndices(idx_dir{dir}(cc), ctr));
                    response_pmax{drug}{dir}(CC, ctr) = sum(OnRes)/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
                end
                CC = CC + 1;
            end
        end
        response_pmax_norm{drug}{dir} = response_pmax{drug}{dir}./repmat(max(response_pmax{1}{dir}, [], 2), 1, size(response_pmax{drug}{dir},2));
%         response_pmax_norm{drug}{dir} = response_pmax{drug}{dir}./repmat(max(response_pmax{drug}{dir}, [], 2), 1, size(response_pmax{drug}{dir},2));
    end
end


ctr_x = [5 10 20 40 80 150 300];
color = 'brgkc';
figure
set(gcf, 'Position', [1 1 900 800])

for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4
        errorbar(ctr_x, mean(response_pmax_norm{drug}{dir}), std(response_pmax_norm{drug}{dir})/sqrt(size(response_pmax_norm{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
%     xlim([3 400])
end

%% fit ds
for dir = 1:4
    for drug = 1:4
%         figure
%         set(gcf, 'Position', [1 1 900 800])
        for cc = 1:size(response_pmax_norm{drug}{dir}, 1)

            ydata = response_pmax_norm{drug}{dir}(cc, :);
            xdata = log10(ctr_x);
            [f, G] = fit_nr(xdata, ydata, 'upper', [1.1, 100, log10(300), 0]);
            fit_all{drug}{dir}{cc} = f;
            G_all{drug}{dir}{cc} = G;

%             x = linspace(min(xdata), max(xdata), 100);
%             y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;
% 
%             subplot(5,6,cc)
%             plot(xdata, ydata)
%             hold on
%             plot(x, y)

            sigma{dir}(drug, cc) = f.sigma;
            ymax{dir}(drug, cc) = f.ymax;
            aa{dir}(drug, cc) = f.a;
            bb{dir}(drug, cc) = f.b;
        end
    end
end



%%
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'AP5+Hex', 'wash'};
step_size = 60;
for drug = 1:4
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_n_sum{drug}{idx_dir{dir}(cc)})
                for ctr = 7:-1:1
                    a = raster_n_sum{drug}{idx_dir{dir}(cc)}{ctr};
                    hist_temp = hist(a, xx);
                    if drug == 1 && ctr == 7
                        [max_p, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
                        Max_i{dir}(CC) = max_i;
                    else
                        max_p = conv(hist_temp, ones(1,step_size), 'valid');
                        max_p = max_p(Max_i{dir}(CC));
                    end
%                     response_pn{drug}{dir}(CC, ctr) = max(max_p - mean_n, 0);
%                     response_pn{drug}{dir}(CC, ctr) = max_p - mean_n;
%                     response_pn{drug}{dir}(CC, ctr) = abs(max_p - mean_n);
                    response_pn{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
                end
                CC = CC + 1;
            end
        end
        response_pn_norm{drug}{dir} = response_pn{drug}{dir}./repmat(max(response_pn{1}{dir}, [], 2), 1, size(response_pn{drug}{dir},2));
    end
end


ctr_x = [5 10 20 40 80 150 300];
color = 'brgkc';
figure
set(gcf, 'Position', [1 1 900 800])

for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4
        errorbar(ctr_x, mean(response_pn_norm{drug}{dir}), std(response_pn_norm{drug}{dir})/sqrt(size(response_pn_norm{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
%     xlim([3 400])
end


