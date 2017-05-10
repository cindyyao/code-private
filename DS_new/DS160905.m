cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);


datadg = load_data('/Volumes/lab/analysis/2016-09-05-0/data004-sorted/data004-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-09-05-0/stimuli/s04.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

datarun = load_data('/Volumes/lab/analysis/2016-09-05-0/data000-003-map/data000-003-map', opt);
time_points = [1900 3800 5700];
datamb(1:4) = split_datarun(datarun, time_points);
datamb{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-09-05-0/stimuli/s00.txt';
datamb{1} = load_stim(datamb{1}, 'user_defined_trigger_set', [1:2:1120]);
datamb{1}.stimulus.triggers = datamb{1}.stimulus.triggers';
datamb{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-09-05-0/stimuli/s01.txt';
datamb{2} = load_stim(datamb{2}, 'user_defined_trigger_set', [1:2:1120]);
datamb{2}.stimulus.triggers = datamb{2}.stimulus.triggers';
datamb{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-09-05-0/stimuli/s02.txt';
datamb{3} = load_stim(datamb{3}, 'user_defined_trigger_set', [1:2:1120]);
datamb{3}.stimulus.triggers = datamb{3}.stimulus.triggers';
datamb{4}.names.stimulus_path = '/Volumes/lab/analysis/2016-09-05-0/stimuli/s03.txt';
datamb{4} = load_stim(datamb{4}, 'user_defined_trigger_set', [1:2:1120]);
datamb{4}.stimulus.triggers = datamb{4}.stimulus.triggers';

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [2 5]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

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
load('DS160905.mat')
n = 4;
duration = [1776.2 1775.1 1775.2 1775.2]; %sec
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
    trial_dur{i} = get_mb_trial_dur(datamb{i}, 400,400,0.5);
end

ctr_p = 7; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_mb = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}, raster_p_sum_all{i}] = get_pdirection_raster(raster_mb{i}, MB{1}.angle{ctr_p});
    MAG_all_norm_mb{i} = normalize_MAG(MB{i});
    rep = datamb{i}.stimulus.repetitions;
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

%% plot cell summary
for cc = 1:length(ds_id)
    plot_mb_raster_ctr(MB, raster_mb, trial_dur, cc, ds_id(cc), 'NDF3', 4, 7, 1)
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
condition = {'control', 'Hex', 'AP5+Hex', 'wash'};
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
%     for dir = 1:2
%         CC = 1;
%         for cc = 1:length(idx_dir_on{dir})
%             if ~isempty(raster_p_sum{drug}{idx_dir_on{dir}(cc)})
%                 spike_temp = cellfun('length', raster_p_sum{drug}{idx_dir_on{dir}(cc)})/datamb{drug}.stimulus.repetitions;
%                 if sum(spike_temp) ~= 0
%                     pd_dir_on_spikes{drug}{dir}(CC,:) = cellfun('length', raster_p_sum{drug}{idx_dir_on{dir}(cc)})/datamb{drug}.stimulus.repetitions;
% %                     pd_dir_on_spikes_nor{drug}{dir}(CC,:) = pd_dir_on_spikes{drug}{dir}(CC,:)/max(pd_dir_on_spikes{1}{dir}(CC,:));
%                     pd_dir_on_spikes_nor{drug}{dir}(CC,:) = pd_dir_on_spikes{drug}{dir}(CC,:)/max(pd_dir_on_spikes{drug}{dir}(CC,:));
%                     CC = CC + 1;
%                 end
%             end
%         end
%         pd_dir_on_spikes_nor{drug}{dir} = nan2empty(pd_dir_on_spikes_nor{drug}{dir});
%     end
% 
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
    legend(condition, 'location', 'northwest')
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
condition = {'control', 'Hex', 'AP5+Hex', 'wash'};
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
%     for dir = 1:2
%         CC = 1;
%         for cc = 1:length(idx_dir_on{dir})
%             if ~isempty(raster_p_sum{drug}{idx_dir_on{dir}(cc)})
%                 a = raster_p_sum{drug}{idx_dir_on{dir}(cc)};
%                 hist_temp = cellfun(@(a) hist(a, xx), a, 'UniformOutput', false);
%                 hist_temp = cell2mat(squeeze(cellfun(@(hist_temp) max(hist_temp), hist_temp, 'UniformOutput', false)))/datamb{drug}.stimulus.repetitions/bin_size;
%                 if sum(hist_temp) ~= 0
%                     pd_dir_on_spikes{drug}{dir}(CC,:) = hist_temp;
%                     pd_dir_on_spikes_nor{drug}{dir}(CC,:) = pd_dir_on_spikes{drug}{dir}(CC,:)/max(pd_dir_on_spikes{1}{dir}(CC,:));
%                     CC = CC + 1;
%                 end
%             end
%         end
%         pd_dir_on_spikes_nor{drug}{dir} = nan2empty(pd_dir_on_spikes_nor{drug}{dir});
%     end

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
bin_size = 0.025;
xx = bin_size/2:bin_size:3.2-bin_size/2;
ctr_x = [5 10 20 40 80 150 300];
bar_t = datamb{1}.stimulus.params.BAR_WIDTH/datamb{1}.stimulus.params.DELTA/60;
bar_bin = ceil(bar_t/bin_size);
color = 'brgkc';
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'Hex', 'AP5+Hex', 'wash'};
for drug = 1:4    
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                a = raster_p_sum{drug}{idx_dir{dir}(cc)};
                hist_temp = cellfun(@(a) hist(a, xx), a, 'UniformOutput', false);
                if drug == 1
                    plot(hist_temp{7})
                    [x,~] = ginput;
                    x = round(x);
                    T{dir}{CC} = x;
                end
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
% 
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
    ylabel('Normalized Response')
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

            sigma{drug}{dir}(cc) = f.sigma;
            ymax{drug}{dir}(cc) = f.ymax;
            aa{drug}{dir}(cc) = f.a;
            bb{drug}{dir}(cc) = f.b;
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
    xlabel('normalized contrast')
    ylabel('normalized response')
    title(dscell_type{dir})
%     xlim([3 400])
end

%%

drug = 1;
for dir = 1:4
    figure
    CC = 1;
    for cc = 1:length(idx_dir{dir})
        if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
            a = raster_p_sum{drug}{idx_dir{dir}(cc)}{7};
            hist_temp = hist(a, xx);
            subplot(3,4,CC)
            plot(xx, hist_temp)
            hold on
            line([T{dir}{CC}*bin_size, T{dir}{CC}*bin_size], [0 max(hist_temp)], 'color', 'r')
            CC = CC + 1;
        end
    end
end

%% get spontaneous activity
for drug = 1:4
    for dir = 1:4
        for cc = 1:length(id_dir_mb{dir})
            idx = get_cell_indices(datamb{drug},id_dir_mb{dir}(cc));
            spikes_temp = datamb{drug}.spikes{idx};
            bgnd_firing{dir}(drug, cc) = length(spikes_temp(spikes_temp > 1780 & spike_temp < 1800))/20;
        end
    end
end

%% prefer - base
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'HEX', 'AP5+Hex', 'wash'};
step_size = 70;
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
%         figure
%         set(gcf, 'Position', [1 1 900 800])
        for cc = 1:size(response_pn_norm{drug}{dir}, 1)

            ydata = response_pn_norm{drug}{dir}(cc, :);
            xdata = log10(ctr_x);
            [f, G] = fit_nr(xdata, ydata, 'upper', [2, 100, log10(300), min(ydata)]);
            fit_all{drug}{dir}{cc} = f;
            G_all{drug}{dir}{cc} = G;

%             x = linspace(min(xdata), max(xdata), 100);
%             y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;
% 
%             subplot(4,4,cc)
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

% 
% figure
% set(gcf, 'Position', [1 1 900 800])
% for dir = 1:4
%     subplot(2,2,dir)
%     for drug = 1:4
%         for cc = 1:size(response_pn_norm{drug}{dir}, 1)
% %             plot(log10(ctr_x), response_pn_norm{drug}{dir}(cc, :), color(drug))
%             plot(log10(ctr_x)-(sigma{1}{dir}(cc)-1), response_pn_norm{drug}{dir}(cc, :), [color(drug) 'o'])
%             hold on
%         end
%     end
% end


%
figure
set(gcf, 'Position', [1 1 900 800])
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4
        ydata = mean(response_pn_norm{drug}{dir});
        yste = std(response_pn_norm{drug}{dir})/sqrt(size(response_pn_norm{drug}{dir}, 1));
        xdata = log10(ctr_x);
        [f, G] = fit_nr(xdata, ydata, 'upper', [2, 100, log10(300), min(ydata)]);

        x = linspace(min(xdata), max(xdata)+0.5, 100);
        y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;
        h(drug) = plot(10.^x,y,color(drug));
        hold on
        errorbar(10.^xdata, ydata, yste, [color(drug) 'o'])
    end
    legend([h(1), h(2), h(3), h(4)], condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
end


