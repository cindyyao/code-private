cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

datadg = load_data('/Volumes/lab/Experiments/Array/Analysis/2017-06-29-0/data010-sorted/data010-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2017-06-29-0/stimuli/s10.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

params_idx = [1 2]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

datamb{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2017-06-29-0/data001-map/data001-map', opt);
datamb{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2017-06-29-0/stimuli/s01.txt';
datamb{1} = load_stim(datamb{1}, 'user_defined_trigger_set', [1:2:1120]); %1120: number of triggers
datamb{1}.stimulus.triggers = datamb{1}.stimulus.triggers';
datamb{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2017-06-29-0/data003-map/data003-map', opt);
datamb{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2017-06-29-0/stimuli/s03.txt';
datamb{2} = load_stim(datamb{2}, 'user_defined_trigger_set', [1:2:1120]);
datamb{2}.stimulus.triggers = datamb{2}.stimulus.triggers';
datamb{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2017-06-29-0/data005-map/data005-map', opt);
datamb{3}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2017-06-29-0/stimuli/s05.txt';
datamb{3} = load_stim(datamb{3}, 'user_defined_trigger_set', [1:2:1120]);
datamb{3}.stimulus.triggers = datamb{3}.stimulus.triggers';
datamb{4} = load_data('/Volumes/lab/Experiments/Array/Analysis/2017-06-29-0/data007-map/data007-map', opt);
datamb{4}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2017-06-29-0/stimuli/s07.txt';
datamb{4} = load_stim(datamb{4}, 'user_defined_trigger_set', [1:2:1120]);
datamb{4}.stimulus.triggers = datamb{4}.stimulus.triggers';
load('DS170629.mat')
%% dg
% ds_id = datadg.cell_ids;
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);

delta_p = 1; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(DG{i});
rep = datadg.stimulus.repetitions;

%% mb
n = 4;
for i = 1:n
    duration(i) = datamb{i}.triggers(end);
end
bin_size = 0.025; %sec
[raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx, raster_mb_all] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},ds_id,duration(i),bin_size);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb{i} = get_mb_raster(datamb{i}, ds_id, duration(i));
    raster_mb_all{i} = combine_repeats(raster_mb{i});
    for j = 1:length(raster_mb{i})
        if(mb_idx(j))
            raster_mb{i}{j} = [];
            raster_mb_all{i}{j} = [];
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

t = 1;
figure
compass(DG{1}.U{t}(idx_sub{1}), DG{1}.V{t}(idx_sub{1}), 'r')
hold on
compass(DG{1}.U{t}(idx_sub{2}), DG{1}.V{t}(idx_sub{2}), 'b')

%%
d = 1;
t = 1;
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
t = 1;
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
for dir = 1:3
    id_dir_on_mb{dir} = id_dir_on{dir}(~mb_idx(idx_dir_on{dir}));
    idx_dir_on_mb{dir} = idx_dir_on{dir}(~mb_idx(idx_dir_on{dir}));
end

%% plot cell summary
for dir = 1:4
    for cc = 2:length(id_dir{dir})
        plot_mb_raster_ctr(MB, raster_mb, trial_dur, idx_dir{dir}(cc), id_dir{dir}(cc), '', 4, 7, 1)
    end
end

% for dir = 1:3
%     for cc = 1:length(id_dir_on{dir})
%         plot_mb_raster_ctr(MB, raster_mb, trial_dur, idx_dir_on{dir}(cc), id_dir_on{dir}(cc), '', 4, 7, 1)
%     end
% end
%% get spontaneous activity
% trial_dur = mean(diff(datamb{1}.stimulus.triggers));
% bin_size = 0.01;
% xx = bin_size/2:bin_size:trial_dur-bin_size/2;
% step_size = 60;
% for drug = 1:4
%     for dir = 1:4
%         for cc = 1:length(id_dir{dir})
%             if(mb_idx(idx_dir{dir}(cc)))
%                 bgnd_firing{dir}(drug, cc) = nan;
%             else
%                 idx = get_cell_indices(datamb{drug},id_dir{dir}(cc));
%                 spikes_temp = datamb{drug}.spikes{idx};
%                 spikes_temp = spikes_temp(spikes_temp > 2050 & spikes_temp < 2100);
%                 hist_spikes_temp = hist(spikes_temp, xx);
%                 [max_s, max_i] = max(conv(hist_spikes_temp, ones(1,step_size), 'valid'));
%                 bgnd_firing{dir}(drug, cc) = max_s/(bin_size * step_size);
%             end
%         end
%     end
% end
% bgnd_firing_combine{1} = exciseRows_empty(bgnd_firing{1}');
% bgnd_firing_combine{2} = exciseRows_empty(cell2mat(bgnd_firing(2:4))');
% 
%% get spontaneous activity

for drug = 1:4
    for dir = 1:4
        for cc = 1:length(id_dir_mb{dir})
            if(~mb_idx(idx_dir_mb{dir}(cc)))
                idx = get_cell_indices(datamb{drug},id_dir_mb{dir}(cc));
                spikes_temp = datamb{drug}.spikes{idx};
                bgnd_firing{dir}(drug, cc) = length(spikes_temp(spikes_temp > 2050 & spikes_temp < 2100))/50;
%                 bgnd_firing{dir}(drug, cc) = length(spikes_temp(spikes_temp < 70))/70;
            end
        end
    end
end
bgnd_firing_combine{1} = exciseRows_empty(bgnd_firing{1}');
bgnd_firing_combine{2} = exciseRows_empty(cell2mat(bgnd_firing(2:4))');


figure
errorbar(mean(bgnd_firing_combine{2}(:, 1)),  mean(bgnd_firing_combine{1}(:, 1)), std(bgnd_firing_combine{1}(:, 1))/sqrt(size(bgnd_firing_combine{1}, 1)), 'bo')
hold on
errorbar(mean(bgnd_firing_combine{2}(:, 2)),  mean(bgnd_firing_combine{1}(:, 2)), std(bgnd_firing_combine{1}(:, 2))/sqrt(size(bgnd_firing_combine{1}, 1)), 'ro')

herrorbar(mean(bgnd_firing_combine{2}(:, 1)),  mean(bgnd_firing_combine{1}(:, 1)), std(bgnd_firing_combine{2}(:, 1))/sqrt(size(bgnd_firing_combine{2}, 1)), 'bo')
herrorbar(mean(bgnd_firing_combine{2}(:, 2)),  mean(bgnd_firing_combine{1}(:, 2)), std(bgnd_firing_combine{2}(:, 2))/sqrt(size(bgnd_firing_combine{2}, 1)), 'ro')
plot([0 5], [0 5], 'k--')
xlim([0 5])
ylim([0 5])
xlabel('others')
ylabel('superior')
legend('control', 'SR', 'location', 'southeast')
%%
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'SR', 'SR+TPMPA', 'wash'};
step_size = 80;
for drug = 1:4
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_n_sum{drug}{idx_dir{dir}(cc)})
                for ctr = 7:-1:1
                    a = raster_n_sum{drug}{idx_dir{dir}(cc)}{ctr};
                    hist_temp = hist(a, xx);
%                     if drug == 1 && ctr == 7
                        [max_p, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
                        Max_i{dir}(CC) = max_i;
%                     else
%                         max_p = conv(hist_temp, ones(1,step_size), 'valid');
%                         max_p = max_p(Max_i{dir}(CC));
%                     end
%                     response_pn{drug}{dir}(CC, ctr) = max(max_p - mean_n, 0);
%                     response_pn{drug}{dir}(CC, ctr) = max_p - mean_n;
%                     response_pn{drug}{dir}(CC, ctr) = abs(max_p - mean_n);
                    response_pn{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
%                     response_pn{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions;
                end
                CC = CC + 1;
            end
        end
        response_pn_norm{drug}{dir} = response_pn{drug}{dir}./repmat(max(response_pn{drug}{dir}, [], 2), 1, size(response_pn{drug}{dir},2));
    end
end

ctr_x = [5 10 20 40 80 150 300];
color = 'brgkc';
figure
set(gcf, 'Position', [1 1 900 800])

for drug = 1:4
    subplot(2,2,drug)
    for dir = 1:4
        errorbar(ctr_x, nanmean(response_pn_norm{drug}{dir}), nanstd(response_pn_norm{drug}{dir})/sqrt(size(response_pn_norm{drug}{dir}, 1)), 'color', color(dir));
        hold on
    end
    legend(dscell_type, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(condition{drug})
    xlim([3 400])
end



figure
for drug = 1:4
    subplot(2,2,drug)
    for dir = 1:4
        response_s{drug} = exciseRows_empty(response_pn_norm{drug}{1});
        errorbar(ctr_x, nanmean(response_s{drug}), nanstd(response_s{drug})/sqrt(size(response_s{drug}, 1)), 'color', color(1));
        hold on
        response_others{drug} = exciseRows_empty(cell2mat(response_pn_norm{drug}(2:4)'));
        errorbar(ctr_x, nanmean(response_others{drug}), nanstd(response_others{drug})/sqrt(size(response_others{drug}, 1)), 'color', color(2));
    end
    legend('superior', 'others', 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(condition{drug})
    xlim([3 400])
end


ctr_x = [5 10 20 40 80 150 300];
color = 'brgkc';
figure
set(gcf, 'Position', [1 1 900 800])

for drug = 1:4
    subplot(2,2,drug)
    for dir = 1:4
        errorbar(ctr_x, nanmean(response_pn{drug}{dir}), nanstd(response_pn{drug}{dir})/sqrt(size(response_pn{drug}{dir}, 1)), 'color', color(dir));
        hold on
    end
    legend(dscell_type, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(condition{drug})
    xlim([3 400])
end


figure
for drug = 1:4
    subplot(2,2,drug)
    response_s{drug} = max(exciseRows_empty(response_pn{drug}{1}), 0);
    errorbar(ctr_x, nanmean(response_s{drug}), nanstd(response_s{drug})/sqrt(size(response_s{drug}, 1)), 'color', color(1));
    hold on
    response_others{drug} = max(exciseRows_empty(cell2mat(response_pn{drug}(2:4)')), 0);
    errorbar(ctr_x, nanmean(response_others{drug}), nanstd(response_others{drug})/sqrt(size(response_others{drug}, 1)), 'color', color(2));
    legend('superior', 'others', 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(condition{drug})
    xlim([3 400])
end

figure
errorbar(ctr_x, nanmean(response_s{1}), nanstd(response_s{1})/sqrt(size(response_s{1}, 1)), 'color', color(1));
hold on
errorbar(ctr_x, nanmean(response_s{2}), nanstd(response_s{2})/sqrt(size(response_s{2}, 1)), 'color', color(2));
errorbar(ctr_x, nanmean(response_others{1}), nanstd(response_others{1})/sqrt(size(response_others{1}, 1)), 'color', color(3));
errorbar(ctr_x, nanmean(response_others{2}), nanstd(response_others{2})/sqrt(size(response_others{2}, 1)), 'color', color(4));
legend('control superior', 'SR superior', 'control others', 'SR others', 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike rate')
title('WT')
xlim([3 400])


figure
errorbar(ctr_x(1:5), nanmean(response_s{1}(:, 1:5)), nanstd(response_s{1}(:, 1:5))/sqrt(size(response_s{1}(:, 1:5), 1)), 'k');
hold on
errorbar(ctr_x(1:5), nanmean(response_others{1}(:, 1:5)), nanstd(response_others{1}(:, 1:5))/sqrt(size(response_others{1}(:, 1:5), 1)), 'r');
errorbar(ctr_x(1:5), nanmean(response_s{2}(:, 1:5)), nanstd(response_s{2}(:, 1:5))/sqrt(size(response_s{2}(:, 1:5), 1)), 'k--');
errorbar(ctr_x(1:5), nanmean(response_others{2}(:, 1:5)), nanstd(response_others{2}(:, 1:5))/sqrt(size(response_others{2}(:, 1:5), 1)), 'r--');
legend('control superior', 'control others', 'SR superior',  'SR others', 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike rate')
title('WT')
xlim([3 100])

%% fit ds
for dir = 1:4
    for drug = 1:4
        figure
        set(gcf, 'Position', [1 1 900 800])
        response_temp = exciseRows_empty(response_pn_norm{drug}{dir});
        for cc = 1:size(response_temp, 1)

            ydata = response_temp(cc, :);
            xdata = log10(ctr_x);
            [f, G] = fit_nr(xdata, ydata, 'upper', [100, 100, log10(1000), 0]);
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


%%
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'SR', 'SR+TPMPA', 'wash'};
step_size = 200;
for drug = 1:4
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                for ctr = 7:-1:1
                    a = raster_p_sum{drug}{idx_dir{dir}(cc)}{ctr};
                    hist_temp = hist(a, xx);
%                     if drug == 1 && ctr == 7
                        [max_p, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
                        Max_i{dir}(CC) = max_i;
%                     else
%                         max_p = conv(hist_temp, ones(1,step_size), 'valid');
%                         max_p = max_p(Max_i{dir}(CC));
%                     end
%                     response_pn{drug}{dir}(CC, ctr) = max(max_p - mean_n, 0);
%                     response_pn{drug}{dir}(CC, ctr) = max_p - mean_n;
%                     response_pn{drug}{dir}(CC, ctr) = abs(max_p - mean_n);
                    response_p{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
%                     response_pn{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions;
                end
                CC = CC + 1;
            end
        end
        response_p_norm{drug}{dir} = response_p{drug}{dir}./repmat(max(response_p{drug}{dir}, [], 2), 1, size(response_p{drug}{dir},2));
    end
end

ctr_x = [5 10 20 40 80 150 300];
color = 'brgkc';
figure
set(gcf, 'Position', [1 1 900 800])

for drug = 1:4
    subplot(2,2,drug)
    for dir = 1:4
        errorbar(ctr_x, nanmean(response_p_norm{drug}{dir}), nanstd(response_p_norm{drug}{dir})/sqrt(size(response_p_norm{drug}{dir}, 1)), 'color', color(dir));
        hold on
    end
    legend(dscell_type, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(condition{drug})
    xlim([3 400])
end



figure
set(gcf, 'Position', [1 1 900 800])

for drug = 1:4
    subplot(2,2,drug)
    for dir = 1:4
        response_s = exciseRows_empty(response_p_norm{drug}{1});
        errorbar(ctr_x, nanmean(response_s), nanstd(response_s)/sqrt(size(response_s, 1)), 'color', color(1));
        hold on
        response_others = exciseRows_empty(cell2mat(response_p_norm{drug}(2:4)'));
        errorbar(ctr_x, nanmean(response_others), nanstd(response_others)/sqrt(size(response_others, 1)), 'color', color(2));
    end
    legend('superior', 'others', 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(condition{drug})
    xlim([3 400])
end


ctr_x = [5 10 20 40 80 150 300];
color = 'brgkc';
figure
set(gcf, 'Position', [1 1 900 800])

for drug = 1:4
    subplot(2,2,drug)
    for dir = 1:4
        errorbar(ctr_x, nanmean(response_p{drug}{dir}), nanstd(response_p{drug}{dir})/sqrt(size(response_p{drug}{dir}, 1)), 'color', color(dir));
        hold on
    end
    legend(dscell_type, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(condition{drug})
    xlim([3 400])
end


figure
set(gcf, 'Position', [1 1 900 800])

for drug = 1:4
    subplot(2,2,drug)
    for dir = 1:4
        response_s_p{dir} = exciseRows_empty(response_p{drug}{1});
        errorbar(ctr_x, nanmean(response_s_p{dir}), nanstd(response_s_p{dir})/sqrt(size(response_s_p{dir}, 1)), 'color', color(1));
        hold on
        response_others_p{dir} = exciseRows_empty(cell2mat(response_p{drug}(2:4)'));
        errorbar(ctr_x, nanmean(response_others_p{dir}), nanstd(response_others_p{dir})/sqrt(size(response_others_p{dir}, 1)), 'color', color(2));
    end
    legend('superior', 'others', 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(condition{drug})
    xlim([3 400])
end

%%

drug = 1;
ctr = 5;
for dir = 1:4
    for cc = 1:length(id_dir_mb{dir})
        null_temp{dir}(cc) = (length(raster_n_sum{drug}{idx_dir_mb{dir}(cc)}{ctr})/trial_dur - bgnd_firing{dir}(drug, cc))/datamb{drug}.stimulus.repetitions;
    end
end
null{1} = null_temp{1};
null{2} = cell2mat(null_temp(2:4));


drug = 1;
ctr = 5;
for dir = 1:4
    for cc = 1:length(id_dir_mb{dir})
        prefer_temp{dir}(cc) = (length(raster_p_sum{drug}{idx_dir_mb{dir}(cc)}{ctr})/trial_dur - bgnd_firing{dir}(drug, cc))/datamb{drug}.stimulus.repetitions;
    end
end
prefer{1} = prefer_temp{1};
prefer{2} = cell2mat(prefer_temp(2:4));

figure
for cc = 1:length(null{1})
    plot(1, null{1}(cc), 'Marker', 'o', 'color', 'b')
    hold on
end
for cc = 1:length(null{2})
    plot(2, null{2}(cc), 'Marker', 'o', 'color', 'r')
end
h1 = errorbar(1, mean(null{1}), std(null{1})/sqrt(length(null{1})), 'Marker', 'd', 'color', 'b', 'markersize', 12);
h2 = errorbar(2, mean(null{2}), std(null{2})/sqrt(length(null{2})), 'Marker', 'd', 'color', 'r', 'markersize', 12);
ylabel('firing rate')
legend([h1, h2], 'superior', 'others')
xlim([0.5 2.5])

%%
ctr = 5;

idx1 = 45;
idx2 = 27;
figure
subplot(3, 2, 1)
plot_raster(squeeze(raster_n_sum_all{1}{idx1}(:, :, ctr, :)), 0, 3.5)
xlim([0 3.5])
subplot(3, 2, 3)
plot_raster(squeeze(raster_n_sum_all{1}{idx2}(:, :, ctr, :)), 0, 3.5)
xlim([0 3.5])
subplot(3, 2, 2)
plot_raster(squeeze(raster_n_sum_all{2}{idx1}(:, :, ctr, :)), 0, 3.5)
xlim([0 3.5])
subplot(3, 2, 4)
plot_raster(squeeze(raster_n_sum_all{2}{idx2}(:, :, ctr, :)), 0, 3.5)
xlim([0 3.5])


% figure
subplot(3, 2, 5)
bin_size = 0.025;
x = bin_size/2:bin_size:trial_dur;
temp = hist(raster_n_sum{1}{idx1}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{1}.stimulus.repetitions/bin_size;
filter = [0.25 0.5 0.25];
temp = conv(temp, filter, 'same');
plot(x,temp, 'b')
hold on
temp = hist(raster_n_sum{1}{idx2}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{1}.stimulus.repetitions/bin_size;
temp = conv(temp, filter, 'same');
plot(x,temp, 'r')
legend('superior', 'anterior')
xlabel('time (second)')
ylabel('firing rate (Hz)')
title('control')
xlim([0 3.5])

subplot(3, 2, 6)
temp = hist(raster_n_sum{2}{idx1}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{2}.stimulus.repetitions/bin_size;
filter = [0.25 0.5 0.25];
temp = conv(temp, filter, 'same');
plot(x,temp, 'b')
hold on
temp = hist(raster_n_sum{2}{idx2}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{2}.stimulus.repetitions/bin_size;
temp = conv(temp, filter, 'same');
plot(x,temp, 'r')
legend('superior', 'anterior')
xlabel('time (second)')
ylabel('firing rate (Hz)')
title('SR')
xlim([0 3.5])

%%
% ctr = 5;
% figure
% subplot(1, 2, 1)
% bin_size = 0.025;
% x = bin_size/2:bin_size:trial_dur;
% temp = hist(raster_n_sum{1}{45}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{1}.stimulus.repetitions/bin_size;
% % filter = [0.25 0.5 0.25];
% % temp = conv(temp, filter, 'same');
% plot(x,temp, 'b')
% hold on
% temp = hist(raster_n_sum{1}{27}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{1}.stimulus.repetitions/bin_size;
% % temp = conv(temp, filter, 'same');
% plot(x,temp, 'r')
% legend('superior', 'anterior')
% xlabel('time (second)')
% ylabel('firing rate (Hz)')
% title('control')
% 
% subplot(1, 2, 2)
% temp = hist(raster_n_sum{2}{45}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{2}.stimulus.repetitions/bin_size;
% % temp = conv(temp, filter, 'same');
% plot(x,temp, 'b')
% hold on
% temp = hist(raster_n_sum{2}{27}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{2}.stimulus.repetitions/bin_size;
% % temp = conv(temp, filter, 'same');
% plot(x,temp, 'r')
% legend('superior', 'anterior')
% xlabel('time (second)')
% ylabel('firing rate (Hz)')
% title('SR')
%%
for drug = 1:2
    response_sn{drug} = response_s{drug}./repmat(response_s_p{drug}(:, 5), 1, 7);
    response_othersn{drug} = response_others{drug}./repmat(response_others_p{drug}(:, 5), 1, 7);
end

figure
errorbar(ctr_x(1:5), nanmean(response_sn{1}(:, 1:5)), nanstd(response_sn{1}(:, 1:5))/sqrt(size(response_sn{1}(:, 1:5), 1)), 'k');
hold on
errorbar(ctr_x(1:5), nanmean(response_othersn{1}(:, 1:5)), nanstd(response_othersn{1}(:, 1:5))/sqrt(size(response_othersn{1}(:, 1:5), 1)), 'r');
errorbar(ctr_x(1:5), nanmean(response_sn{2}(:, 1:5)), nanstd(response_sn{2}(:, 1:5))/sqrt(size(response_sn{2}(:, 1:5), 1)), 'k--');
errorbar(ctr_x(1:5), nanmean(response_othersn{2}(:, 1:5)), nanstd(response_othersn{2}(:, 1:5))/sqrt(size(response_othersn{2}(:, 1:5), 1)), 'r--');
legend('control superior', 'control others', 'SR superior',  'SR others', 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike rate')
title('WT')
xlim([3 100])

for drug = 1:2
    response_sn_ko{drug} = response_s_ko{drug}./repmat(response_s_p_ko{drug}(:, 5), 1, 7);
    response_othersn_ko{drug} = response_others_ko{drug}./repmat(response_others_p_ko{drug}(:, 5), 1, 7);
end

figure
errorbar(ctr_x(1:5), nanmean(response_sn_ko{1}(:, 1:5)), nanstd(response_sn_ko{1}(:, 1:5))/sqrt(size(response_sn_ko{1}(:, 1:5), 1)), 'k');
hold on
errorbar(ctr_x(1:5), nanmean(response_othersn_ko{1}(:, 1:5)), nanstd(response_othersn_ko{1}(:, 1:5))/sqrt(size(response_othersn_ko{1}(:, 1:5), 1)), 'r');
errorbar(ctr_x(1:5), nanmean(response_sn_ko{2}(:, 1:5)), nanstd(response_sn_ko{2}(:, 1:5))/sqrt(size(response_sn_ko{2}(:, 1:5), 1)), 'k--');
errorbar(ctr_x(1:5), nanmean(response_othersn_ko{2}(:, 1:5)), nanstd(response_othersn_ko{2}(:, 1:5))/sqrt(size(response_othersn_ko{2}(:, 1:5), 1)), 'r--');
legend('control superior', 'control others', 'SR superior',  'SR others', 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike rate')
title('WT')
xlim([3 100])


%%
% ratio_s = (response_s{2} - response_s{1});%./response_s{1};
ratio_s = response_s{2}./response_s{1};
ratio_s_mean = mean(ratio_s);
ratio_s_ste = std(ratio_s)/sqrt(size(ratio_s, 1));
% ratio_others = (response_others{2} - response_others{1});%./response_others{1};
ratio_others = response_others{2}./response_others{1};
ratio_others(ratio_others == Inf) = nan;
ratio_others_mean = nanmean(ratio_others);
ratio_others_ste = nanstd(ratio_others)./sqrt(sum(~isnan(ratio_others)));

ratio_mean = [ratio_s_mean(1) ratio_s_mean(5); ratio_others_mean(1) ratio_others_mean(5)];
ratio_ste = [ratio_s_ste(1) ratio_s_ste(5); ratio_others_ste(1) ratio_others_ste(5)];
Ticks = {'low contrast', 'high contrast'};
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = Ticks;
model_series = ratio_mean';
model_error = ratio_ste';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('SR - control (spike #)')
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end
xlim([0.5 2.5])


%%
ctr = 5;
figure
for i = 1:size(response_s{1}, 1)
    plot([1 2], [response_s{1}(i, ctr) response_s{2}(i, ctr)], '-ko')
    hold on
end

for i = 1:size(response_others{1}, 1)
    plot([3 4], [response_others{1}(i, ctr) response_others{2}(i, ctr)], '-ko')
    hold on
end
xlim([0.5 4.5])

%% direction tuning curve (sliding window)

id_dir{2} = cell2mat(id_dir(2:4));
idx_dir{2} = cell2mat(idx_dir(2:4));
id_dir = id_dir(1:2);
idx_dir = idx_dir(1:2);

for dir = 1:2
    id_dir_mb{dir} = id_dir{dir}(~mb_idx(idx_dir{dir}));
    idx_dir_mb{dir} = idx_dir{dir}(~mb_idx(idx_dir{dir}));
end


clear Max_i
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'SR', 'SR+TPMPA', 'wash'};
step_size = 80;
for drug = 1:4
    for dir = 1:2
        for theta = 1:8
            CC = 1;
            for cc = 1:length(idx_dir{dir})
                if ~isempty(raster_mb_all{drug}{idx_dir{dir}(cc)})
                    for ctr = 7:-1:1
                        a = raster_mb_all{drug}{idx_dir{dir}(cc)}{1,1,theta,ctr};
                        hist_temp = hist(a, xx);
%                         if drug == 1 && ctr == 7
                            [max_fr, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
                            Max_i{dir}(theta, CC) = max_i;
%                         else
%                             max_fr = conv(hist_temp, ones(1,step_size), 'valid');
%                             max_fr = max_fr(Max_i{dir}(theta, CC));
%                         end
                        response_max{drug}{dir}(CC, theta, ctr) = max_fr/datamb{drug}.stimulus.repetitions - bgnd_firing_combine{dir}(CC, drug)*bin_size*step_size;
                    end
                    CC = CC + 1;
                end
            end
            response_max_norm{drug}{dir}(:, theta, :) = squeeze(response_max{drug}{dir}(:, theta, :))./repmat(squeeze(max(response_max{1}{dir}(:, theta, :), [], 1)), 1, size(response_max{drug}{dir},1))';
        end
    end
end

for drug = 1:4
    response_max_all{drug} = cell2mat(response_max{drug}');
end


for dir = 1:2
    for cc = 1:length(idx_dir_mb{dir})
        pindex{dir}(cc) = get_pindex(response_max{1}{dir}(cc, :, 5));
        for drug = 1:4
            for ctr = 1:7
                response_max{drug}{dir}(cc, :, ctr) = circshift(response_max{drug}{dir}(cc, :, ctr), [0, 4-pindex{dir}(cc)]);
                response_max_norm{drug}{dir}(cc, :, ctr) = response_max{drug}{dir}(cc, :, ctr)/max(response_max{drug}{dir}(cc, :, ctr));
            end
        end
    end
end

for drug = 1:4
    for dir = 1:2
        for ctr = 1:7
            response_max_mean{drug}{dir}{ctr} = mean(response_max{drug}{dir}(:, :, ctr));
            response_max_ste{drug}{dir}{ctr} = std(response_max{drug}{dir}(:, :, ctr))./sqrt(size(response_max{drug}{dir}, 1));
            response_max_norm_mean{drug}{dir}{ctr} = nanmean(response_max_norm{drug}{dir}(:, :, ctr));
            response_max_norm_ste{drug}{dir}{ctr} = nanstd(response_max_norm{drug}{dir}(:, :, ctr))./sqrt(size(response_max_norm{drug}{dir}, 1));
        end
    end
end

theta = linspace(-pi, pi, 9);
theta = theta(1:8);
color = 'brgkcmy';
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
contrast = {'5%', '10%', '20%', '40%', '80%', '150%', '300%'};

figure
for dir = 1:2
    for ctr = 1:7
        subplot(2, 7, 7 * (dir - 1) + ctr);
        for drug = 1:2
            errorbar(theta, response_max_norm_mean{drug}{dir}{ctr}, response_max_norm_ste{drug}{dir}{ctr}, 'color', color(drug));
            hold on
        end
        xlim([-3.5 3])
        if dir == 1
            title(contrast{ctr})
            if ctr == 7
                legend('control', 'SR')
            end
        end
        if (ctr == 1)
            ylabel(dscell_type{dir});
        end
        if (dir == 4)
            xlabel('direction')
        end
    end
end

figure
for dir = 1:2
    for ctr = 1:7
        subplot(2, 7, 7 * (dir - 1) + ctr);
        for drug = 1:2
            errorbar(theta, response_max_mean{drug}{dir}{ctr}, response_max_ste{drug}{dir}{ctr}, 'color', color(drug));
            hold on
        end
        xlim([-3.5 3])
        if dir == 1
            title(contrast{ctr})
            if ctr == 7
                legend('control', 'SR')
            end
        end
        if (ctr == 1)
            ylabel(dscell_type{dir});
        end
        if (dir == 4)
            xlabel('direction')
        end
    end
end

figure
for drug = 1:4
    for ctr = 1:7
        subplot(4, 7, 7 * (drug - 1) + ctr);
        for dir = 1:2
            errorbar(theta, response_max_norm_mean{drug}{dir}{ctr}, response_max_norm_ste{drug}{dir}{ctr}, 'color', color(dir));
            hold on
        end
        xlim([-3.5 3])
        if drug == 1
            title(contrast{ctr})
            if ctr == 7
                legend('superior', 'others')
            end
        end
        if (ctr == 1)
            ylabel(condition{drug});
        end
        if (drug == 4)
            xlabel('direction')
        end
        ylim([0 1])
    end
end

for dir = 1:2
figure
for drug = 1:4
    for ctr = 1:7
        subplot(4, 7, 7 * (drug - 1) + ctr);
        
        for cc = 1:size(response_max_norm{1}{1}, 1)
            if max(response_max_norm{drug}{dir}(cc, :, ctr)) == 1
                plot(theta, response_max_norm{drug}{dir}(cc, :, ctr), 'color', color(dir));
                hold on
            end
        end
        xlim([-3.5 3])
        if drug == 1
            title(contrast{ctr})
        end
        if (ctr == 1)
            ylabel(condition{drug});
        end
        if (drug == 4)
            xlabel('direction')
        end
        ylim([0 1])
    end
end
end


%% unnormalized tuning curves
figure
for drug = 1:4
    for ctr = 1:7
        subplot(4, 7, 7 * (drug - 1) + ctr);
        for dir = 1:2
            errorbar(theta, response_max_mean{drug}{dir}{ctr}, response_max_ste{drug}{dir}{ctr}, 'color', color(dir));
            hold on
        end
        xlim([-3.5 3])
        if drug == 1
            title(contrast{ctr})
            if ctr == 7
                legend('superior', 'others')
            end
        end
        if (ctr == 1)
            ylabel(condition{drug});
        end
        if (drug == 4)
            xlabel('direction')
        end
%         ylim([0 1])
    end
end

for dir = 1:2
figure
for drug = 1:4
    for ctr = 1:7
        subplot(4, 7, 7 * (drug - 1) + ctr);
        
        for cc = 1:size(response_max{1}{1}, 1)
            plot(theta, response_max{drug}{dir}(cc, :, ctr), 'color', color(dir));
            hold on
        end
        xlim([-3.5 3])
        if drug == 1
            title(contrast{ctr})
        end
        if (ctr == 1)
            ylabel(condition{drug});
        end
        if (drug == 4)
            xlabel('direction')
        end
%         ylim([0 1])
    end
end
end

%% SNR at SR condition
color = 'brgk';
drug = 2;
rep = datamb{drug}.stimulus.repetitions;

for ctr = 1:5
for ct = 1:4
    [SpikeN_temp{ct}, SpikeN_mean_temp{ct}, SpikeN_std_temp{ct}] = deal(cell(length(id_dir{ct}), 1));
    for cc = 1:length(id_dir{ct})
        raster_temp = raster_mb{drug}{idx_dir{ct}(cc)};
        if ~isempty(raster_temp)
            SpikeN_temp{ct}{cc} = squeeze(cellfun('length',raster_temp(1,1,:,ctr,:)));
            SpikeN_mean_temp{ct}{cc} = mean(SpikeN_temp{ct}{cc},2);
            SpikeN_std_temp{ct}{cc} = std(SpikeN_temp{ct}{cc},[],2);
        end
    end
end
SpikeN = SpikeN_temp;
SpikeN_mean = SpikeN_mean_temp;
SpikeN_std = SpikeN_std_temp;

clear SpikeN_mean_all SpikeN_std_all
SpikeN_mean_all{1} = SpikeN_mean{1};
SpikeN_std_all{1} = SpikeN_std{1};
SpikeN_mean_all{2} = cat(1,SpikeN_mean{2:4});
SpikeN_std_all{2} = cat(1,SpikeN_std{2:4});
for i = 1:2
    SpikeN_mean_all{i}(any(cellfun(@isempty, SpikeN_mean_all{i}), 2), :) = [];
    SpikeN_std_all{i}(any(cellfun(@isempty, SpikeN_std_all{i}), 2), :) = [];
end

for i = 1:2
    [rmax{i}, index] = cellfun(@max, SpikeN_mean_all{i});
    for cc = 1:size(SpikeN_mean_all{i}, 1)
        stdmax{i}(cc) = SpikeN_std_all{i}{cc}(index(cc));
    end
    SNRmax{i}(:, ctr) = rmax{i}./stdmax{i}';
%     SNRmax{i} = exciseRows_empty(SNRmax{i});
end

end


ctr_x = [5 10 20 40 80];
figure
for i = 1:2
    errorbar(ctr_x, mean(SNRmax{i}), std(SNRmax{i})/sqrt(size(SNRmax{i}, 1)))
    hold on
end
set(gca,'XScale', 'log')
legend('superior','others')
ylim([0 5])
xlim([3 100])
ylabel('SNR')
xlabel('% contrast')
title('SR PD')