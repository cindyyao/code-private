%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load moving bar data
datarun = load_data('/Volumes/lab/analysis/2016-10-17-0/data003-007-map/data003-007-map', opt);
time_points = [1900 3800 5700];
datamb(1:4) = split_datarun(datarun, time_points);

datamb{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-10-17-0/stimuli/s03.txt';
datamb{1} = load_stim(datamb{1}, 'user_defined_trigger_set', [1:2:1120]);
datamb{1}.stimulus.triggers = datamb{1}.stimulus.triggers';
datamb{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-10-17-0/stimuli/s05.txt';
datamb{2} = load_stim(datamb{2}, 'user_defined_trigger_set', [1:2:1120]);
datamb{2}.stimulus.triggers = datamb{2}.stimulus.triggers';
datamb{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-10-17-0/stimuli/s06.txt';
datamb{3} = load_stim(datamb{3}, 'user_defined_trigger_set', [1:2:1120]);
datamb{3}.stimulus.triggers = datamb{3}.stimulus.triggers';
datamb{4}.names.stimulus_path = '/Volumes/lab/analysis/2016-10-17-0/stimuli/s07.txt';
datamb{4} = load_stim(datamb{4}, 'user_defined_trigger_set', [1:2:1120]);
datamb{4}.stimulus.triggers = datamb{4}.stimulus.triggers';

load('DS161017.mat')

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

%% plot single cell

for ct = 2:2
    for cc = 5:5 %length(id_dir{ct});
        plot_mb_raster_ctr_one(MB(1), raster_mb(1), trial_dur(1), idx_dir{ct}(cc), id_dir{ct}(cc), '', 1, 1, 0)
    end
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

%% get spontaneous activity
for drug = 1:4
    for dir = 1:4
        for cc = 1:length(id_dir_mb{dir})
            idx = get_cell_indices(datamb{drug},id_dir_mb{dir}(cc));
            spikes_temp = datamb{drug}.spikes{idx};
            bgnd_firing{dir}(drug, cc) = length(spikes_temp(spikes_temp < 1830 & spikes_temp > 1780))/50;
        end
    end
    for cc = 1:length(id_sub{1})
        if ~mb_idx(idx_sub{1}(cc))
            idx = get_cell_indices(datamb{drug},id_sub{1}(cc));
            spikes_temp = datamb{drug}.spikes{idx};
            bgnd_firing_on(drug, cc) = length(spikes_temp(spikes_temp < 1830 & spikes_temp > 1780))/50;
        end
    end
    for cc = 1:length(id_ot_mb)
        idx = get_cell_indices(datamb{drug},id_ot_mb(cc));
        spikes_temp = datamb{drug}.spikes{idx};
        bgnd_firing_ot(drug, cc) = length(spikes_temp(spikes_temp < 1830 & spikes_temp > 1780))/50;
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
                    response_pmax{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
                end
                CC = CC + 1;
            end
        end
        response_pmax_norm{drug}{dir} = response_pmax{drug}{dir}./repmat(max(response_pmax{1}{dir}, [], 2), 1, size(response_pmax{drug}{dir},2));
%         response_pmax_norm{drug}{dir} = response_pmax{drug}{dir}./repmat(max(response_pmax{drug}{dir}, [], 2), 1, size(response_pmax{drug}{dir},2));
    end
    response_pmax_all{drug} = cell2mat(response_pmax{drug}');
    response_pmax_norm_all{drug} = cell2mat(response_pmax_norm{drug}');
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

%% direction tuning (sliding window)

clear Max_i
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'wash'};
% step_size = 10;
for drug = 1:4
    for dir = 1:4
        for theta = 1:8
            CC = 1;
            for cc = 1:length(idx_dir{dir})
                if ~isempty(raster_mb_all{drug}{idx_dir{dir}(cc)})
                    for ctr = 7:-1:1
                        a = raster_mb_all{drug}{idx_dir{dir}(cc)}{1,1,theta,ctr};
                        hist_temp = hist(a, xx);
                        if drug == 1 && ctr == 7
                            [max_fr, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
                            Max_i{dir}(theta, CC) = max_i;
                        else
                            max_fr = conv(hist_temp, ones(1,step_size), 'valid');
                            max_fr = max_fr(Max_i{dir}(theta, CC));
                        end
                        response_max{drug}{dir}(CC, theta, ctr) = max_fr/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
                    end
                    CC = CC + 1;
                end
            end
            response_max_norm{drug}{dir}(:, theta, :) = squeeze(response_max{drug}{dir}(:, theta, :))./repmat(squeeze(max(response_max{1}{dir}(:, theta, :), [], 1)), 1, size(response_max{drug}{dir},1))';
        end
    end
end


for dir = 1:4
    for cc = 1:length(idx_dir_mb{dir})
        pindex{dir}(cc) = get_pindex(response_max{1}{dir}(cc, :, 7));
        for drug = 1:4
            for ctr = 1:7
                response_max{drug}{dir}(cc, :, ctr) = circshift(response_max{drug}{dir}(cc, :, ctr), [0, 4-pindex{dir}(cc)]);
                response_max_norm{drug}{dir}(cc, :, ctr) = response_max{drug}{dir}(cc, :, ctr)/max(response_max{drug}{dir}(cc, :, ctr));
            end
        end
    end
end

for drug = 1:4
    for dir = 1:4
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
for dir = 1:4
    for ctr = 1:7
        subplot(4, 7, 7 * (dir - 1) + ctr);
        for drug = [1 2 4]
            errorbar(theta, response_max_norm_mean{drug}{dir}{ctr}, response_max_norm_ste{drug}{dir}{ctr}, 'color', color(drug));
            hold on
        end
        xlim([-3.5 3])
        if dir == 1
            title(contrast{ctr})
            if ctr == 7
                legend('control', 'AP5', 'wash')
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

%% paper figure
% direction tuning
theta = linspace(-180, 180, 9);
theta = theta(2:9);
ctr = 6; % contrast = 80%
figure(4)
subplot(3, 2, 1)
for drug = 1:2
    response_max_all{drug} = cell2mat(response_max{drug}');
    tuning_mean(drug, :) = squeeze(mean(response_max_all{drug}(:, :, ctr)));
    tuning_ste(drug, :) = squeeze(std(response_max_all{drug}(:, :, ctr))./sqrt(size(response_max_all{drug}, 1)));
    errorbar(theta, tuning_mean(drug, :), tuning_ste(drug, :), [color(drug) 'o'])
    hold on
    [f, g] = fit_cos(theta/180*pi, tuning_mean(drug, :));
    xfit = linspace(-180, 180, 100);
    yfit = f.ymax * (0.5 + 0.5 * cos(xfit/180 *pi + f.phi)).^f.alpha + f.b;
    h(drug) = plot(xfit, yfit, color(drug));
end
xlabel('direction')
ylabel('spike number')
legend([h(1), h(2)], 'control', 'AP5')

ctr = 3;
subplot(3, 2, 2)
for drug = 1:2
    tuning_mean(drug, :) = squeeze(mean(response_max_all{drug}(:, :, ctr)));
    tuning_ste(drug, :) = squeeze(std(response_max_all{drug}(:, :, ctr))./sqrt(size(response_max_all{drug}, 1)));
    errorbar(theta, tuning_mean(drug, :), tuning_ste(drug, :), [color(drug) 'o'])
    hold on
    [f, g] = fit_cos(theta/180*pi, tuning_mean(drug, :));
    xfit = linspace(-180, 180, 100);
    yfit = f.ymax * (0.5 + 0.5 * cos(xfit/180 *pi + f.phi)).^f.alpha + f.b;
    h(drug) = plot(xfit, yfit, color(drug));
end
xlabel('direction')
ylabel('spike number')
legend([h(1), h(2)], 'control', 'AP5')


% heat map
conditions = {'control', 'AP5'};
[xx, yy] = meshgrid(theta, ctr_x);
for drug = 1:2
    figure
    response_max_all_mean{drug} = squeeze(mean(response_max_all{drug}));
    s = surf(xx, yy, response_max_all_mean{drug}'/max(response_max_all_mean{1}(:)));
    s.EdgeColor = 'none';
    set(gca, 'yscale', 'log')
%     xlabel('direction (degree)')
%     ylabel('contrast')
%     zlabel('spike number')
%     title(conditions{drug})
    xlim([min(xx(:)) max(xx(:))])
    ylim([min(yy(:)) max(yy(:))])
    caxis([0 1])
    axis off
end

% CRF 
figure(4)
subplot(3, 2, 3)
dir = 3;
for drug = 1:2
    crf_mean(drug, :) = mean(response_max_all{drug}(:, dir, :));
    crf_ste(drug, :) = std(response_max_all{drug}(:, dir, :))./sqrt(size(response_max_all{drug}, 1));
    errorbar(ctr_x, crf_mean(drug, :), crf_ste(drug, :), [color(drug) 'o'])
    hold on
    [f, g] = fit_nr(log10(ctr_x), crf_mean(drug, :), 'Upper', [Inf, Inf, max(log10(ctr_x)), min(crf_mean(drug, :))]);
    xfit = linspace(log10(4), log10(1000), 100);
    yfit = f.ymax * xfit.^f.a./(xfit.^f.a + f.sigma^f.a) + f.b;
    h(drug) = plot(10.^xfit, yfit, color(drug));
end
set(gca, 'xscale', 'log')
xlabel('contrast')
ylabel('firing rate')
legend([h(1), h(2)], 'control', 'AP5')

dir = 7;
subplot(3, 2, 4)
for drug = 1:2
    crf_mean(drug, :) = mean(response_max_all{drug}(:, dir, :));
    crf_ste(drug, :) = std(response_max_all{drug}(:, dir, :))./sqrt(size(response_max_all{drug}, 1));
    errorbar(ctr_x, crf_mean(drug, :), crf_ste(drug, :), [color(drug) 'o'])
    hold on
    [f, g] = fit_nr(log10(ctr_x), crf_mean(drug, :), 'Upper', [Inf, Inf, max(log10(ctr_x)), max(min(crf_mean(drug, :)), 0)]);
    xfit = linspace(log10(4), log10(1000), 100);
    yfit = f.ymax * xfit.^f.a./(xfit.^f.a + f.sigma^f.a) + f.b;
    h(drug) = plot(10.^xfit, yfit, color(drug));
end
set(gca, 'xscale', 'log')
xlabel('contrast')
ylabel('firing rate')
legend([h(1), h(2)], 'control', 'AP5')
ylim([0 2])

a = (response_max_all{1} - response_max_all{2})./response_max_all{1};
for i = 1:8
    for j = 1:7
        temp = a(:, i, j);
        temp(isinf(temp)) = [];
        temp(isnan(temp)) = [];
        temp(temp < 0) = 0;
        temp(temp > 1) = 1;
        blockMean(i, j) = mean(temp);
        blockSte(i, j) = std(temp)/sqrt(length(temp));
        sizetemp(i, j) = length(temp);
    end
end
        
% a(isinf(a)) = nan;
% a(a<0) = 0;
% a(a>1) = 1;
% b = squeeze(nanmean(a));
% c = squeeze(nanstd(a)./sqrt(size(a, 1)));
figure(3)
s = surf(xx, yy, blockMean');
s.EdgeColor = 'none';
set(gca, 'yscale', 'log')
% xlabel('direction (degree)')
% ylabel('contrast')
% zlabel('spike number')
xlim([min(min(xx)) max(max(xx))])
ylim([min(yy(:)) max(yy(:))])
caxis([0 1])
axis off
colorbar



i = 4;
figure
errorbar(theta, blockMean(:, i), blockSte(:, i))
ylim([0 1])
% xlim([-0.1 1.1])
legend('150%', '20%')

j = 4;
figure
errorbar(ctr_x, blockMean(j, :), blockSte(j, :))
ylim([0 1])
xlim([4 400])
set(gca, 'xscale', 'log')



%% paper figure (revert)
theta = linspace(-180, 180, 9);
theta = theta(2:9);
ctr = 4; % contrast = 40%
figure(4)
subplot(3, 2, 1)
for drug = 1:2
    response_max_all{drug} = cell2mat(response_max{drug}');
    tuning_mean(drug, :) = squeeze(mean(response_max_all{drug}(:, :, ctr)));
    tuning_ste(drug, :) = squeeze(std(response_max_all{drug}(:, :, ctr))./sqrt(size(response_max_all{drug}, 1)));
    errorbar(theta, tuning_mean(drug, :), tuning_ste(drug, :), [color(drug) 'o'])
    hold on
    [f, g] = fit_cos(theta/180*pi, tuning_mean(drug, :));
    xfit = linspace(-180, 180, 100);
    yfit = f.ymax * (0.5 + 0.5 * cos(xfit/180 *pi + f.phi)).^f.alpha + f.b;
    h(drug) = plot(xfit, yfit, color(drug));
end
xlabel('direction')
ylabel('spike number')
legend([h(1), h(2)], 'control', 'AP5')


% CRF 
figure(4)
subplot(3, 2, 3)
dir = 4;
for drug = 1:2
    crf_mean(drug, :) = mean(response_max_all{drug}(:, dir, :));
    crf_ste(drug, :) = std(response_max_all{drug}(:, dir, :))./sqrt(size(response_max_all{drug}, 1));
    errorbar(ctr_x, crf_mean(drug, :), crf_ste(drug, :), [color(drug) 'o'])
    hold on
    [f, g] = fit_nr(log10(ctr_x), crf_mean(drug, :), 'Upper', [Inf, Inf, max(log10(ctr_x)), min(crf_mean(drug, :))]);
    xfit = linspace(log10(4), log10(1000), 100);
    yfit = f.ymax * xfit.^f.a./(xfit.^f.a + f.sigma^f.a) + f.b;
    h(drug) = plot(10.^xfit, yfit, color(drug));
end
set(gca, 'xscale', 'log')
xlabel('contrast')
ylabel('firing rate')
legend([h(1), h(2)], 'control', 'AP5')

a = (response_max_all{1} - response_max_all{2})./response_max_all{1};
for i = 1:8
    for j = 1:7
        temp = a(:, i, j);
        temp(isinf(temp)) = [];
        temp(isnan(temp)) = [];
        temp(temp < 0) = 0;
        temp(temp > 1) = 1;
        blockMean(i, j) = mean(temp);
        blockSte(i, j) = std(temp)/sqrt(length(temp));
        sizetemp(i, j) = length(temp);
    end
end
        
% a(isinf(a)) = nan;
% a(a<0) = 0;
% a(a>1) = 1;
% b = squeeze(nanmean(a));
% c = squeeze(nanstd(a)./sqrt(size(a, 1)));
[xx, yy] = meshgrid(theta, ctr_x);
figure(3)
s = surf(xx, yy, blockMean');
s.EdgeColor = 'none';
set(gca, 'yscale', 'log')
% xlabel('direction (degree)')
% ylabel('contrast')
% zlabel('spike number')
xlim([min(min(xx)) max(max(xx))])
ylim([min(yy(:)) max(yy(:))])
caxis([0 1])
axis off
colorbar

%
i = 4;
figure
errorbar(theta, blockMean(:, i), blockSte(:, i))
ylim([0 1])
% xlim([-0.1 1.1])
legend('150%', '20%')

j = 4;
figure
errorbar(ctr_x, blockMean(j, :), blockSte(j, :))
ylim([0 1])
xlim([4 400])
set(gca, 'xscale', 'log')

% heat map
conditions = {'control', 'AP5'};
[xx, yy] = meshgrid(theta, ctr_x);
for drug = 1:2
    figure(drug)
    response_max_all_mean{drug} = squeeze(mean(response_max_all{drug}));
    s = surf(xx, yy, response_max_all_mean{drug}'/max(response_max_all_mean{1}(:)));
    s.EdgeColor = 'none';
    set(gca, 'yscale', 'log')
%     xlabel('direction (degree)')
%     ylabel('contrast')
%     zlabel('spike number')
%     title(conditions{drug})
    xlim([min(xx(:)) max(xx(:))])
    ylim([min(yy(:)) max(yy(:))])
    caxis([0 1])
    axis off
end

%%
a_null = squeeze(a(:, 1, :));
a_pref = squeeze(a(:, 4, :));
figure
subplot(1, 2, 1)
for cc = 1:size(a_null, 1)
plot(ctr_x, a_null(cc, :))
hold on
end
set(gca, 'xscale', 'log')
xlabel('contrast')
ylabel('% blocked')
xlim([4, 400])
title('null direction')

subplot(1, 2, 2)
for cc = 1:size(a_pref, 1)
plot(ctr_x, a_pref(cc, :))
hold on
end
set(gca, 'xscale', 'log')
xlabel('contrast')
ylabel('% blocked')
xlim([4, 400])
title('prefered direction')

for c = 1:7
    a_null_temp = a(:, 1, c);
    a_null_temp(isnan(a_null_temp)) = [];
    a_pref_temp = a(:, 4, c);
    a_pref_temp(isnan(a_pref_temp)) = [];
    a_null_mean(c) = mean(a_null_temp);
    a_null_ste(c) = std(a_null_temp)/sqrt(length(a_null_temp));
    a_pref_mean(c) = mean(a_pref_temp);
    a_pref_ste(c) = std(a_pref_temp)/sqrt(length(a_pref_temp));
end


figure
errorbar(ctr_x, a_pref_mean, a_pref_ste, 'b')
hold on
errorbar(ctr_x, a_null_mean, a_null_ste, 'r')
set(gca, 'xscale', 'log')
xlabel('contrast')
ylabel('% blocked')
legend('pref', 'null')


%%

for i = 1:8
    spike = response_max_all{1}(:, i, :);
    spike = spike(:);
    ratio = a(:, i, :);
    ratio = ratio(:);
    spike(ratio == inf | isnan(ratio) | ratio < 0 | ratio > 1) = [];
    ratio(ratio == inf | isnan(ratio) | ratio < 0 | ratio > 1) = [];
    ratio_spike_dir{i}(1, :) = spike;
    ratio_spike_dir{i}(2, :) = ratio;
end

for i = 1:7
    spike = response_max_all{1}(:, :, i);
    spike = spike(:);
    ratio = a(:, :, i);
    ratio = ratio(:);
    spike(ratio == inf | isnan(ratio) | ratio < 0 | ratio > 1) = [];
    ratio(ratio == inf | isnan(ratio) | ratio < 0 | ratio > 1) = [];
    ratio_spike_ctr{i}(1, :) = spike;
    ratio_spike_ctr{i}(2, :) = ratio;
end

clear spikeCountAvg blockAvg disError
i = 6;j = 3;
[spikeCountAvg{1}(1, :), blockAvg{1}(1, :), disError{1}(1, :)] = curve_from_binning(ratio_spike_ctr{i}(1, :), ratio_spike_ctr{i}(2, :), 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:2:40);
[spikeCountAvg{1}(2, :), blockAvg{1}(2, :), disError{1}(2, :)] = curve_from_binning(ratio_spike_dir{j}(1, :), ratio_spike_dir{j}(2, :), 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:2:40);


i = 3;j = 7;
[spikeCountAvg{2}(1, :), blockAvg{2}(1, :), disError{2}(1, :)] = curve_from_binning(ratio_spike_ctr{i}(1, :), ratio_spike_ctr{i}(2, :), 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:0.5:10);
[spikeCountAvg{2}(2, :), blockAvg{2}(2, :), disError{2}(2, :)] = curve_from_binning(ratio_spike_dir{j}(1, :), ratio_spike_dir{j}(2, :), 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:0.5:10);

figure(5)
subplot(1, 2, 1)
errorbar(spikeCountAvg{1}(1, :), blockAvg{1}(1, :), disError{1}(1, :), 'b');
hold on
errorbar(spikeCountAvg{1}(2, :), blockAvg{1}(2, :), disError{1}(1, :), 'r');
xlabel('spike number')
ylabel('% block')
legend('150%', 'PD')


subplot(1, 2, 2)
errorbar(spikeCountAvg{2}(1, :), blockAvg{2}(1, :), disError{2}(1, :), 'b');
hold on
errorbar(spikeCountAvg{2}(2, :), blockAvg{2}(2, :), disError{2}(1, :), 'r');
xlabel('spike number')
ylabel('% block')
legend('20%', 'ND')

%%


theta = linspace(-pi, pi, 9);
theta = theta(2:9);
color = 'brgkcmy';
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
contrast = {'5%', '10%', '20%', '40%', '80%', '150%', '300%'};

figure
for ctr = 1:7
    subplot(1, 7, ctr);
    for drug = 1:2
        response_max_all_norm{drug}(:, :, ctr) = response_max_all{drug}(:, :, ctr)./repmat(max(response_max_all{drug}(:, :, ctr), [], 2), 1, 8);
        errorbar(theta, nanmean(response_max_all_norm{drug}(:, :, ctr)), nanstd(response_max_all_norm{drug}(:, :, ctr))/sqrt(size(response_max_all_norm{drug}, 1)), 'color', color(drug));
        hold on
    end
    xlim([-3.5 3])
    title(contrast{ctr})
    if ctr == 7
        legend('control', 'AP5')
    end
    xlabel('direction')
end


%%
color = 'brgkcmy';
ct = 2; cc = 5; theta = linspace(-pi, pi, 9); theta = theta(2:9);
raster_temp = squeeze(raster_mb{1}{idx_dir{ct}(cc)});
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;

for dir = 1:8
    for ctr = 1:7
        for repeat = 1:10
            a = raster_mb{1}{idx_dir{ct}(cc)}{1,1,dir,ctr,repeat};
            hist_temp = hist(a, xx);
            max_fr = conv(hist_temp, ones(1, step_size), 'valid');
            max_fr = max_fr(Max_i{ct}(dir, cc));
            raster_mb_max(dir, ctr, repeat) = max_fr;
        end
    end
end
raster_mb_max_mean = mean(raster_mb_max, 3);
raster_mb_max_mean = circshift(raster_mb_max_mean, [3, 0]);
raster_mb_max_ste = std(raster_mb_max, [], 3)/sqrt(size(raster_mb_max, 3));
raster_mb_max_ste = circshift(raster_mb_max_ste, [3, 0]);

figure
for ctr = 1:7
    errorbar(theta/pi*180, raster_mb_max_mean(:, ctr), raster_mb_max_ste(:, ctr), [color(8-ctr) 'o'])
    hold on 
    [f, g] = fit_cos(theta, raster_mb_max_mean(:, ctr)');
    xfit = linspace(-180, 180, 100);
    yfit = f.ymax * (0.5 + 0.5 * cos(xfit/180 *pi + f.phi)).^f.alpha + f.b;
    plot(xfit, yfit, color(8-ctr));
end
    
