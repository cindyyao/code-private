cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

datadg = load_data('/Volumes/lab/analysis/2017-07-24-0/data012-sorted/data012-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2017-07-24-0/stimuli/s12.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);
% datadg.stimulus.triggers = datadg.stimulus.triggers';
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

params_idx = [1 3]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

datarun = load_data('/Volumes/lab/analysis/2017-07-24-0/datamb-map/datamb-map', opt);
time_points = [2100 3200 4273];
datamb(1:4) = split_datarun(datarun, time_points);
datamb{1}.names.stimulus_path = '/Volumes/lab/analysis/2017-07-24-0/stimuli/s00.txt';
datamb{1} = load_stim(datamb{1}, 'user_defined_trigger_set', [1:2:1120]);
datamb{1}.stimulus.triggers = datamb{1}.stimulus.triggers';
datamb{2}.names.stimulus_path = '/Volumes/lab/analysis/2017-07-24-0/stimuli/s02.txt';
datamb{2} = load_stim(datamb{2}, 'user_defined_trigger_set', [1:2:560]);
datamb{2}.stimulus.triggers = datamb{2}.stimulus.triggers';
datamb{3}.names.stimulus_path = '/Volumes/lab/analysis/2017-07-24-0/stimuli/s05.txt';
datamb{3} = load_stim(datamb{3}, 'user_defined_trigger_set', [1:2:560]);
datamb{3}.stimulus.triggers = datamb{3}.stimulus.triggers';
datamb{4}.names.stimulus_path = '/Volumes/lab/analysis/2017-07-24-0/stimuli/s07.txt';
datamb{4} = load_stim(datamb{4}, 'user_defined_trigger_set', [1:2:560]);
datamb{4}.stimulus.triggers = datamb{4}.stimulus.triggers';
load DS170724.mat
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

%%
d = 1;
t = 1;
h = figure;
dirn = 2;
set(h, 'Position', [1 1 1080 500])
compass(DG{d}.U{t}, DG{d}.V{t});
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}, DG{d}.V{t}, x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = I;
    id_dir{i} = ds_id(idx_dir{i});
end



%% 
for dir = 1:2
    id_dir_mb{dir} = id_dir{dir}(~mb_idx(idx_dir{dir}));
    idx_dir_mb{dir} = idx_dir{dir}(~mb_idx(idx_dir{dir}));
end
for dir = 1:3
    id_dir_on_mb{dir} = id_dir_on{dir}(~mb_idx(idx_dir_on{dir}));
    idx_dir_on_mb{dir} = idx_dir_on{dir}(~mb_idx(idx_dir_on{dir}));
end

%% plot cell summary
for dir = 1:4
    for cc = 1:length(id_dir{dir})
        plot_mb_raster_ctr(MB, raster_mb, trial_dur, idx_dir{dir}(cc), id_dir{dir}(cc), '', 4, 7, 1)
    end
end

% for dir = 1:3
%     for cc = 1:length(id_dir_on{dir})
%         plot_mb_raster_ctr(MB, raster_mb, trial_dur, idx_dir_on{dir}(cc), id_dir_on{dir}(cc), '', 4, 7, 1)
%     end
% end


%% get spontaneous activity

for drug = 1:4
    for dir = 1:2
        for cc = 1:length(id_dir_mb{dir})
            if(~mb_idx(idx_dir_mb{dir}(cc)))
                idx = get_cell_indices(datamb{drug},id_dir_mb{dir}(cc));
                spikes_temp = datamb{drug}.spikes{idx};
                if drug == 2
                    bgnd_firing{dir}(drug, cc) = length(spikes_temp(spikes_temp > 1020 & spikes_temp < 1040))/20;
                else
                    bgnd_firing{dir}(drug, cc) = length(spikes_temp(spikes_temp > datamb{drug}.duration - 20 & spikes_temp < datamb{drug}.duration))/20;
%                 bgnd_firing{dir}(drug, cc) = length(spikes_temp(spikes_temp > 40 & spikes_temp < 60))/20;
            
                end
            end
        end
    end
end
bgnd_firing_combine{1} = exciseRows_empty(bgnd_firing{1}');
bgnd_firing_combine{2} = exciseRows_empty(bgnd_firing{2}');


figure
errorbar(mean(bgnd_firing_combine{2}(:, 1)),  mean(bgnd_firing_combine{1}(:, 1)), std(bgnd_firing_combine{1}(:, 1))/sqrt(size(bgnd_firing_combine{1}, 1)), 'bo')
hold on
errorbar(mean(bgnd_firing_combine{2}(:, 2)),  mean(bgnd_firing_combine{1}(:, 2)), std(bgnd_firing_combine{1}(:, 2))/sqrt(size(bgnd_firing_combine{1}, 1)), 'ro')

herrorbar(mean(bgnd_firing_combine{2}(:, 1)),  mean(bgnd_firing_combine{1}(:, 1)), std(bgnd_firing_combine{2}(:, 1))/sqrt(size(bgnd_firing_combine{2}, 1)), 'bo')
herrorbar(mean(bgnd_firing_combine{2}(:, 2)),  mean(bgnd_firing_combine{1}(:, 2)), std(bgnd_firing_combine{2}(:, 2))/sqrt(size(bgnd_firing_combine{2}, 1)), 'ro')
plot([0 4], [0 4], 'k--')
xlim([0 4])
ylim([0 4])
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
    for dir = 1:2
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
    for dir = 1:2
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
    for dir = 1:2
        response_s = exciseRows_empty(response_pn_norm{drug}{1});
        errorbar(ctr_x, nanmean(response_s), nanstd(response_s)/sqrt(size(response_s, 1)), 'color', color(1));
        hold on
        response_others = exciseRows_empty(response_pn_norm{drug}{2});
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
    for dir = 1:2
        if dir == 2
            errorbar(ctr_x, nanmean(response_pn{drug}{dir}), nanstd(response_pn{drug}{dir})/sqrt(size(response_pn{drug}{dir}, 1)), 'color', color(dir));
        else
            errorbar(ctr_x, nanmean(response_pn{drug}{dir}), nanstd(response_pn{drug}{dir})/sqrt(size(response_pn{drug}{dir}, 1)), 'color', color(dir));
        end
        hold on
    end
    legend('superior', 'others', 'location', 'northwest')
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
    for dir = 1:2
        for cc = 1:size(response_pn{drug}{dir}, 1)
            plot(ctr_x, response_pn{drug}{dir}(cc, :), 'color', color(dir));
            hold on
        end
    end
    legend('superior', 'others', 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(condition{drug})
    xlim([3 400])
end


figure
for drug = 1:4
    subplot(2,2,drug)
    response_s{drug} = exciseRows_empty(response_pn{drug}{1});
    errorbar(ctr_x, nanmean(response_s{drug}), nanstd(response_s{drug})/sqrt(size(response_s{drug}, 1)), 'color', color(1));
    hold on
    response_others{drug} = exciseRows_empty(response_pn{drug}{2});
    errorbar(ctr_x, nanmean(response_others{drug}), nanstd(response_others{drug})/sqrt(size(response_others{drug}, 1)), 'color', color(2));
    legend('superior', 'others', 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(condition{drug})
    xlim([3 400])
end

%%
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'SR', 'SR+TPMPA', 'wash'};
step_size = 80;
for drug = 1:4
    for dir = 1:2
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
for drug = 1:4
    subplot(2,2,drug)
    response_s_p_ko_0724{drug} = exciseRows_empty(response_p{drug}{1});
    errorbar(ctr_x, nanmean(response_s_p_ko_0724{drug}), nanstd(response_s_p_ko_0724{drug})/sqrt(size(response_s_p_ko_0724{drug}, 1)), 'color', color(1));
    hold on
    response_others_p_ko_0724{drug} = exciseRows_empty(response_p{drug}{2});
    errorbar(ctr_x, nanmean(response_others_p_ko_0724{drug}), nanstd(response_others_p_ko_0724{drug})/sqrt(size(response_others_p_ko_0724{drug}, 1)), 'color', color(2));
    legend('superior', 'others', 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(condition{drug})
    xlim([3 400])
end


%%
load('DS170724.mat', 'response_s_0724', 'response_others_0724');
load('DS170802.mat', 'response_s_0802', 'response_others_0802');

for drug = 1:4
    response_s_ko{drug} = max([response_s_0724{drug}; response_s_0802{drug}], 0);
    response_others_ko{drug} = max([response_others_0724{drug}; response_others_0802{drug}], 0);
%     response_s_p_ko{drug} = max([response_s_p_ko_0724{drug}; response_s_p_ko_0802{drug}], 0);
%     response_others_p_ko{drug} = max([response_others_p_ko_0724{drug}; response_others_p_ko_0802{drug}], 0);
%     response_s_ko{drug} = abs(response_s_0724{drug});
%     response_others_ko{drug} = abs(response_others_0724{drug});

end

load('DS180413.mat', 'response_s', 'response_others')
response_s_0413 = response_s;
response_others_0413 = response_others;
load('DS170629.mat', 'response_s', 'response_others')
response_s_0629 = response_s;
response_others_0629 = response_others;
response_s_0629(3) = [];
response_others_0629(3) = [];
for drug = 1:3
    response_s_all{drug} = [response_s_0629{drug}; response_s_0413{drug}];
    response_others_all{drug} = [response_others_0629{drug}; response_others_0413{drug}];
end

figure
for drug = 1:4
    subplot(2,2,drug)
    errorbar(ctr_x, mean(response_s_ko{drug}), std(response_s_ko{drug})/sqrt(size(response_s_ko{drug}, 1)), 'color', color(1));
    hold on
    errorbar(ctr_x, mean(response_others_ko{drug}), std(response_others_ko{drug})/sqrt(size(response_others_ko{drug}, 1)), 'color', color(2));
    legend('superior', 'others', 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(condition{drug})
    xlim([3 400])
end

figure
errorbar(ctr_x, mean(response_s_ko{1}), std(response_s_ko{1})/sqrt(size(response_s_ko{1}, 1)), 'color', color(1));
hold on
errorbar(ctr_x, mean(response_s_ko{2}), std(response_s_ko{2})/sqrt(size(response_s_ko{2}, 1)), 'color', color(2));
errorbar(ctr_x, mean(response_others_ko{1}), std(response_others_ko{1})/sqrt(size(response_others_ko{1}, 1)), 'color', color(3));
errorbar(ctr_x, mean(response_others_ko{2}), std(response_others_ko{2})/sqrt(size(response_others_ko{2}, 1)), 'color', color(4));
legend('control superior', 'SR superior', 'control others', 'SR others', 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike rate')
title('KO')
xlim([3 400])

figure
temp = response_s_ko{1};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'k');
hold on
temp = response_others_ko{1};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'r');
temp = response_s_ko{2};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'k--');
temp = response_others_ko{2};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'r--');
legend('control superior', 'control others', 'SR superior', 'SR others', 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike rate')
title('KO')
xlim([3 100])


figure
temp = response_s_all{1};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'k');
hold on
temp = response_s_ko{1};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'r');
temp = response_s_all{2};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'k--');
temp = response_s_ko{2};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'r--');
legend('WT control', 'KO control', 'WT SR', 'KO SR', 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike rate')
title('superior')
xlim([3 100])


figure
temp = response_others_all{1};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'k');
hold on
temp = response_others_ko{1};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'r');
temp = response_others_all{2};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'k--');
temp = response_others_ko{2};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'r--');
legend('WT control', 'KO control', 'WT SR', 'KO SR', 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike rate')
title('others')
xlim([3 100])


figure
temp = response_s_all{1};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'k');
hold on
temp = response_s_ko{1};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'k--');
temp = response_others_all{1};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'r');
temp = response_others_ko{1};
errorbar(ctr_x(1:5), mean(temp(:, 1:5)), std(temp(:, 1:5))/sqrt(size(temp(:, 1:5), 1)), 'r--');
legend('WT superior', 'KO superior', 'WT others', 'KO others', 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike rate')
xlim([4 100])


%% bar graph
% WT
l_ctr_i = 1; h_ctr_i = 5;
r_mean_wt(1,1) = mean(response_others{1}(:, l_ctr_i));
r_mean_wt(2,1) = mean(response_s{1}(:, l_ctr_i));
r_mean_wt(3,1) = mean(response_s{2}(:, l_ctr_i));
r_mean_wt(4,1) = mean(response_others{2}(:, l_ctr_i));

r_mean_wt(1,2) = mean(response_others{1}(:, h_ctr_i));
r_mean_wt(2,2) = mean(response_s{1}(:, h_ctr_i));
r_mean_wt(3,2) = mean(response_s{2}(:, h_ctr_i));
r_mean_wt(4,2) = mean(response_others{2}(:, h_ctr_i));

r_ste_wt(1,1) = std(response_others{1}(:, l_ctr_i))/sqrt(length(response_others{1}(:, l_ctr_i)));
r_ste_wt(2,1) = std(response_s{1}(:, l_ctr_i))/sqrt(length(response_s{1}(:, l_ctr_i)));
r_ste_wt(3,1) = std(response_s{2}(:, l_ctr_i))/sqrt(length(response_s{2}(:, l_ctr_i)));
r_ste_wt(4,1) = std(response_others{2}(:, l_ctr_i))/sqrt(length(response_others{2}(:, l_ctr_i)));

r_ste_wt(1,2) = std(response_others{1}(:, h_ctr_i))/sqrt(length(response_others{1}(:, h_ctr_i)));
r_ste_wt(2,2) = std(response_s{1}(:, h_ctr_i))/sqrt(length(response_s{1}(:, h_ctr_i)));
r_ste_wt(3,2) = std(response_s{2}(:, h_ctr_i))/sqrt(length(response_s{2}(:, h_ctr_i)));
r_ste_wt(4,2) = std(response_others{2}(:, h_ctr_i))/sqrt(length(response_others{2}(:, h_ctr_i)));

r_mean_wt(5, :) = r_mean_wt(1, :);
r_ste_wt(5, :) = r_ste_wt(1, :);
r_mean_wt(1, :) = [];
r_ste_wt(1, :) = [];

Ticks = {'low contrast', 'high contrast'};

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = Ticks;
model_series = r_mean_wt';
model_error = r_ste_wt';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('spike rate (Hz)')
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
% KO
l_ctr_i = 1; h_ctr_i = 5;
r_mean_ko(1,1) = mean(response_others_ko{1}(:, l_ctr_i));
r_mean_ko(2,1) = mean(response_s_ko{1}(:, l_ctr_i));
r_mean_ko(3,1) = mean(response_s_ko{2}(:, l_ctr_i));
r_mean_ko(4,1) = mean(response_others_ko{2}(:, l_ctr_i));

r_mean_ko(1,2) = mean(response_others_ko{1}(:, h_ctr_i));
r_mean_ko(2,2) = mean(response_s_ko{1}(:, h_ctr_i));
r_mean_ko(3,2) = mean(response_s_ko{2}(:, h_ctr_i));
r_mean_ko(4,2) = mean(response_others_ko{2}(:, h_ctr_i));

r_ste_ko(1,1) = std(response_others_ko{1}(:, l_ctr_i))/sqrt(length(response_others_ko{1}(:, l_ctr_i)));
r_ste_ko(2,1) = std(response_s_ko{1}(:, l_ctr_i))/sqrt(length(response_s_ko{1}(:, l_ctr_i)));
r_ste_ko(3,1) = std(response_s_ko{2}(:, l_ctr_i))/sqrt(length(response_s_ko{2}(:, l_ctr_i)));
r_ste_ko(4,1) = std(response_others_ko{2}(:, l_ctr_i))/sqrt(length(response_others_ko{2}(:, l_ctr_i)));

r_ste_ko(1,2) = std(response_others_ko{1}(:, h_ctr_i))/sqrt(length(response_others_ko{1}(:, h_ctr_i)));
r_ste_ko(2,2) = std(response_s_ko{1}(:, h_ctr_i))/sqrt(length(response_s_ko{1}(:, h_ctr_i)));
r_ste_ko(3,2) = std(response_s_ko{2}(:, h_ctr_i))/sqrt(length(response_s_ko{2}(:, h_ctr_i)));
r_ste_ko(4,2) = std(response_others_ko{2}(:, h_ctr_i))/sqrt(length(response_others_ko{2}(:, h_ctr_i)));

r_mean_ko(5, :) = r_mean_ko(1, :);
r_ste_ko(5, :) = r_ste_ko(1, :);
r_mean_ko(1, :) = [];
r_ste_ko(1, :) = [];
Ticks = {'low contrast', 'high contrast'};

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = Ticks;
model_series = r_mean_ko';
model_error = r_ste_ko';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('spike rate (Hz)')
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
Ticks = {'WT', 'KO'};

r_mean_l = [r_mean_wt(:, 1) r_mean_ko(:, 1)];
r_mean_h = [r_mean_wt(:, 2) r_mean_ko(:, 2)];
r_ste_l = [r_ste_wt(:, 1) r_ste_ko(:, 1)];
r_ste_h = [r_ste_wt(:, 2) r_ste_ko(:, 2)];

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = Ticks;
model_series = r_mean_l';
model_error = r_ste_l';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('spike rate (Hz)')
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

Ticks = {'WT', 'KO'};
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = Ticks;
model_series = r_mean_h';
model_error = r_ste_h';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('spike rate (Hz)')
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
ctr = 5;
figure
subplot(1, 2, 1)
bin_size = 0.2;
% trial_dur = trial_dur{1};
x = bin_size/2:bin_size:trial_dur;
plot(x,hist(raster_n_sum{1}{52}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{1}.stimulus.repetitions/bin_size, 'b')
hold on
plot(x,hist(raster_n_sum{1}{42}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{1}.stimulus.repetitions/bin_size, 'r')
legend('superior', 'anterior')
xlabel('time (second)')
ylabel('firing rate (Hz)')
title('control')

subplot(1, 2, 2)
bin_n = 10;
plot(x,hist(raster_n_sum{2}{52}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{2}.stimulus.repetitions/bin_size, 'b')
hold on
plot(x,hist(raster_n_sum{2}{42}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{2}.stimulus.repetitions/bin_size, 'r')
legend('superior', 'anterior')
xlabel('time (second)')
ylabel('firing rate (Hz)')
title('SR')


%% 
% trial_dur = trial_dur{1};
load('DS170629.mat', 'raster_n_sum', 'raster_n_sum_all')
ctr = 5;

idx1 = 45;
idx2 = 27;
figure
subplot(2, 3, 4)
plot_raster(squeeze(raster_n_sum_all{1}{idx1}(:, :, ctr, :)), 0, 3.5)
xlim([0 3.5])
xlabel('time (second)')

subplot(2, 3, 1)
bin_size = 0.025;
x = bin_size/2:bin_size:trial_dur;
temp = hist(raster_n_sum{1}{idx1}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{1}.stimulus.repetitions/bin_size;
filter = [0.25 0.5 0.25];
temp = conv(temp, filter, 'same');
plot(x,temp, 'k')
ylim([0 40])
title('control superior')
xlim([0 3.5])

subplot(2, 3, 6)
plot_raster(squeeze(raster_n_sum_all{1}{idx2}(:, :, ctr, :)), 0, 3.5)
xlim([0 3.5])
xlabel('time (second)')

subplot(2, 3, 3)
temp = hist(raster_n_sum{1}{idx2}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{1}.stimulus.repetitions/bin_size;
temp = conv(temp, filter, 'same');
plot(x,temp, 'k')
ylabel('firing rate (Hz)')
title('control')
xlim([0 3.5])
ylim([0 40])
title('control posterior')


load('DS170724.mat', 'raster_n_sum', 'raster_n_sum_all')
idx3 = 52;
subplot(2, 3, 5)
plot_raster(squeeze(raster_n_sum_all{1}{idx3}(:, :, ctr, :)), 0, 3.5)
xlim([0 3.5])
xlabel('time (second)')

subplot(2, 3, 2)
temp = hist(raster_n_sum{1}{idx3}{ctr}, bin_size/2:bin_size:trial_dur)/datamb{1}.stimulus.repetitions/bin_size;
temp = conv(temp, filter, 'same');
plot(x,temp, 'k')
ylabel('firing rate (Hz)')
title('KO superior')
xlim([0 3.5])
ylim([0 40])


%% preferred direction CRF
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'SR', 'SR+TPMPA', 'wash'};
step_size = 80;
for drug = 1:4
    for dir = 1:2
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
                    response_pn{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
%                     response_pn{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions;
                end
                CC = CC + 1;
            end
        end
        response_pn_norm{drug}{dir} = response_pn{drug}{dir}./repmat(max(response_pn{drug}{dir}, [], 2), 1, size(response_pn{drug}{dir},2));
    end
end

figure
for drug = 1:4
    subplot(2,2,drug)
    response_s{drug} = max(exciseRows_empty(response_pn{drug}{1}), 0);
    errorbar(ctr_x, nanmean(response_s{drug}), nanstd(response_s{drug})/sqrt(size(response_s{drug}, 1)), 'color', color(1));
    hold on
    response_others{drug} = max(exciseRows_empty(response_pn{drug}{2}), 0);
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
title('FACx')
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
title('FACx')
xlim([3 100])




%% direction tuning curves

% n = 4;
% for i = 1:n
%     duration(i) = datamb{i}.triggers(end);
% end
% bin_size = 0.08; %sec
% [raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx, raster_mb_all] = deal(cell(n, 1));
% for i = 1:n    
%     [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},ds_id,duration(i),bin_size);
%     MB{i} = sort_direction(mbcellanalysis(MaxRate, StimComb, datamb{i}));
%     raster_mb{i} = get_mb_raster(datamb{i}, ds_id, duration(i));
%     raster_mb_all{i} = combine_repeats(raster_mb{i});
%     for j = 1:length(raster_mb{i})
%         if(mb_idx(j))
%             raster_mb{i}{j} = [];
%             raster_mb_all{i}{j} = [];
%         end
%     end
%     trial_dur{i} = get_mb_trial_dur(datamb{i}, 400, 400, 0.5);
% end
% 
% ctr_p = 7; % choose which params to use to calculate prefer direction indices 
% MAG_all_norm_mb = cell(n, 1);
% 
% for i = 1:n
%     [raster_p_sum{i}, p_idx{i}, raster_p_sum_all{i}] = get_pdirection_raster(raster_mb{i}, MB{1}.angle{ctr_p});
%     [raster_n_sum{i}, n_idx{i}, raster_n_sum_all{i}] = get_ndirection_raster(raster_mb{i}, MB{1}.angle{ctr_p});
%     MAG_all_norm_mb{i} = normalize_MAG(MB{i});
%     rep = datamb{i}.stimulus.repetitions;
% end
% close all
% 
% all ds cells
ctr = {'5%', '10%', '20%', '40%', '80%', '150%', '300%'};
color = 'brgkcmy';
dirn = 2;
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
clear rho_mb_all RHO_mb_all rho_mb_mean rho_mb_ste dsi_mb_mean dsi_mb_ste rho_mb dsi_mb dsi_mb_mean_all dsi_mb_ste_all dsi_mb_mean_all_on dsi_mb_ste_all_on

for drug = 1:4
    for cl = 1:7
%         subplot(3, 7, (drug-1)*7+cl)
        rho_mb_all{drug}{cl} = [];
        RHO_mb_all{drug}{cl} = [];
        for i = 1:dirn
            rho_mb{drug}{i}{cl} = [];
            RHO_mb{drug}{i}{cl} = [];
            dsi_mb{drug}{i}{cl} = [];
            for cc = 1:length(idx_dir{i})
                if ~mb_idx(idx_dir{i}(cc))% && sum(MB{drug}.RHO{T, BW,cl}(idx_dir{i}(cc), :))>0
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
                rho_mb_all{drug}{cl} = [rho_mb_all{drug}{cl}; y_temp(seq)];
                RHO_mb_all{drug}{cl} = [RHO_mb_all{drug}{cl}; Y_temp(seq)];
                end
            end
            if ~isempty(rho_mb{drug}{i}{cl})
                rho_mb_mean{cl}{drug}(i, :) = mean(rho_mb{drug}{i}{cl}, 1);
                rho_mb_ste{cl}{drug}(i, :) = std(rho_mb{drug}{i}{cl}, [], 1)/sqrt(size(rho_mb{drug}{i}{cl}, 1));
                RHO_mb_mean{cl}{drug}(i, :) = mean(RHO_mb{drug}{i}{cl}, 1);
                RHO_mb_ste{cl}{drug}(i, :) = std(RHO_mb{drug}{i}{cl}, [], 1)/sqrt(size(RHO_mb{drug}{i}{cl}, 1));
                dsi_mb_mean{cl}{drug}(i) = mean(dsi_mb{drug}{i}{cl});
                dsi_mb_ste{cl}{drug}(i) = std(dsi_mb{drug}{i}{cl})/sqrt(length(dsi_mb{drug}{i}{cl}));
                rho_mb_all_mean{cl}{drug} = mean(rho_mb_all{drug}{cl});
                rho_mb_all_ste{cl}{drug} = std(rho_mb_all{drug}{cl}, [], 1)/sqrt(size(rho_mb_all{drug}{cl}, 1));
                RHO_mb_all_mean{cl}{drug} = mean(RHO_mb_all{drug}{cl});
                RHO_mb_all_ste{cl}{drug} = std(RHO_mb_all{drug}{cl}, [], 1)/sqrt(size(RHO_mb_all{drug}{cl}, 1));
            else
                rho_mb_mean{cl}{drug}(i, :) = nan;
                rho_mb_ste{cl}{drug}(i, :) = nan;
                RHO_mb_mean{cl}{drug}(i, :) = nan;
                RHO_mb_ste{cl}{drug}(i, :) = nan;
                dsi_mb_mean{cl}{drug}(i) = nan;
                dsi_mb_ste{cl}{drug}(i) = nan;
                rho_mb_all_mean{cl}{drug} = nan;
                rho_mb_all_ste{cl}{drug} = nan;
                RHO_mb_all_mean{cl}{drug} = nan;
                RHO_mb_all_ste{cl}{drug} = nan;
            end
        end
%         xlabel('direction (rad)')
%         ylabel('spike number')
%             title(ll{d})
%         xlim([-pi pi])
        dsi_mb_mean_all{cl} = cell2mat(dsi_mb_mean{cl}');
        dsi_mb_ste_all{cl} = cell2mat(dsi_mb_ste{cl}'); 
    end
end
dsi_mb_mean_all = reshape(cell2mat(dsi_mb_mean_all), size(dsi_mb_mean_all{1}, 1), size(dsi_mb_mean_all{1}, 2), length(dsi_mb_mean_all));
dsi_mb_ste_all = reshape(cell2mat(dsi_mb_ste_all), size(dsi_mb_ste_all{1}, 1), size(dsi_mb_ste_all{1}, 2), length(dsi_mb_ste_all));

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


% plot average (cell type)
% Drug = {'control', 'AP5', 'AP5+HEX', 'wash'};
Drug = {'control', 'SR', 'SR+TPMPA', 'wash'};

h = figure;
set(h, 'Position', [1 1 1520,1080])
for i = 1:2
    for cl = 1:7
        subplot(2, 7, (i-1)*7+cl)
        for drug = 1:4
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
        ylim([0 1])
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
% ctr = {'5%', '10%', '20%', '40%', '80%', '150%', '300%'};
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


%% sliding window direction tuning curve
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
                        response_max{drug}{dir}(CC, theta, ctr) = max_fr/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
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
        for drug = 1:4
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
    end
end
end
