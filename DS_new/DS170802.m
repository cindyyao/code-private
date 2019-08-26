cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

datadg = load_data('/Volumes/lab/analysis/2017-08-02-0/data008-sorted/data008-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2017-08-02-0/stimuli/s08.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

params_idx = [2 4]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);


datarun = load_data('/Volumes/lab/analysis/2017-08-02-0/datamb-map/datamb-map', opt);
time_points = [2100 3200 4300];
datamb(1:4) = split_datarun(datarun, time_points);
datamb{1}.names.stimulus_path = '/Volumes/lab/analysis/2017-08-02-0/stimuli/s00.txt';
datamb{1} = load_stim(datamb{1}, 'user_defined_trigger_set', [1:2:1120]);
datamb{1}.stimulus.triggers = datamb{1}.stimulus.triggers';
datamb{2}.names.stimulus_path = '/Volumes/lab/analysis/2017-08-02-0/stimuli/s02.txt';
datamb{2} = load_stim(datamb{2}, 'user_defined_trigger_set', [1:2:560]);
datamb{2}.stimulus.triggers = datamb{2}.stimulus.triggers';
datamb{3}.names.stimulus_path = '/Volumes/lab/analysis/2017-08-02-0/stimuli/s04.txt';
datamb{3} = load_stim(datamb{3}, 'user_defined_trigger_set', [1:2:560]);
datamb{3}.stimulus.triggers = datamb{3}.stimulus.triggers';
datamb{4}.names.stimulus_path = '/Volumes/lab/analysis/2017-08-02-0/stimuli/s06.txt';
datamb{4} = load_stim(datamb{4}, 'user_defined_trigger_set', [1:2:1120]);
datamb{4}.stimulus.triggers = datamb{4}.stimulus.triggers';
load DS170802.mat
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

%% plot cell summary
for dir = 1:4
    for cc = 1:length(id_dir{dir})
        plot_mb_raster_ctr(MB, raster_mb, trial_dur, idx_dir{dir}(cc), id_dir{dir}(cc), '', 4, 7, 1)
    end
end

%% get spontaneous activity

for drug = 1:4
    for dir = 1:2
        for cc = 1:length(id_dir_mb{dir})
            if(~mb_idx(idx_dir_mb{dir}(cc)))
                idx = get_cell_indices(datamb{drug},id_dir_mb{dir}(cc));
                spikes_temp = datamb{drug}.spikes{idx};
                bgnd_firing{dir}(drug, cc) = length(spikes_temp(spikes_temp > datamb{drug}.duration - 60 & spikes_temp < datamb{drug}.duration))/60;
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
% xlim([0 4])
% ylim([0 4])
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
    for dir = 1:4
        response_s{drug} = exciseRows_empty(response_pn_norm{drug}{1});
        errorbar(ctr_x, nanmean(response_s{drug}), nanstd(response_s{drug})/sqrt(size(response_s{drug}, 1)), 'color', color(1));
        hold on
        response_others{drug} = exciseRows_empty(response_pn_norm{drug}{2});
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
    for dir = 1:2
        errorbar(ctr_x, nanmean(response_pn{drug}{dir}), nanstd(response_pn{drug}{dir})/sqrt(size(response_pn{drug}{dir}, 1)), 'color', color(dir));
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
    response_s_p_ko_0802{drug} = exciseRows_empty(response_p{drug}{1});
    errorbar(ctr_x, nanmean(response_s_p_ko_0802{drug}), nanstd(response_s_p_ko_0802{drug})/sqrt(size(response_s_p_ko_0802{drug}, 1)), 'color', color(1));
    hold on
    response_others_p_ko_0802{drug} = exciseRows_empty(response_p{drug}{2});
    errorbar(ctr_x, nanmean(response_others_p_ko_0802{drug}), nanstd(response_others_p_ko_0802{drug})/sqrt(size(response_others_p_ko_0802{drug}, 1)), 'color', color(2));
    legend('superior', 'others', 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(condition{drug})
    xlim([3 400])
end

