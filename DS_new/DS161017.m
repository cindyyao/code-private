%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load drifting grating data
datadg = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-10-17-0/data002/data002', opt);
datadg.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-10-17-0/stimuli/s02.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);
datadg = load_ei(datadg, id_dir{1});
% identify DS cells
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [3 4]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

% load moving bar data
datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-10-17-0/data003-007-map/data003-007-map', opt);
time_points = [1900 3800 5700];
datamb(1:4) = split_datarun(datarun, time_points);

datamb{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-10-17-0/stimuli/s03.txt';
datamb{1} = load_stim(datamb{1}, 'user_defined_trigger_set', [1:2:1120]);
datamb{1}.stimulus.triggers = datamb{1}.stimulus.triggers';
datamb{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-10-17-0/stimuli/s05.txt';
datamb{2} = load_stim(datamb{2}, 'user_defined_trigger_set', [1:2:1120]);
datamb{2}.stimulus.triggers = datamb{2}.stimulus.triggers';
datamb{3}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-10-17-0/stimuli/s06.txt';
datamb{3} = load_stim(datamb{3}, 'user_defined_trigger_set', [1:2:1120]);
datamb{3}.stimulus.triggers = datamb{3}.stimulus.triggers';
datamb{4}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-10-17-0/stimuli/s07.txt';
datamb{4} = load_stim(datamb{4}, 'user_defined_trigger_set', [1:2:1120]);
datamb{4}.stimulus.triggers = datamb{4}.stimulus.triggers';

load('DS161017.mat')
% load white noise data
datawn = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-10-17-0/data001-map/data001-map', opt);
% datawn = load_sta(datawn);
datawn = load_ei(datawn, cell2mat(id_dir));

% id_ot = get_cell_ids(datawn, 'ON type 1');
% id_ot_mb = intersect(id_ot, datamb{1}.cell_ids);
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load drifting grating data
datadg = load_data('/Volumes/janacek/Analysis/2016-10-17-0/data002/data002', opt);
datadg.names.stimulus_path = '/Volumes/janacek/Analysis/2016-10-17-0/stimuli/s02.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

% identify DS cells
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [3 4]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

% load moving bar data
datarun = load_data('/Volumes/janacek/Analysis/2016-10-17-0/data003-007-map/data003-007-map', opt);
time_points = [1900 3800 5700];
datamb(1:4) = split_datarun(datarun, time_points);

datamb{1}.names.stimulus_path = '/Volumes/janacek/Analysis/2016-10-17-0/stimuli/s03.txt';
datamb{1} = load_stim(datamb{1}, 'user_defined_trigger_set', [1:2:1120]);
datamb{1}.stimulus.triggers = datamb{1}.stimulus.triggers';
datamb{2}.names.stimulus_path = '/Volumes/janacek/Analysis/2016-10-17-0/stimuli/s05.txt';
datamb{2} = load_stim(datamb{2}, 'user_defined_trigger_set', [1:2:1120]);
datamb{2}.stimulus.triggers = datamb{2}.stimulus.triggers';
datamb{3}.names.stimulus_path = '/Volumes/janacek/Analysis/2016-10-17-0/stimuli/s06.txt';
datamb{3} = load_stim(datamb{3}, 'user_defined_trigger_set', [1:2:1120]);
datamb{3}.stimulus.triggers = datamb{3}.stimulus.triggers';
datamb{4}.names.stimulus_path = '/Volumes/janacek/Analysis/2016-10-17-0/stimuli/s07.txt';
datamb{4} = load_stim(datamb{4}, 'user_defined_trigger_set', [1:2:1120]);
datamb{4}.stimulus.triggers = datamb{4}.stimulus.triggers';

% load white noise data
datawn = load_data('/Volumes/janacek/Analysis/2016-10-17-0/data001-map/data001-map', opt);
datawn = load_sta(datawn);

load('DS161017.mat')
% id_ot = get_cell_ids(datawn, 'ON type 1');
% id_ot_mb = intersect(id_ot, datamb{1}.cell_ids);
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

%% dg
% ds_id = datadg.cell_ids;
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);

[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,nonds_id,0,1);
DG_nds{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg_nds{i} = get_ds_raster(datadg, nonds_id);

ll = 1; t = 3;
ds_vs = DG{ll}.mag{t};
nds_vs = DG_nds{ll}.mag{t};
max_vs = max([ds_vs nds_vs]);
min_vs = min([ds_vs nds_vs]);

X = 0.05:0.1:2.15;
ds_vs_hist = hist(ds_vs, X);
nds_vs_hist = hist(nds_vs, X);
figure
bar(X, [ds_vs_hist' nds_vs_hist'], 1, 'stacked')
xlim([0 2.2])

delta_p = 4; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(DG{i});
rep = datadg.stimulus.repetitions;
%% plot individual cells
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [3]);
for cc = 2:5 %length(ds_id)
    plot_ds_raster(DG_cut, raster_dg_cut, cc, ds_id(cc), '', 1, 1, 0)
end


ll = 1; t = 3; cc = 4;
[idx, xx, yy] = subplot_idx(1, 1);
tt = DG{1}.theta{1}(1, :);
FigHandle = figure;
set(FigHandle, 'Position', [1 1 400 400])
h = subplot(xx, yy, idx(1)); 
u_temp = DG{ll}.U{t}(cc);
v_temp = DG{ll}.V{t}(cc);
P = polar(0, 1);
set(P, 'Visible', 'off')
hold on
polar(tt, DG{ll}.rho{t}(cc, :), 'b');
polar_theta_off(h)
for i = 2:9
    subplot(xx, yy, idx(i)); plot_raster(squeeze(raster_dg{ll}{cc}(1, t, i-1, :)), 0, 8);
    if idx(i) == 7
        ylabel('trials')
        xlabel('time (s)')
    end
end

%%
ll = 1; t = 3; cc = 8;
[idx, xx, yy] = subplot_idx(1, 1);
tt = DG_nds{1}.theta{1}(1, :);
FigHandle = figure;
set(FigHandle, 'Position', [1 1 400 400])
h = subplot(xx, yy, idx(1)); 
u_temp = DG_nds{ll}.U{t}(cc);
v_temp = DG_nds{ll}.V{t}(cc);
P = polar(0, 1);
set(P, 'Visible', 'off')
hold on
polar(tt, DG_nds{ll}.rho{t}(cc, :), 'b');
polar_theta_off(h)
for i = 2:9
    subplot(xx, yy, idx(i)); plot_raster(squeeze(raster_dg_nds{ll}{cc}(1, t, i-1, :)), 0, 8);
    if idx(i) == 7
        ylabel('trials')
        xlabel('time (s)')
    end
end

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
mag_pca = MAG_all_norm_dg{L}(2:9,:);
mag_pca = mag_pca';
[id_sub, idx_sub] = deal(cell(2, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(mag_pca);
pc1 = 1; pc2 = 3;
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

t = 4;
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
for dir = 2:2
    for cc = 1:1%length(id_dir{dir})
        plot_mb_raster_ctr(MB, raster_mb, trial_dur, idx_dir{dir}(cc), id_dir{dir}(cc), '', 4, 7, 0)
    end
end

for dir = 1:3
    for cc = 1:length(id_dir_on{dir})
        plot_mb_raster_ctr(MB, raster_mb, trial_dur, idx_dir_on{dir}(cc), id_dir_on{dir}(cc), '', 4, 7, 1)
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
%                     if drug == 1 && ctr == 7
                        [max_p, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
                        Max_i{dir}(CC) = max_i;
%                     else
%                         max_p = conv(hist_temp, ones(1,step_size), 'valid');
%                         max_p = max_p(Max_i{dir}(CC));
%                     end
                    response_pmax{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
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
        figure
        set(gcf, 'Position', [1 1 900 800])
        for cc = 1:size(response_pmax{drug}{dir}, 1)

            ydata = response_pmax{drug}{dir}(cc, :);
            xdata = log10(ctr_x);
            [f, G] = fit_nr(xdata, ydata, 'upper', [100, 100, log10(300), 0]);
            fit_all{drug}{dir}{cc} = f;
            G_all{drug}{dir}{cc} = G;

            x = linspace(min(xdata), max(xdata), 100);
            y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;

            subplot(5,6,cc)
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

fitting = struct();
fitting.sigma = sigma{1};
fitting.ymax = ymax{1};
fitting.aa = aa{1};
fitting.bb = bb{1};
%%

for dir = 1:4
    for drug = 1:4
        for cc = 1:size(response_pmax_norm{drug}{dir}, 1)
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

figure
set(gcf, 'Position', [1 1 900 800])
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4
        ydata = mean(response_pmax_norm{drug}{dir});
        yste = std(response_pmax_norm{drug}{dir})/sqrt(size(response_pmax_norm{drug}{dir}, 1));
        xdata = log10(ctr_x);
        [f, G] = fit_nr(xdata, ydata, 'upper', [2, 100, log10(300), min(ydata)]);
        sigma_avg(dir, drug) = f.sigma;
        x = linspace(min(xdata), max(xdata)+0.5, 100);
        y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;
        h(drug) = plot(10.^x,y,'color', color(drug));
        hold on
        errorbar(10.^xdata, ydata, yste, 'color', color(drug), 'marker', 'o', 'linestyle', 'none')
    end
    legend([h(1), h(2), h(3), h(4)], condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
end

% ON DSGC (combine direction)
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'AP5+Hex', 'wash'};
step_size = 60;
for drug = 1:4
    CC = 1;
    for cc = 1:length(id_sub{1})
        if ~isempty(raster_p_sum{drug}{idx_sub{1}(cc)})
            for ctr = 7:-1:1
                a = raster_p_sum{drug}{idx_sub{1}(cc)}{ctr};
                hist_temp = hist(a, xx);
%                 if drug == 1 && ctr == 7
                    [max_p, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
%                     Max_i{dir}(CC) = max_i;
%                 else
%                     max_p = conv(hist_temp, ones(1,step_size), 'valid');
%                     max_p = max_p(Max_i{dir}(CC));
%                 end
                response_pmax_on{drug}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions - bgnd_firing_on(drug, CC)*bin_size*step_size;
            end
            CC = CC + 1;
        end
    end
    response_pmax_on_norm{drug} = response_pmax_on{drug}./repmat(max(response_pmax_on{1}, [], 2), 1, size(response_pmax_on{drug},2));
end


ctr_x = [5 10 20 40 80 150 300];
color = 'brgkc';
figure
set(gcf, 'Position', [1 1 900 800])

for drug = 1:4
    errorbar(ctr_x, mean(response_pmax_on_norm{drug}, 1), std(response_pmax_on_norm{drug}, [], 1)/sqrt(size(response_pmax_on_norm{drug}, 1)), 'color', color(drug));
    hold on
end
legend(condition, 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike rate')

% fit ds
clear fit_all G_all
for drug = 1:4
    figure
    set(gcf, 'Position', [1 1 900 800])
    for cc = 1:size(response_pmax_on_norm{drug}, 1)

        ydata = response_pmax_on_norm{drug}(cc, :);
        xdata = log10(ctr_x);
        [f, G] = fit_nr(xdata, ydata, 'upper', [2, 100, log10(300), 0]);
        fit_all{drug}{cc} = f;
        G_all{drug}{cc} = G;

        x = linspace(min(xdata), max(xdata), 100);
        y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;

        subplot(5,6,cc)
        plot(xdata, ydata)
        hold on
        plot(x, y)

        sigma{6}(drug, cc) = f.sigma;
        ymax{6}(drug, cc) = f.ymax;
        aa{6}(drug, cc) = f.a;
        bb{6}(drug, cc) = f.b;
    end
end


for dir = 1:3
    for drug = 1:4
        for cc = 1:size(response_pmax_on_norm{drug}{dir}, 1)
            x = linspace(sigma{dir}(1,cc)-0.9, sigma{dir}(1, cc)+1.1, 10);
            y = ymax{dir}(drug, cc)*x.^aa{dir}(drug, cc)./(x.^aa{dir}(drug, cc) + sigma{dir}(drug, cc)^aa{dir}(drug, cc))+bb{dir}(drug, cc);
            pmax_y_on{drug}{dir}(cc, :) = y;
%             pause
        end
    end
end
        
figure
set(gcf, 'Position', [1 1 900 800])
for dir = 1:3
    subplot(2,2,dir)
    for drug = 1:4 
        errorbar(linspace(0.1,2.1,10), mean(pmax_y_on{drug}{dir}, 1), std(pmax_y_on{drug}{dir}, [], 1)/sqrt(size(pmax_y_on{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
%     set(gca, 'Xscale', 'log')
    xlabel('normalized contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
%     xlim([3 400])
end

% color = {[0,0,0]/255,[0,191,255]/255, [255,20,60]/255, [0.5 0.5 0.5]};
figure
set(gcf, 'Position', [1 1 900 800])
for drug = 1:4
    ydata = mean(response_pmax_on_norm{drug}, 1);
    yste = std(response_pmax_on_norm{drug}, [], 1)/sqrt(size(response_pmax_on_norm{drug}, 1));
    xdata = log10(ctr_x);
    [f, G] = fit_nr(xdata, ydata, 'upper', [2, 100, log10(300), min(ydata)]);

    x = linspace(min(xdata), max(xdata)+0.5, 100);
    y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;
    h(drug) = plot(10.^x,y,'color', color(drug));
    hold on
    errorbar(10.^xdata, ydata, yste, 'color', color(drug), 'marker', 'o', 'linestyle', 'none')
end
legend([h(1), h(2), h(3), h(4)], condition, 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('contrast')
ylabel('spike rate')
title('ON DSGC')

%% figure
color = {[0,0,0]/255,[0,191,255]/255, [255,20,60]/255, [0.5 0.5 0.5]};

figure(2)
% set(gcf, 'Position', [1 1 900 800])
for dir = 1:4
    subplot(2,4,dir)
    for drug = 1:3
        ydata = mean(response_pmax_norm{drug}{dir});
        yste = std(response_pmax_norm{drug}{dir})/sqrt(size(response_pmax_norm{drug}{dir}, 1));
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
    xlim([3 1000])
    ylim([-0.1 1.2])

end
subplot(2,4,5)
for drug = 1:3
    ydata = mean(response_pmax_on_norm{drug}, 1);
    yste = std(response_pmax_on_norm{drug}, [], 1)/sqrt(size(response_pmax_on_norm{drug}, 1));
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
title('ON DSGC')
xlim([3 1000])
ylim([-0.1 1.2])

subplot(2,4,6)
for drug = 1:3
    ydata = mean(response_pmax_ot_norm{drug}, 1);
    yste = std(response_pmax_ot_norm{drug}, [], 1)/sqrt(size(response_pmax_ot_norm{drug}, 1));
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
title('ON transient')
xlim([3 1000])
ylim([-0.1 1.2])

%% mb
n = 4;
bin_size = 0.025; %sec

% id_ot_mb = get_cell_ids(datawn, 'OFF type2');
% id_ot_mb = id_ot_mb([1 3 4 5 10 12 15]);
clear MB_ot raster_mb_ot raster_p_sum_ot
% 
[raster_mb_ot, MB_ot, trial_dur, raster_p_sum_ot, p_idx_ot] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},id_ot_mb,duration(i),bin_size);
    MB_ot{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb_ot{i} = get_mb_raster(datamb{i}, id_ot_mb, duration(i));
    trial_dur{i} = get_mb_trial_dur(datamb{i}, 400, 400, 0.5);
end

ctr_p = 7; % choose which params to use to calculate prefer direction indices 

for i = 1:n
    [raster_p_sum_ot{i}, p_idx_ot{i}, raster_p_sum_all{i}] = get_pdirection_raster(raster_mb_ot{i}, MB_ot{1}.angle{ctr_p});
    rep = datamb{i}.stimulus.repetitions;
end
close all
%% plot cell summary
drug = 1;
for cc = 1:length(id_ot_mb)
    plot_mb_raster_ctr(MB_ot(drug), raster_mb_ot(drug), trial_dur, cc, id_ot_mb(cc), '', 1, 7, 0)
    pause
    close all
end

%% ON transient cell CRF
clear Max_i response_pmax_ot_norm response_pmax_ot
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
condition = {'control', 'AP5', 'AP5+Hex', 'wash'};
step_size = 60;
for drug = 1:4
    for cc = 1:length(id_ot_mb)
        for ctr = 7:-1:1
            a = raster_p_sum_ot{drug}{cc}{ctr};
            hist_temp = hist(a, xx);
%             if drug == 1 && ctr == 7
                [max_p, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
%                 Max_i(cc) = max_i;
%             else
%                 max_p = conv(hist_temp, ones(1,step_size), 'valid');
%                 max_p = max_p(Max_i(cc));
%             end
            response_pmax_ot{drug}(cc, ctr) = max_p/datamb{drug}.stimulus.repetitions - bgnd_firing_ot(drug, cc)*bin_size*step_size;
        end
    end
    response_pmax_ot_norm{drug} = response_pmax_ot{drug}./repmat(max(response_pmax_ot{1}, [], 2), 1, size(response_pmax_ot{drug},2));
end


ctr_x = [5 10 20 40 80 150 300];
color = 'brgkc';
figure
set(gcf, 'Position', [1 1 900 800])

for drug = 1:4
    errorbar(ctr_x, mean(response_pmax_ot_norm{drug}, 1), std(response_pmax_ot_norm{drug}, [], 1)/sqrt(size(response_pmax_ot_norm{drug}, 1)), 'color', color(drug));
    hold on
end
legend(condition, 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike rate')
title('ON transient')
xlim([3 400])

% fit ds
clear fit_all G_all
i = 5;
for drug = 1:4
    figure
    set(gcf, 'Position', [1 1 900 800])
    for cc = 1:size(response_pmax_ot_norm{drug}, 1)

        ydata = response_pmax_ot{drug}(cc, :);
        xdata = log10(ctr_x);
        [f, G] = fit_nr(xdata, ydata, 'upper', [100, 100, log10(300), 0]);
        fit_all{drug}{cc} = f;
        G_all{drug}{cc} = G;

        x = linspace(min(xdata), max(xdata), 100);
        y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;

        subplot(5,6,cc)
        plot(xdata, ydata)
        hold on
        plot(x, y)

        sigma{i}(drug, cc) = f.sigma;
        ymax{i}(drug, cc) = f.ymax;
        aa{i}(drug, cc) = f.a;
        bb{i}(drug, cc) = f.b;
    end
end
 

for dir = 1:6
    sigma_mean{dir} = mean(10.^sigma{dir}');
    sigma_ste{dir} = std(10.^sigma{dir}')/sqrt(size(sigma{dir}, 2));
end

sigma_oods = 10.^cell2mat(sigma(1:4))';
sigma_oods_mean = mean(sigma_oods);
sigma_oods_ste = std(sigma_oods)/sqrt(size(sigma_oods, 2));
%% sigma
cell_type = {'superior', 'anterior', 'inferior', 'posterior', 'ON DS', 'ON transient'};
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = cell_type;
model_series = cell2mat(sigma_mean');
model_error = cell2mat(sigma_ste');
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('C50')
title('NDF 0')
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

% 
% for dir = 1:3
%     for drug = 1:4
%         for cc = 1:size(response_pmax_on_norm{drug}{dir}, 1)
%             x = linspace(sigma{dir}(1,cc)-0.9, sigma{dir}(1, cc)+1.1, 10);
%             y = ymax{dir}(drug, cc)*x.^aa{dir}(drug, cc)./(x.^aa{dir}(drug, cc) + sigma{dir}(drug, cc)^aa{dir}(drug, cc))+bb{dir}(drug, cc);
%             pmax_y_on{drug}{dir}(cc, :) = y;
% %             pause
%         end
%     end
% end
%         
% figure
% set(gcf, 'Position', [1 1 900 800])
% for dir = 1:3
%     subplot(2,2,dir)
%     for drug = 1:4 
%         errorbar(linspace(0.1,2.1,10), mean(pmax_y_on{drug}{dir}, 1), std(pmax_y_on{drug}{dir}, [], 1)/sqrt(size(pmax_y_on{drug}{dir}, 1)), 'color', color(drug));
%         hold on
%     end
%     legend(condition, 'location', 'northwest')
% %     set(gca, 'Xscale', 'log')
%     xlabel('normalized contrast')
%     ylabel('spike rate')
%     title(dscell_type{dir})
% %     xlim([3 400])
% end

% color = {[0,0,0]/255,[0,191,255]/255, [255,20,60]/255, [0.5 0.5 0.5]};
figure
set(gcf, 'Position', [1 1 900 800])
for drug = 1:4
    ydata = mean(response_pmax_ot_norm{drug}, 1);
    yste = std(response_pmax_ot_norm{drug}, [], 1)/sqrt(size(response_pmax_ot_norm{drug}, 1));
    xdata = log10(ctr_x);
    [f, G] = fit_nr(xdata, ydata, 'upper', [2, 100, log10(300), min(ydata)]);

    x = linspace(min(xdata), max(xdata)+0.5, 100);
    y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;
    h(drug) = plot(10.^x,y,'color', color(drug));
    hold on
    errorbar(10.^xdata, ydata, yste, 'color', color(drug), 'marker', 'o', 'linestyle', 'none')
end
legend([h(1), h(2), h(3), h(4)], condition, 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('contrast')
ylabel('spike rate')
title('ON transient')

%% c20

for dir = 1:6
    for drug = 1:4
        for cc = 1:size(sigma{dir},2)
            c20 = sigma{dir}(1,cc)/1^(1/aa{dir}(drug,cc));
            c20r{dir}(drug, cc) = ymax{dir}(drug,cc)*c20^aa{dir}(drug,cc)/(c20^aa{dir}(drug,cc)+sigma{dir}(drug,cc)^aa{dir}(drug,cc))+bb{dir}(drug,cc);
        end
    end
    c20r_mean{dir} = mean(c20r{dir}, 2);
    c20r_ste{dir} = std(c20r{dir}, [], 2)/sqrt(size(sigma{dir},2));
end

cell_type = {'superior', 'anterior', 'inferior', 'posterior', 'ON DS', 'ON transient'};
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = cell_type;
model_series = cell2mat(c20r_mean)';
model_error = cell2mat(c20r_ste)';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('c20r')
title('NDF 0')
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
load('DS160904_test.mat', 'c20r_combined_ndf3_mean', 'c20r_combined_ndf3_ste')
c20r_combined{1} = cell2mat(c20r(1:5));
c20r_combined{2} = c20r{6};
for ct = 1:2
    c20r_combined_mean(:, ct) = mean(c20r_combined{ct}, 2);
    c20r_combined_ste(:, ct) = std(c20r_combined{ct}, [], 2)/sqrt(size(c20r_combined{ct}, 2));
end

cell_type = {'DSGC scotopic', 'DSGC photopic', 'ON transient photopic'};
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = cell_type;
model_series = [c20r_combined_ndf3_mean c20r_combined_mean]';
model_error = [c20r_combined_ndf3_ste c20r_combined_ste]';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('c50r')
% title('NDF 0')
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
%% combine anterior, inferior, posterior
ct = {'superior', 'other'};
id_ct{1} = id_dir{1};
id_ct{2} = cell2mat(id_dir(2:4));
clear spike_count_var_p spike_count_mean_p spike_timing_var_p
for i = 1:4
    for dir = 1:2
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
        fitcoeff_count{i}{dir} = polyfit(spike_count_mean_p{i}{dir}, spike_count_var_p{i}{dir}, 1);
    end
end


figure
for dir = 1:2
    subplot(1,2,dir)
    for drug = 1:3
        h(drug) = plot(spike_count_mean_p{drug}{dir}, spike_count_var_p{drug}{dir}, [color(drug) 'o']);
        hold on
        x = linspace(min(spike_count_mean_p{drug}{dir})*0.8, max(spike_count_mean_p{drug}{dir})*1.2, 100);
        xx = [x; ones(1,100)];
        yy = fitcoeff_count{drug}{dir}*xx;
        plot(xx(1,:), yy, color(drug))

    end
    xlabel('spike count mean')
    ylabel('spike count variance')
    legend([h(1), h(2), h(3)], 'control', 'AP5', 'AP5+HEX')
    title(ct{dir})
end

% 1st spike timing variance
for i = 1:4
    for dir = 1:2
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
        fitcoeff_timing{i}{dir} = polyfit(spike_count_mean_p{i}{dir}, log10(spike_timing_var_p{i}{dir}), 1);
    end
end


figure
for dir = 1:2
    subplot(1,2,dir)
    for drug = 1:3
        h(drug) = plot(spike_count_mean_p{drug}{dir}, log10(spike_timing_var_p{drug}{dir}), [color(drug) 'o']);
        hold on
        x = linspace(min(spike_count_mean_p{drug}{dir})*0.8, max(spike_count_mean_p{drug}{dir})*1.2, 100);
        xx = [x; ones(1,100)];
        yy = fitcoeff_timing{drug}{dir}*xx;
        plot(xx(1,:), yy, color(drug))
    end
    xlabel('spike count mean')
    ylabel('log(first spike timing variance)')

%     set(gca, 'YScale', 'log')
    legend([h(1), h(2), h(3)], 'control', 'AP5', 'AP5+HEX')
    title(ct{dir})
end

%% figure
% example cells:
% Superior: 618 Anterior: 1066 Inferior: 347 Posterior: 6183
color = {[0,0,0]/255,[0,191,255]/255, [255,20,60]/255, [0.5 0.5 0.5]};
cell_id = [3634 6798 4983 4054];
ctr_i = [3 4 4 4];
% direction = [7 8 3 5];
% bin_size = 0.1; % second
% XX = bin_size/2:bin_size:3-bin_size/2; 

figure;
for cc = 1:length(cell_id)
    cell_idx = find(ds_id == cell_id(cc));
    for drug = 1:4
        h = subplot(4,4,4*(drug-1)+cc);
        plot_raster(squeeze(raster_p_sum_all{drug}{cell_idx}(:,:, ctr_i(cc),:)), 0,3, 'color', color{drug})
        set(h, 'xtick', [])
        set(h, 'ytick', [])
        if drug == 1
            title(dscell_type{cc})
        end
    end
end


%% direction tuning curves
% all ds cells
ctr = {'10%', '20%', '40%', '80%', '150%', '300%', '400%'};
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
clear rho_mb_all RHO_mb_all rho_mb_mean rho_mb_ste dsi_mb_mean dsi_mb_ste rho_mb dsi_mb dsi_mb_mean_all dsi_mb_ste_all dsi_mb_mean_all_on dsi_mb_ste_all_on

for drug = 1:3
    for cl = 1:7
%         subplot(3, 7, (drug-1)*7+cl)
        rho_mb_all{drug}{cl} = [];
        RHO_mb_all{drug}{cl} = [];
        for i = 1:4
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
for drug = 1:3
    for cl = 1:7
        subplot(3, 7, (drug-1)*7+cl)
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
for drug = 1:3
    for cl = 1:7
        subplot(3, 7, (drug-1)*7+cl)
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
        for drug = 1:3
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

%% separate ON and OFF responses for CRF
ctr = 7;
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
for drug = 1:1
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
%             [x,~] = ginput;
%             BreakIndices{drug}(cc, :) = round(x'/bin_size);
%             close(1)
%         else
%             BreakIndices{drug}(cc, :) = nan(1, ctr);

pause
close(1)
        
        end
    end
end

for dir = 1:1
    for cc = 1:length(id_dir{dir})
        if ~mb_idx(idx_dir{dir}(cc))
            h = figure(1);
            set(h, 'position', [1 1 1000 1000])
            for i = 7:-1:1
                subplot(7,1,i)
                temp = hist(raster_p_sum{drug}{idx_dir{dir}(cc)}{ctr-7+i}, xx);
                plot(xx, temp)
                hold on
            end
%             [x,~] = ginput;
%             BreakIndices(idx_dir{dir}(cc), :) = round(x'/bin_size);
%             close(1)
%         else
%             BreakIndices(idx_dir{dir}(cc), :) = nan(1, ctr);
pause
close(1)
        end
    end
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
            clear raster_temp
            for ctr = 1:7
                raster_temp{ctr} = raster_p_sum{drug}{cc}{ctr};
                raster_temp{ctr}(raster_temp{ctr} > BreakIndices(cc, ctr)*bin_size) = [];
            end
            pd_spikes{drug}(cc,:) = cellfun('length', raster_temp)/datamb{drug}.stimulus.repetitions;
        end
    end
    
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                clear raster_temp
                for ctr = 1:7
                    raster_temp{ctr} = raster_p_sum{drug}{idx_dir{dir}(cc)}{ctr};
                    raster_temp{ctr}(raster_temp{ctr} > BreakIndices(idx_dir{dir}(cc), ctr)*bin_size) = [];
                end

                spike_temp = cellfun('length', raster_temp)/datamb{drug}.stimulus.repetitions;
                if sum(spike_temp) ~= 0
                    pd_dir_spikes{drug}{dir}(CC,:) = spike_temp;
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

%% max window
% on-off DSGC
clear Max_i
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'AP5+Hex', 'wash'};
step_size = 30;
for drug = 1:4
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                for ctr = 7:-1:1
                    a = raster_p_sum{drug}{idx_dir{dir}(cc)}{ctr};
                    hist_temp = hist(a, xx);
                    hist_temp(BreakIndices(idx_dir{dir}(cc), ctr):end) = zeros(1, length(hist_temp)-BreakIndices(idx_dir{dir}(cc), ctr)+1);
%                     if drug == 1 && ctr == 7
                        [max_p, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
%                         Max_i{dir}(CC) = max_i;
%                     else
%                         max_p = conv(hist_temp, ones(1,step_size), 'valid');
%                         max_p = max_p(Max_i{dir}(CC));
%                     end
                    response_pmax{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
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
%                     response_pmax{drug}{dir}(CC, ctr) = sum(OnRes)/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
                    response_pmax{drug}{dir}(CC, ctr) = sum(OnRes)/datamb{drug}.stimulus.repetitions;
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
    subplot(1,4,dir)
    for drug = 1:4
        errorbar(ctr_x, nanmean(response_pmax_norm{drug}{dir}), nanstd(response_pmax_norm{drug}{dir})/sqrt(size(response_pmax_norm{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
    xlim([3 400])
end

for dir = 1:4
    ratio = abs(response_pmax{2}{dir}./response_pmax{1}{dir});
    ratio(isnan(ratio)) = 0;
    ratio(isinf(ratio)) = 0;
    
    RatioCRF{dir} = ratio;
    RatioCRFMean(dir, :) = robust_mean(ratio);
    RatioCRFSte(dir,:) = robust_std(ratio)/sqrt(size(ratio, 1));
    
end

figure
for dir = 1:4
    subplot(2,2,dir)
    plot(RatioCRF{dir}')
end

figure
for dir = 1:4
    errorbar(ctr_x, RatioCRFMean(dir, :), RatioCRFSte(dir, :));
    hold on
end
set(gca, 'xscale', 'log')
ylabel('AP5/control')
xlabel('contrast')
legend(dscell_type)

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

%% paper figure
% A
clear MbHist
color = {[0,0,0]/255,[255,20,60]/255,[0,191,255]/255,  [0.5 0.5 0.5]};
cell_id = 438;
cell_idx = find(ds_id == cell_id);
ctr_i = 5;
bin_size = 0.03;
xx = bin_size/2:bin_size:trial_dur{1}-bin_size/2;

for drug = 1:4
    for dir = 1:8
        MbHist{drug}(dir, :) = hist(cell2mat(squeeze(raster_mb{drug}{cell_idx}(1,1,dir,ctr_i,:))), xx)/datamb{drug}.stimulus.repetitions/bin_size;
    end
end

subplot_i = [6 3 2 1 4 7 8 9];
figure
subplot(3,3,5)
polar(MB{1}.theta{1}(1, :), MB{drug}.rho{ctr_i}(cell_idx, :));
for dir = 1:8
    subplot(3,3,subplot_i(dir))
    for drug = [1 2]
        plot(xx, MbHist{drug}(dir,:), 'color', color{drug})
        hold on
    end
    xlim([0.5 2.1])
    ylim([0 100])
end

% B
Drug = {'control', 'AP5', 'wash'};

h = figure;
set(h, 'Position', [1 1 1520,1080])
for cl = 1:7
    subplot(1, 7, cl)
    for drug = 1:3
        errorbar(xsort/pi*180, rho_mb_all_mean{cl}{drug}, rho_mb_all_ste{cl}{drug}, color(drug));
        hold on
    end
    xlabel('degrees')

    if cl == 1
        ylabel(ct{i})
    end

end
legend(Drug)

% C
ctr_i = 5;
x = 0:45:315;
temp = RHO_mb_all{1}{ctr_i} - RHO_mb_all{2}{ctr_i};
temp_mean = mean(temp, 1);
temp_ste = std(temp, [], 1)/sqrt(size(temp, 1));
figure
errorbar(x,temp_mean, temp_ste)
xlabel('degree')
ylabel('absolute suppression (spike #)')

figure
plot(repmat(x', 1, size(temp, 1)), temp', 'color', [.5 .5 .5]) 
xlabel('degree')
ylabel('absolute suppression (spike #)')


temp = (RHO_mb_all{1}{ctr_i} - RHO_mb_all{2}{ctr_i})./RHO_mb_all{1}{ctr_i};
temp(isinf(temp)) = 0;
temp(isnan(temp)) = 0;
temp_mean = nanmean(temp, 1);
temp_ste = nanstd(temp, [], 1)/sqrt(size(temp, 1));
figure
errorbar(x, temp_mean, temp_ste)
ylim([0 1])
xlabel('degree')
ylabel('percentage suppression')

figure
plot(repmat(x', 1, size(temp, 1)), temp', 'color', [.5 .5 .5]) 
xlabel('degree')
ylabel('percentage suppression')

% G
ctr_x = [5 10 20 40 80 150 300];
temp = response_pmax_all{1} - response_pmax_all{2};
temp_mean = mean(temp, 1);
temp_ste = std(temp, [], 1)/sqrt(size(temp, 1));
figure
errorbar(ctr_x, temp_mean, temp_ste)
set(gca, 'xscal', 'log')
xlabel('% contrast')
ylabel('absolute suppression (spike #)')

figure
plot(repmat(ctr_x', 1, size(temp, 1)), temp', 'color', [.5 .5 .5]) 
set(gca, 'xscal', 'log')
xlabel('% contrast')
ylabel('absolute suppression (spike #)')


temp = (response_pmax_all{1} - response_pmax_all{2})./response_pmax_all{1};
temp(isinf(temp)) = nan;
temp_mean = nanmean(temp, 1);
temp_ste = nanstd(temp, [], 1)/sqrt(size(temp, 1));
figure
errorbar(ctr_x, temp_mean, temp_ste)
ylim([0 1])
set(gca, 'xscal', 'log')
xlabel('% contrast')
ylabel('percentage suppression')

figure
plot(repmat(ctr_x', 1, size(temp, 1)), temp', 'color', [.5 .5 .5]) 
set(gca, 'xscal', 'log')
xlabel('% contrast')
ylabel('percentage suppression')

%% direction tuning (sliding window)

clear Max_i
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'wash'};
step_size = 50;
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

for drug = 1:2
    response_max_all{drug} = cell2mat(response_max{drug}');
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


clear RHO_thresh RHO_thresh_mean RHO_thresh_ste
for drug = 1:3
    for dir = 1:4
        CC = 1;
        RHO_thresh{drug}{dir} = [];
        for cc = 1:length(id_dir_mb{dir})
            if ~mb_idx(idx_dir_mb{dir}(cc))
                for ctr = 7:-1:1
                    spikeNum = length(raster_p_sum{drug}{idx_dir_mb{dir}(cc)}{ctr});
                    if spikeNum > 0 && spikeNum < 50
                        RHO_thresh{drug}{dir}(CC,:) = RHO_mb{drug}{dir}{ctr}(cc,:);
                        CC = CC + 1;
                        break;
                    end
                end
            end
        end
        if (~isempty(RHO_thresh{drug}{dir}))
            RHO_thresh_mean{drug}(dir, :) = mean(RHO_thresh{drug}{dir}, 1);
            RHO_thresh_ste{drug}(dir, :) = std(RHO_thresh{drug}{dir}, [], 1)/sqrt(size(RHO_thresh{drug}{dir}, 1));
        end
    end
end

figure
for dir = 1:4
    subplot(2,2,dir)
    errorbar(theta, RHO_thresh_mean{1}(dir, :), RHO_thresh_ste{1}(dir, :), 'b');
    hold on
    errorbar(theta, RHO_thresh_mean{2}(dir, :), RHO_thresh_ste{2}(dir, :), 'r');
    xlabel('direction')
    ylabel('spike #')
    title(dscell_type{dir})
    legend('control', 'AP5')
end
    
figure
RHO_thresh_temp = cell2mat(RHO_thresh{1}');
for cc = 1:size(RHO_thresh_temp, 1)
    subplot(8,8,cc)
    plot(theta, RHO_thresh_temp(cc, :))
end

figure
RHO_thresh_temp = cell2mat(RHO_thresh{2}');
for cc = 1:size(RHO_thresh_temp, 1)
    subplot(8,8,cc)
    plot(theta, RHO_thresh_temp(cc, :))
end

%% 
for drug = 1:2
    for cc = 1:length(ds_id)
        if ~mb_idx(cc)
            for ctr = 1:7
                dsiVector{drug}(cc, ctr) = MB{drug}.dsindex{ctr}(cc);
                spikeCount{drug}(cc, ctr) = length(raster_p_sum{drug}{cc}{ctr})/rep;
            end
        end
    end
end

figure
scatter(spikeCount{1}(:), dsiVector{1}(:), 'b');
hold on
scatter(spikeCount{2}(:), dsiVector{2}(:), 'r');
xlabel('spike# at PD')
ylabel('DSI')
title('2016-10-17-0')
legend('control', 'AP5', 'location', 'southeast')
ylim([0 1])

xlim([0 10])

%%
for drug = 1:2
    CC = 1;
    for dir = 2:4
    for cc = 1:length(idx_dir{dir})
        if ~mb_idx(idx_dir{dir}(cc))
            for ctr = 1:7
                dsiVector{drug}(CC, ctr) = MB{drug}.dsindex{ctr}(idx_dir{dir}(cc));
                spikeCount{drug}(CC, ctr) = length(raster_p_sum{drug}{idx_dir{dir}(cc)}{ctr})/rep;
            end
            CC = CC + 1;
        end
    end
    end
end


figure
[spikeCountAvg{1}, dsiVectorAvg{1}, disError{1}] = curve_from_binning(spikeCount{1}(:), dsiVector{1}(:), 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:2:40);
h1 = errorbar(spikeCountAvg{1}, dsiVectorAvg{1}, disError{1}, 'b');
hold on
scatter(spikeCount{1}(:), dsiVector{1}(:), [], [0.7 0.7 1]);
[spikeCountAvg{2}, dsiVectorAvg{2}, disError{2}] = curve_from_binning(spikeCount{2}(:), dsiVector{2}(:), 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:2:40);
h2 = errorbar(spikeCountAvg{2}, dsiVectorAvg{2}, disError{2}, 'r');
scatter(spikeCount{2}(:), dsiVector{2}(:), [], [1 0.7 0.7]);


xlabel('spike# at PD')
ylabel('DSI')
title('2016-10-17-0')
legend([h1, h2],'control', 'AP5', 'location', 'southeast')
ylim([0 1])


%% paper figures
% surface plot
color = 'brgk';
for ctr = 1:7
    for drug = 1:2
        RHO_surf{drug}(ctr, :) = RHO_mb_all_mean{ctr}{drug};
    end
end

theta = linspace(-180, 180, 9);
theta = theta(2:9);
ctr_x = [5 10 20 40 80 150 300];

[xx, yy] = meshgrid(theta, ctr_x);

conditions = {'control', 'AP5'};
figure
for drug = 1:2
    subplot(1, 2, drug)
    surf(xx, yy, RHO_surf{drug})
    set(gca, 'yscale', 'log')
    xlabel('direction (degree)')
    ylabel('contrast')
    zlabel('spike number')
    title(conditions{drug})
end

% direction tuning
ctr = 5; % contrast = 80%
figure
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
title('contrast = 80%')

% CRF 
figure
for drug = 1:2
    crf_mean(drug, :) = mean(response_pmax_all{drug});
    crf_ste(drug, :) = std(response_pmax_all{drug})./sqrt(size(response_pmax_all{drug}, 1));
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

% heat map
conditions = {'control', 'AP5'};
[xx, yy] = meshgrid(theta, ctr_x);
figure
for drug = 1:2
    response_max_all_mean{drug} = squeeze(mean(response_max_all{drug}));
    subplot(1, 2, drug)
    surf(xx, yy, response_max_all_mean{drug}'/max(response_max_all_mean{1}(:)))
    set(gca, 'yscale', 'log')
    xlabel('direction (degree)')
    ylabel('contrast')
    zlabel('spike number')
    title(conditions{drug})
    xlim([min(xx(:)) max(xx(:))])
    ylim([min(yy(:)) max(yy(:))])
    caxis([0 1])
end

%
a = (RHO_surf{1} -  RHO_surf{2})./RHO_surf{1};
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
[xx, yy] = meshgrid(theta, ctr_x);
s = surf(xx, yy, a);
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

%% 
figure
subplot(1, 2, 1)
hist(Sigma1, 10)
subplot(1, 2, 2)
hist(Sigma2, 10)

%% fit ds
clear Max_i
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'wash'};
step_size = 60;
for drug = 1:4
    for dir = 1:4
        for theta = 1:8
            CC = 1;
            for cc = 1:length(idx_dir{dir})
                if ~isempty(raster_mb_all{drug}{idx_dir{dir}(cc)})
                    for repeat = 1:10
                        for ctr = 7:-1:1
                            a = raster_mb{drug}{idx_dir{dir}(cc)}{1,1,theta,ctr,repeat};
                            hist_temp = hist(a, xx);
%                             if drug == 1 && ctr == 7
                                [max_fr, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
%                                 Max_i{dir}(theta, CC) = max_i;
%                             else
%                                 max_fr = conv(hist_temp, ones(1,step_size), 'valid');
%                                 max_fr = max_fr(Max_i{dir}(theta, CC));
%                             end
                            response_max{drug}{dir}(CC, theta, ctr, repeat) = max_fr - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
                        end
                    end
                    CC = CC + 1;
                end
            end
%             response_max{drug}{dir} = mean(response_max{drug}{dir}, 4);
%             response_max_norm{drug}{dir}(:, theta, :) = squeeze(response_max{drug}{dir}(:, theta, :))./repmat(squeeze(max(response_max{1}{dir}(:, theta, :), [], 1)), 1, size(response_max{drug}{dir},1))';
        end
        response_max{drug}{dir} = mean(response_max{drug}{dir}, 4);

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

for drug = 1:2
    response_max_all{drug} = cell2mat(response_max{drug}');
end

dir = 4;
CC = 1;
x = log10([80 150 300]);
for cc = 1:size(response_max_all{drug}, 1)
%     figure
%     set(gcf, 'Position', [1 1 900 800])
    
    for drug = 1:2

        ydata = squeeze(response_max_all{drug}(cc, dir, :));
        xdata = log10(ctr_x);
        [f, G] = fit_nr(xdata, ydata', 'upper', [100, 100, log10(300), 0]);
        fit_all{drug}{cc} = f;
        G_all{drug}{cc} = G;

        y(drug, :) = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;
    end
    if y(1, 2) > (y(1, 3) + y(1, 1))/2 && y(2, 2) > (y(2, 3) + y(2, 1))/2 && fit_all{1}{cc}.sigma < log10(150) && fit_all{2}{cc}.sigma < log10(150)
        for drug = 1:2    
            sigma(drug, CC) = fit_all{drug}{cc}.sigma;
            ymax(drug, CC) = fit_all{drug}{cc}.ymax;
            aa(drug, CC) = fit_all{drug}{cc}.a;
            bb(drug, CC) = fit_all{drug}{cc}.b;
        end
        CC = CC + 1;
    end

%         subplot(5,6,cc)
%         plot(xdata, ydata)
%         hold on
%         plot(x, y)


end
p = signrank(10.^sigma(1, :), 10.^sigma(2, :))
[~, p] = ttest(sigma(1, :), sigma(2, :))
%
figure
subplot(1, 2, 1)
hist(sigma(1, :))
xlabel('c50')
ylabel('cell #')
title('control')
subplot(1, 2, 2)
hist(sigma(2, :))
title('AP5')
xlabel('c50')
ylabel('cell #')

%%
color_temp = 'br';
figure
i = 1;
for cc = 1:size(response_max_all{drug}, 1)
    if i > 21
        i = 1;
        figure
    end
    subplot(3, 7, i)
    for drug = 1:2
        ydata = squeeze(response_max_all{drug}(cc, dir, :));
        xdata = log10(ctr_x);
        [f, G] = fit_nr(xdata, ydata', 'upper', [100, 100, log10(300), 0]);
        x = linspace(min(xdata), max(xdata), 100);
        y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;
        sigma(drug, cc) = f.sigma;
        plot(xdata, ydata, color_temp(drug))
        hold on
        plot(x, y, color_temp(drug))
    end
    xlabel('log(contrast)')
    ylabel('response')
    title(['s1: ' num2str(sigma(1, cc)) '  s2: ' num2str(sigma(2, cc))])
    i = i + 1;
end

%% 
figure
scatter(10.^sigma(1, :), 10.^sigma(2, :))
hold on
plot([0 150], [0 150])
xlim([0 150])
ylim([0 150])
xlabel('C50 (%) at control')
ylabel('C50 (%) at AP5')

%% white noise cross correlation
duration = 2700;
bin_size = 0.00025;
max_lag = 40;
ct = 1;
N = 10000;
xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
corr_cells_test = [];
dis = [];
area_ccf = [];
area_ccf_min = [];
pos = datawn.ei.position;
mode = 'neg';

for c1 = 1:length(id_dir{ct})-1
%     FigHandle = figure;
%     set(FigHandle, 'Position', [1 1 1980 1080])
    for c2 = c1+1:length(id_dir{ct})
        if c1 ~= c2
            id1 = id_dir{ct}(c1);
            id2 = id_dir{ct}(c2);
            idx1 = get_cell_indices(datawn, id1);
            idx2 = get_cell_indices(datawn, id2);
            spikes1 = datawn.spikes{idx1};
            spikes1_TF= ceil(spikes1/bin_size);
            spikes1 = zeros(duration/bin_size, 1);
            spikes1(spikes1_TF) = 1;
            
            spikes2 = datawn.spikes{idx2};
            spikes2_TF= ceil(spikes2/bin_size);
            spikes2 = zeros(duration/bin_size, 1);
            spikes2(spikes2_TF) = 1;
            
            A = xcorr(spikes1, spikes2, max_lag, 'coeff');
%             [maxv(c1, c2), maxi(c1, c2)] = max(A);
%             a = round(0.001/bin_size)+max_lag;
%             b = conv(A, ones(1, 11), 'valid');
%             ratio(c1, c2) = (sum(A(a:a+10)) + sum(A(max_lag*2-a-10:max_lag*2-a)) - min(A)*22)/(min(b)*2 - min(A)*22);
% %             ratio(c1, c2) = (sum(A(a:a+10)) + sum(A(max_lag*2-a-10:max_lag*2-a)) - min(A)*20)/ ...
% %                 (sum(A(1:11)) + sum(A(end-10:end)) - min(A)*20);
% %             A = xcov(spikes1, spikes2, max_lag, 'coeff');
% %             A_count = xcorr(spikes1, spikes2, max_lag);
%             [h, filteredA] = find_smallest_h(A);
% %             h_hat = zeros(1, N);
% %             for i = 1:N
% %                 samples = sample_dist(filteredA, round(sum(A_count)));
% %                 A2 = hist(samples, max_lag*2+1);
% %                 h_hat(i) = find_smallest_h(A2');
% %             end
% %             
%             [bootstat,bootsam] = bootstrp(N,@find_smallest_h_hist,rude(round(filteredA*sum(A)/sum(filteredA)), 1:max_lag*2+1), max_lag);
% %             p(c1, c2) = sum(h_hat > h)/N;
%             p(c1, c2) = sum(bootstat > h)/N;
%             subplot(5, 5, c2)
%             if p(c1, c2) < 0.05 && ratio(c1, c2) > 2 && maxi(c1, c2) > 0.75*max_lag && maxi(c1, c2) < 1.25*max_lag+1
% %                 bar(xx, A, 'r')
%                corr_cells_test = [corr_cells_test; id1 id2];
%             else
% %                 bar(xx, A, 'b')
%             end
% % %             subplot(5, 5, c2)
% % %             bar(xx, A, 'b')
% %             title([num2str(id1) '  ' num2str(id2) '  ' num2str(p(c1, c2))])
% %             xlim([-0.01 0.01])
% % 
% % %             pause
%             ei1 = datawn.ei.eis{idx1};
%             com1 = ei_com_xy(ei1, pos, 30*3, mode);
%             ei2 = datawn.ei.eis{idx2};
%             com2 = ei_com_xy(ei2, pos, 30*3, mode);
%             dis = [dis pdist([com1;com2])];
            dp_temp = sum(A(21:61) - min(A(21:61)));
%             area_ccf = [area_ccf dp_temp];
            area_ccf_min = [area_ccf_min dp_temp];
        end
    end
    c1
%     print_close(1, [24 12], num2str(id1));
end

[area_ccf_sort, i] = sort(area_ccf);
dis_sort = dis(i);


figure
plot(dis_sort(4:end), area_ccf_sort(4:end), 'o')
xlabel('distance (um)')
ylabel('correlated activity (area under CCF)')

figure
plot(dis, area_ccf_min, 'o')
xlabel('distance (um)')
ylabel('correlated activity (area under CCF)')
%% neighboring pairs
pos = datawn.ei.position;
mode = 'neg';
neighbors = [];
ct = 1;
for cc1 = 1:length(id_dir{ct})
    for cc2 = cc1+1:length(id_dir{ct})
        id1 = id_dir{ct}(cc1);
        idx1 = get_cell_indices(datawn, id1);
        ei1 = datawn.ei.eis{idx1};
        com1 = ei_com_xy(ei1, pos, 30*3, mode);
        id2 = id_dir{ct}(cc2);
        idx2 = get_cell_indices(datawn, id2);
        ei2 = datawn.ei.eis{idx2};
        com2 = ei_com_xy(ei2, pos, 30*3, mode);
        if pdist([com1;com2]) < 150
            neighbors = [neighbors; id1 id2];
        end
    end
end

cn = 0;
ct = 1;
celltype = {'superior', 'anterior', 'inferior', 'posterior'};
coms = [];
for cc = 1:length(id_dir{ct})
    id = id_dir{ct}(cc);
    idx = get_cell_indices(datawn, id);
    ei = datawn.ei.eis{idx};
    com = ei_com_xy(ei, pos, 30*3, mode);
    coms = [coms; com];
end

corner_i = [4 126 195 264 386 455 4];
corner_position = datawn.ei.position(corner_i, :);
figure
for cc = 1:length(id_dir{ct})
    plot(coms(cc, 1), coms(cc, 2),'ko')
    hold on
    text(coms(cc, 1)+5, coms(cc, 2)+5, num2str(id_dir{ct}(cc)), 'FontSize', 10)
    
end

% for cp = 1:size(corr_cells)
%     idx1 = find(id_dir{1} == corr_cells(cp, 1));
%     idx2 = find(id_dir{1} == corr_cells(cp, 2));
%     plot([coms(idx1, 1), coms(idx2, 1)], [coms(idx1, 2), coms(idx2, 2)], 'k');
% end
% 
plot(corner_position(:, 1), corner_position(:, 2), 'color', [.5 .5 .5])
axis off
title(celltype{ct})

% ct = 1;
cp_i = [];
for c1 = 1:length(id_dir{ct})-1
    for c2 = c1+1:length(id_dir{ct})
        if norm([coms(c1, :) - coms(c2, :)]) < 150
            cn = cn + 1;
            cp_i = [cp_i; c1 c2];
        end
    end
end
%%
duration = 2700;
bin_size = 0.00025;
max_lag = 40;
ct = 2;
N = 10000;

xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
FigHandle = figure;
set(FigHandle, 'Position', [1 1 2000 2000])
A_all = [];
for cp = 1:size(neighbors, 1)
%     id1 = corr_cells(cp, 1);
%     id2 = corr_cells(cp, 2);
    
    id1 = neighbors(cp, 1);
    id2 = neighbors(cp, 2);

    idx1 = get_cell_indices(datawn, id1);
    idx2 = get_cell_indices(datawn, id2);
    spikes1 = datawn.spikes{idx1};
    spikes1_TF= ceil(spikes1/bin_size);
    spikes1 = zeros(duration/bin_size, 1);
    spikes1(spikes1_TF) = 1;

    spikes2 = datawn.spikes{idx2};
    spikes2_TF= ceil(spikes2/bin_size);
    spikes2 = zeros(duration/bin_size, 1);
    spikes2(spikes2_TF) = 1;

    A = xcorr(spikes1, spikes2, max_lag, 'coeff');
    A_all = [A_all A];
% %     A = xcorr(spikes1, spikes2, max_lag);
%     [h, filteredA] = find_smallest_h(A);
%     [bootstat,bootsam] = bootstrp(N,@find_smallest_h_hist,rude(round(filteredA), 1:max_lag*2+1), max_lag);
%     p = sum(bootstat > h)/N;
    subplot(7, 6, cp)
%     if p(c1, c2) < 0.05 && ratio(c1, c2) > 2 && maxi(c1, c2) > 0.75*max_lag && maxi(c1, c2) < 1.25*max_lag+1 && maxv(c1, c2) > 10
       bar(xx, A, 'k')
%     else
%        bar(xx, A, 'r')
%     end
    title([num2str(id1) '  ' num2str(id2) '  ' num2str(p(find(id_dir{1} == id1), find(id_dir{1} == id2)))])
    xlim([-0.01 0.01])

end

%%
A_superior_mean = mean(A_all_superior, 2);
A_superior_ste = std(A_all_superior, [], 2)/sqrt(size(A_all_superior, 2));
figure
patch([xx fliplr(xx)], [A_superior_mean + A_superior_ste; flipud(A_superior_mean - A_superior_ste)]', [1 1 1]*0.8)
hold on
plot(xx, A_superior_mean, 'k')
% axis off

A_anterior_mean = mean(A_all_anterior, 2);
A_anterior_ste = std(A_all_anterior, [], 2)/sqrt(size(A_all_anterior, 2));
figure
patch([xx fliplr(xx)], [A_anterior_mean + A_anterior_ste; flipud(A_anterior_mean - A_anterior_ste)]', [1 1 1]*0.8)
hold on
plot(xx, A_anterior_mean, 'k')
% axis off
ylim([0 0.015])

%% test for multimodality 
% find smallest h that makes filtered distribution unimodal
[h, filteredA] = find_smallest_h(A);
[bootstat,bootsam] = bootstrp(500,@find_smallest_h_hist,rude(round(filteredA), 1:max_lag*2+1), max_lag);
figure
bar(A)
hold on
plot(filteredA)

%% null direction responses

trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'AP5+HEX', 'wash'};
step_size = 30;
for drug = 1:4
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_n_sum{drug}{idx_dir{dir}(cc)})
                for ctr = 5:-1:1
                    a = raster_n_sum{drug}{idx_dir{dir}(cc)}{ctr};
                    hist_temp = hist(a, xx);
                    if drug == 1 && ctr == 5
                        [max_p, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
                        Max_i{dir}(CC) = max_i;
                    else
                        max_p = conv(hist_temp, ones(1,step_size), 'valid');
                        max_p = max_p(Max_i{dir}(CC));
                    end
%                     response_pn{drug}{dir}(CC, ctr) = max(max_p - mean_n, 0);
%                     response_pn{drug}{dir}(CC, ctr) = max_p - mean_n;
%                     response_pn{drug}{dir}(CC, ctr) = abs(max_p - mean_n);
                    response_pn{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions;
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
for drug = 1:4
    subplot(2,2,drug)
%     response_s{drug} = max(exciseRows_empty(response_pn{drug}{1}), 0);
    response_s{drug} = exciseRows_empty(response_pn{drug}{1});
    mean_temp = nanmean(response_s{drug});
    ste_temp = nanstd(response_s{drug})/sqrt(size(response_s{drug}, 1));
    errorbar(ctr_x(1:5), mean_temp(1:5), ste_temp(1:5), 'color', color(1));
    hold on
%     response_others{drug} = max(exciseRows_empty(cell2mat(response_pn{drug}(2:4)')), 0);
    response_others{drug} = exciseRows_empty(cell2mat(response_pn{drug}(2:4)'));
    mean_temp = nanmean(response_others{drug});
    ste_temp = nanstd(response_others{drug})/sqrt(size(response_others{drug}, 1));
    errorbar(ctr_x(1:5), mean_temp(1:5), ste_temp(1:5), 'color', color(2));
    legend('superior', 'others', 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(condition{drug})
    xlim([3 400])
end

%% white noise cross correlation
duration = 2700;
bin_size = 0.00025;
max_lag = 40;
ct = 1;
N = 10000;

xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
FigHandle = figure;
set(FigHandle, 'Position', [1 1 2000 2000])

for cp = 1:size(cp_i, 1)
    id1 = id_dir{1}(cp_i(cp, 1));
    id2 = id_dir{1}(cp_i(cp, 2));
    idx1 = get_cell_indices(datawn, id1);
    idx2 = get_cell_indices(datawn, id2);
    spikes1 = datawn.spikes{idx1};
    spikes1_TF= ceil(spikes1/bin_size);
    spikes1 = zeros(duration/bin_size, 1);
    spikes1(spikes1_TF) = 1;

    spikes2 = datawn.spikes{idx2};
    spikes2_TF= ceil(spikes2/bin_size);
    spikes2 = zeros(duration/bin_size, 1);
    spikes2(spikes2_TF) = 1;

    A = xcorr(spikes1, spikes2, max_lag, 'coeff');
%     subplot(5, 6, cp)
%     bar(xx, A)
    cp_index(cp) = (max(A(8:10)) + max(A(12:14)) - 2*min(A)) / (2*(A(11) - min(A)));

end

load('DS161017.mat', 'cp_index_wt')
figure
subplot(2, 1, 1)
hist(cp_index_wt)
subplot(2, 1, 2)
hist(cp_index_ko)

%%
duration = 2700/2;
bin_size = 0.00025/2;
max_lag = 40*2;
ct = 1;
N = 10000; 
xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
FigHandle = figure;
set(FigHandle, 'Position', [1 1 2000 2000])
A_all = [];
for cp = 1:size(corr_cells, 1)
    id1 = corr_cells(cp, 1);
    id2 = corr_cells(cp, 2);
    
%     id1 = neighbors(cp, 1);
%     id2 = neighbors(cp, 2);

    idx1 = get_cell_indices(datawn, id1);
    idx2 = get_cell_indices(datawn, id2);
    spikes1 = datawn.spikes{idx1};
    spikes1(spikes1 > duration) = [];
    spikes1_TF= ceil(spikes1/bin_size);
    spikes1 = zeros(duration/bin_size, 1);
    spikes1(spikes1_TF) = 1;

    spikes2 = datawn.spikes{idx2};
    spikes2(spikes2 > duration) = [];
    spikes2_TF= ceil(spikes2/bin_size);
    spikes2 = zeros(duration/bin_size, 1);
    spikes2(spikes2_TF) = 1;

%     A = xcorr(spikes1, spikes2, max_lag, 'coeff');
%     A_all = [A_all A];
    A = xcorr(spikes1, spikes2, max_lag);
    [maxv, maxi] = max(A);
    a = round(0.001/bin_size)+max_lag;
    b = conv(A, ones(1, 11), 'valid');
    ratio = (sum(A(a:a+10)) + sum(A(max_lag*2-a-10:max_lag*2-a)) - min(A)*22)/(min(b)*2 - min(A)*22);

    [h, filteredA] = find_smallest_h(A);
    [bootstat,bootsam] = bootstrp(N,@find_smallest_h_hist,rude(round(filteredA), 1:max_lag*2+1), max_lag);
    p = sum(bootstat > h)/N;
    subplot(5, 6, cp)
%     if p(c1, c2) < 0.05 && ratio(c1, c2) > 2 && maxi(c1, c2) > 0.75*max_lag && maxi(c1, c2) < 1.25*max_lag+1 && maxv(c1, c2) > 10
    if p < 0.05 && ratio > 2 && maxi > 0.75*max_lag && maxi < 1.25*max_lag+1 && maxv > 10
       bar(xx, A, 'k')
    else
       bar(xx, A, 'r')
    end
%     title([num2str(id1) '  ' num2str(id2) '  ' num2str(p(find(id_dir{1} == id1), find(id_dir{1} == id2)))])
    xlim([-0.01 0.01])
    title([num2str(id1) '  ' num2str(id2)])

end

%%
duration = 2700/2;
bin_size = 0.00025;
max_lag = 40;
ct = 1;
N = 10000; 
xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
FigHandle = figure;
set(FigHandle, 'Position', [1 1 2000 2000])
A_all = [];
for i = 1:3
    indices = [1 12 10];
    cp = indices(i);
%     id1 = corr_cells(cp, 1);
%     id2 = corr_cells(cp, 2);
    
    id1 = neighbors(cp, 1);
    id2 = neighbors(cp, 2);

    idx1 = get_cell_indices(datawn, id1);
    idx2 = get_cell_indices(datawn, id2);
    spikes1 = datawn.spikes{idx1};
    spikes1(spikes1 > duration) = [];
    spikes1_TF= ceil(spikes1/bin_size);
    spikes1 = zeros(duration/bin_size, 1);
    spikes1(spikes1_TF) = 1;

    spikes2 = datawn.spikes{idx2};
    spikes2(spikes2 > duration) = [];
    spikes2_TF= ceil(spikes2/bin_size);
    spikes2 = zeros(duration/bin_size, 1);
    spikes2(spikes2_TF) = 1;

    A = xcorr(spikes1, spikes2, max_lag, 'coeff');
%     A_all = [A_all A];
%     A = xcorr(spikes1, spikes2, max_lag);
%     [maxv(cp), maxi(cp)] = max(A);
%     a = round(0.001/bin_size)+max_lag;
%     b = conv(A, ones(1, 11), 'valid');
%     ratio(cp) = (sum(A(a:a+10)) + sum(A(max_lag*2-a-10:max_lag*2-a)) - min(A)*22)/(min(b)*2 - min(A)*22);
% 
%     [h, filteredA] = find_smallest_h(A);
%     [bootstat,bootsam] = bootstrp(N,@find_smallest_h_hist,rude(round(filteredA), 1:max_lag*2+1), max_lag);
%     p(cp) = sum(bootstat > h)/N;
    subplot(1, 3, i)
% %     if p(c1, c2) < 0.05 && ratio(c1, c2) > 2 && maxi(c1, c2) > 0.75*max_lag && maxi(c1, c2) < 1.25*max_lag+1 && maxv(c1, c2) > 10
%     if p(cp) < 0.05 && ratio(cp) > 2 && maxi(cp) > 0.75*max_lag && maxi(cp) < 1.25*max_lag+1
       bar(xx, A, 'k')
%     else
%        bar(xx, A, 'r')
%     end
%     title([num2str(id1) '  ' num2str(id2) '  ' num2str(p(find(id_dir{1} == id1), find(id_dir{1} == id2)))])
    xlim([-0.01 0.01])
    title([num2str(id1) '  ' num2str(id2)])

end


%%
duration = 2700;
bin_size = 0.00025;
max_lag = 40;
ct = 1;
N = 10000; 
xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
FigHandle = figure;
set(FigHandle, 'Position', [1 1 2000 2000])
A_all = [];
for cp = 1:length(corr_cells)
    id1 = corr_cells(cp, 1);
    id2 = corr_cells(cp, 2);
    
%     id1 = neighbors(cp, 1);
%     id2 = neighbors(cp, 2);

    idx1 = get_cell_indices(datawn, id1);
    idx2 = get_cell_indices(datawn, id2);
    spikes1 = datawn.spikes{idx1};
    spike_n = min(5500, length(spikes1));
    spikes1 = spikes1(1:spike_n);
    spikes1_TF= ceil(spikes1/bin_size);
    spikes1 = zeros(duration/bin_size, 1);
    spikes1(spikes1_TF) = 1;

    spikes2 = datawn.spikes{idx2};
    spike_n = min(5500, length(spikes2));
    spikes2 = spikes2(1:spike_n);
    spikes2_TF= ceil(spikes2/bin_size);
    spikes2 = zeros(duration/bin_size, 1);
    spikes2(spikes2_TF) = 1;

%     A = xcorr(spikes1, spikes2, max_lag, 'coeff');
%     A_all = [A_all A];
    A = xcorr(spikes1, spikes2, max_lag);
    [maxv(cp), maxi(cp)] = max(A);
    a = round(0.001/bin_size)+max_lag;
    b = conv(A, ones(1, 11), 'valid');
    ratio(cp) = (sum(A(a:a+10)) + sum(A(max_lag*2-a-10:max_lag*2-a)) - min(A)*22)/(min(b)*2 - min(A)*22);

    [h, filteredA] = find_smallest_h(A);
    [bootstat,bootsam] = bootstrp(N,@find_smallest_h_hist,rude(round(filteredA), 1:max_lag*2+1), max_lag);
    p(cp) = sum(bootstat > h)/N;
    subplot(5, 6, cp)
% %     if p(c1, c2) < 0.05 && ratio(c1, c2) > 2 && maxi(c1, c2) > 0.75*max_lag && maxi(c1, c2) < 1.25*max_lag+1 && maxv(c1, c2) > 10
    if p(cp) < 0.05 && ratio(cp) > 2 && maxi(cp) > 0.75*max_lag && maxi(cp) < 1.25*max_lag+1
       bar(xx, A, 'k')
    else
       bar(xx, A, 'r')
    end
%     title([num2str(id1) '  ' num2str(id2) '  ' num2str(p(find(id_dir{1} == id1), find(id_dir{1} == id2)))])
    xlim([-0.01 0.01])
    title([num2str(id1) '  ' num2str(id2)])

end
