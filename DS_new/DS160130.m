cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);


datadg{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data000-map/data000-map', opt);
datadg{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/stimuli/s00.mat';
datadg{1} = load_stim_matlab(datadg{1}, 'user_defined_trigger_interval', 10);
datadg{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data004-map/data004-map', opt);
datadg{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/stimuli/s04.mat';
datadg{2} = load_stim_matlab(datadg{2}, 'user_defined_trigger_interval', 10);
datadg{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data007-map/data007-map', opt);
datadg{3}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/stimuli/s07.mat';
datadg{3} = load_stim_matlab(datadg{3}, 'user_defined_trigger_interval', 10);
datadg{4} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data010-map/data010-map', opt);
datadg{4}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/stimuli/s10.mat';
datadg{4} = load_stim_matlab(datadg{4}, 'user_defined_trigger_interval', 10);
datadg{5} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data013-map/data013-map', opt);
datadg{5}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/stimuli/s13.mat';
datadg{5} = load_stim_matlab(datadg{5}, 'user_defined_trigger_interval', 10);

datamb{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data001-map/data001-map', opt);
datamb{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/stimuli/s01.mat';
datamb{1} = load_stim_matlab(datamb{1});
datamb{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data005-map/data005-map', opt);
datamb{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/stimuli/s05.mat';
datamb{2} = load_stim_matlab(datamb{2});
datamb{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data008-map/data008-map', opt);
datamb{3}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/stimuli/s08.mat';
datamb{3} = load_stim_matlab(datamb{3});
datamb{4} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data011-map/data011-map', opt);
datamb{4}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/stimuli/s11.mat';
datamb{4} = load_stim_matlab(datamb{4});
datamb{5} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data014-map/data014-map', opt);
datamb{5}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/stimuli/s14.mat';
datamb{5} = load_stim_matlab(datamb{5});

dataffp{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data002-map/data002-map', opt);
dataffp{1}.triggers = dataffp{1}.triggers(2:end);
dataffp{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data006-map/data006-map', opt);
dataffp{2}.triggers = dataffp{2}.triggers(2:end);
dataffp{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data009-map/data009-map', opt);
dataffp{3}.triggers = dataffp{3}.triggers(2:end);
dataffp{4} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data012-map/data012-map', opt);
dataffp{4}.triggers = dataffp{4}.triggers(2:end);
dataffp{5} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data015-map/data015-map', opt);
dataffp{5}.triggers = dataffp{5}.triggers(2:end);

i = 5;
[NumSpikesCell,MaxRate, StimComb] = get_spikescellstim(datadg{i},datadg{i}.cell_ids,0,0.025);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg{i});
params_idx = [3 4]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg{i}, ds_struct, params_idx);
title('Normalized Vector Sum Amplitude')
xlabel('High speed - 480 micron/s')
ylabel('Low speed - 240 micron/s')
% ds_id(20) = [];
load('DS160130.mat')

%% 
% spike count
n = 5;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg{i},ds_id,0,1);
    DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg{i}));
    raster_dg{i} = get_ds_raster(datadg{i}, ds_id);
    for j = 1:length(raster_dg{i})
        if(dg_idx(j, i))
            raster_dg{i}{j} = [];
        end
    end
end

% max firing rate
n = 5;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim(datadg{i},ds_id,0,0.05);
    DG{i} = sort_direction(dscellanalysis(MaxRate, StimComb,datadg{i}));
    raster_dg{i} = get_ds_raster(datadg{i}, ds_id);
    for j = 1:length(raster_dg{i})
        if(dg_idx(j, i))
            raster_dg{i}{j} = [];
        end
    end
end

[raster_dg_nds, DG_nds] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg{i},nonds_id,0,1);
    DG_nds{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg{i}));
    raster_dg_nds{i} = get_ds_raster(datadg{i}, nonds_id);
end

delta_p = 4; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

for i = 1:n
    if ismember(i, [2 5])
        [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
    else
        [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p-3});
    end
    MAG_all_norm_dg{i} = normalize_MAG(DG{i});
    rep = datadg{i}.stimulus.repetitions;
end

ll = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};

%%
ll = 5; t = 5;
ds_vs = DG{ll}.mag{t};
nds_vs = DG_nds{ll}.mag{t};
max_vs = max([ds_vs nds_vs]);
min_vs = min([ds_vs nds_vs]);

X = 0.05:0.1:2.55;
ds_vs_hist = hist(ds_vs, X);
nds_vs_hist = hist(nds_vs, X);
bar(X, [ds_vs_hist' nds_vs_hist'], 1, 'stacked')
xlim([0 2.2])

%% plot cell summary
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [4 5]);

for cc = 1:length(ds_id)
    plot_ds_raster(DG_cut, raster_dg_cut, cc, ds_id(cc), ll, 2, 3, 1)
end

%% plot single cell tuning curve

% use unnormalized vector sum as response
figure
v = datadg{2}.stimulus.params.SPATIAL_PERIOD./datadg{2}.stimulus.params.TEMPORAL_PERIOD*4;
subplot(1, 2, 1)
semilogx(v, exciseColumn(MAG_all_norm_dg{2}), 'b')
xlabel('micron/second')
ylabel('Response')
title(ll{2})
xlim([v(end) v(1)])

subplot(1, 2, 2)
semilogx(v, exciseColumn(MAG_all_norm_dg{5}), 'b')
xlabel('micron/second')
ylabel('Response')
title(ll{5})
xlim([v(end) v(1)])

%% classification based on speed tunning
L = 2;
ds_id_L = intersect(ds_id, datadg{L}.cell_ids);
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
ylabel('3rd Principal Component')
title('NDF 0')

% % automatic classification
% % pca
% mag_pca = MAG_all_norm_dg{L};
% mag_pca = mag_pca';
% mag_pca = nan2empty(mag_pca);
% 
% [id_sub, idx_sub] = deal(cell(2, 1));
% 
% FigHandle = figure;
% set(FigHandle, 'Position', [1 1 380 400])
% 
% [~,scores,~,~] = princomp(mag_pca);
% pc1 = 1; pc2 = 3;
% plot(scores(:, pc1), scores(:, pc2), 'o')
% hold on
% [x, y] = ginput;
% plot(x, y)
% xlabel('1st Principal Component')
% ylabel('3rd Principal Component')
% title('NDF 0')
% 
% IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
% [~, I] = find(IN' == 1);
% id_init = ds_id_L(I);
% 
% [C ia ib] = intersect(id_init, ds_id_L);
% vc = ones(length(ds_id_L),1);
% vc(ib) = 2; %initializing ds cells to cluster 2, everything else cluster 1
% close all;
% X = [];
% N = [];
% p = [];
% [idx N p] = clustering_analysis_plots(scores(:, [pc1 pc2]), 0,1, 2, 0, 1, 0, 0, 0,0, vc);
% idx_sub{1} = find(idx == 1);
% idx_sub{2} = find(idx == 2);
% id_sub{1} = ds_id_L(idx_sub{1});
% id_sub{2} = ds_id_L(idx_sub{2});

figure
plot(scores(idx_sub{1}, pc1), scores(idx_sub{1}, pc2), 'bo', scores(idx_sub{2}, pc1), scores(idx_sub{2}, pc2), 'ro')
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
title('NDF 0')
% 
% several cells are missing during mapping between NDF0 grating and NDF3
% grating, simply put them into ON-OFF DSGC class here
id_sub{2} = setdiff(ds_id, id_sub{1});
[~, idx_sub{1}] = intersect(ds_id, id_sub{1});
[~, idx_sub{2}] = intersect(ds_id, id_sub{2});



% figure
% v = 4*datadg{2}.stimulus.params.SPATIAL_PERIOD./datadg{2}.stimulus.params.TEMPORAL_PERIOD;
% subplot(1, 2, 1)
% semilogx(v, exciseColumn(MAG_all_norm_dg{2}(:, idx_sub{1})), 'r')
% hold on
% semilogx(v, exciseColumn(MAG_all_norm_dg{2}(:, idx_sub{2})), 'b')
% xlabel('micron/second')
% ylabel('Response')
% % title(ll{2})
% xlim([v(end) v(1)])

% subplot(1, 2, 2)
figure
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
% title(ll{5})
xlim([v(end) v(1)])

figure
subplot(1,2,1)
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{1})), 'b')
xlabel('micron/second')
ylabel('Response')
% title(ll{5})
xlim([v(end) v(1)])
subplot(1,2,2)
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{2})), 'r')
xlabel('micron/second')
ylabel('Response')
% title(ll{5})
xlim([v(end) v(1)])

t = 4;
figure
compass(DG{5}.U{t}(idx_sub{1}), DG{5}.V{t}(idx_sub{1}), 'r')
hold on
compass(DG{5}.U{t}(idx_sub{2}), DG{5}.V{t}(idx_sub{2}), 'b')

figure
subplot(1,2,1)
compass(DG{5}.U{t}(idx_sub{1}), DG{5}.V{t}(idx_sub{1}), 'b')
subplot(1,2,2)
compass(DG{5}.U{t}(idx_sub{2}), DG{5}.V{t}(idx_sub{2}), 'r')


%% plot average tunning curve
color = 'brk';
figure
idx_temp = [2 5];
for i = 1:2
    v = 4*datadg{idx_temp(i)}.stimulus.params.SPATIAL_PERIOD./datadg{idx_temp(i)}.stimulus.params.TEMPORAL_PERIOD;
    subplot(1, 2, i)
    for ct = 1:2
        mag_temp = exciseColumn(MAG_all_norm_dg{idx_temp(i)}(:, idx_sub{ct}));
        tuning_avg{i}(:, ct) = mean(mag_temp, 2);
        tuning_ste{i}(:, ct) = std(mag_temp, [], 2)/sqrt(size(mag_temp, 2));
        errorbar(v, tuning_avg{i}(:, ct), tuning_ste{i}(:, ct), color(ct), 'LineWidth', 2)
        hold on
    end
    set(gca, 'XScale', 'log')
    title(ll{idx_temp(i)})
    xlim([min(v) max(v)])
    xlabel('speed')
    ylabel('response')
end
legend('on-off DSGC', 'on DSGC', 'location', 'southeast')


%% compare across light level
CT = {'ON DSGC', 'ON-OFF DSGC'};
figure
v = 4*datadg{2}.stimulus.params.SPATIAL_PERIOD./datadg{2}.stimulus.params.TEMPORAL_PERIOD;
for j = 1:2
    subplot(1, 2, j)
    errorbar(v, tuning_avg{1}(:, j), tuning_ste{1}(:, j), 'b', 'LineWidth', 2)
    hold on
    errorbar(v, tuning_avg{2}(:, j), tuning_ste{2}(:, j), 'r', 'LineWidth', 2)
    title(CT{j})
    set(gca, 'XScale', 'log')
    xlim([min(v) max(v)])
    xlabel('speed')
    ylabel('response')

end
legend('NDF 3', 'NDF 0', 'location', 'northeast')

%% classify DSGC into subtypes (directions)
d = 5;
t = 3;
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

d = 5;
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

%% DS tuning curves (drifting grating)
ct = {'superior', 'anterior', 'inferior', 'posterior'};
% all ds cells
color = 'brgkc';

% t = 2;
dirn = 4;
D = 5;
T = 2;
figure
for d = 1:n
    p_direction = DG_cut{D}.angle{T}';
    xx = 0:pi/4:7*pi/4;
    xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
    xx(xx>pi) = xx(xx>pi)-2*pi;
    xx(xx<-pi) = xx(xx<-pi)+2*pi;


    subplot(2, 3, d)
    for i = 1:dirn
        for cc = 1:length(ds_id)
            if ~dg_idx(cc, d)
            [xsort, seq] = sort(xx(cc, :));
            y_temp = DG_cut{d}.rho{T}(cc, :);
            plot(xsort, y_temp(seq), 'b')
            hold on
            end
        end
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    title(ll{d})
    xlim([-pi pi])
end

%subtypes
figure
clear rho_dg_mean rho_dg_ste dsi_dg_mean dsi_dg_ste
for d = 1:5
    subplot(2, 3, d)
    for i = 1:dirn
        rho_dg{d}{i} = [];
        dsi_dg{d}{i} = [];
        for cc = 1:length(idx_dir{i})
            if ~dg_idx(idx_dir{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx_dir{i}(cc), :))>0
            [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
            y_temp = DG_cut{d}.rho{T}(idx_dir{i}(cc), :);
            plot(xsort, y_temp(seq), color(i))
            ylim([0 1])
%             pause
            hold on
            rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
            dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx_dir{i}(cc))];
            end
        end
        rho_dg_mean{d}(i, :) = mean(rho_dg{d}{i});
        rho_dg_ste{d}(i, :) = std(rho_dg{d}{i})/sqrt(size(rho_dg{d}{i}, 1));
        dsi_dg_mean{d}(i) = mean(dsi_dg{d}{i});
        dsi_dg_ste{d}(i) = std(dsi_dg{d}{i})/sqrt(length(dsi_dg{d}{i}));
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    title(ll{d})
    xlim([-pi pi])
end
dsi_dg_mean = cell2mat(dsi_dg_mean');
dsi_dg_ste = cell2mat(dsi_dg_ste');
% plot average (cell type)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for d = 1:5
    subplot(2, 3, d)
    for i = 1:dirn
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(i));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
    title(ll{d});
end
legend(ct)

% plot average (light level)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:dirn
    subplot(2, 2, i)
    for d = 1:5
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(d));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
    title(ct{i})
end
legend(ll)
% DSI
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = ct;
model_series = dsi_dg_mean';
model_error = dsi_dg_ste';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('DSI')
legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

%% on-off vs on
%subtypes
clear rho_dg_mean rho_dg_ste dsi_dg_mean dsi_dg_ste
for d = 5:5
%     subplot(2, 3, d)
    for i = 1:2
        rho_dg{d}{i} = [];
        dsi_dg{d}{i} = [];
        for cc = 1:length(idx_sub{i})
            if ~dg_idx(idx_sub{i}(cc), d) && sum(DG_cut{d}.RHO{T}(idx_sub{i}(cc), :))>0
            [xsort, seq] = sort(xx(idx_sub{i}(cc), :));
            y_temp = DG_cut{d}.RHO{T}(idx_sub{i}(cc), :);
            plot(xsort, y_temp(seq), color(i))
%             ylim([0 1])
%             pause
            hold on
            rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
            dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx_sub{i}(cc))];
            end
        end
        rho_dg_mean{d}(i, :) = mean(rho_dg{d}{i});
        rho_dg_ste{d}(i, :) = std(rho_dg{d}{i})/sqrt(size(rho_dg{d}{i}, 1));
        dsi_dg_mean{d}(i) = mean(dsi_dg{d}{i});
        dsi_dg_ste{d}(i) = std(dsi_dg{d}{i})/sqrt(length(dsi_dg{d}{i}));
    end
    xlabel('direction (rad)')
    ylabel('max rate (Hz)')
    title(ll{d})
    xlim([-pi pi])
end
dsi_dg_mean = cell2mat(dsi_dg_mean');
dsi_dg_ste = cell2mat(dsi_dg_ste');
% plot average (cell type)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for d = 5:5
%     subplot(2, 3, d)
    for i = 1:2
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(i));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('max rate (Hz)')
    title(ll{d});
end
legend('on DSGC', 'on-off DSGC')

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
d = 5;
for i = 1:2
    errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(i));
    hold on
end
xlabel('direction (rad)')
ylabel('max rate (Hz)')
title(ll{d});
legend('on DSGC', 'on-off DSGC')

%%
duration = 573;
[raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, ~, StimComb] = get_spikescellstim_mb(datamb{i},ds_id,duration,1);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb{i} = get_mb_raster(datamb{i}, ds_id, duration);
%     for j = 1:length(raster_mb{i})
%         if(mb_idx(j, i))
%             raster_mb{i}{j} = [];
%         end
%     end
    trial_dur{i} = get_mb_trial_dur(datamb{i}, 600, 800, 0);
end

delta_p = 2; % choose which params to use to calculate prefer direction indices 
ll_p = 5;
MAG_all_norm_mb = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

for i = 1:n
    [raster_p_sum_mb{i}, p_idx{i}] = get_pdirection_raster(raster_mb{i}, MB{ll_p}.angle{delta_p});
    MAG_all_norm_mb{i} = normalize_MAG(MB{i});
    rep = datamb{i}.stimulus.repetitions;
end

ll = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};
ct = {'superior', 'posterior', 'inferior', 'anterior'};
close all

%% plot cell summary
Title = {'BarWidth 120, white', 'BarWidth 240, white', 'BarWidth 360, white'; 'BarWidth 120, black', 'BarWidth 240, black', 'BarWidth 360, black'};
for cc = 2:length(ds_id)
    plot_mb_raster_one(MB, raster_mb, trial_dur, cc, ds_id(cc), num2str(ds_id(cc)), 1, 1, 1)
end

%% frequency analysis
clear ratio_avg ratio_ste ratio_dir_mean ratio_dir_ste ratio
duration = 8;
bin_rate = 10000;
hist_spikes = cell(n, 1);
signal_length = duration*bin_rate;                
NFFT = 2^nextpow2(signal_length);
f = bin_rate/2*linspace(0,1,NFFT/2+1);
fft_spikes = cell(n, 1);
[DC, F1, F2] = deal(cell(n, 1));
for i = [2 5]
%     tp = datadg{1}.stimulus.params.TEMPORAL_PERIOD;
    tp = datadg{2}.stimulus.params.TEMPORAL_PERIOD;
    [DC{i}, F1{i}, F2{i}] = deal(zeros(length(ds_id), length(tp)));
    for rgc = 1:length(ds_id)
%         if ~isempty(raster_dg_cut{i}{rgc}) && ~dg_idx(rgc, i)
        if ~isempty(raster_dg{i}{rgc}) && ~dg_idx(rgc, i)
            for time = 1:length(tp)
    %             spikes = floor(raster_p_sum_cut{i}{rgc}{time}*bin_rate);
                spikes = floor(raster_p_sum{i}{rgc}{time}*bin_rate);
                tmp_binned_spikes = zeros(1, signal_length);
                tmp_binned_spikes(spikes) = 1;
                hist_spikes{i}{rgc}(time, :) = tmp_binned_spikes;

                f1 = 1/tp(time); %Hz
                f2 = f1*2;
                f_diff1 = f - f1;
                f_diff2 = f - f2;
                [~,f1_index] = min(abs(f_diff1));
                [~,f2_index] = min(abs(f_diff2));
                tmp_fft = fft(tmp_binned_spikes, NFFT)./ signal_length;
                fft_spikes{i}{rgc}(time, :) = 2*abs(tmp_fft(1:NFFT/2+1));
                if f1_index > 1
                    fund_power(time) = sum(fft_spikes{i}{rgc}(time, f1_index-1:f1_index+1)); % f1_index+2???
                    sec_power(time) = sum(fft_spikes{i}{rgc}(time, f2_index-1:f2_index+1));
                else
                    fund_power(time) = sum(fft_spikes{i}{rgc}(time, f1_index:f1_index+2)); % f1_index+2???
                    sec_power(time) = sum(fft_spikes{i}{rgc}(time, f2_index:f2_index+2));
                end
                DC_power(time) = fft_spikes{i}{rgc}(time, 1);
            end
            % stores info for this cell into the matrix tuning curves
            F1{i}(rgc,:) = fund_power ./ max(DC_power);
            F2{i}(rgc,:) = sec_power ./ max(DC_power);
            DC{i}(rgc,:) = DC_power ./ max(DC_power);

            clear fund_power sec_power DC_power

        end
        
    end
    ratio{i} = F2{i}./F1{i};
    for ct = 1:4
        ratio_dir{ct}{i} = ratio{i}(idx_dir{ct}, :);
        ratio_dir{ct}{i} = exciseRows_empty(ratio_dir{ct}{i});
        ratio_dir_mean(i, ct, :) = mean(ratio_dir{ct}{i});
        ratio_dir_ste(i, ct, :) = std(ratio_dir{ct}{i})/sqrt(size(ratio_dir{ct}{i}, 1));
    end
    for ct = 1:3
        ratio_dir_on{ct}{i} = ratio{i}(idx_dir_on{ct}, :);
        ratio_dir_on{ct}{i} = exciseRows_empty(ratio_dir_on{ct}{i});
        ratio_dir_mean_on(i, ct, :) = mean(ratio_dir_on{ct}{i}, 1);
        ratio_dir_ste_on(i, ct, :) = std(ratio_dir_on{ct}{i}, [], 1)/sqrt(size(ratio_dir_on{ct}{i}, 1));
    end
    for ct = 1:2
        ratio_onoff{ct}{i} = ratio{i}(idx_sub{ct}, :);
        ratio_onoff{ct}{i} = exciseRows_empty(ratio_onoff{ct}{i});
        ratio_onoff_mean(i, ct, :) = mean(ratio_onoff{ct}{i}, 1);
        ratio_onoff_ste(i, ct, :) = std(ratio_onoff{ct}{i}, [], 1)/sqrt(size(ratio_onoff{ct}{i}, 1));
    end
    ratio{i} = exciseRows_empty(ratio{i});
    ratio_mean(i, :) = mean(ratio{i});
    ratio_ste(i, :) = std(ratio{i})/sqrt(size(ratio{i}, 1));
end

for i = 1:n
    tp = datadg{1}.stimulus.params.TEMPORAL_PERIOD;
    [DC{i}, F1{i}, F2{i}] = deal(zeros(length(ds_id), length(tp)));
    for rgc = 1:length(ds_id)
        if ~isempty(raster_dg_cut{i}{rgc}) && ~dg_idx(rgc, i)
        for time = 1:length(tp)
            spikes = floor(raster_p_sum_cut{i}{rgc}{time}*bin_rate);
            tmp_binned_spikes = zeros(1, signal_length);
            tmp_binned_spikes(spikes) = 1;
            hist_spikes{i}{rgc}(time, :) = tmp_binned_spikes;
            
            f1 = 1/tp(time); %Hz
            f2 = f1*2;
            f_diff1 = f - f1;
            f_diff2 = f - f2;
            [~,f1_index] = min(abs(f_diff1));
            [~,f2_index] = min(abs(f_diff2));
            tmp_fft = fft(tmp_binned_spikes, NFFT)./ signal_length;
            fft_spikes{i}{rgc}(time, :) = 2*abs(tmp_fft(1:NFFT/2+1));
            if f1_index > 1
                fund_power(time) = sum(fft_spikes{i}{rgc}(time, f1_index-1:f1_index+1)); % f1_index+2???
                sec_power(time) = sum(fft_spikes{i}{rgc}(time, f2_index-1:f2_index+1));
            else
                fund_power(time) = sum(fft_spikes{i}{rgc}(time, f1_index:f1_index+2)); % f1_index+2???
                sec_power(time) = sum(fft_spikes{i}{rgc}(time, f2_index:f2_index+2));
            end
            DC_power(time) = fft_spikes{i}{rgc}(time, 1);
        end
    % stores info for this cell into the matrix tuning curves
        F1{i}(rgc,:) = fund_power ./ max(DC_power);
        F2{i}(rgc,:) = sec_power ./ max(DC_power);
        DC{i}(rgc,:) = DC_power ./ max(DC_power);
        
        clear fund_power sec_power DC_power

        end
        
    end
    ratio{i} = F2{i}./F1{i};
    for ct = 1:4
        ratio_dir{ct}{i} = ratio{i}(idx_dir{ct}, :);
        ratio_dir{ct}{i} = exciseRows_empty(ratio_dir{ct}{i});
        ratio_dir_mean(i, ct, :) = mean(ratio_dir{ct}{i});
        ratio_dir_ste(i, ct, :) = std(ratio_dir{ct}{i})/sqrt(size(ratio_dir{ct}{i}, 1));
    end
    ratio_oo{i} = ratio{i}(idx_sub{2}, :);
    ratio_oo{i} = exciseRows_empty(ratio_oo{i});
    ratio_oo_mean(i, :) = mean(ratio_oo{i});
    ratio_oo_ste(i, :) = std(ratio_oo{i})/sqrt(size(ratio_oo{i}, 1));
end

% plot average f2/f1
ct = {'posterior', 'inferior', 'anterior', 'superior'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = {'high speed'; 'low speed'};
model_series = ratio_oo_mean';
model_error = ratio_oo_ste';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('F2/F1')
legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

% plot F2/F1 tuning curve
v3 = [4800 2400 960 480 240 120 80 60 40];
figure
for i = 1:2
    subplot(1, 2, i)
    errorbar(v3, ratio_onoff_mean(2, i, :), ratio_onoff_ste(2, i, :), 'b')
    hold on
    errorbar(v3, ratio_onoff_mean(5, i, :), ratio_onoff_ste(5, i, :), 'r')
    set(gca, 'XScale', 'log')
    xlim([100 inf])
end

% F2/F1 tuning for individual direction
celltype = {'posterior', 'inferior', 'anterior', 'superior'};
figure
for ct = 1:4
    subplot(2, 2, ct)
    errorbar(v3, ratio_dir_mean(2, ct, :), ratio_dir_ste(2, ct, :), 'b')
    hold on
    errorbar(v3, ratio_dir_mean(5, ct, :), ratio_dir_ste(5, ct, :), 'r')
    set(gca, 'XScale', 'log')
    title(celltype{ct})
end

%% plot cell type specific f2/f1 
for tp = 1:2
    figure
    set(gcf, 'DefaultLineLineWidth', 1.5)
    xtick = ct;
    model_series = ratio_dir_mean(:,:,tp)';
    model_error = ratio_dir_ste(:,:,tp)';
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('F2/F1')
    legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
    if tp == 1
        title('high speed')
    else
        title('low speed')
    end
    hold on;

    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 

    groupwidth = min(0.8, numbars/(numbars+1.5));

    for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end
end


%% full field pulses

n_ffp = 5;
bin_size = 0.1;
XX = bin_size/2:bin_size:6-bin_size/2;

[raster_ff, raster_ff_all] = deal(cell(n_ffp, 1));
for d = 1:n_ffp
    [raster_ff{d}, raster_ff_all{d}] = get_ffp_raster(dataffp{d}, ds_id, 3);
    for j = 1:length(raster_ff{d})
%         if(ffp_idx(j, d))
%             raster_ff{d}{j} = [];
%             raster_ff_all{d}{j} = [];
%         end
        step_raster = get_raster(raster_ff_all{d}{j}, [0 6], 'plot', 0);
        hist_on = hist(step_raster{1}, XX);
        hist_off = hist(step_raster{2}, XX);
        ratio_ffp{d}(j) = max(hist_off)/max(hist_on);
    end
    ratio_temp = ratio_ffp{d};
    ratio_temp(isnan(ratio_temp)) = [];
    ratio_ffp_avg(d) = mean(ratio_temp);
    ratio_ffp_ste(d) = std(ratio_temp)/sqrt(length(ratio_temp));
end

figure
errorbar([0:4], ratio_ffp_avg, ratio_ffp_ste)
title('ffp off on ratio')
xlabel('light level')
ylabel('off on ratio')
xlim([-0.5 4.5])

for i = 1:length(ds_id) 
%     if ~isempty(raster_ff{1}{i}) || ~isempty(raster_ff{2}{i})
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 1800 800])
        for d = 1:n_ffp
        subplot(1, n_ffp, d)
        plot_ffp(raster_ff{d}, raster_ff_all{d}, i, 3)
        title([num2str(ds_id(i)) ' ' ll{d}])
        end
        
        print_close(1, [24, 12], num2str(ds_id(i)))
%     end
end

%%
ct = 1;
pos = datadg{5}.ei.position;

figure
for cc = 1:1 %length(id_dir{ct})
    cell_idx = get_cell_indices(datadg{5}, id_dir_on{ct}(cc));
    ei = datadg{5}.ei.eis{cell_idx};
    plot_ei_(ei, pos)
    hold on
end
%% find nn of on/on-off pair of similar direction
load('DS160130.mat', 'nni')
clear rho_dg_mean rho_dg_ste dsi_dg_mean dsi_dg_ste
d = 5;
T = 1;
p_direction = DG_cut{d}.angle{T}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;
for i = 1:3
    figure
    for p = 1:size(nni{i}, 1)
        subplot(2, 3, p)
        idx = nni{i}(p, :);
        id_on = id_dir_on{i}(idx(1));
        id_oo = id_dir{i}(idx(2));
        [xsort, seq] = sort(xx(idx_dir_on{i}(idx(1)), :));
        y_temp = DG_cut{d}.RHO{T}(idx_dir_on{i}(idx(1)), :);
        plot(xsort, y_temp(seq), 'b')
        hold on
        [xsort, seq] = sort(xx(idx_dir{i}(idx(2)), :));
        y_temp = DG_cut{d}.RHO{T}(idx_dir{i}(idx(2)), :);
        plot(xsort, y_temp(seq), 'r')
        xlabel('direction (rad)')
        ylabel('max rate (Hz)')
        title(['cell-id: ' num2str(id_on) '-' num2str(id_oo)])
        if p == 1
            legend('on DSGC', 'on-off DSGC')
        end
    end
end


%% max firing rate speed tuning

max_rate_prefer = cell(5,1);
max_rate_null = cell(5,1);
max_rate_prefer_sub = cell(5, 1);
max_rate_null_sub = cell(5, 1);
n_idx = ceil(mod(p_idx{5}+4-0.1,8));
for i = [2 5]
    clear max_rate_prefer_temp max_rate_null_temp
    for tp = 1:length(DG{2}.RHO)
        for cc = 1:length(ds_id)
            if ~dg_idx(cc, i)
                max_rate_prefer_temp{cc}(tp) = DG{i}.RHO{tp}(cc, p_idx{5}(cc));
                max_rate_null_temp{cc}(tp) = DG{i}.RHO{tp}(cc, n_idx(cc));
            end
        end
    end
    max_rate_prefer{i} = max_rate_prefer_temp';
    max_rate_null{i} = max_rate_null_temp';
    for sub = 1:2
        max_rate_prefer_sub{i}{sub} = max_rate_prefer{i}(idx_sub{sub});
        max_rate_null_sub{i}{sub} = max_rate_null{i}(idx_sub{sub});
    end
end

% unnormalized individual tuning curve
color = 'br';
for i = [2 5]
    figure
    subplot(1, 2, 1)
    for sub = 1:2
        tuning_p = cell2mat(max_rate_prefer_sub{i}{sub});
        semilogx(repmat(v', 1, size(tuning_p, 1)), tuning_p', 'color', color(sub))
        hold on
    end
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')
    ylabel('Hz')
    subplot(1, 2, 2)
    for sub = 1:2
        tuning_n = cell2mat(max_rate_null_sub{i}{sub});
        semilogx(repmat(v', 1, size(tuning_n, 1)), tuning_n', 'color', color(sub))
        hold on
    end 
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')
    ylabel('Hz')

end

% unnormalized average tuning curve
for i = [2 5]
    figure
    subplot(1, 2, 1)
    for sub = 1:2
        tuning_p = cell2mat(max_rate_prefer_sub{i}{sub});
        errorbar(v, mean(tuning_p),std(tuning_p)/sqrt(size(tuning_p, 1)), 'color', color(sub))
        hold on
    end
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')
    ylabel('Hz')
    legend('ON', 'ON-OFF', 'Location', 'northwest')
    title('preferred direction')
    subplot(1, 2, 2)
    for sub = 1:2
        tuning_n = cell2mat(max_rate_null_sub{i}{sub});
        errorbar(v, mean(tuning_n), std(tuning_n)/sqrt(size(tuning_n, 1)), 'color', color(sub))
        hold on
    end
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')
    ylabel('Hz')
    legend('ON', 'ON-OFF', 'Location', 'northwest')
    title('null direction')
end

% normalized individual tuning curve
color = 'br';
for i = [2 5]
    figure
    subplot(1, 2, 1)
    for sub = 1:2
        tuning_p = cell2mat(max_rate_prefer_sub{i}{sub});
        tuning_p_n = tuning_p./repmat(max(tuning_p, [], 2), 1,9);
        tuning_p_n = nan2empty(tuning_p_n);
        semilogx(repmat(v', 1, size(tuning_p_n, 1)), tuning_p_n', 'color', color(sub))
        hold on
    end
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')
    subplot(1, 2, 2)
    for sub = 1:2
        tuning_n = cell2mat(max_rate_null_sub{i}{sub});
        tuning_n_n = tuning_n./repmat(max(tuning_n, [], 2), 1,9);
        tuning_n_n = nan2empty(tuning_n_n);
        semilogx(repmat(v', 1, size(tuning_n_n, 1)), tuning_n_n', 'color', color(sub))
        hold on
    end 
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')

end

% normalized average tuning curve
for i = [2 5]
    figure
    subplot(1, 2, 1)
    for sub = 1:2
        tuning_p = cell2mat(max_rate_prefer_sub{i}{sub});
        tuning_p_n = tuning_p./repmat(max(tuning_p, [], 2), 1,9);
        tuning_p_n = nan2empty(tuning_p_n);
        errorbar(v, mean(tuning_p_n),std(tuning_p_n)/sqrt(size(tuning_p_n, 1)), 'color', color(sub))
        hold on
    end
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')
    legend('ON', 'ON-OFF', 'Location', 'northwest')
    title('preferred direction')
    subplot(1, 2, 2)
    for sub = 1:2
        tuning_n = cell2mat(max_rate_null_sub{i}{sub});
        tuning_n_n = tuning_n./repmat(max(tuning_n, [], 2), 1,9);
        tuning_n_n = nan2empty(tuning_n_n);
        errorbar(v, mean(tuning_n_n), std(tuning_n_n)/sqrt(size(tuning_n_n, 1)), 'color', color(sub))
        hold on
    end
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')
    legend('ON', 'ON-OFF', 'Location', 'northwest')
    title('null direction')
end

%% fit direction tuning curves
T = 2;
for dir = 1:4
    figure
    for i = 1:5
        CC = 1;
        for cc = 1:length(id_dir{dir})
            if ~(dg_idx(idx_dir{dir}(cc), i))
                xdata = DG_cut{i}.theta{1}(idx_dir{dir}(cc), :);
                ydata = DG_cut{i}.rho{T}(idx_dir{dir}(cc), :);
                [f, g] = fit_cos(xdata, ydata);
                Ymax{dir}{i}(CC) = f.ymax;
                Phi{dir}{i}(CC) = f.phi;
                Alpha{dir}{i}(CC) = f.alpha;
                Width{dir}{i}(CC) = acos(2*0.5^(1/f.alpha) - 1)/pi*360;
                B{dir}{i}(CC) = f.b;
                CC = CC + 1;
                
                xfit = linspace(0, 2*pi, 100);
                yfit = f.ymax.*(0.5+0.5*cos(xfit+f.phi)).^f.alpha+f.b;
                subplot(4, 6, cc)
                plot(xdata, ydata, 'b')
                hold on
                plot(xfit, yfit, 'r')
                ylim([0 1.1])
                width = acos(2 * (0.5.^(1/f.alpha) - 0.5))/pi*360;
                title(['width = ', num2str(width)])


            end
        end
        Ymax_mean(dir, i) = robust_mean(Ymax{dir}{i});
        Phi_mean(dir, i) = robust_mean(Phi{dir}{i});
        Alpha_mean(dir, i) = robust_mean(Alpha{dir}{i});
        B_mean(dir, i) = robust_mean(B{dir}{i});
    end
end


for dir = 1:4
    figure

    for i = 1:5
        subplot(2, 3, i)
        hist(Width{dir}{i})
        title(num2str(mean(Width{dir}{i})))
    end
end

for i = 1:5
    for dir = 1:4
        WidthMean(dir, i) = mean(Width{dir}{i});
        WidthSte(dir, i) = std(Width{dir}{i})/sqrt(length(Width{dir}{i}));
    end
end

marker = 'xosd';
figure
for dir = 1:4
    errorbar(0:4, WidthMean(dir, :), WidthSte(dir, :), 'Color', 'k', 'Marker', marker(dir), 'MarkerSize', 10)
    hold on
end
legend('superior', 'anterior', 'inferior', 'posterior')
ylim([0 300])
xlabel('light level')
ylabel('tuning width (degree)')
xlim([-0.5 4.5])

%%
for ll = 1:5
    for ct = 1:4
        Width_all{ct}{ll} = [Width{ct}{ll} Width0130{ct}{ll}];
    end
end

for dir = 1:4
    figure

    for i = 1:5
        subplot(2, 3, i)
        hist(Width_all{dir}{i})
        title(num2str(mean(Width_all{dir}{i})))
    end
end

for i = 1:5
    for dir = 1:4
        WidthMean(dir, i) = mean(Width_all{dir}{i});
        WidthSte(dir, i) = std(Width_all{dir}{i})/sqrt(length(Width_all{dir}{i}));
    end
end

marker = 'xosd';
figure
for dir = 1:4
    errorbar(0:4, WidthMean(dir, :), WidthSte(dir, :), 'Color', 'k', 'Marker', marker(dir), 'MarkerSize', 10)
    hold on
end
legend('superior', 'anterior', 'inferior', 'posterior')
ylim([0 300])
xlabel('light level')
ylabel('tuning width (degree)')
xlim([-0.5 4.5])

%% mb preferred direction hist
color = 'brgkm';
% get spike histgram of preferred direction
bin_size = 0.2;
for i = 1:n
    for t = 1:2
        dur = max(trial_dur{i}(t));
        XX = bin_size/2:bin_size:dur;
        hist_mb_p{i}{t} = zeros(length(ds_id), length(XX));
        for cc = 1:length(ds_id)
            if ~isempty(raster_p_sum_mb{i}{cc})
                hist_mb_p{i}{t}(cc, :) = hist(raster_p_sum_mb{i}{cc}{t}, XX);
            end
        end

    end
end

%plot
% x = [4 4 4 4]; y = [4 4 5 5];
x = [2 2 2]; y = [3 2 2];

for t = 1:1

    for ct = 1:3
        figure
        set(gcf, 'Position', [1,1,1080,1080])
        for cc = 1:length(idx_dir_on{ct})
            subplot(x(ct), y(ct), cc)
    %         subplot(3, 4, cc-3)
            dur = max(trial_dur{1}(t,bw));
            XX = bin_size/2:bin_size:dur;
    %         for i = 4:4
            for i = 3:3
                plot(XX, 1.25*hist_mb_p{i}{t}(idx_dir_on{ct}(cc), :), color(i))
                hold on
            end
            title(num2str(ds_id(idx_dir_on{ct}(cc))))
            xlim([0 dur])
        end
        legend('NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0')
    end
end


x = [4 4 2 2]; y = [5 4 3 3];

for t = 1:1

    for ct = 1:4
        figure
        set(gcf, 'Position', [1,1,1080,1080])
        for cc = 1:length(idx_dir{ct})
            subplot(x(ct), y(ct), cc)
    %         subplot(3, 4, cc-3)
            dur = max(trial_dur{1}(t,bw));
            XX = bin_size/2:bin_size:dur;
    %         for i = 4:4
            for i = 3:3
                plot(XX, 1.25*hist_mb_p{i}{t}(idx_dir{ct}(cc), :), color(i))
                hold on
            end
            title(num2str(ds_id(idx_dir{ct}(cc))))
            xlim([0 dur])
        end
        legend('NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0')
    end
end


%% subtract bgnd firing
% estimate background activity
n = 5;
interval = 1;
raster_interval_dg = deal(cell(n, 1));
for i = 1:n    
    raster_interval_dg{i} = get_ds_interval_raster(datadg{i}, ds_id, interval);
    bgfr{i} = zeros(1, length(ds_id));
    for j = 1:length(raster_interval_dg{i})
        if iscell(raster_interval_dg{i}{j})
            bgfr{i}(j) = mean(cellfun(@length, raster_interval_dg{i}{j}))/interval;
        end
    end
    for dir = 1:4
        bgfr_ct{i}{dir} = [];
        for cc = 1:length(idx_dir{dir})
            if ~dg_idx(idx_dir{dir}(cc), i) && sum(DG_cut{i}.rho{1}(idx_dir{dir}(cc), :))>0
                bgfr_ct{i}{dir} = [bgfr_ct{i}{dir} bgfr{i}(idx_dir{dir}(cc))];
            end
        end
        bgfr_ct_mean(i, dir) = mean(bgfr_ct{i}{dir});
        bgfr_ct_ste(i, dir) = std(bgfr_ct{i}{dir})/sqrt(length(bgfr_ct{i}{dir}));
    end
end

celltype = {'superior', 'anterior', 'inferior', 'posterior'};
marker = 'xosd';
figure
for i = 1:4
    errorbar([0:4], bgfr_ct_mean(:, i), bgfr_ct_ste(:, i), 'Color', 'k', 'Marker', marker(i), 'MarkerSize', 10);
    hold on
end
xlim([-0.5 4.5])
xlabel('R*/rod/s')
ylabel('background firing (Hz)')
legend(celltype)

n = 5;
% spike count
DG_bgnd = cell(n, 1);
for i = 1:n    
    [NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg{i},ds_id,0, 0.05);
    DG_bgnd{i} = sort_direction(dscellanalysis_bgnd_subtract(NumSpikesCell, StimComb, datadg{i}, bgfr{i}));
end
[DG_bgnd_cut, ~, ~] = cut_dg(DG_bgnd, raster_dg, raster_p_sum, 2, [4 5]);

%% DS tuning curves (drifting grating)
% all ds cells

% t = 2;
dirn = 4;
D = 5;
T = 1;
color = 'bkrgc';
ct = {'superior', 'anterior', 'inferior', 'posterior'};
ll = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};
DG_cut = DG_bgnd_cut;

for d = 1:n
    p_direction = DG{D}.angle{1}';
    xx = 0:pi/4:7*pi/4;
    xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
    xx(xx>pi) = xx(xx>pi)-2*pi;
    xx(xx<-pi) = xx(xx<-pi)+2*pi;


    subplot(2, 3, d)
    for i = 1:dirn
        for cc = 1:length(ds_id)
            if ~dg_idx(cc, d)
            [xsort, seq] = sort(xx(cc, :));
%             xsort = xsort/pi*180;
            y_temp = DG_cut{d}.rho{T}(cc, :);
            plot(xsort, y_temp(seq), 'b')
            hold on
            end
        end
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    title(ll{d})
    xlim([-pi pi])
end

%subtypes
clear rho_dg_mean rho_dg_ste dsi_dg_mean dsi_dg_ste
for d = 1:5
%     subplot(2, 3, d)
    for i = 1:dirn
        rho_dg{d}{i} = [];
        RHO_dg{d}{i} = [];
        dsi_dg{d}{i} = [];
        for cc = 1:length(idx_dir{i})
            if ~dg_idx(idx_dir{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx_dir{i}(cc), :))>0
            [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
            xsort = xsort/pi*180;
            Y_TEMP = DG_cut{d}.RHO{T}(idx_dir{i}(cc), :);
            y_temp = DG_cut{d}.rho{T}(idx_dir{i}(cc), :);
%             if i == 4
%                 y_temp = circshift(y_temp, -1, 2);
%             end
%             y_temp = y_temp/max(y_temp);
%             plot(xsort, y_temp(seq), color(i))
%             ylim([0 1])
% %             pause
%             hold on
            rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
            RHO_dg{d}{i} = [RHO_dg{d}{i}; Y_TEMP(seq)];
            dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx_dir{i}(cc))];
            end
        end
        rho_dg_mean{d}(i, :) = mean(rho_dg{d}{i});
        rho_dg_ste{d}(i, :) = std(rho_dg{d}{i})/sqrt(size(rho_dg{d}{i}, 1));
        RHO_dg_mean{d}(i, :) = mean(RHO_dg{d}{i});
        RHO_dg_ste{d}(i, :) = std(RHO_dg{d}{i})/sqrt(size(RHO_dg{d}{i}, 1));
        dsi_dg_mean{d}(i) = mean(dsi_dg{d}{i});
        dsi_dg_ste{d}(i) = std(dsi_dg{d}{i})/sqrt(length(dsi_dg{d}{i}));
    end
%     xlabel('direction (rad)')
%     ylabel('normalized response')
%     title(ll{d})
%     xlim([-pi pi])
end
dsi_dg_mean = cell2mat(dsi_dg_mean');
dsi_dg_ste = cell2mat(dsi_dg_ste');
% plot average (cell type)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for d = 1:5
    subplot(2, 3, d)
    for i = 1:dirn
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(i));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
    title(ll{d});
end
legend(ct)

% plot average (light level)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:4
    subplot(2, 2, i)
    for d = 1:5
%         [f, g] = fit_cos(xsort/180*pi, rho_dg_mean{d}(i, :));
%         fitting{i, d} = f;
%         xfit = linspace(-160, 190, 100);
%         yfit = f.ymax * (0.5 + 0.5 * cos(xfit/180 *pi + f.phi)).^f.alpha + f.b;
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(d));
        hold on
%         plot(xfit, yfit, color(d))
    end
    xlabel('direction (degree)')
    ylabel('normalized average response')
    title(ct{i})
    ylim([0 1.1])
end
% legend(ll)
% DSI
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = ct;
model_series = dsi_dg_mean';
model_error = dsi_dg_ste';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('DSI')
legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

% dsi curve
marker = 'xosd';
figure
for i = 1:4
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
errorbar(0:4, model_series(i,:), model_error(i,:), 'Color', 'k', 'Marker', marker(i), 'MarkerSize', 10);
hold on
end
xlim([-0.5 4.5])
ylim([0 1.1])
legend(ct)

color = 'bkrgc';
figure
for i = 1:4
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
errorbar(0:4, model_series(i,:), model_error(i,:), 'Color', color(i));
hold on
end
xlim([-0.5 4.5])
ylim([0 1.1])
legend(ct)
xlabel('log(light intensity)')
ylabel('DSI')

% tuning curve of all directions
figure
j = 1;
for d = [2 5]
    xsort_all = xsort;
    subplot(2, 1, j)
    for i = 1:4
        [f, g] = fit_cos(xsort_all/180*pi, rho_dg_mean{d}(i, :));
        xfit = linspace(min(xsort_all), max(xsort_all), 100);
        yfit = f.ymax * (0.5 + 0.5 * cos(xfit/180 *pi + f.phi)).^f.alpha + f.b;

        errorbar(xsort_all, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), 'k');
        hold on
%         plot(xfit, yfit, 'k')
        xsort_all = xsort_all + 90;
    end
    xlim([-150 460])
    j = j + 1;
end


rho_dg_bgnd = rho_dg;
rho_dg_mean_bgnd = rho_dg_mean;
rho_dg_ste_bgnd = rho_dg_ste;
dsi_dg_wt_bgnd = dsi_dg;