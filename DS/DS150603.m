% addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/DS' cell analysis'/
cd /Users/xyao/matlab/code-private/

opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load data
datadg{1} = load_data('/Volumes/lab/analysis/2015-06-03-0/data000-map/data000-map', opt);
datadg{1}.names.stimulus_path = '/Volumes/lab/analysis/2015-06-03-0/stimuli/s00.mat';
datadg{1} = load_stim_matlab(datadg{1}, 'user_defined_trigger_interval', 10);
datadg{2} = load_data('/Volumes/lab/analysis/2015-06-03-0/data003-map/data003-map', opt);
datadg{2}.names.stimulus_path = '/Volumes/lab/analysis/2015-06-03-0/stimuli/s03.mat';
datadg{2} = load_stim_matlab(datadg{2}, 'user_defined_trigger_interval', 10);
datadg{3} = load_data('/Volumes/lab/analysis/2015-06-03-0/data006-map/data006-map', opt);
datadg{3}.names.stimulus_path = '/Volumes/lab/analysis/2015-06-03-0/stimuli/s06.mat';
datadg{3} = load_stim_matlab(datadg{3}, 'user_defined_trigger_interval', 10);
datadg{4} = load_data('/Volumes/lab/analysis/2015-06-03-0/data009-map/data009-map', opt);
datadg{4}.names.stimulus_path = '/Volumes/lab/analysis/2015-06-03-0/stimuli/s09.mat';
datadg{4} = load_stim_matlab(datadg{4}, 'user_defined_trigger_interval', 10, 'user_defined_trigger_interval_error', 0.6);
datadg{5} = load_data('/Volumes/lab/analysis/2015-06-03-0/data012-map/data012-map', opt);
datadg{5}.names.stimulus_path = '/Volumes/lab/analysis/2015-06-03-0/stimuli/s12.mat';
datadg{5} = load_stim_matlab(datadg{5}, 'user_defined_trigger_interval', 10);

datamb{1} = load_data('/Volumes/lab/analysis/2015-06-03-0/data001-map/data001-map', opt);
datamb{1}.names.stimulus_path = '/Volumes/lab/analysis/2015-06-03-0/stimuli/s01.mat';
datamb{1} = load_stim_matlab(datamb{1});
datamb{2} = load_data('/Volumes/lab/analysis/2015-06-03-0/data004-map/data004-map', opt);
datamb{2}.names.stimulus_path = '/Volumes/lab/analysis/2015-06-03-0/stimuli/s04.mat';
datamb{2} = load_stim_matlab(datamb{2});
datamb{2} = delete_last_repeat(datamb{2});
datamb{3} = load_data('/Volumes/lab/analysis/2015-06-03-0/data007-map/data007-map', opt);
datamb{3}.names.stimulus_path = '/Volumes/lab/analysis/2015-06-03-0/stimuli/s07.mat';
datamb{3} = load_stim_matlab(datamb{3});
datamb{4} = load_data('/Volumes/lab/analysis/2015-06-03-0/data010-map/data010-map', opt);
datamb{4}.names.stimulus_path = '/Volumes/lab/analysis/2015-06-03-0/stimuli/s10.mat';
datamb{4} = load_stim_matlab(datamb{4});
datamb{5} = load_data('/Volumes/lab/analysis/2015-06-03-0/data013-map/data013-map', opt);
datamb{5}.names.stimulus_path = '/Volumes/lab/analysis/2015-06-03-0/stimuli/s13.mat';
datamb{5} = load_stim_matlab(datamb{5});

dataffp{1} = load_data('/Volumes/lab/analysis/2015-06-03-0/data002-map/data002-map', opt);
dataffp{1}.triggers = dataffp{1}.triggers(2:end);
dataffp{2} = load_data('/Volumes/lab/analysis/2015-06-03-0/data005-map/data005-map', opt);
dataffp{2}.triggers = dataffp{2}.triggers(2:end);
dataffp{3} = load_data('/Volumes/lab/analysis/2015-06-03-0/data008-map/data008-map', opt);
dataffp{3}.triggers = dataffp{3}.triggers(2:end);
dataffp{4} = load_data('/Volumes/lab/analysis/2015-06-03-0/data011-map/data011-map', opt);
dataffp{4}.triggers = dataffp{4}.triggers(2:end);
dataffp{5} = load_data('/Volumes/lab/analysis/2015-06-03-0/data014-map/data014-map', opt);
dataffp{5}.triggers = dataffp{5}.triggers(2:end);

%%

n = 5;
i = 5;
duration = 573;
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg{i},datadg{i}.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);

params_idx = [3 4]; % which parameters to use for classification
[ds_id, nonds_id] = classify_ds(datadg{i}, ds_struct, params_idx);

%% drifting grating

load('DS150603.mat')
% spike count
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim(datadg{i},ds_id,0, 0.05);
    DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster_dg{i} = get_ds_raster(datadg{i}, ds_id);
    for j = 1:length(raster_dg{i})
        if(dg_idx(j, i))
            raster_dg{i}{j} = [];
        end
    end
end

% max firing rate
% [raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
% for i = 1:n    
%     [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim(datadg{i},ds_id,0, 0.05);
%     DG{i} = sort_direction(dscellanalysis(MaxRate, StimComb));
%     raster_dg{i} = get_ds_raster(datadg{i}, ds_id);
%     for j = 1:length(raster_dg{i})
%         if(dg_idx(j, i))
%             raster_dg{i}{j} = [];
%         end
%     end
% end

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

%% plot cell summary
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [4 5]);

for cc = 3:3 %length(ds_id)
    plot_ds_raster(DG_cut, raster_dg_cut, cc, ds_id(cc), ll, 2, 3, 0)
end

% plot all light levels in one figure
for cc = 52:52 %length(ds_id)
    plot_ds_raster_one(DG_cut, raster_dg_cut, cc, ds_id(cc), {num2str(ds_id(cc))}, 1, 1, 0)
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
%pca
mag_pca = MAG_all_norm_dg{2};
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
title('NDF 3')

% several cells are missing during mapping between NDF0 grating and NDF3
% grating, simply put them into ON-OFF DSGC class here
idx_all = 1:length(ds_id);
idx_sub{2} = setdiff(idx_all, idx_sub{1});
id_sub{2} = ds_id(idx_sub{2});

figure
plot(scores(idx_sub{1}, pc1), scores(idx_sub{1}, pc2), 'ro')
hold on
plot(scores(idx_sub{2}, pc1), scores(idx_sub{2}, pc2), 'bo')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
title('NDF 3')


figure
v = 4*datadg{2}.stimulus.params.SPATIAL_PERIOD./datadg{2}.stimulus.params.TEMPORAL_PERIOD;
subplot(1, 2, 1)
semilogx(v, exciseColumn(MAG_all_norm_dg{2}(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg{2}(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
title(ll{2})
xlim([v(end) v(1)])

subplot(1, 2, 2)
semilogx(v, exciseColumn(MAG_all_norm_dg{5}(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg{5}(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
title(ll{5})
xlim([v(end) v(1)])

 
figure
compass(DG{5}.U{3}(idx_sub{1}), DG{5}.V{3}(idx_sub{1}), 'r')
hold on
compass(DG{5}.U{3}(idx_sub{2}), DG{5}.V{3}(idx_sub{2}), 'b')

%% plot average tunning curve
color = 'brk';
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
idx_temp = [2 5];
for i = 1:2
    v = 4*datadg{idx_temp(i)}.stimulus.params.SPATIAL_PERIOD./datadg{idx_temp(i)}.stimulus.params.TEMPORAL_PERIOD;
    subplot(1, 2, i)
    for ct = 1:2
        mag_temp = exciseColumn(MAG_all_norm_dg{idx_temp(i)}(:, idx_sub{ct}));
        tuning_avg{i}(:, ct) = mean(mag_temp, 2);
        tuning_ste{i}(:, ct) = std(mag_temp, [], 2)/sqrt(size(mag_temp, 2));
        errorbar(v, tuning_avg{i}(:, ct), tuning_ste{i}(:, ct), color(ct))
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
set(gcf, 'DefaultLineLineWidth', 1.5)
v = 4*datadg{2}.stimulus.params.SPATIAL_PERIOD./datadg{2}.stimulus.params.TEMPORAL_PERIOD;
for j = 1:2
    subplot(1, 2, j)
    errorbar(v, tuning_avg{1}(:, j), tuning_ste{1}(:, j), 'b')
    hold on
    errorbar(v, tuning_avg{2}(:, j), tuning_ste{2}(:, j), 'r')
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
color = 'bkrgc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}(idx_sub{2}), DG{d}.V{t}(idx_sub{2}), x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = idx_sub{2}(I);
    id_dir{i} = ds_id(idx_dir{i});
end

h = figure;
dirn = 2;
set(h, 'Position', [1 1 1080 500])
compass(DG{d}.U{t}(idx_sub{1}), DG{d}.V{t}(idx_sub{1}));
color = 'bkrgc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}(idx_sub{1}), DG{d}.V{t}(idx_sub{1}), x, y);
    [~, I] = find(IN == 1);
    idx_dir_on{i} = idx_sub{1}(I);
    id_dir_on{i} = ds_id(idx_dir_on{i});
end

%% classification based on MB response polarity

[raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},ds_id,duration,1);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb));
    raster_mb{i} = get_mb_raster(datamb{i}, ds_id, duration);
    for j = 1:length(raster_mb{i})
        if(mb_idx(j, i))
            raster_mb{i}{j} = [];
        end
    end
    trial_dur{i} = get_mb_trial_dur(datamb{i});
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
ct = {'superior', 'anterior', 'inferior', 'posterior'};

%% plot cell summary

for cc = 2:2 %length(ds_id)
    plot_mb_raster_one(MB, raster_mb, trial_dur, cc, ds_id(cc), ll, 1, 1, 0)
end


%% ffp

n_ffp = 5;
bin_size = 0.1;
XX = bin_size/2:bin_size:6-bin_size/2;

[raster_ff, raster_ff_all] = deal(cell(n_ffp, 1));
for d = 1:n_ffp
    [raster_ff{d}, raster_ff_all{d}] = get_ffp_raster(dataffp{d}, ds_id, 3);
    for j = 1:length(raster_ff{d})
        if(ffp_idx(j, d))
            raster_ff{d}{j} = [];
            raster_ff_all{d}{j} = [];
        end
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

%% Direction tuning (moving bar)
% all ds cells

% t = 2;
dirn = 4;
D = 5;
T = 2;
color = 'bkrgc';
ct = {'superior', 'anterior', 'inferior', 'posterior'};

figure
for d = 1:n
    p_direction = MB{D}.angle{T}';
    xx = 0:pi/4:7*pi/4;
    xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
    xx(xx>pi) = xx(xx>pi)-2*pi;
    xx(xx<-pi) = xx(xx<-pi)+2*pi;


    subplot(2, 3, d)
    for i = 1:dirn
        for cc = 1:length(ds_id)
            if ~mb_idx(cc, d)
            [xsort, seq] = sort(xx(cc, :));
            y_temp = MB{d}.rho{T}(cc, :);
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
clear rho_mb_mean rho_mb_ste dsi_mb_mean dsi_mb_ste
figure
for d = 1:5
    subplot(2, 3, d)
    for i = 1:dirn
        rho_mb{d}{i} = [];
        dsi_mb{d}{i} = [];
        for cc = 1:length(idx_dir{i})
            if ~mb_idx(idx_dir{i}(cc), d) && sum(MB{d}.rho{T}(idx_dir{i}(cc), :))>0
            [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
            y_temp = MB{d}.rho{T}(idx_dir{i}(cc), :);
            plot(xsort, y_temp(seq), color(i))
            ylim([0 1])
%             pause
            hold on
            rho_mb{d}{i} = [rho_mb{d}{i}; y_temp(seq)];
            dsi_mb{d}{i} = [dsi_mb{d}{i}; MB{d}.dsindex{T}(idx_dir{i}(cc))];
            end
        end
        rho_mb_mean{d}(i, :) = mean(rho_mb{d}{i});
        rho_mb_ste{d}(i, :) = std(rho_mb{d}{i})/sqrt(size(rho_mb{d}{i}, 1));
        dsi_mb_mean{d}(i) = mean(dsi_mb{d}{i});
        dsi_mb_ste{d}(i) = std(dsi_mb{d}{i})/sqrt(length(dsi_mb{d}{i}));
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    title(ll{d})
    xlim([-pi pi])
end
dsi_mb_mean = cell2mat(dsi_mb_mean');
dsi_mb_ste = cell2mat(dsi_mb_ste');
% plot average (cell type)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for d = 1:5
    subplot(2, 3, d)
    for i = 1:4
        errorbar(xsort/pi*180, rho_mb_mean{d}(i, :), rho_mb_ste{d}(i, :), color(i));
        hold on
    end
    xlabel('degrees')
    ylabel('normalized mean spike number')
    title(ll{d});
end
legend(ct)

% plot average (light level)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:4
    subplot(2, 2, i)
    for d = 1:5
        errorbar(xsort/pi*180, rho_mb_mean{d}(i, :), rho_mb_ste{d}(i, :), color(d));
        hold on
    end
    ylim([0 1])
    xlabel('degrees')
    ylabel('normalized average response')
    title(ct{i})
end
legend(ll)
% DSI
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = ct;
model_series = [dsi_mb_mean(1,1) dsi_mb_mean(2,1) dsi_mb_mean(3,1) dsi_mb_mean(4,1) dsi_mb_mean(5,1); dsi_mb_mean(1,2) dsi_mb_mean(2,2) dsi_mb_mean(3,2) dsi_mb_mean(4,2) dsi_mb_mean(5,2);dsi_mb_mean(1,3) dsi_mb_mean(2,3) dsi_mb_mean(3,3) dsi_mb_mean(4,3) dsi_mb_mean(5,3); dsi_mb_mean(1,4) dsi_mb_mean(2,4) dsi_mb_mean(3,4) dsi_mb_mean(4,4) dsi_mb_mean(5,4)];   
model_error = [dsi_mb_ste(1,1) dsi_mb_ste(2,1) dsi_mb_ste(3,1) dsi_mb_ste(4,1) dsi_mb_ste(5,1); dsi_mb_ste(1,2) dsi_mb_ste(2,2) dsi_mb_ste(3,2) dsi_mb_ste(4,2) dsi_mb_ste(5,2);dsi_mb_ste(1,3) dsi_mb_ste(2,3) dsi_mb_ste(3,3) dsi_mb_ste(4,3) dsi_mb_ste(5,3); dsi_mb_ste(1,4) dsi_mb_ste(2,4) dsi_mb_ste(3,4) dsi_mb_ste(4,4) dsi_mb_ste(5,4)];
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

%% DS tuning curves (drifting grating)
% all ds cells

% t = 2;
dirn = 4;
D = 5;
T = 1;
color = 'bkrgc';
ct = {'superior', 'anterior', 'inferior', 'posterior'};

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
    subplot(2, 3, d)
    for i = 1:dirn
        rho_dg{d}{i} = [];
        dsi_dg{d}{i} = [];
        for cc = 1:length(idx_dir{i})
            if ~dg_idx(idx_dir{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx_dir{i}(cc), :))>0
            [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
            xsort = xsort/pi*180;
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
    for i = 1:4
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
        [f, g] = fit_cos(xsort/180*pi, rho_dg_mean{d}(i, :));
        xfit = linspace(-160, 190, 100);
        yfit = f.ymax * (0.5 + 0.5 * cos(xfit/180 *pi + f.phi)).^f.alpha + f.b;
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), [color(d) 'o']);
        hold on
        plot(xfit, yfit, color(d))
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

color = 'rbgk';
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
for d = [2 5]
    xsort_all = xsort;
    figure
    for i = 1:4
        [f, g] = fit_cos(xsort_all/180*pi, rho_dg_mean{d}(i, :));
        xfit = linspace(min(xsort_all), max(xsort_all), 100);
        yfit = f.ymax * (0.5 + 0.5 * cos(xfit/180 *pi + f.phi)).^f.alpha + f.b;

        errorbar(xsort_all, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), 'ko');
        hold on
        plot(xfit, yfit, 'k')
        xsort_all = xsort_all + 90;
    end
    xlim([-150 460])
end
%% frequency analysis
duration = 8;
bin_rate = 10000;
hist_spikes = cell(n, 1);
signal_length = duration*bin_rate;                
NFFT = 2^nextpow2(signal_length);
f = bin_rate/2*linspace(0,1,NFFT/2+1);
fft_spikes = cell(n, 1);
[DC, F1, F2] = deal(cell(n, 1));

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
        F1{i}(rgc,:) = fund_power;
        F2{i}(rgc,:) = sec_power;
        DC{i}(rgc,:) = DC_power;

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

% p value
for time = 1:length(tp)
    for i = 1:n
        for j = 1:n
            [~, p] = ttest2(ratio_oo{i}(:, time), ratio_oo{j}(:, time));
            P_all(i,j,time) = p;
            for ct = 1:4
                [~, p] = ttest2(ratio_dir{ct}{i}(:, time), ratio_dir{ct}{j}(:, time));
                P_dir{ct}(i,j,time) = p;
            end
        end
    end
end
        
% plot average f2/f1
ct = {'superior', 'anterior', 'inferior', 'posterior'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = {'high speed'; 'low speed'};
model_series = [ratio_oo_mean(1,1) ratio_oo_mean(2,1) ratio_oo_mean(3,1) ratio_oo_mean(4,1) ratio_oo_mean(5,1); ratio_oo_mean(1,2) ratio_oo_mean(2,2) ratio_oo_mean(3,2) ratio_oo_mean(4,2) ratio_oo_mean(5,2)];   
model_error = [ratio_oo_ste(1,1) ratio_oo_ste(2,1) ratio_oo_ste(3,1) ratio_oo_ste(4,1) ratio_oo_ste(5,1); ratio_oo_ste(1,2) ratio_oo_ste(2,2) ratio_oo_ste(3,2) ratio_oo_ste(4,2) ratio_oo_ste(5,2)];
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

% plot f2/f1 curve with errorbar
% low speed (tp = 4)
figure
for ct = 1:4
    errorbar(ratio_dir_mean(:, ct, 2), ratio_dir_ste(:, ct, 2), '--')
    hold on 
end
errorbar(ratio_oo_mean(:, 2), ratio_oo_ste(:, 2))

legend('S', 'A', 'I', 'P', 'all')
ylim([0 3])
xlim([0.5 5.5])

%% frequency analysis across speed
clear ratio ratio_dir ratio_dir_mean ratio_dir_ste ratio_oo ratio_oo_mean ratio_oo_ste
duration = 8;
bin_rate = 10000;
hist_spikes = cell(n, 1);
signal_length = duration*bin_rate;                
NFFT = 2^nextpow2(signal_length);
f = bin_rate/2*linspace(0,1,NFFT/2+1);
fft_spikes = cell(n, 1);
[DC, F1, F2] = deal(cell(n, 1));

for i = [2 5]
    tp = datadg{i}.stimulus.params.TEMPORAL_PERIOD(1:6);
    [DC{i}, F1{i}, F2{i}] = deal(zeros(length(ds_id), 6));
    for rgc = 1:length(ds_id)
        if ~isempty(raster_dg{i}{rgc}) && ~dg_idx(rgc, i)
        for time = 1:length(tp)
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
        F1{i}(rgc,:) = fund_power;
        F2{i}(rgc,:) = sec_power;
        DC{i}(rgc,:) = DC_power;

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

% p value

for time = 1:length(tp)
    [~, p] = ttest2(ratio_oo{2}(:, time), ratio_oo{5}(:, time));
    P_all(time) = p;
    for ct = 1:4
        [~, p] = ttest2(ratio_dir{ct}{2}(:, time), ratio_dir{ct}{5}(:, time));
        P_dir{ct}(i,j,time) = p;
    end
end
        
% plot f2/f1 curve with errorbar
% low speed (tp = 4)
sp = datadg{5}.stimulus.params.SPATIAL_PERIOD;
speed = sp./tp*4;
figure
errorbar(speed, ratio_oo_mean(2, :), ratio_oo_ste(2, :))
hold on
errorbar(speed, ratio_oo_mean(5, :), ratio_oo_ste(5, :))
set(gca, 'XScale', 'log')
legend('rod', 'cone')

% ylim([0 3])
xlim([100 6000])

%% plot cell type specific f2/f1 

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
ct = {'superior', 'posterior', 'inferior', 'anterior'};
xtick = ['all', ct];
model_series = [ratio_oo_mean(:, 1) ratio_dir_mean(:, :, 1)]';   
model_error = [ratio_oo_ste(:, 1) ratio_dir_ste(:, :, 1)]';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('F2/F1')
legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
title('high speed')
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
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = ['all', ct];
model_series = [ratio_oo_mean(:, 2) ratio_dir_mean(:, :, 2)]';   
model_error = [ratio_oo_ste(:, 2) ratio_dir_ste(:, :, 2)]';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('F2/F1')
legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
hold on;
title('low speed')

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

%% frequency doubling of moving bar response

% get spike histgram of preferred direction
bin_size = 0.2;
for i = 1:n
% for i = 4:4
    for t = 1:2
        dur = max(trial_dur{i}(t));
        XX = bin_size/2:bin_size:dur;
        hist_mb_p{i}{t} = zeros(length(ds_id), length(XX));
        for cc = 1:length(ds_id)
            if ~mb_idx(cc, i) && ~isempty(raster_p_sum_mb{i}{cc})
                hist_mb_p{i}{t}(cc, :) = hist(raster_p_sum_mb{i}{cc}{t}, XX);
            end
        end
    end
end

x = [4 5 3 3]; y = [5 5 3 4];
t = 2;
for ct = 1:4
    figure(ct)
    for cc = 1:length(idx_dir{ct})
        subplot(x(ct), y(ct), cc)
%         subplot(3, 4, cc-3)
        dur = max(trial_dur{1}(t));
        XX = bin_size/2:bin_size:dur;
%         for i = 4:4
        for i = 1:n
            plot(XX, 1.25*hist_mb_p{i}{t}(idx_dir{ct}(cc), :), color(i))
            hold on
        end
        title(num2str(ds_id(idx_dir{ct}(cc))))
    end
    legend('NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0')
end
            
    
%%c F2/F1 ratio
delta = datamb{1}.stimulus.params.DELTA(t);
bar_width = datamb{1}.stimulus.params.BAR_WIDTH;
refresh_rate = 60.3578;
bar_time = bar_width/delta/refresh_rate;
range = 1/3;
for ct = 1:4
    for cc = 1:length(idx_dir{ct})
        figure(1)
        dur = max(trial_dur{1}(t));
        XX = bin_size/2:bin_size:dur;
        for i = 1:n
            plot(XX, hist_mb_p{i}{t}(idx_dir{ct}(cc), :), color(i))
            hold on
        end
        title(num2str(ds_id(idx_dir{ct}(cc))))
        [x, ~] = ginput;
        on_time{ct}(cc, t) = x;
        left_on = max(round((on_time{ct}(cc, t)-bar_time*range)/bin_size), 1);
        right_on = min(round((on_time{ct}(cc, t)+bar_time*range)/bin_size), size(hist_mb_p{1}{t}, 2));
        left_off = max(round((on_time{ct}(cc, t)+bar_time*(1-range))/bin_size), 1);
        right_off = min(round((on_time{ct}(cc, t)+bar_time*(1+range))/bin_size), size(hist_mb_p{1}{t}, 2));
        for i = 1:n
            max_onoff{i}{ct}{t}(cc, 1) = max(hist_mb_p{i}{t}(idx_dir{ct}(cc), left_on:right_on));
            max_onoff{i}{ct}{t}(cc, 2) = max(hist_mb_p{i}{t}(idx_dir{ct}(cc), left_off:right_off));
        end
        hold off
    end
    legend('NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0')
end

T = t;
% clear ratio_ste
for i = 1:5
    for ct = 1:4
        for t = T:T
            for cc = 1:size(max_onoff{i}{ct}{t}, 1)
                ratio_onoff{ct}{t}{i}(cc) = max_onoff{i}{ct}{t}(cc, 2)/max_onoff{i}{ct}{t}(cc, 1);
            end
            index = logical(isnan(ratio_onoff{ct}{t}{i}) + isinf(ratio_onoff{ct}{t}{i}));
            ratio_onoff{ct}{t}{i}(index) = [];
            ratio_avg{ct}{t}(i) = mean(ratio_onoff{ct}{t}{i});
            ratio_oo_ste{ct}{t}(i) = std(ratio_onoff{ct}{t}{i})/sqrt(length(ratio_onoff{ct}{t}{i}));
        end
    end
end

%%
t = 2;
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = {'superior'; 'posterior'; 'inferior'; 'anterior'};
model_series = [ratio_avg{1}{t}(1) ratio_avg{1}{t}(2) ratio_avg{1}{t}(3) ratio_avg{1}{t}(4) ratio_avg{1}{t}(5); ratio_avg{2}{t}(1) ratio_avg{2}{t}(2) ratio_avg{2}{t}(3) ratio_avg{2}{t}(4) ratio_avg{2}{t}(5); ratio_avg{3}{t}(1) ratio_avg{3}{t}(2) ratio_avg{3}{t}(3) ratio_avg{3}{t}(4) ratio_avg{3}{t}(5); ratio_avg{4}{t}(1) ratio_avg{4}{t}(2) ratio_avg{4}{t}(3) ratio_avg{4}{t}(4) ratio_avg{4}{t}(5)];   
model_error = [ratio_oo_ste{1}{t}(1) ratio_oo_ste{1}{t}(2) ratio_oo_ste{1}{t}(3) ratio_oo_ste{1}{t}(4) ratio_oo_ste{1}{t}(5); ratio_oo_ste{2}{t}(1) ratio_oo_ste{2}{t}(2) ratio_oo_ste{2}{t}(3) ratio_oo_ste{2}{t}(4) ratio_oo_ste{2}{t}(5); ratio_oo_ste{3}{t}(1) ratio_oo_ste{3}{t}(2) ratio_oo_ste{3}{t}(3) ratio_oo_ste{3}{t}(4) ratio_oo_ste{3}{t}(5); ratio_oo_ste{4}{t}(1) ratio_oo_ste{4}{t}(2) ratio_oo_ste{4}{t}(3) ratio_oo_ste{4}{t}(4) ratio_oo_ste{4}{t}(5)];   
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('OFF/ON')
legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
hold on;
title('speed 2')

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

%% fit mb direction tuning

xdata = MB{1}.theta{1}(1, :);
for i = 1:length(rho_mb)
    for j = 1:length(rho_mb{i})
        for cc = 1:size(rho_mb{i}{j}, 1)
            [dst_fit{i}{j}{cc}{1}, dst_fit{i}{j}{cc}{2}] = fit_gaussian(xdata, rho_mb{i}{j}(cc, :));
        end
    end
end

for i = 1:length(rho_mb)
    rmse = [];
    for j = 1:length(rho_mb{i})
        for cc = 1:size(rho_mb{i}{j}, 1)
            rmse = [rmse dst_fit{i}{j}{cc}{2}.rmse];
        end
    end
    RMSE{i} = rmse;
end

figure
for i = 1:5
    subplot(2, 3, i)
    hist(RMSE{i})
    xlim([0 .4])
end

x = 0:0.01:6;
figure
for i = 1:n
    figure
    celln = sum(cellfun(@length, dst_fit{i}(:)));
    temp = cell2mat(rho_mb{i}');
    temp_fit = [dst_fit{i}{:}];
    for cc = 1:celln
        subplot(7, 8, cc)
        plot(xdata, temp(cc, :), 'o')
        hold on
        f = temp_fit{cc}{1};
        plot(x, f.a*exp(-((x-f.b)/f.c).^2)+f.d)
        ylim([0 1.2])
        xlim([0 6])
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
    h = figure;
    set(h, 'Position', [1,1,900 450])
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
    h = figure;
    set(h, 'Position', [1,1,900 450])
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
    h = figure;
    set(h, 'Position', [1,1,900 450])
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
    h = figure;
    set(h, 'Position', [1,1,900 400])
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

%% get spontaneous activity

for ll = 1:length(datadg)
    for cc = 1:length(ds_id)
        if ~dg_idx(cc,ll)
            idx = get_cell_indices(datadg{ll},ds_id(cc));
            spikes_temp = datadg{ll}.spikes{idx};
            if ismember(ll, [2,5])
                bgnd_firing{ll}{cc} = spikes_temp(spikes_temp > 2888)-2888;
            else
                bgnd_firing{ll}{cc} = spikes_temp(spikes_temp > 648)-648;
            end
        end
    end
end

%% fit tuning curves


for dir = 1:4
    for i = 1:5
        CC = 1;
        for cc = 1:length(id_dir{dir})
            if ~(dg_idx(idx_dir{dir}(cc), i))
                xdata = DG_cut{i}.theta{1}(idx_dir{dir}(cc), :);
                ydata = DG_cut{i}.RHO{1}(idx_dir{dir}(cc), :);
                [f, g] = fit_cos(xdata, ydata);
                Ymax{dir}{i}(CC) = f.ymax;
                Phi{dir}{i}(CC) = f.phi;
                Alpha{dir}{i}(CC) = f.alpha;
                Width{dir}{i}(CC) = acos(2*0.5^(1/f.alpha) - 1)/pi*360;
                B{dir}{i}(CC) = f.b;
                CC = CC + 1;
                
%                 xfit = linspace(0, 2*pi, 100);
%                 yfit = f.ymax.*(0.5+0.5*cos(xfit+f.phi)).^f.alpha+f.b;
%                 subplot(5, length(id_dir{dir}), length(id_dir{dir})*(i-1) + cc)
%                 plot(xdata, ydata, 'b')
%                 hold on
%                 plot(xfit, yfit, 'r')


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
ylim([0 250])
xlabel('light level')
ylabel('tuning width (degree)')
xlim([-0.5 4.5])

