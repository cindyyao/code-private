clear all
close all
cd /Users/xyao/matlab/code-private/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load data
datadg{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-06-03-0/data000-map/data000-map', opt);
datadg{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-06-03-0/stimuli/s00.mat';
datadg{1} = load_stim_matlab(datadg{1}, 'user_defined_trigger_interval', 10);
datadg{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-06-03-0/data003-map/data003-map', opt);
datadg{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-06-03-0/stimuli/s03.mat';
datadg{2} = load_stim_matlab(datadg{2}, 'user_defined_trigger_interval', 10);
datadg{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-06-03-0/data006-map/data006-map', opt);
datadg{3}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-06-03-0/stimuli/s06.mat';
datadg{3} = load_stim_matlab(datadg{3}, 'user_defined_trigger_interval', 10);
datadg{4} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-06-03-0/data009-map/data009-map', opt);
datadg{4}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-06-03-0/stimuli/s09.mat';
datadg{4} = load_stim_matlab(datadg{4}, 'user_defined_trigger_interval', 10, 'user_defined_trigger_interval_error', 0.6);
datadg{5} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-06-03-0/data012-map/data012-map', opt);
datadg{5}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-06-03-0/stimuli/s12.mat';
datadg{5} = load_stim_matlab(datadg{5}, 'user_defined_trigger_interval', 10);

%% drifting grating
n = 5;
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
close all
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [4 5]);


%% estimate background activity
interval = 1;
raster_interval_dg = deal(cell(n, 1));
for i = 1:n    
    raster_interval_dg{i} = get_ds_interval_raster(datadg{i}, ds_id, interval);
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

n = 5;
% spike count
DG_bgnd = cell(n, 1);
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim(datadg{i},ds_id,0, 0.05);
    DG_bgnd{i} = sort_direction(dscellanalysis_bgnd_subtract(NumSpikesCell, StimComb, bgfr{i}));
end
[DG_bgnd_cut, ~, ~] = cut_dg(DG_bgnd, raster_dg, raster_p_sum, 2, [4 5]);
close all

%% DS tuning curves (drifting grating)
dirn = 4;
T = 1;
color = 'bkrgc';
ct = {'superior', 'anterior', 'inferior', 'posterior'};
ll = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};
DG_cut = DG_bgnd_cut;

p_direction = DG{5}.angle{1}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

%subtypes
clear rho_dg_mean rho_dg_ste dsi_dg_mean dsi_dg_ste
for d = 1:5
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
end
dsi_dg_mean = cell2mat(dsi_dg_mean');
dsi_dg_ste = cell2mat(dsi_dg_ste');

%% fit bgnd subtracted tuning curve

for dir = 1:4
    figure
    for i = 1:5
        CC = 1;
        for cc = 1:length(id_dir{dir})
            xdata = DG_bgnd_cut{i}.theta{1}(idx_dir{dir}(cc), :);
            ydata = DG_bgnd_cut{i}.rho{1}(idx_dir{dir}(cc), :);
            if ~(dg_idx(idx_dir{dir}(cc), i)) && sum(ydata) ~= 0
                [f, g] = fit_cos_(xdata, ydata);
                Ymax{dir}{i}(CC) = f.ymax;
                Phi{dir}{i}(CC) = f.phi;
                Alpha{dir}{i}(CC) = f.alpha;
                Width{dir}{i}(CC) = acos(2*0.5^(1/f.alpha) - 1)/pi*360;
                B{dir}{i}(CC) = f.b;
                CC = CC + 1;
                
                xfit = linspace(0, 2*pi, 100);
                yfit = f.ymax.*((0.5+0.5*cos(xfit+f.phi)).^f.alpha.*f.b+1-f.b);
                subplot(4, 6, cc)
                plot(xdata, ydata, 'b')
                hold on
                plot(xfit, yfit, 'r')
                ylim([0 1.1])
                width = acos(2 * (0.5.^(1/f.alpha) - 0.5))/pi*360;
                title(['width = ', num2str(width)])


            end
        end
        Ymax_mean(dir, i) = mean(Ymax{dir}{i});
        Phi_mean(dir, i) = mean(Phi{dir}{i});
        Alpha_mean(dir, i) = mean(Alpha{dir}{i});
        B_mean(dir, i) = mean(B{dir}{i});
        B_ste(dir, i) = std(B{dir}{i})/sqrt(length(B{dir}{i}));
    end
end
close all

%% fig2A-D
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [4 5]);
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG_cut, raster_dg_cut, raster_p_sum_cut, 1, [1]);

cell_ids = [4352 5941 6197 7128];
[~, cell_idxs] = intersect(ds_id, cell_ids);
for cc = 1:4
    plot_ds_raster_one_old(DG_cut, raster_dg_cut, cell_idxs(cc), cell_ids(cc), 0)
end

%% fig2E-H
% plot average (light level)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:4
    subplot(2, 2, i)
    for d = 1:5
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(d));
        hold on
    end
    xlabel('direction (degree)')
    ylabel('normalized average response')
    title(ct{i})
    ylim([0 1.1])
end
legend(ll)

%% fig2K
% tuning curve of all directions
figure
j = 1;
for d = [2 5]
    xsort_all = xsort;
    subplot(2, 1, j)
    for i = 1:4
        errorbar(xsort_all, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), 'k');
        hold on
        xsort_all = xsort_all + 90;
    end
    xlim([-150 460])
    j = j + 1;
end


rho_dg_bgnd = rho_dg;
rho_dg_mean_bgnd = rho_dg_mean;
rho_dg_ste_bgnd = rho_dg_ste;

%% fig2I-J
load('DS150603.mat', 'dsi_dg')
dsi_150603 = dsi_dg;
load('DS160130.mat', 'dsi_dg')
dsi_160130 = dsi_dg;
ct = {'superior', 'anterior', 'inferior', 'posterior'};

for ll = 1:5
    for ct = 1:4
        dsi_all{ll}{ct} = [dsi_150603{ll}{ct}; dsi_160130{ll}{ct}];
        dsi_mean{ll}(ct) = mean(dsi_all{ll}{ct});
        dsi_ste{ll}(ct) = std(dsi_all{ll}{ct})/sqrt(length(dsi_all{ll}{ct}));
    end
end
dsi_mean = cell2mat(dsi_mean');
dsi_ste = cell2mat(dsi_ste');


xtick = ct;
model_series = dsi_mean';
model_error = dsi_ste';

% dsi curve
marker = 'xosd';
figure
subplot(2,1,1)
for i = 1:4
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
errorbar(0:4, model_series(i,:), model_error(i,:), 'Color', 'k', 'Marker', marker(i), 'MarkerSize', 10);
hold on
end
xlim([-0.5 4.5])
ylim([0 1.1])
legend('superior', 'anterior', 'inferior', 'posterior')


% width combine
load('DS150603.mat', 'Width')
Width_150603 = Width;
load('DS160130.mat', 'Width')
Width_160130 = Width;
for ll = 1:5
    for ct = 1:4
        Width_all{ct}{ll} = [Width_150603{ct}{ll} Width_160130{ct}{ll}];
        WidthMean(ct, ll) = mean(Width_all{ct}{ll});
        WidthSte(ct, ll) = std(Width_all{ct}{ll})/sqrt(length(Width_all{ct}{ll}));
    end
end

marker = 'xosd';
subplot(2,1,2)
for dir = 1:4
    errorbar(0:4, WidthMean(dir, :), WidthSte(dir, :), 'Color', 'k', 'Marker', marker(dir), 'MarkerSize', 10)
    hold on
end
legend('superior', 'anterior', 'inferior', 'posterior')
ylim([0 250])
xlabel('light level')
ylabel('tuning width (degree)')
xlim([-0.5 4.5])


%% figS2E
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
%% figS2F
ct = {'superior', 'anterior', 'inferior', 'posterior'};
marker = 'xosd';
figure
for i = 1:4
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
errorbar(0:4, B_mean(i, :), B_ste(i, :), 'Color', 'k', 'Marker', marker(i), 'MarkerSize', 10);
hold on
end
xlim([-0.5 4.5])
ylim([0 1.1])
legend(ct)
xlabel('light level')
ylabel('beta')
title('drifting grating fitting (fig2)')

%% figS2A-D
% 2015-06-03-0
load('DS150603.mat')
id = {[1607 3902 4998 6737], [6347], []};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);

dirn = 3;
D = 5;
T = 2;
p_direction = DG_cut{D}.angle{T}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

%subtypes
clear rho_dg dsi_dg rho_dg_mean rho_dg_ste dsi_dg_mean dsi_dg_ste
for d = 1:5
    for i = 1:dirn
        rho_dg{d}{i} = [];
        dsi_dg{d}{i} = [];
        if ~isempty(idx{i})
            for cc = 1:length(idx{i})
                if ~dg_idx(idx{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx{i}(cc), :))>0
                [xsort, seq] = sort(xx(idx{i}(cc), :));
                y_temp = DG_cut{d}.rho{T}(idx{i}(cc), :);
                rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
                dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx{i}(cc))];
                end
            end
        end
    end
end

clearvars -EXCEPT rho_dg dsi_dg dirn D T color
% 2015-06-18-0
load('DS150618.mat')
id = {[6766], [], [1952]};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);

p_direction = DG_cut{D}.angle{T}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

%subtypes
for d = 1:5
    for i = 1:dirn
        if ~isempty(idx{i})
            for cc = 1:length(idx{i})
                if ~dg_idx(idx{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx{i}(cc), :))>0
                [xsort, seq] = sort(xx(idx{i}(cc), :));
                y_temp = DG_cut{d}.rho{T}(idx{i}(cc), :);
                rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
                dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx{i}(cc))];
                end
            end
        end
    end
end

clearvars -EXCEPT rho_dg dsi_dg dirn D T color
% 2015-07-03-0
dsi_dg_temp = dsi_dg;
load('DS150703-1.mat')
dsi_dg = dsi_dg_temp;
id = {[261 4576 4923 5087 5239 6017], [2596 3391 4352 7397], [287 1653 6123]};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);

p_direction = DG_cut{D}.angle{T}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

%subtypes
for d = 1:5
    for i = 1:dirn
        if ~isempty(idx{i})
            for cc = 1:length(idx{i})
                if ~dg_idx(idx{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx{i}(cc), :))>0
                [xsort, seq] = sort(xx(idx{i}(cc), :));
                y_temp = DG_cut{d}.rho{T}(idx{i}(cc), :);
                rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
                dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx{i}(cc))];
                end
            end
        end
    end
end

clearvars -EXCEPT rho_dg dsi_dg dirn D T color

% 2016-01-30-0
dsi_dg_temp = dsi_dg;
load('DS160130.mat')
dsi_dg = dsi_dg_temp;
id = {[2582 3123 6917 7306], [4502], [392 4863]};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);

p_direction = DG_cut{D}.angle{T}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

%subtypes
for d = 1:5
    for i = 1:dirn
        if ~isempty(idx{i})
            for cc = 1:length(idx{i})
                if ~dg_idx(idx{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx{i}(cc), :))>0
                [xsort, seq] = sort(xx(idx{i}(cc), :));
                y_temp = DG_cut{d}.rho{T}(idx{i}(cc), :);
                rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
                dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx{i}(cc))];
                end
            end
        end
    end
end

clearvars -EXCEPT rho_dg dsi_dg dirn D T color xsort

% average
ct = {'superior', 'anterior', 'inferior'};
ll = {'NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0'};
for d = 1:5
    for i = 1:dirn
        rho_dg_mean{d}(i, :) = mean(rho_dg{d}{i});
        rho_dg_ste{d}(i, :) = std(rho_dg{d}{i})/sqrt(size(rho_dg{d}{i}, 1));
        dsi_dg_mean{d}(i) = mean(dsi_dg{d}{i});
        dsi_dg_ste{d}(i) = std(dsi_dg{d}{i})/sqrt(length(dsi_dg{d}{i}));
    end
end
dsi_dg_mean = cell2mat(dsi_dg_mean');
dsi_dg_ste = cell2mat(dsi_dg_ste');

% plot average (cell type)
for i = 1:5
    for dir = 1:3
        dsi_mean(dir, i) = mean(dsi_dg{i}{dir});
        dsi_ste(dir, i) = std(dsi_dg{i}{dir})/sqrt(length(dsi_dg{i}{dir}));
    end
end

% plot average (light level)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:dirn
    subplot(2, 2, i)
    for d = 1:5
        errorbar(xsort/pi*180, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(d));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
    title(ct{i})
    xlim([-200 200])
end
legend(ll)

marker = 'xo*d';
subplot(2, 2, 4)
for dir = 1:3
    errorbar([0:4], dsi_mean(dir, :), dsi_ste(dir, :), ['k-' marker(dir)], 'MarkerSize', 10);
    hold on
end
ylim([0 1.1])
xlim([-0.5 4.5])
legend(ct)
