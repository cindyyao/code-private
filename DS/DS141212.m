%% load data
addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/Classification/
addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/DS' cell analysis'/

opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load data
datadg{1} = load_data('/Analysis/xyao/2014-12-12-0/data001-map/data001-map', opt);
datadg{1}.names.stimulus_path = '/Analysis/xyao/2014-12-12-0/stimuli/s01.mat';
datadg{1} = load_stim_matlab(datadg{1}, 'user_defined_trigger_interval', 10);
datadg{2} = load_data('/Analysis/xyao/2014-12-12-0/data004-map/data004-map', opt);
datadg{2}.names.stimulus_path = '/Analysis/xyao/2014-12-12-0/stimuli/s04.mat';
datadg{2} = load_stim_matlab(datadg{2}, 'user_defined_trigger_interval', 10);
datadg{3} = load_data('/Analysis/xyao/2014-12-12-0/data007-map/data007-map', opt);
datadg{3}.names.stimulus_path = '/Analysis/xyao/2014-12-12-0/stimuli/s07.mat';
datadg{3} = load_stim_matlab(datadg{3}, 'user_defined_trigger_interval', 10);
datadg{4} = load_data('/Analysis/xyao/2014-12-12-0/data010-map/data010-map', opt);
datadg{4}.names.stimulus_path = '/Analysis/xyao/2014-12-12-0/stimuli/s10.mat';
datadg{4} = load_stim_matlab(datadg{4}, 'user_defined_trigger_interval', 10);
datadg{5} = load_data('/Analysis/xyao/2014-12-12-0/data013-map/data013-map', opt);
datadg{5}.names.stimulus_path = '/Analysis/xyao/2014-12-12-0/stimuli/s13.mat';
datadg{5} = load_stim_matlab(datadg{5}, 'user_defined_trigger_interval', 10);

datamb{1} = load_data('/Analysis/xyao/2014-12-12-0/data000-map/data000-map', opt);
datamb{1}.names.stimulus_path = '/Analysis/xyao/2014-12-12-0/stimuli/s00.mat';
datamb{1} = load_stim_matlab(datamb{1});
datamb{2} = load_data('/Analysis/xyao/2014-12-12-0/data003-map/data003-map', opt);
datamb{2}.names.stimulus_path = '/Analysis/xyao/2014-12-12-0/stimuli/s03.mat';
datamb{2} = load_stim_matlab(datamb{2});
datamb{3} = load_data('/Analysis/xyao/2014-12-12-0/data006-map/data006-map', opt);
datamb{3}.names.stimulus_path = '/Analysis/xyao/2014-12-12-0/stimuli/s06.mat';
datamb{3} = load_stim_matlab(datamb{3});
datamb{4} = load_data('/Analysis/xyao/2014-12-12-0/data009-map/data009-map', opt);
datamb{4}.names.stimulus_path = '/Analysis/xyao/2014-12-12-0/stimuli/s09.mat';
datamb{4} = load_stim_matlab(datamb{4});
datamb{5} = load_data('/Analysis/xyao/2014-12-12-0/data012-map/data012-map', opt);
datamb{5}.names.stimulus_path = '/Analysis/xyao/2014-12-12-0/stimuli/s12.mat';
datamb{5} = load_stim_matlab(datamb{5});

dataffp{1} = load_data('/Analysis/xyao/2014-12-12-0/data002-map/data002-map', opt);
dataffp{1}.triggers = dataffp{1}.triggers(2:end);
dataffp{2} = load_data('/Analysis/xyao/2014-12-12-0/data005-map/data005-map', opt);
dataffp{2}.triggers = dataffp{2}.triggers(2:end);
dataffp{3} = load_data('/Analysis/xyao/2014-12-12-0/data008-map/data008-map', opt);
dataffp{3}.triggers = dataffp{3}.triggers(2:end);
dataffp{4} = load_data('/Analysis/xyao/2014-12-12-0/data011-map/data011-map', opt);
dataffp{4}.triggers = dataffp{4}.triggers(2:end);
dataffp{5} = load_data('/Analysis/xyao/2014-12-12-0/data014-map/data014-map', opt);
dataffp{5}.triggers = dataffp{5}.triggers(2:end);

%% Moving Bar

% bin_size = [1];
% for i = 1:length(bin_size)
%     [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datarun{2},datarun{2}.cell_ids,0, bin_size(i));
%     mb_struct = mbcellanalysis(MaxRate, StimComb);
%     title(['bin size ' num2str(bin_size(i))])
% end
i = 4;
[NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},datamb{i}.cell_ids,0, 1);
mb_struct = mbcellanalysis(NumSpikesCell, StimComb);

n = 5;
duration = 871;
% pull out DS cells

figure
plot(mb_struct.mag{1, 1}, mb_struct.mag{2, 1}, 'o')
title('ndf1-map')
xlabel('TP 1')
ylabel('TP 2')
hold on
[x, y] = ginput;
plot(x, y);
IN = inpolygon(mb_struct.mag{1, 1}, mb_struct.mag{2, 1}, x, y);
[~, I] = find(IN == 1);
id_init = datamb{i}.cell_ids(I);

[C ia ib] = intersect(id_init, datamb{i}.cell_ids);
vc = ones(length(datamb{i}.cell_ids),1);
vc(ib) = 2; %initializing ds cells to cluster 2, everything else cluster 1

close all;
X = [];
N = [];
p = [];
X(:,1) = log(mb_struct.mag{1,1})';
X(:,2) = log(mb_struct.mag{2,1})';
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 0, 1, 0, 0, 0,0, vc);

ds_id = [];
ds_id = datamb{i}.cell_ids(idx==2);
nonds_id = datamb{i}.cell_ids(idx==1);

%%
load('DS141212.mat')
[raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},ds_id,0,1);
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
ll_p = 4;
MAG_all_norm_mb = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

for i = 1:n
    [raster_p_sum_mb{i}, p_idx{i}] = get_pdirection_raster(raster_mb{i}, MB{ll_p}.angle{delta_p});
    MAG_all_norm_mb{i} = normalize_MAG(MB{i});
    rep = datamb{i}.stimulus.repetitions;
end

ll = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};
ct = {'posterior', 'inferior', 'anterior', 'superior'};

%% plot cell summary

for d = 4:4
for cc = 6:length(ds_id)
    plot_mb_raster(MB, raster_mb, trial_dur, idx_dir{d}(cc), ds_id(idx_dir{d}(cc)), ll, 2, 3, 1)
end
end

%% classify DSGC into subtypes (directions)
d = 4;
t = 2;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(MB{d}.U{t}, MB{d}.V{t})
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(MB{d}.U{t}, MB{d}.V{t}, x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = I;
end

%% DS tuning curves (moving bar)
% all ds cells

% t = 2;
dirn = 4;
D = 4;
T = 2;

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
        if i == 3 && d == 1
            rho_mb{1}{3}(1, :) = [rho_mb{1}{3}(1, 2:8) 0];
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

xtick = {'posterior'; 'inferior'; 'anterior'; 'superior'};
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

%% drifting grating

[raster_dg, DG, trial_dur, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, StimComb] = get_spikescellstim(datadg{i},ds_id,0);
    DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster_dg{i} = get_ds_raster(datadg{i}, ds_id);
    for j = 1:length(raster_dg{i})
        if(dg_idx(j, i))
            raster_dg{i}{j} = [];
        end
    end
end

delta_p = 2; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
    MAG_all_norm_dg{i} = normalize_MAG(DG{i});
    rep = datadg{i}.stimulus.repetitions;
end

ll = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};

%% DS tuning curves (drifting grating)
% all ds cells

% t = 2;
dirn = 4;
D = 4;
T = 2;

for d = 1:n
    p_direction = DG{D}.angle{T}';
    xx = 0:pi/9:17*pi/9;
    xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 18);
    xx(xx>pi) = xx(xx>pi)-2*pi;
    xx(xx<-pi) = xx(xx<-pi)+2*pi;


    subplot(2, 3, d)
    for i = 1:dirn
        for cc = 1:length(ds_id)
            if ~dg_idx(cc, d)
            [xsort, seq] = sort(xx(cc, :));
            y_temp = DG{d}.rho{T}(cc, :);
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
            if ~dg_idx(idx_dir{i}(cc), d) && sum(DG{d}.rho{T}(idx_dir{i}(cc), :))>0
            [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
            y_temp = DG{d}.rho{T}(idx_dir{i}(cc), :);
            plot(xsort, y_temp(seq), color(i))
            ylim([0 1])
%             pause
            hold on
            rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
            dsi_dg{d}{i} = [dsi_dg{d}{i}; DG{d}.dsindex{T}(idx_dir{i}(cc))];
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
for i = 1:4
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

xtick = {'posterior'; 'inferior'; 'anterior'; 'superior'};
model_series = [dsi_dg_mean(1,1) dsi_dg_mean(2,1) dsi_dg_mean(3,1) dsi_dg_mean(4,1) dsi_dg_mean(5,1); dsi_dg_mean(1,2) dsi_dg_mean(2,2) dsi_dg_mean(3,2) dsi_dg_mean(4,2) dsi_dg_mean(5,2);dsi_dg_mean(1,3) dsi_dg_mean(2,3) dsi_dg_mean(3,3) dsi_dg_mean(4,3) dsi_dg_mean(5,3); dsi_dg_mean(1,4) dsi_dg_mean(2,4) dsi_dg_mean(3,4) dsi_dg_mean(4,4) dsi_dg_mean(5,4)];   
model_error = [dsi_dg_ste(1,1) dsi_dg_ste(2,1) dsi_dg_ste(3,1) dsi_dg_ste(4,1) dsi_dg_ste(5,1); dsi_dg_ste(1,2) dsi_dg_ste(2,2) dsi_dg_ste(3,2) dsi_dg_ste(4,2) dsi_dg_ste(5,2);dsi_dg_ste(1,3) dsi_dg_ste(2,3) dsi_dg_ste(3,3) dsi_dg_ste(4,3) dsi_dg_ste(5,3); dsi_dg_ste(1,4) dsi_dg_ste(2,4) dsi_dg_ste(3,4) dsi_dg_ste(4,4) dsi_dg_ste(5,4)];
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

%% plot dg cell summary
DG_8dir = DG;
for d = 1:5
    DG_8dir{d}.rho = [];
    DG_8dir{d}.theta = [];
    for t = 1:2
        i = [1 3 5 8 10 12 14 17];
        DG_8dir{d}.rho{t} = DG{d}.rho{t}(:, i);
        DG_8dir{d}.theta{t} = DG{d}.theta{t}(:, i);
    end
    for cc = 1:length(ds_id)
        if ~isempty(raster_dg{d}{cc})
            raster_dg_8dir{d, 1}{cc, 1} = raster_dg{d}{cc}(:, :, i, :);
        else
            raster_dg_8dir{d, 1}{cc, 1} = [];
        end
    end
end

for cc = 1:length(ds_id)
    plot_ds_raster(DG_8dir, raster_dg_8dir, cc, ds_id(cc), ll, 2, 3, 1)
end

%% ffp

n_ffp = 5;

[raster_ff, raster_ff_all] = deal(cell(n_ffp, 1));
for d = 1:n_ffp
    [raster_ff{d}, raster_ff_all{d}] = get_ffp_raster(dataffp{d}, ds_id, 3);
    for j = 1:length(raster_ff{d})
        if(ffp_idx(j, d))
            raster_ff{d}{j} = [];
            raster_ff_all{d}{j} = [];
        end
    end
end

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

%% Classify ON and ON-OFF
% bin_size = 0.5 second
MB_maxrate = cell(n, 1);
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},ds_id,0,0.5);
    MB_maxrate{i} = sort_direction(mbcellanalysis(MaxRate, StimComb));
end

MAG_all_norm_mb_maxrate = cell(n, 1);

for i = 1:n
    MAG_all_norm_mb_maxrate{i} = normalize_MAG(MB{i});
    for cc = 1:length(ds_id)
        if mb_idx(cc, i)
            MAG_all_norm_mb_maxrate{i}(:, cc) = nan;
        end
    end
    MAG_all_norm_mb_maxrate_avg(:, i) = mean(exciseColumn(MAG_all_norm_mb_maxrate{i}), 2);
    MAG_all_norm_mb_maxrate_ste(:, i) = std(exciseColumn(MAG_all_norm_mb_maxrate{i}), 0, 2)/sqrt(size(exciseColumn(MAG_all_norm_mb_maxrate{i}), 2));
end

figure
x = [1 2 4 20] * 240;
for i = 1:5
    errorbar(x, MAG_all_norm_mb_maxrate_avg(:, i), MAG_all_norm_mb_maxrate_ste(:, i), color(i))
    hold on
end
set(gca, 'XScale', 'log')
legend(ll)
xlim([0.8 30] * 240)
xlabel('micron/s')
ylabel('normalized response')
title('speed tuning (spike number)')
%% pca
mag_pca = MAG_all_norm_mb_maxrate{5};
% mag_pca = MAG_all_norm{1}(:, idx4);
mag_pca = mag_pca';
[id_sub, idx_sub] = deal(cell(3, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 800 800])
subplot(2, 2, 1)
[~,scores,~,~] = princomp(mag_pca);
pc1 = 1; pc2 = 2;
plot(scores(:, pc1), scores(:, pc2), 'o')
% subplot(2, 2, 2)
% [~,scores,~,~] = princomp(mag_pca);
% pc1 = 1; pc2 = 3;
% plot(scores(:, pc1), scores(:, pc2), 'o')
% subplot(2, 2, 3)
% [~,scores,~,~] = princomp(mag_pca);
% pc1 = 2; pc2 = 3;
% plot(scores(:, pc1), scores(:, pc2), 'o')
% 
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

% frequency analysis
    
duration = 8;
bin_rate = 10000;
hist_spikes = cell(n, 1);
signal_length = duration*bin_rate;                
NFFT = 2^nextpow2(signal_length);
f = bin_rate/2*linspace(0,1,NFFT/2+1);
fft_spikes = cell(n, 1);
[DC, F1, F2] = deal(cell(n, 1));

for i = 1:n
    tp = datadg{i}.stimulus.params.TEMPORAL_PERIOD;
    [DC{i}, F1{i}, F2{i}] = deal(zeros(length(ds_id), length(tp)));
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
        F1{i}(rgc,:) = fund_power ./ max(DC_power);
        F2{i}(rgc,:) = sec_power ./ max(DC_power);
        DC{i}(rgc,:) = DC_power ./ max(DC_power);
        
        clear fund_power sec_power DC_power

        end
        
    end
    ratio{i} = F2{i}./F1{i};
    for ct = 1:4
        ratio_dir{ct}{i} = ratio{i}(idx_dir{ct}, :);
        ratio_dir{ct}{i} = exciseRows(ratio_dir{ct}{i});
        ratio_dir_mean(i, ct, :) = mean(ratio_dir{ct}{i});
        ratio_dir_ste(i, ct, :) = std(ratio_dir{ct}{i})/sqrt(size(ratio_dir{ct}{i}, 1));
    end
    ratio{i} = exciseRows(ratio{i});
    ratio_mean(i, :) = mean(ratio{i});
    ratio_ste(i, :) = std(ratio{i})/sqrt(size(ratio{i}, 1));
end

% plot average f2/f1

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = {'high speed'; 'low speed'};
model_series = [ratio_mean(1,1) ratio_mean(2,1) ratio_mean(3,1) ratio_mean(4,1) ratio_mean(5,1); ratio_mean(1,2) ratio_mean(2,2) ratio_mean(3,2) ratio_mean(4,2) ratio_mean(5,2)];   
model_error = [ratio_ste(1,1) ratio_ste(2,1) ratio_ste(3,1) ratio_ste(4,1) ratio_ste(5,1); ratio_ste(1,2) ratio_ste(2,2) ratio_ste(3,2) ratio_ste(4,2) ratio_ste(5,2)];
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

%% plot cell type specific f2/f1 

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = {'posterior'; 'inferior'; 'anterior'; 'superior'};
model_series = [ratio_dir_mean(1,1,1) ratio_dir_mean(2,1,1) ratio_dir_mean(3,1,1) ratio_dir_mean(4,1,1) ratio_dir_mean(5,1,1); ratio_dir_mean(1,2,1) ratio_dir_mean(2,2,1) ratio_dir_mean(3,2,1) ratio_dir_mean(4,2,1) ratio_dir_mean(5,2,1); ratio_dir_mean(1,3,1) ratio_dir_mean(2,3,1) ratio_dir_mean(3,3,1) ratio_dir_mean(4,3,1) ratio_dir_mean(5,3,1); ratio_dir_mean(1,4,1) ratio_dir_mean(2,4,1) ratio_dir_mean(3,4,1) ratio_dir_mean(4,4,1) ratio_dir_mean(5,4,1)];   
model_error = [ratio_dir_ste(1,1,1) ratio_dir_ste(2,1,1) ratio_dir_ste(3,1,1) ratio_dir_ste(4,1,1) ratio_dir_ste(5,1,1); ratio_dir_ste(1,2,1) ratio_dir_ste(2,2,1) ratio_dir_ste(3,2,1) ratio_dir_ste(4,2,1) ratio_dir_ste(5,2,1); ratio_dir_ste(1,3,1) ratio_dir_ste(2,3,1) ratio_dir_ste(3,3,1) ratio_dir_ste(4,3,1) ratio_dir_ste(5,3,1); ratio_dir_ste(1,4,1) ratio_dir_ste(2,4,1) ratio_dir_ste(3,4,1) ratio_dir_ste(4,4,1) ratio_dir_ste(5,4,1)];
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
xtick = {'posterior'; 'inferior'; 'anterior'; 'superior'};
model_series = [ratio_dir_mean(1,1,2) ratio_dir_mean(2,1,2) ratio_dir_mean(3,1,2) ratio_dir_mean(4,1,2) ratio_dir_mean(5,1,2); ratio_dir_mean(1,2,2) ratio_dir_mean(2,2,2) ratio_dir_mean(3,2,2) ratio_dir_mean(4,2,2) ratio_dir_mean(5,2,2); ratio_dir_mean(1,3,2) ratio_dir_mean(2,3,2) ratio_dir_mean(3,3,2) ratio_dir_mean(4,3,2) ratio_dir_mean(5,3,2); ratio_dir_mean(1,4,2) ratio_dir_mean(2,4,2) ratio_dir_mean(3,4,2) ratio_dir_mean(4,4,2) ratio_dir_mean(5,4,2)];   
model_error = [ratio_dir_ste(1,1,2) ratio_dir_ste(2,1,2) ratio_dir_ste(3,1,2) ratio_dir_ste(4,1,2) ratio_dir_ste(5,1,2); ratio_dir_ste(1,2,2) ratio_dir_ste(2,2,2) ratio_dir_ste(3,2,2) ratio_dir_ste(4,2,2) ratio_dir_ste(5,2,2); ratio_dir_ste(1,3,2) ratio_dir_ste(2,3,2) ratio_dir_ste(3,3,2) ratio_dir_ste(4,3,2) ratio_dir_ste(5,3,2); ratio_dir_ste(1,4,2) ratio_dir_ste(2,4,2) ratio_dir_ste(3,4,2) ratio_dir_ste(4,4,2) ratio_dir_ste(5,4,2)];
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
% for i = 1:n
for i = 4:4
    for t = 1:4
        dur = max(trial_dur{i}(t, :));
        XX = bin_size/2:bin_size:dur;
        hist_mb_p{i}{t} = zeros(length(ds_id), length(XX));
        for cc = 1:length(ds_id)
            if ~mb_idx(cc, i) && ~isempty(raster_p_sum_mb{i}{cc})
                hist_mb_p{i}{t}(cc, :) = hist(raster_p_sum_mb{i}{cc}{t}, XX);
            end
        end
    end
end

x = [5 3 3 6]; y = [5 3 1 9];
t = 2;
for ct = 1:1
    figure(ct)
    for cc = 4:15 %length(idx_dir{ct})
%         subplot(x(ct), y(ct), cc)
        subplot(3, 4, cc-3)
        dur = max(trial_dur{1}(t, :));
        XX = bin_size/2:bin_size:dur;
        for i = 4:4
%         for i = 1:n
            plot(XX, 1.25*hist_mb_p{i}{t}(idx_dir{ct}(cc), :), color(i))
            hold on
        end
%         title(num2str(ds_id(idx_dir{ct}(cc))))
    end
    legend('NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0')
end
            
    
% F2/F1 ratio
delta = datamb{1}.stimulus.params.DELTA(t);
bar_width = datamb{1}.stimulus.params.BAR_WIDTH;
refresh_rate = 60.3578;
bar_time = bar_width/delta/refresh_rate;
range = 1/3;
for ct = 1:4
    for cc = 1:length(idx_dir{ct})
        figure(1)
        dur = max(trial_dur{1}(t, :));
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
for i = 1:5
    for ct = 1:4
        for t = T:T
            for cc = 1:size(max_onoff{i}{ct}{t}, 1)
                ratio_onoff{ct}{t}{i}(cc) = max_onoff{i}{ct}{t}(cc, 2)/max_onoff{i}{ct}{t}(cc, 1);
            end
            index = logical(isnan(ratio_onoff{ct}{t}{i}) + isinf(ratio_onoff{ct}{t}{i}));
            ratio_onoff{ct}{t}{i}(index) = [];
            ratio_avg{ct}{t}(i) = mean(ratio_onoff{ct}{t}{i});
            ratio_ste{ct}{t}(i) = std(ratio_onoff{ct}{t}{i})/sqrt(length(ratio_onoff{ct}{t}{i}));
        end
    end
end

%%
t = 2;
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = {'posterior'; 'inferior'; 'anterior'; 'superior'};
model_series = [ratio_avg{1}{t}(1) ratio_avg{1}{t}(2) ratio_avg{1}{t}(3) ratio_avg{1}{t}(4) ratio_avg{1}{t}(5); ratio_avg{2}{t}(1) ratio_avg{2}{t}(2) ratio_avg{2}{t}(3) ratio_avg{2}{t}(4) ratio_avg{2}{t}(5); ratio_avg{3}{t}(1) ratio_avg{3}{t}(2) ratio_avg{3}{t}(3) ratio_avg{3}{t}(4) ratio_avg{3}{t}(5); ratio_avg{4}{t}(1) ratio_avg{4}{t}(2) ratio_avg{4}{t}(3) ratio_avg{4}{t}(4) ratio_avg{4}{t}(5)];   
model_error = [ratio_ste{1}{t}(1) ratio_ste{1}{t}(2) ratio_ste{1}{t}(3) ratio_ste{1}{t}(4) ratio_ste{1}{t}(5); ratio_ste{2}{t}(1) ratio_ste{2}{t}(2) ratio_ste{2}{t}(3) ratio_ste{2}{t}(4) ratio_ste{2}{t}(5); ratio_ste{3}{t}(1) ratio_ste{3}{t}(2) ratio_ste{3}{t}(3) ratio_ste{3}{t}(4) ratio_ste{3}{t}(5); ratio_ste{4}{t}(1) ratio_ste{4}{t}(2) ratio_ste{4}{t}(3) ratio_ste{4}{t}(4) ratio_ste{4}{t}(5)];   
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

%%
d = 4;
t = 2;
id = [3452 4608 4866 6121 7653 4490 4178 4862];
h = figure(1);
set(h, 'Position', [1 1 500 500])
compass(MB{d}.U{t}, MB{d}.V{t})
hold on
for i = 1:length(id)
    compass(MB{d}.U{t}(ds_id == id(i)), MB{d}.V{t}(ds_id == id(i)), 'r')
end

