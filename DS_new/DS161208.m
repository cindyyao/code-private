%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load drifting grating data
datadg = load_data('/Volumes/lab/analysis/2016-12-08-0/data009/data009', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-12-08-0/stimuli/s09.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

% identify DS cells
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [2 3]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

clear datadg
dataDG{1} = load_data('/Volumes/lab/analysis/2016-12-08-0/data002-map/data002-map', opt);
dataDG{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-12-08-0/stimuli/s02.txt';
dataDG{1} = load_stim(dataDG{1}, 'user_defined_trigger_interval', 10);
dataDG{2} = load_data('/Volumes/lab/analysis/2016-12-08-0/data003-map/data003-map', opt);
dataDG{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-12-08-0/stimuli/s03.txt';
dataDG{2} = load_stim(dataDG{2}, 'user_defined_trigger_interval', 10);
dataDG{3} = load_data('/Volumes/lab/analysis/2016-12-08-0/data004-map/data004-map', opt);
dataDG{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-12-08-0/stimuli/s04.txt';
dataDG{3} = load_stim(dataDG{3}, 'user_defined_trigger_interval', 10);
dataDG{4} = load_data('/Volumes/lab/analysis/2016-12-08-0/data005-map/data005-map', opt);
dataDG{4}.names.stimulus_path = '/Volumes/lab/analysis/2016-12-08-0/stimuli/s05.txt';
dataDG{4} = load_stim(dataDG{4}, 'user_defined_trigger_interval', 10);
dataDG{5} = load_data('/Volumes/lab/analysis/2016-12-08-0/data009/data009', opt);
dataDG{5}.names.stimulus_path = '/Volumes/lab/analysis/2016-12-08-0/stimuli/s09.txt';
dataDG{5} = load_stim(dataDG{5}, 'user_defined_trigger_interval', 10);
dataDG{5} = load_ei(dataDG{5}, 'all');

datadg = load_data('/Volumes/lab/analysis/2016-12-08-0/data007-map/data007-map', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-12-08-0/stimuli/s07.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

datafs = load_data('/Volumes/lab/analysis/2016-12-08-0/data010-map/data010-map', opt);
datafs.names.stimulus_path = '/Volumes/lab/analysis/2016-12-08-0/stimuli/s10.mat';
datafs = load_stim_mfs(datafs);

%
dataffp = load_data('/Volumes/lab/analysis/2016-12-08-0/data008-map/data008-map', opt);
dataffp.triggers([1:50]*3) = [];

% 
datarun1 = load_data('/Volumes/lab/analysis/2016-12-08-0/data000/data000', opt);
datarun = load_data('/Volumes/lab/analysis/2016-12-08-0/data000-001-map/data000-001-map', opt);
dataflash(1:2) = split_datarun(datarun, 1081); 
dataflash{2}.DfParams.NDF =   [5,5,5,5,4,4,4,4,4,3,3,3,2,2] ; % on filter turret 
dataflash{2}.DfParams.Ftime = [1,4,2,8,2,3,4,6,8,2,4,8,2,8] ; % ms
dataflash{2}.DfParams.interFlashInt = [3] ; % sec

load('DS161208.mat')
%% dg
n = 1;
i = 1;
[raster_dg, dg, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
dg{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);


delta_p = 4; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, dg{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(dg{i});
rep = datadg.stimulus.repetitions;

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

t = 5;
figure
compass(dg{1}.U{t}(idx_sub{1}), dg{1}.V{t}(idx_sub{1}), 'r')
hold on
compass(dg{1}.U{t}(idx_sub{2}), dg{1}.V{t}(idx_sub{2}), 'b')


%%
d = 1;
t = 5;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(dg{d}.U{t}(idx_sub{2}), dg{d}.V{t}(idx_sub{2}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(dg{d}.U{t}(idx_sub{2}), dg{d}.V{t}(idx_sub{2}), x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = idx_sub{2}(I);
    id_dir{i} = ds_id(idx_dir{i});
end

d = 1;
t = 3;
h = figure;
dirn = 3;
set(h, 'Position', [1 1 1080 500])
compass(dg{d}.U{t}(idx_sub{1}), dg{d}.V{t}(idx_sub{1}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(dg{d}.U{t}(idx_sub{1}), dg{d}.V{t}(idx_sub{1}), x, y);
    [~, I] = find(IN == 1);
    idx_dir_on{i} = idx_sub{1}(I);
    id_dir_on{i} = ds_id(idx_dir_on{i});
end

%% DG
n = 5;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, ~, StimComb] = get_spikescellstim(dataDG{i},ds_id,0,1);
    DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,dataDG{i}));
    raster_dg{i} = get_ds_raster(dataDG{i}, ds_id);
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
    rep = dataDG{i}.stimulus.repetitions;
end

ll = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};

%% plot cell summary
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [4 5]);

for dir = 1:3
    for cc = 1:length(id_dir_on{dir})
        plot_ds_raster(DG, raster_dg, idx_dir_on{dir}(cc), id_dir_on{dir}(cc), ll, 2, 3, 1)
    end
end

%% DS tuning curves (drifting grating)
% all ds cells
color = 'brgkc';
LL = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};
ct = {'superior', 'anterior', 'inferior', 'posterior'};
% t = 2;
dirn = 4;
D = 5;
T = 3;

p_direction = DG{D}.angle{T}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

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
%     title(ll{d})
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
%     title(ll{d});
end
legend(ct)

% plot average (light level)
color = 'rbgkc';
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
legend(LL)
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
marker = 'xo*d';
figure
for i = 1:4
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
errorbar(0:4, model_series(i,:), model_error(i,:), 'color', color(i));
hold on
end
xlim([-0.5 4.5])
ylim([0 1.1])
legend(ct)
xlabel('log(light intensity)')
ylabel('DSI')


for dir = 1:4
    figure
    for ll = 1:5
        subplot(2, 3, ll)
        hist(dsi_dg{ll}{dir}, 0:0.05:1)
        xlim([0 1])
        title(LL{ll})
        xlabel('DSI')
        ylabel('cell#')
    end
end

%% compare with WT
load('DS150603.mat', 'dsi_dg_wt')
figure
X = 0.05:0.1:0.95;
dsi_wt = hist(dsi_dg_wt{1}{1}, X);
dsi_ko = hist(dsi_dg{1}{1}, X);
bar(X, [dsi_wt' dsi_ko'], 1, 'stacked')
xlim([0 1])
xlabel('DSI')
ylabel('cell#')
legend('WT', 'KO')
title('NDF 4')
%% plot rfs
cell_type = {'anterior', 'inferior', 'posterior', 'superior'};
field_width = 13; field_height = 13;
field_width_sta = 40; field_height_sta = 40;
subregion = 0;
stop = 0.5; %second
%     fs_raster{i} = get_fs_raster(datafs{i}, ds_id, 'stop', 0.5);
fs_raster = get_fs_raster(datafs, ds_id);
for cc = 1:length(ds_id)
    if fs_idx(cc)
        fs_raster{cc} = [];
    end
end
fs_spike = fs_get_spike(fs_raster);
[rf_all, rf_std] = get_fs_rf(fs_spike, field_width, field_height,subregion);

%%
for cc = 1:length(ds_id)
    figure(1)
    set(gcf, 'Position', [1 1 900 300])
    id = ds_id(cc);
    if ~isempty(rf_all{cc})
%             rf = padarray(rf_all{i}{cc},[7,7]);
        rf = rf_all{cc};

        subplot(1,3,1)
        imagesc(sum(rf,3))
        colormap gray
        axis image
        axis off

        subplot(1,3,2)
        imagesc(rf(:,:,1))
        colormap gray
        axis image
        axis off

        subplot(1,3,3)
        imagesc(rf(:,:,2))
        colormap gray
        axis image
        axis off

    end
    print_close(1,[15 5],num2str(id))
end

%
% fs_spike_temp = cellfun(@(fs_spike) sum(fs_spike,2), fs_spike, 'UniformOutput', false);
% fs_spike_temp = cellfun(@(fs_spike_temp) sum(fs_spike_temp,3), fs_spike_temp, 'UniformOutput', false);
% [~, I] = cellfun(@(fs_spike_temp) sort(fs_spike_temp, 'descend'), fs_spike_temp, 'UniformOutput', false);trigger = datafs.triggers(2:2:end);
trigger = datafs.triggers(1:2:end);
list = datafs.stimulus.trial_list;
repeat = datafs.stimulus.repetitions;
for i = 1:max(list)
    index(i, :) = find(list == i);
end
raster_onoff = cell(length(ds_id), 1);
for i = 1:length(ds_id)
    if ~fs_idx(i)
        idx = get_cell_indices(datafs, ds_id(i));
        raster = get_raster(datafs.spikes{idx}, trigger, 'plot', false);
        raster_onoff{i} = raster(index);
    end
end

%
LL = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};
bin_size = 0.02;
xx = [bin_size/2:bin_size:2-bin_size/2];
XX = [bin_size/2:bin_size:1-bin_size/2];
for cc = 1:length(raster_onoff)
    if ~isempty(raster_onoff{cc})
        for p = 1:169%size(raster_onoff{ll}{1}, 1)
            raster_all_onoff{cc}{p} = sort(cell2mat(raster_onoff{cc}(p,:)'));
            histg_onoff{cc}{p} = hist(raster_all_onoff{cc}{p}, xx);
            for onoff = 1:2
                raster_all{cc}{onoff}{p} = sort(cell2mat(fs_raster{cc}(p,:,onoff)'));
                histg_all{cc}{onoff}{p} = hist(raster_all{cc}{onoff}{p}, XX);
            end
        end
    end
end

h = subplot(field_height, field_width, 1); p1 = get(h, 'pos');
h = subplot(field_height, field_width, 2); p2 = get(h, 'pos');
width = (p2(1) - p1(1))*0.8;

h = subplot(field_height, field_width, 1); p1 = get(h, 'pos');
h = subplot(field_height, field_width, field_width+1); p2 = get(h, 'pos');
height = p1(2) - p2(2);
for cc = 1:length(ds_id)
    for ll = 1:5
        H = figure(1);
        set(H, 'Position', [1 1 1080 1080])
        if ~fs_idx(cc,ll)
            for p = 1:289
                h = subplot(field_width, field_height, p);
                l = get(h, 'pos');
                l(3) = width; l(4) = height;
                set(h, 'pos', l);
                plot(xx, histg_onoff{ll}{cc}{p})
                ylim([0 max(max(cell2mat(histg_onoff{ll}{cc}')))])
                axis off
            end
        end
%     pause 
%     close all
        name = [num2str(ds_id(cc)) '_' LL{ll}];
        print_close(1, [14 14], name);
    end
end

%% gaussian filter
PixelArea = (30*4)^2/10^6;
threshold = 0.3;

clear rf_wt rf_wt_area rf_wt_area_mean rf_wt_area_ste rf_area rf_area_clean
tau = 0.05;
tt = -3*tau:bin_size:3*tau;
filter = exp(-tt.^2/(2*tau^2));

filter = filter/norm(filter);
npixel = 5;
for dir = 1:4

    for cc = 1:length(id_dir{dir})
        if ~fs_idx(idx_dir{dir}(cc))
            for onoff = 1:2
                % use the npixel brightest pixels to calculate a psth template
                [~, idx] = sort(cellfun(@max,histg_all{idx_dir{dir}(cc)}{onoff}), 'descend');
                idx = idx(1:npixel);
                temp = histg_all{idx_dir{dir}(cc)}{onoff}(idx);
                for i = 1:npixel
                    temp{i} = temp{i}/norm(temp{i});
                end
                temp_mean = mean(cell2mat(temp'));
                temp_mean_norm = temp_mean/norm(temp_mean);
                
%                 figure
%                 subplot(1, 2, 1); plot(temp_mean_norm)
                
                temp_mean_norm = conv(temp_mean_norm, filter, 'same');

%                 subplot(1, 2, 2); plot(temp_mean_norm)
%                 pause
                
                for p = 1:169
                    rf_wt{dir}{cc}{onoff}(p) = histg_all{idx_dir{dir}(cc)}{onoff}{p}*temp_mean_norm';
                end
                rf_wt{dir}{cc}{onoff} = reshape(rf_wt{dir}{cc}{onoff}, field_width, field_height);
                rf_wt_area{dir}{onoff}{cc} = sum(sum(rf_wt{dir}{cc}{onoff} > max(rf_wt{dir}{cc}{onoff}(:))*threshold))*PixelArea;
            end
        end
    end
    for onoff = 1:2
        rf_wt_area{dir}{onoff} = cell2mat(rf_wt_area{dir}{onoff});
        rf_wt_area_mean(dir, onoff) = mean(rf_wt_area{dir}{onoff});
        rf_wt_area_ste(dir, onoff) = std(rf_wt_area{dir}{onoff})/sqrt(length(rf_wt_area{dir}{onoff}));
    end
end


%% fit RF with Gaussian                
PixelArea = (30*4)^2/10^6;
% fit and compute rf area
for dir = 1:4
    clear rf_area_temp
    for onoff = 1:2
        rf_area_temp = [];
        for cc = 1:length(id_dir{dir})
            if ~fs_idx(idx_dir{dir}(cc))
                data = rf_all{idx_dir{dir}(cc)}(:, :, onoff);
%                     data = rf_wt{dir}{cc}{onoff};
                if sum(sum(data > mean(data(:))+5*std(data(:))))>0
%                         figure(100)
%                         imagesc(data)
%                         colormap gray
%                         pause

                    params = fit_2d_gaussian(data);
%                     Gaussian_params{ll}{dir}{cc}{onoff} = params;
                    rf_area_temp = [rf_area_temp params.xd * params.yd * pi * PixelArea];
                end
            end
        end
        rf_area{dir}{onoff} = rf_area_temp;
    end
end


% exclude outliers
stdn = 100;
for dir = 1:4
    for onoff = 1:2
        notdone = 1;
        rf_area_temp = rf_area{dir}{onoff};
        while notdone
            a = length(rf_area_temp);
            rf_area_temp(rf_area_temp > std(rf_area_temp)*stdn + mean(rf_area_temp)) = [];
            b = length(rf_area_temp);
            if a == b
                notdone = 0;
                rf_area_clean{dir}{onoff} = rf_area_temp;
            end
        end
        rf_area_clean_mean{onoff}(dir) = mean(rf_area_clean{dir}{onoff});
        rf_area_clean_ste{onoff}(dir) = std(rf_area_clean{dir}{onoff})/sqrt(length(rf_area_clean{dir}{onoff}));
    end
end


% plot 
ll = 1;
color = 'brgkc';
figure
for onoff = 1:2
    subplot(1,2,onoff)
    for dir = 1:4
        n = size(rf_area_clean{dir}{onoff}, 1);
        h = plot((dir-1)*5+ll*ones(n,1), rf_area_clean{dir}{onoff}, [color(ll) 'o']);
%             n = length(rf_area{ll}{dir}{onoff});
%             h{ll} = plot((dir-1)*5+ll*ones(n,1), rf_area{ll}{dir}{onoff}, [color(ll) 'o']);
        hold on
    end
%     set(gca, 'yscale', 'log')
%     legend([h{1}(1), h{2}(1), h{3}(1), h{4}(1), h{5}(1)], 'NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0')
    if onoff == 1
        title('ON')
    else
        title('OFF')
    end
    ylabel('RF area (mm^2)')
    set(gca, 'xtick', [])
%         ylim([0 0.3])

end


%% full field pulses

n_ffp = 1;

[raster_ff, raster_ff_all] = deal(cell(n_ffp, 1));
for d = 1:n_ffp
    [raster_ff{d}, raster_ff_all{d}] = get_ffp_raster(dataffp, ds_id, 3);
    for j = 1:length(raster_ff{d})
        if(ffp_idx(j))
            raster_ff{d}{j} = [];
            raster_ff_all{d}{j} = [];
        end
    end
end

for i = 1:length(ds_id) 
    if ~isempty(raster_ff{1}{i})
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 400 800])
        for d = 1:n_ffp
        subplot(1, n_ffp, d)
        plot_ffp(raster_ff{d}, raster_ff_all{d}, i, 3)
        title(num2str(ds_id(i)))
        end
        
        print_close(1, [6, 12], num2str(ds_id(i)))
    end
end
%% dim flashes
ds_id_flash = ds_id(~flash_idx);
ds_idx_flash = get_cell_indices(dataflash{2}, ds_id_flash);

for ct = 1:4
    [id_dir_flash{ct}, idx_dir_flash{ct}] = intersect(ds_id_flash, id_dir{ct});
end
for ct = 1:3
    [id_dir_on_flash{ct}, idx_dir_on_flash{ct}] = intersect(ds_id_flash, id_dir_on{ct});
end

%% parameters
interFlashIntVar = 0.005; % (sec) expected precision of trigger intervals

%% find sets of flash stimuli
ts=1 ; 
trigger_set_i{1} = [] ;
for a=1:length(dataflash{2}.triggers) ; % for each trigger
    trigger_set_i{ts} = [trigger_set_i{ts},a] ; % put it in a set   
    if a<length(dataflash{2}.triggers) ; % if its not the last trigger
        if sum(abs(dataflash{2}.triggers(a+1)-dataflash{2}.triggers(a)-dataflash{2}.DfParams.interFlashInt)<interFlashIntVar)>0 ; % next trigger is has the right interval
            ts = ts ; % keep it in the set
        else
            ts = ts+1 ; % put it in a new set
            trigger_set_i{ts} = [] ;
        end
    end
end

i = 1;
while(i<=length(trigger_set_i))
    if length(trigger_set_i{i}) < 60
        trigger_set_i(i) = [];
    else
        i = i+1;
    end
end

trigger_set_i(2) = [];

%%
bin_size = 0.02; 
start = 0;
XX = start+bin_size/2:bin_size:dataflash{2}.DfParams.interFlashInt-bin_size/2;
% load /Volumes/lab/Experiments/Calibration/NdfCalibration
Irel = (dataflash{2}.DfParams.Ftime/1000).*NdfCalibration(2,dataflash{2}.DfParams.NDF+1) ;
[Irel, i] = sort(Irel);

ds_flash_raster = cell(length(ds_id_flash), 1);
ds_flash_hist_trial = cell(length(ds_id_flash), 1);
ds_flash_hist_mean = cell(length(ds_id_flash), 1);
for cc = 1:length(ds_id_flash)
    raster_1cell = cell(length(trigger_set_i), 1);
    raster_1cell_all = cell(length(trigger_set_i), 1);
    hist_1cell = cell(length(trigger_set_i), 1);
    hist_1cell_all = cell(length(trigger_set_i), 1);
    for ts = 1:length(trigger_set_i)
        raster_1cell{ts} = get_raster(dataflash{2}.spikes{ds_idx_flash(cc)}, ...
            dataflash{2}.triggers(trigger_set_i{ts}), 'start', start, 'stop', ...
            dataflash{2}.DfParams.interFlashInt, 'plot', 0);
        for t = 1:length(raster_1cell{ts})
            hist_1cell{ts}{t} = hist(raster_1cell{ts}{t}, XX);
        end
        raster_1cell_all{ts} = sort(cell2mat(raster_1cell{ts}));
        hist_1cell_all{ts} = hist(raster_1cell_all{ts}, XX)/length(trigger_set_i{ts});
    end
    ds_flash_raster{cc} = raster_1cell(i, :);
    ds_flash_hist_trial{cc} = hist_1cell(i, :);
    ds_flash_hist_mean{cc} = hist_1cell_all(i, :);
end

%% get raster and psth for dark trials
% use 400-580 second (roughly) as dark trials

trigger = 1:3:300;
ds_dark_raster = cell(length(ds_id_flash), 1);
ds_dark_raster_all = cell(length(ds_id_flash), 1);
ds_dark_hist_trial = cell(length(ds_id_flash), 1);
ds_dark_hist_mean = cell(length(ds_id_flash), 1);
for cc = 1:length(ds_id_flash)
    ds_dark_raster{cc} = get_raster(dataflash{1}.spikes{ds_idx_flash(cc)}, ...
        trigger, 'stop', dataflash{2}.DfParams.interFlashInt, ...
        'plot', 0);
    for t = 1:length(ds_dark_raster{cc})
        ds_dark_hist_trial{cc}{t} = hist(ds_dark_raster{cc}{t}, XX);
    end
    ds_dark_raster_all{cc} = sort(cell2mat(ds_dark_raster{cc}));
    ds_dark_hist_mean{cc} = hist(ds_dark_raster_all{cc}, XX)/length(trigger);
end

%% 2 alternative forced choice test
% use 400-580 second (roughly) as dark trials
window = 1;
bin_n = window/bin_size;

Pc = zeros(length(ds_id_flash), length(trigger_set_i));
for cc = 1:length(ds_id_flash)
    for ts = 1:length(trigger_set_i)
%         trial_n = length(ds_flash_raster{1}{ts});
        trial_n = 60;
        corr_flash = zeros(trial_n, 1);
        corr_dark = zeros(trial_n, 1);
        temp = ds_dark_hist_trial{cc}(1:trial_n);
        ds_dark_hist_sum = sum(cell2mat(temp'));
        ds_flash_hist_sum = sum(cell2mat(ds_flash_hist_trial{cc}{ts}'));
        for t = 1:trial_n
            template_flash = ds_flash_hist_sum - ds_flash_hist_trial{cc}{ts}{t};
            template_flash = template_flash/norm(template_flash);
            template_flash(isnan(template_flash)) = 0;
            template_dark = ds_dark_hist_sum - ds_dark_hist_trial{cc}{t};
            template_dark = template_dark/norm(template_dark);
            template_dark(isnan(template_dark)) = 0;
            DV = template_flash - template_dark;
            corr_flash(t) = ds_flash_hist_trial{cc}{ts}{t}(1:bin_n) * DV(1:bin_n)';
            corr_dark(t) = ds_dark_hist_trial{cc}{t}(1:bin_n) * DV(1:bin_n)';
        end
        Pc(cc, ts) = (sum(corr_flash > 0) + sum(corr_dark <= 0))/(trial_n*2);
    end
end


color = 'rbgk';
figure
for ct = 1:4
    plot(log10(Irel), Pc(idx_dir_flash{ct}, :)', 'color', color(ct))
    hold on
end
xlabel('log(R*/rod)')
ylabel('probability')

% compare across cell type
for ct = 1:4
    Pc_dir{ct} = Pc(idx_dir_flash{ct}, :);
    Pc_dir_mean(ct, :) = mean(Pc_dir{ct}, 1);
    Pc_dir_ste(ct, :) = std(Pc_dir{ct}, [], 1)/sqrt(length(idx_dir_flash{ct}));
end

figure
for ct = 1:4
    errorbar(log10(Irel), Pc_dir_mean(ct, :), Pc_dir_ste(ct, :), 'color', color(ct));
    hold on
end
legend('superior', 'Anterior', 'inferior')
xlabel('log(R*/rod)')
ylabel('probability')
%% intensity-response curve
window = 1;
for cc = 1:length(ds_id_flash)
    for ts = 1:length(trigger_set_i)
%         response(cc, ts) = sum(ds_flash_hist_mean{cc}{ts}(1:window/bin_size)) - sum(ds_flash_hist_mean{cc}{ts}(end-window/bin_size+1:end));
        response(cc, ts) = sum(ds_flash_hist_mean{cc}{ts}(1:window/bin_size)) - sum(ds_dark_hist_mean{cc}(1:window/bin_size));
    end
end
 
response_norm = response./repmat(max(response, [], 2), 1, length(trigger_set_i));

figure
for ct = 1:4
    plot(log10(Irel), response_norm(idx_dir{ct}, :), 'color', color(ct));
    hold on
end
%% raster plot

a = ceil((length(trigger_set_i)+1)/2);
for cc =2:length(ds_id_flash);
    ts = 1;
    b = 1;
    trial_n = length(ds_dark_raster{1});
    FigHandle = figure;
    set(FigHandle, 'Position', [0, 0, 1920, 1080]);
    while ts <= length(trigger_set_i)+1
        if  ts == a+1
            b = 3;
        end
        subplot(a, 4, b)
        if ts > 1
            trial_n = length(ds_flash_raster{cc}{ts-1});
        end
        for j = 1:trial_n
            if ts == 1
                SpikeTime = ds_dark_raster{cc}{j};
            else
                SpikeTime = ds_flash_raster{cc}{ts-1}{j};
            end
            SpikeTime = SpikeTime';
            X = [SpikeTime; SpikeTime];
            Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
            line(X, Y, 'color', 'b');
            axis([0, 3, 0, trial_n]);
            hold on
        end
        if b == 1
            title(num2str(ds_id_flash(cc)))
        end
        b = b+1;
        subplot(a, 4, b)
        if ts > 1
            bar(XX, ds_flash_hist_mean{cc}{ts-1}, 1)
        else
            bar(XX, ds_dark_hist_mean{cc}, 1)
        end
        xlim([0 3])

        b = b+3;
        ts = ts+1;
    end

    print_close(1, [24 12], num2str(ds_id_flash(cc)))

end


%% fit
Pc_temp = Pc;
Irel_temp = Irel;

% superior
ct = 1;
idx = idx_dir_flash{ct};
for cc = 1:length(idx)

    ydata = Pc_temp(idx(cc), :)-0.5;
    xdata = log10(Irel_temp)+4;
    [f, G] = fit_mm(xdata, ydata);
    fit_all{ct}{cc} = f;
    G_all{ct}{cc} = G;
end

% anterior
ct = 2;
idx = idx_dir_flash{ct};
for cc = 1:length(idx)

    ydata = Pc_temp(idx(cc), 1:11)-0.5;
    xdata = log10(Irel_temp(1:11))+4;
    [f, G] = fit_mm(xdata, ydata);
    fit_all{ct}{cc} = f;
    G_all{ct}{cc} = G;
end

% inferior
ct = 3;
idx = idx_dir_flash{ct};
for cc = 1:length(idx)

    ydata = Pc_temp(idx(cc), 1:10)-0.5;
    xdata = log10(Irel_temp(1:10))+4;
    [f, G] = fit_mm(xdata, ydata);
    fit_all{ct}{cc} = f;
    G_all{ct}{cc} = G;
end

% posterior
ct = 4;
idx = idx_dir_flash{ct};
for cc = 1:length(idx)

    ydata = Pc_temp(idx(cc), 1:14)-0.5;
    xdata = log10(Irel_temp(1:14))+4;
    [f, G] = fit_mm(xdata, ydata);
    fit_all{ct}{cc} = f;
    G_all{ct}{cc} = G;
end
xdata = log10(Irel_temp)+4;
x = linspace(min(xdata), max(xdata), 1000);

threshold = 0.34;
for ct = 1:4
    idx = idx_dir_flash{ct};
    for cc = 1:length(idx)
        f = fit_all{ct}{cc};
        y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a);
        [~, i] = min(abs(y-threshold));
%         xthreshold{ct}(cc) = 10^(x(i)-4);
        xthreshold{ct}(cc) = x(i)-4;
        rmse{ct}(cc) = G_all{ct}{cc}.rmse;
    end
end

color = 'brgk';
figure
for ct = 1:4
    plot(ct*ones(1, length(xthreshold{ct})), xthreshold{ct}, [color(ct) 'o'])
    hold on
    errorbar(ct+0.2, mean(xthreshold{ct}), std(xthreshold{ct})/sqrt(length(xthreshold{ct})), [color(ct) 'd']);
end
xlim([0.5 4.5])
ylabel('log(R*/rod)')
title('Pc = 0.83')
xtick = {'superior'; 'anterior'; 'inferior'; 'posterior'};
set(gca,'XTicklabel',xtick)

% fit

for ct = 1:4
    Pc_dir_mean(ct, :) = mean(Pc_dir{ct}, 1);
    Pc_dir_ste(ct, :) = std(Pc_dir{ct}, [], 1)/sqrt(size(Pc_dir{ct}, 1));
end

color = 'brgkc';
Pc_temp = Pc_dir_mean;
Pc_ste_temp = Pc_dir_ste;
Irel_temp = Irel;
a = 13;
figure
for ct = 1:4
    ydata = Pc_temp(ct, :)-0.5;
    xdata = log10(Irel_temp)+4;
    [f, G] = fit_mm(xdata, ydata, 'Startpoints', [0.5 2 0]);
    fit_avg{ct} = f;
    G_avg{ct} = G;

    x = linspace(min(xdata(1:a)), max(xdata(1:a)), 100);
    y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a);

    errorbar(log10(Irel(1:a)), Pc_dir_mean(ct, 1:a), Pc_dir_ste(ct, 1:a), [color(ct) 'o']);
    hold on
    h(ct) = plot(x-4, y+0.5, 'color', color(ct));
end
xlim([-inf 0.5])
legend([h(1), h(2), h(3), h(4)], 'superior', 'Anterior', 'inferior', 'posterior')
xlabel('Log(intensity)(R*/rod)')
ylabel('% correct')


%% map ei with fluorescent image

im_s = imread('/Volumes/lab/Experiments/Array/Images/2016-12-08-0/WN.jpg'); % load stimulus picture taken by camera
im_array = imread('/Volumes/lab/Experiments/Array/Images/2016-12-08-0/10X2.jpg'); % load array image taken by camera
stixel_size = 30; % frame shown in WN.jpg
movie_path = '/Volumes/lab/acquisition/movie-xml/BW-30-6-0.48-11111-20x20-60.35.xml';

mov = get_movie(movie_path, 0, 1);
mov_frame = matrix_scaled_up(squeeze(mov(:,:,1)), stixel_size);
clear movingPoints fixedPoints
cpselect(im_s, mov_frame) % select 4 control points

%% register two images
tform = fitgeotrans(movingPoints, fixedPoints, 'projective');
registered = imwarp(im_s, tform,'OutputView',imref2d(size(mov_frame)));
figure 
imshow(registered);
figure
imshowpair(mov_frame,registered,'blend');

% transform array image into display coordinates
registered_array = imwarp(im_array, tform, 'OutputView', imref2d(size(mov_frame)));
figure
imshow(im_array);

% get array location in display coordinates

%                 EI                               DISPLAY
%
%               195 (1)                         386(5)  264(6)
%                 / \                               ______
%               /     \                            /      \
%   264 (6)    |       |    126 (2)               /        \
%   386 (5)    |       |    4   (3)       455(4)  \        / 195(1)
%               \     /                            \      /
%                 \ /                               ------
%                455 (4)                          4(3)  126(2)
array_location_display = ginput;

% get array location in ei coordinates
elec_corner = [195 126 4 455];
array_location_ei = dataDG{5}.ei.position(elec_corner,:);
Tform = maketform('projective', array_location_ei, array_location_display);
test = tformfwd(Tform, array_location_ei)-array_location_display % should be equal or close to zeros

ds_idx = get_cell_indices(dataDG{5}, ds_id);
frame = dataDG{5}.ei.nrPoints + dataDG{5}.ei.nlPoints + 1;
elec = size(dataDG{5}.ei.position, 1);
distance = 1;
center_ei_ds = zeros(length(ds_id),2);
for i = 1:length(ds_id)
    distance_temp = distance;
    full_elec_n = sum(1:distance_temp)*6+1; % Assume all arrays have hexagonal-arranged electrodes
    ei = dataDG{5}.ei.eis{ds_idx(i)};
    ei = ei';
    [~,I] = max(abs(ei(:)));
    elec_n = ceil(I/frame);
    elecs_n = get_ei_neighbors(elec_n, elec, distance_temp);
    while(length(elecs_n) < full_elec_n)
        distance_temp = distance_temp - 1;
        elecs_n = get_ei_neighbors(elec_n, elec, distance_temp);
        full_elec_n = sum(1:distance_temp)*6+1;
    end
    points = dataDG{5}.ei.position(elecs_n,:);
    mass = max(abs(ei(:,elecs_n)));
    center_ei_ds(i,:) = centroid(points, mass);
end
center_ei_ds = tformfwd(Tform, center_ei_ds);


%% 
center_ei_ds_s = center_ei_ds(idx_dir{1}, :);
figure(2)
for i = 1:length(id_dir{1})
    imshow(im_array);
    axis off
    hold on
    plot(center_ei_ds_s(i, 1), center_ei_ds_s(i, 2), 'ro')
    pause
    hold off
end

%%
pos = dataDG{5}.ei.position;
mode = 'neg';

figure
imshow(im_array);
array_location_image = ginput;
elec_corner = [195 126 4 455];
array_location_ei = pos(elec_corner,:);
Tform = maketform('projective', array_location_image, array_location_ei);
test = tformfwd(Tform, array_location_image)-array_location_ei % should be equal or close to zeros

% get coordinates of GFP cell bodies
soma_location_image = ginput_label('r');
soma_location_ei = tformfwd(Tform, soma_location_image);
% get array location in display coordinates

%                 EI                               DISPLAY
%
%               195 (1)                         386(5)  264(6)
%                 / \                               ______
%               /     \                            /      \
%   264 (6)    |       |    126 (2)               /        \
%   386 (5)    |       |    4   (3)       455(4)  \        / 195(1)
%               \     /                            \      /
%                 \ /                               ------
%                455 (4)                          4(3)  126(2)

%% distance distribution

for ct = 1:length(id_dir)
    for cc = 1:length(id_dir{ct})
        id = id_dir{ct}(cc);
        idx = get_cell_indices(dataDG{5}, id);
        ei = dataDG{5}.ei.eis{idx};
        com = ei_com_xy(ei, pos, 30*3, mode);
        com_oo{ct}(cc, :) = com;
        com_oo_image{ct}(cc, :) = tforminv(Tform, com);
        dis = sqrt(sum((soma_location_ei - repmat(com, size(soma_location_ei, 1), 1)).^2, 2));
        [disMin_oo{ct}(cc), matchI{ct}(cc)] = min(dis);
    end
end

%% plot com of ei on microscopy image
figure
imshow(im_array);
hold on
plot(com_oo_image{1}(:, 1), com_oo_image{1}(:, 2), 'o', 'color', 'r', 'MarkerSize', 10)
axis off

XX = 0:10:120;
figure
a = hist(disMin_oo{1}, XX);
bar(XX, a, 1)
xlabel('um')
ylabel('number of cells')
xlim([0 150])
ylim([0 8])

kop = sum(disMin_oo{1} < 50)/length(disMin_oo{1});
figure
bar([1 2], [kop 1-kop]*100, 1)
xlim([0 3])
