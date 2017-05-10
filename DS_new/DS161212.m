%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load drifting grating data
datadg = load_data('/Volumes/lab/analysis/2016-12-12-0/data005-sorted/data005-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-12-12-0/stimuli/s05.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

% identify DS cells
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [1 2]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

% clear datadg
dataDG{1} = load_data('/Volumes/lab/analysis/2016-12-12-0/data000-map/data000-map', opt);
dataDG{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-12-12-0/stimuli/s00.txt';
dataDG{1} = load_stim(dataDG{1}, 'user_defined_trigger_interval', 10);
dataDG{2} = load_data('/Volumes/lab/analysis/2016-12-12-0/data002-map/data002-map', opt);
dataDG{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-12-12-0/stimuli/s02.txt';
dataDG{2} = load_stim(dataDG{2}, 'user_defined_trigger_interval', 10);
dataDG{3} = load_data('/Volumes/lab/analysis/2016-12-12-0/data004-map/data004-map', opt);
dataDG{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-12-12-0/stimuli/s04.txt';
dataDG{3} = load_stim(dataDG{3}, 'user_defined_trigger_interval', 10);
dataDG{4} = load_data('/Volumes/lab/analysis/2016-12-12-0/data005-sorted/data005-sorted', opt);
dataDG{4}.names.stimulus_path = '/Volumes/lab/analysis/2016-12-12-0/stimuli/s05.txt';
dataDG{4} = load_stim(dataDG{4}, 'user_defined_trigger_interval', 10);
% dataDG{4} = load_ei(dataDG{4}, 'all');
dataDG{5} = load_data('/Volumes/lab/analysis/2016-12-12-0/data009-map/data009-map', opt);
dataDG{5}.names.stimulus_path = '/Volumes/lab/analysis/2016-12-12-0/stimuli/s09.txt';
dataDG{5} = load_stim(dataDG{5}, 'user_defined_trigger_interval', 10);

% datadg = load_data('/Volumes/lab/analysis/2016-12-12-0/data007-map/data007-map', opt);
% datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-12-12-0/stimuli/s07.txt';
% datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

datafs{1} = load_data('/Volumes/lab/analysis/2016-12-12-0/data001-map/data001-map', opt);
datafs{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-12-12-0/stimuli/s01.mat';
datafs{1} = load_stim_mfs(datafs{1});
datafs{2} = load_data('/Volumes/lab/analysis/2016-12-12-0/data003-map/data003-map', opt);
datafs{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-12-12-0/stimuli/s03.mat';
datafs{2} = load_stim_mfs(datafs{2});
datafs{3} = load_data('/Volumes/lab/analysis/2016-12-12-0/data008-map/data008-map', opt);
datafs{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-12-12-0/stimuli/s08.mat';
datafs{3} = load_stim_mfs(datafs{3});

load('DS161212.mat')

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

t = 5;
figure
compass(dg{1}.U{t}(idx_sub{1}), dg{1}.V{t}(idx_sub{1}), 'r')
hold on
compass(dg{1}.U{t}(idx_sub{2}), dg{1}.V{t}(idx_sub{2}), 'b')


%%
d = 1;
t = 1;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(dg{d}.U{t}, dg{d}.V{t});
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(dg{d}.U{t}, dg{d}.V{t}, x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = I;
    id_dir{i} = ds_id(idx_dir{i});
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

for dir = 1:4
    for cc = 1:length(id_dir{dir})
        plot_ds_raster(DG, raster_dg, idx_dir{dir}(cc), id_dir{dir}(cc), ll, 2, 3, 1)
    end
end

%% DS tuning curves (drifting grating)
% all ds cells
ll = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};
ct = {'superior', 'anterior', 'inferior', 'posterior'};
% t = 2;
dirn = 4;
D = 4;
T = 1;

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
        rho_dg_mean{d}(i, :) = mean(rho_dg{d}{i}, 1);
        rho_dg_ste{d}(i, :) = std(rho_dg{d}{i}, [], 1)/sqrt(size(rho_dg{d}{i}, 1));
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
legend(ll)
% DSI
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = ct;
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
    hist(dsi_dg{ll}{dir})
end
end

%% compare with WT
load('DS150603.mat', 'dsi_dg_wt')
figure
X = 0.025:0.05:0.975;
dsi_wt = hist(dsi_dg_wt{1}{1}, X);
dsi_ko = hist(dsi_dg{1}{1}, X);
bar(X, [dsi_wt' dsi_ko'], 1, 'stacked')
xlim([0 1])
xlabel('DSI')
ylabel('cell#')
legend('WT', 'KO')
title('NDF 4')

%% plot rfs
cell_type = {'superior', 'anterior', 'inferior', 'posterior'};
field_width = 17; field_height = 17;
field_width_sta = 40; field_height_sta = 40;
subregion = 0;
stop = 0.5; %second
for i = 1:3
%     fs_raster{i} = get_fs_raster(datafs{i}, ds_id, 'stop', 0.5);
    fs_raster{i} = get_fs_raster(datafs{i}, ds_id);
    for cc = 1:length(ds_id)
        if fs_idx(cc,i)
            fs_raster{i}{cc} = [];
        end
    end
    fs_spike{i} = fs_get_spike(fs_raster{i});
    [rf_all{i}, rf_std{i}] = get_fs_rf(fs_spike{i}, field_width, field_height,subregion);
end
%%
for cc = 1:length(ds_id)
    figure(1)
    set(gcf, 'Position', [1 1 1000 1000])
    id = ds_id(cc);
    for i = 1:3
        if ~isempty(rf_all{i}{cc})
%             rf = padarray(rf_all{i}{cc},[7,7]);
            rf = rf_all{i}{cc};
    
            subplot(3,3,3*(i-1)+1)
            imagesc(sum(rf,3))
            colormap gray
            axis image
            axis off

            subplot(3,3,3*(i-1)+2)
            imagesc(rf(:,:,1))
            colormap gray
            axis image
            axis off

            subplot(3,3,3*i)
            imagesc(rf(:,:,2))
            colormap gray
            axis image
            axis off

        end
    end
    print_close(1,[24 12],num2str(id))
%     pause
%     close(1)
end

%% fit RF with Gaussian                
PixelArea = (30*4)^2/10^6;
% fit and compute rf area
for ll = 1:3
    for dir = 1:4
        clear rf_area_temp
        for onoff = 1:2
            rf_area_temp = [];
            for cc = 1:length(id_dir{dir})
                if ~fs_idx(idx_dir{dir}(cc), ll)
                    data = rf_all{ll}{idx_dir{dir}(cc)}(:, :, onoff);
%                     data = rf_wt{ll}{dir}{cc}{onoff};
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
            rf_area{ll}{dir}{onoff} = rf_area_temp;
        end
    end
end


% exclude outliers
stdn = 2;
for ll = 1:3
    for dir = 1:4
        for onoff = 1:2
            notdone = 1;
            rf_area_temp = rf_area{ll}{dir}{onoff};
            while notdone
                a = length(rf_area_temp);
                rf_area_temp(rf_area_temp > std(rf_area_temp)*stdn + mean(rf_area_temp)) = [];
                b = length(rf_area_temp);
                if a == b
                    notdone = 0;
                    rf_area_clean{ll}{dir}{onoff} = rf_area_temp;
                end
            end
            rf_area_clean_mean{onoff}(ll, dir) = mean(rf_area_clean{ll}{dir}{onoff});
            rf_area_clean_ste{onoff}(ll, dir) = std(rf_area_clean{ll}{dir}{onoff})/sqrt(length(rf_area_clean{ll}{dir}{onoff}));
        end
    end
end


% plot 
color = 'brgkc';
figure
for onoff = 1:2
    subplot(1,2,onoff)
    for dir = 4:-1:1
        for ll = 1:3
            if ~isempty(rf_area_clean{ll}{dir}{onoff})
                n = size(rf_area_clean{ll}{dir}{onoff}, 1);
                h{ll} = plot((dir-1)*5+ll*ones(n,1), rf_area_clean{ll}{dir}{onoff}, [color(ll) 'o']);
    %             n = length(rf_area{ll}{dir}{onoff});
    %             h{ll} = plot((dir-1)*5+ll*ones(n,1), rf_area{ll}{dir}{onoff}, [color(ll) 'o']);
                hold on
            end
        end
    end
%     set(gca, 'yscale', 'log')
    legend([h{1}(1), h{2}(1), h{3}(1)], 'NDF 4', 'NDF 2', 'NDF 0')
    if onoff == 1
        title('ON')
    else
        title('OFF')
    end
    ylabel('RF area (mm^2)')
    set(gca, 'xtick', [])
%         ylim([0 0.3])

end

figure
for onoff = 1:2
    subplot(1,2,onoff)
    for dir = 1:4
        errorbar(0:2:4, rf_area_clean_mean{onoff}(:, dir), rf_area_clean_ste{onoff}(:, dir), 'color', color(dir))
        hold on
    end
    ylim([0 0.14])
    xlim([-1 5])
    xlabel('log(background intensity)')
    ylabel('RF area (mm^2)')
    legend(cell_type)
    if onoff == 1
        title('ON')
    else
        title('OFF')
    end

end
%%
stdn = 2;
for ll = 1:3
    for onoff = 1:2
        rf_area_temp = [];
        for dir = 2:4
            rf_area_temp = [rf_area_temp rf_area{ll}{dir}{onoff}];
        end
        notdone = 1;
        while notdone
            a = length(rf_area_temp);
            rf_area_temp(rf_area_temp > mean(rf_area_temp) + std(rf_area_temp)*stdn) = [];
            rf_area_temp(rf_area_temp < mean(rf_area_temp) - std(rf_area_temp)*stdn) = [];
            b = length(rf_area_temp);
            if a == b
                notdone = 0;
                rf_area_clean_all{ll}{onoff} = rf_area_temp;
            end
        end
        rf_area_clean_all_mean{onoff}(ll) = mean(rf_area_clean_all{ll}{onoff});
        rf_area_clean_all_ste{onoff}(ll) = std(rf_area_clean_all{ll}{onoff})/sqrt(length(rf_area_clean_all{ll}{onoff}));
    end
end

figure
for onoff = 1:2
    errorbar(0:2:4, rf_area_clean_all_mean{onoff}, rf_area_clean_all_ste{onoff})
    hold on
%     ylim([0 0.05])
    xlim([-1 5])
    xlabel('log(background intensity)')
    ylabel('RF area (mm^2)')
    legend('ON', 'OFF')
end

%% compare with WT
load('DS150603.mat', 'dsi_dg_wt')
dsi_0603 = dsi_dg_wt{1}{1};
% load('DS150703-1.mat', 'dsi_dg')
% dsi_0703 = dsi_dg{1}{4};
load('DS160130.mat', 'dsi_dg')
dsi_0130 = dsi_dg{1}{1};
load('DS161208.mat', 'dsi_dg')
dsi_1208 = dsi_dg{1}{1};
load('DS161212.mat', 'dsi_dg')
dsi_1212 = dsi_dg{1}{1};

figure
X = 0.025:0.05:0.975;
dsi_wt1 = hist(dsi_0603, X);
dsi_wt2 = hist(dsi_0130, X);

dsi_ko1 = hist(dsi_1208, X);
dsi_ko2 = hist(dsi_1212, X);

bar(X, [dsi_wt1' dsi_wt2' dsi_ko1' dsi_ko2'], 1, 'stacked')
xlim([0 1])
xlabel('DSI')
ylabel('cell#')
legend('WT', 'KO')
title('NDF 4')

figure
X = 0.025:0.05:0.975;
dsi_wt = hist(dsi_0603, X);

dsi_ko = hist(dsi_1208, X);

bar(X, [dsi_wt' dsi_ko'], 1, 'stacked')
xlim([0 1])
xlabel('DSI')
ylabel('cell#')
legend('WT', 'KO')
title('NDF 4')
