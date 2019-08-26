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
dataDG{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/data000-map/data000-map', opt);
dataDG{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/stimuli/s00.txt';
dataDG{1} = load_stim(dataDG{1}, 'user_defined_trigger_interval', 10);
dataDG{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/data002-map/data002-map', opt);
dataDG{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/stimuli/s02.txt';
dataDG{2} = load_stim(dataDG{2}, 'user_defined_trigger_interval', 10);
dataDG{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/data004-map/data004-map', opt);
dataDG{3}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/stimuli/s04.txt';
dataDG{3} = load_stim(dataDG{3}, 'user_defined_trigger_interval', 10);
dataDG{4} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/data005-sorted/data005-sorted', opt);
dataDG{4}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/stimuli/s05.txt';
dataDG{4} = load_stim(dataDG{4}, 'user_defined_trigger_interval', 10);
dataDG{4} = load_ei(dataDG{4}, 'all', 'array_type', 519);
dataDG{5} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/data009-map/data009-map', opt);
dataDG{5}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/stimuli/s09.txt';
dataDG{5} = load_stim(dataDG{5}, 'user_defined_trigger_interval', 10);

datadg = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/data007-map/data007-map', opt);
datadg.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/stimuli/s07.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

datafs{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/data001-map/data001-map', opt);
datafs{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/stimuli/s01.mat';
datafs{1} = load_stim_mfs(datafs{1});
datafs{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/data003-map/data003-map', opt);
datafs{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/stimuli/s03.mat';
datafs{2} = load_stim_mfs(datafs{2});
datafs{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/data008-map/data008-map', opt);
datafs{3}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/stimuli/s08.mat';
datafs{3} = load_stim_mfs(datafs{3});

datawn = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-12-12-0/data006-map/data006-map', opt);
datawn = load_ei(datawn, ds_id);
load('DS161212.mat')

%% dg
n = 1;
i = 1;
[raster_dg, dg, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
dg{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);


delta_p = 1; % choose which params to use to calculate prefer direction indices 
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
color = 'brgkc';

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
%                     if sum(sum(data > mean(data(:))+3*std(data(:))))>0
%                         figure(100)
%                         imagesc(data)
%                         colormap gray
%                         pause
                        
                        params = fit_2d_gaussian(data);
    %                     Gaussian_params{ll}{dir}{cc}{onoff} = params;
                        rf_area_temp = [rf_area_temp params.xd * params.yd * pi * PixelArea];
%                     end
                end
            end
            rf_area{ll}{dir}{onoff} = rf_area_temp;
        end
    end
end


% exclude outliers
% stdn = 2;
stdn = 100;

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


%% map ei with fluorescent image

im_s = imread('/Volumes/lab/Experiments/Array/Images/2016-12-12-0/WN.jpg'); % load stimulus picture taken by camera
im_array = imread('/Volumes/lab/Experiments/Array/Images/2016-12-12-0/10X2.jpg'); % load array image taken by camera
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
array_location_ei = dataDG{4}.ei.position(elec_corner,:);
Tform = maketform('projective', array_location_ei, array_location_display);
test = tformfwd(Tform, array_location_ei)-array_location_display % should be equal or close to zeros

ds_idx = get_cell_indices(dataDG{4}, ds_id);
frame = dataDG{4}.ei.nrPoints + dataDG{4}.ei.nlPoints + 1;
elec = size(dataDG{4}.ei.position, 1);
distance = 1;
center_ei_ds = zeros(length(ds_id),2);
for i = 1:length(ds_id)
    distance_temp = distance;
    full_elec_n = sum(1:distance_temp)*6+1; % Assume all arrays have hexagonal-arranged electrodes
    ei = dataDG{4}.ei.eis{ds_idx(i)};
    ei = ei';
    [~,I] = max(abs(ei(:)));
    elec_n = ceil(I/frame);
    elecs_n = get_ei_neighbors(elec_n, elec, distance_temp);
    while(length(elecs_n) < full_elec_n)
        distance_temp = distance_temp - 1;
        elecs_n = get_ei_neighbors(elec_n, elec, distance_temp);
        full_elec_n = sum(1:distance_temp)*6+1;
    end
    points = dataDG{4}.ei.position(elecs_n,:);
    mass = max(abs(ei(:,elecs_n)));
    center_ei_ds(i,:) = centroid(points, mass);
end
center_ei_ds = tformfwd(Tform, center_ei_ds);

%% 
center_ei_ds_s = center_ei_ds(idx_dir{1}, :);
figure(2)
imshow(im_array);
axis off
hold on
for i = 1:length(id_dir{1})
    plot(center_ei_ds_s(i, 1), center_ei_ds_s(i, 2), 'ro')
    hold on
%     pause
%     hold off
end

%%
pos = dataDG{4}.ei.position;
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
        idx = get_cell_indices(dataDG{4}, id);
        ei = dataDG{4}.ei.eis{idx};
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

%%
xx = [0:pi/4:7*pi/4] - pi;
xx = xx/pi*180+22.5;

gfp = disMin_oo{1} < 50;
gfp_rho = rho_dg{1}{1}(gfp, :);
gfp_n_rho = rho_dg{1}{1}(~gfp, :);
for ll = 1:5
    if ll == 5
        gfp_temp = gfp(2:end);
    else
        gfp_temp = gfp;
    end
    gfp_dsi{ll} = dsi_dg{ll}{1}(gfp_temp);
    gfp_n_dsi{ll} = dsi_dg{ll}{1}(~gfp_temp);
end
figure
errorbar(xx, mean(gfp_rho), std(gfp_rho)/sqrt(size(gfp_rho, 1)), 'b')
hold on
errorbar(xx, mean(gfp_n_rho), std(gfp_n_rho)/sqrt(size(gfp_n_rho, 1)), 'r')

%%
load DS161212.mat
gfp_rho_1212 = gfp_rho;
gfp_n_rho_1212 = gfp_n_rho;
gfp_dsi_1212 = gfp_dsi{1};
gfp_n_dsi_1212 = gfp_n_dsi{1};
load DS161208.mat
gfp_rho_1208 = gfp_rho;
gfp_n_rho_1208 = gfp_n_rho;
gfp_dsi_1208 = gfp_dsi{1};
gfp_n_dsi_1208 = gfp_n_dsi{1};

gfp_rho = [gfp_rho_1208; gfp_rho_1212];
gfp_n_rho = [gfp_n_rho_1208; gfp_n_rho_1212];
gfp_dsi = [gfp_dsi_1208; gfp_dsi_1212];
gfp_n_dsi = [gfp_n_dsi_1208; gfp_n_dsi_1212];

figure
errorbar(xx, mean(gfp_rho), std(gfp_rho)/sqrt(size(gfp_rho, 1)), 'b')
hold on
errorbar(xx, mean(gfp_n_rho), std(gfp_n_rho)/sqrt(size(gfp_n_rho, 1)), 'r')


for ll = 1:5
    gfp_dsi_mean(ll) = mean(gfp_dsi{ll});
    gfp_dsi_ste(ll) = std(gfp_dsi{ll})/sqrt(size(gfp_dsi{ll}, 1));
    gfp_n_dsi_mean(ll) = mean(gfp_n_dsi{ll});
    gfp_n_dsi_ste(ll) = std(gfp_n_dsi{ll})/sqrt(size(gfp_n_dsi{ll}, 1));
end
figure
errorbar([0:4], gfp_dsi_mean, gfp_dsi_ste, 'b')
hold on
errorbar([0:4], gfp_n_dsi_mean, gfp_n_dsi_ste, 'r')
ylim([0 1])


%% cross correlation
duration = 2700;
bin_size = 0.00025;
max_lag = 40;
xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
N = 10000;
corr_cells_test = [];

for  c1 = 1:length(id_dir{1})-1
%     FigHandle = figure;
%     set(FigHandle, 'Position', [1 1 2000 2000])
    for c2 = c1+1:length(id_dir{1})
        id1 = id_dir{1}(c1);
        id2 = id_dir{1}(c2);
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

%         A = xcorr(spikes1, spikes2, max_lag, 'coeff');
        A = xcorr(spikes1, spikes2, max_lag);
        [maxv(c1, c2), maxi(c1, c2)] = max(A);
        a = round(0.001/bin_size)+max_lag;
        b = conv(A, ones(1, 11), 'valid');
        ratio(c1, c2) = (sum(A(a:a+10)) + sum(A(max_lag*2-a-10:max_lag*2-a)) - min(A)*22)/(min(b)*2 - min(A)*22);
        [h, filteredA] = find_smallest_h(A);
        if sum(round(filteredA*sum(A)/sum(filteredA))) > 0
            [bootstat,bootsam] = bootstrp(N,@find_smallest_h_hist,rude(round(filteredA*sum(A)/sum(filteredA)), 1:max_lag*2+1), max_lag);
            p(c1, c2) = sum(bootstat > h)/N;
        else
            p(c1, c2) = 1;
        end
%         subplot(4, 5, c2)
        if p(c1, c2) < 0.05 && ratio(c1, c2) > 2 && maxi(c1, c2) > 0.75*max_lag && maxi(c1, c2) < 1.25*max_lag+1 && maxv(c1, c2) > 10
%            bar(xx, A, 'r')
           corr_cells_test = [corr_cells_test; id1 id2];
        else
%            bar(xx, A, 'b')
        end

%         title([num2str(id1) '  ' num2str(id2)])
%         xlim([-0.01 0.01])

    end
%     print_close(1, [24 12], num2str(c1));
    c1
end

%% neighboring pairs
% id_dir{1}(9) = [];
pos = datawn.ei.position;
mode = 'neg';
neighbors = [];
corner_i = [4 126 195 264 386 455 4];
corner_position = datawn.ei.position(corner_i, :);

for cc1 = 1:length(id_dir{1})
    for cc2 = cc1+1:length(id_dir{1})
        id1 = id_dir{1}(cc1);
        idx1 = get_cell_indices(datawn, id1);
        ei1 = datawn.ei.eis{idx1};
        com1 = ei_com_xy(ei1, pos, 30*3, mode);
        id2 = id_dir{1}(cc2);
        idx2 = get_cell_indices(datawn, id2);
        ei2 = datawn.ei.eis{idx2};
        com2 = ei_com_xy(ei2, pos, 30*3, mode);
        if pdist([com1;com2]) < 150
            neighbors = [neighbors; id1 id2];
        end
    end
end

coms = [];
for cc = 1:length(id_dir{1})
    id = id_dir{1}(cc);
    idx = get_cell_indices(datawn, id);
    ei = datawn.ei.eis{idx};
    com = ei_com_xy(ei, pos, 30*3, mode);
    coms = [coms; com];
end

figure
for cc = 1:length(id_dir{1})
    plot(coms(cc, 1), coms(cc, 2),'ko')
    hold on
    text(coms(cc, 1)+5, coms(cc, 2)+5, num2str(id_dir{1}(cc)), 'FontSize', 10)
    
end


for cp = 1:size(corr_cells_test)
    idx1 = find(id_dir{1} == corr_cells_test(cp, 1));
    idx2 = find(id_dir{1} == corr_cells_test(cp, 2));
    plot([coms(idx1, 1), coms(idx2, 1)], [coms(idx1, 2), coms(idx2, 2)], 'k');
end
plot(corner_position(:, 1), corner_position(:, 2), 'color', [.5 .5 .5])
axis off
%%
ct = 1;
cp_i = [];
cn = 0;
for c1 = 1:length(id_dir{ct})-1
    for c2 = c1+1:length(id_dir{ct})
        if norm([coms(c1, :) - coms(c2, :)]) < 150
            cn = cn + 1;
            cp_i = [cp_i; c1 c2];
        end
    end
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
A_all = [];

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
    A_all = [A_all A];
    subplot(4, 4, cp)
    c1 = cp_i(cp, 1);
    c2 = cp_i(cp, 2);
    if p(c1, c2) < 0.05 && ratio(c1, c2) > 2 && maxi(c1, c2) > 0.75*max_lag && maxi(c1, c2) < 1.25*max_lag+1 && maxv(c1, c2) > 10
       bar(xx, A, 'k')
    else
       bar(xx, A, 'r')
    end
    xlim([-0.01 0.01])
    title([num2str(id1) '  ' num2str(id2)])
%     cp_index_ko(cp) = (max(A(8:10)) + max(A(12:14)) - 2*min(A)) / (2*(A(11) - min(A)));

end

load('DS161017.mat', 'cp_index_wt')
figure
subplot(2, 1, 1)
hist(cp_index_wt,20)
xlabel('ratio')
ylabel('# of cell pairs')
title('WT')
subplot(2, 1, 2)
hist(cp_index_ko, 20)
xlabel('ratio')
ylabel('# of cell pairs')
title('FACx')


%%
A_facx_mean = mean(A_all_facx, 2);
A_facx_ste = std(A_all_facx, [], 2)/sqrt(size(A_all_facx, 2));
figure
patch([xx fliplr(xx)], [A_facx_mean + A_facx_ste; flipud(A_facx_mean - A_facx_ste)]', [1 1 1]*0.8)
hold on
plot(xx, A_facx_mean, 'k')
% axis off
% ylim([0 0.015])

%% estimate background activity
n = 5;
interval = 1;
raster_interval_dg = deal(cell(n, 1));
for i = 1:n    
    raster_interval_dg{i} = get_ds_interval_raster(dataDG{i}, ds_id, interval);
    bgfr{i} = zeros(1,length(ds_id));
    for j = 1:length(raster_interval_dg{i})
        if iscell(raster_interval_dg{i}{j})
            bgfr{i}(j) = mean(cellfun(@length, raster_interval_dg{i}{j}))/interval;
        end
    end
    for dir = 1:4
        bgfr_ct{i}{dir} = [];
        for cc = 1:length(idx_dir{dir})
            if ~dg_idx(idx_dir{dir}(cc), i) && sum(DG{i}.rho{1}(idx_dir{dir}(cc), :))>0
                bgfr_ct{i}{dir} = [bgfr_ct{i}{dir} bgfr{i}(idx_dir{dir}(cc))];
            end
        end
    end
    bgfr_ct_mean(i, 1) = mean(bgfr_ct{i}{1});
    bgfr_ct_ste(i, 1) = std(bgfr_ct{i}{1})/sqrt(length(bgfr_ct{i}{1}));
    temp = cell2mat(bgfr_ct{i}(2:4));
    bgfr_ct_mean(i, 2) = mean(temp);
    bgfr_ct_ste(i, 2) = std(temp)/sqrt(length(temp));
end

celltype = {'superior', 'anterior', 'inferior', 'posterior'};
figure
for i = 1:2
    errorbar(10.^[0:4], bgfr_ct_mean(:, i), bgfr_ct_ste(:, i));
    hold on
end
set(gca, 'XScale', 'log')
xlim([0.1 10^5])
xlabel('R*/rod/s')
ylabel('background firing (Hz)')
legend('superior', 'others')
title('FACx')

n = 5;
% spike count
DG_bgnd = cell(n, 1);
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim(dataDG{i},ds_id,0, 0.05);
    DG_bgnd{i} = sort_direction(dscellanalysis_bgnd_subtract(NumSpikesCell, StimComb, dataDG{i}, bgfr{i}));
end

%% DS tuning curves (drifting grating)
% all ds cells

% t = 2;
dirn = 4;
D = 5;
T = 1;
color = 'bkrgc';
ct = {'superior', 'anterior', 'inferior', 'posterior'};
ll = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};
DG_cut = DG_bgnd;

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
        dsi_dg{d}{i} = [];
        for cc = 1:length(idx_dir{i})
            if ~dg_idx(idx_dir{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx_dir{i}(cc), :))>0
            [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
            xsort = xsort/pi*180;
            y_temp = DG_cut{d}.rho{T}(idx_dir{i}(cc), :);
%             plot(xsort, y_temp(seq), color(i))
%             ylim([0 1])
% %             pause
%             hold on
            rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
            dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx_dir{i}(cc))];
            end
        end
        rho_dg_mean{d}(i, :) = mean(rho_dg{d}{i});
        rho_dg_ste{d}(i, :) = std(rho_dg{d}{i})/sqrt(size(rho_dg{d}{i}, 1));
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


%%
xx = [0:pi/4:7*pi/4] - pi;
xx = xx/pi*180+22.5;

gfp = disMin_oo{1} < 50;
gfp_rho = rho_dg{1}{1}(gfp, :);
gfp_n_rho = rho_dg{1}{1}(~gfp, :);
for ll = 1:5
    if ll == 5
        gfp_temp = gfp(2:end);
    else
        gfp_temp = gfp;
    end
    gfp_dsi{ll} = dsi_dg{ll}{1}(gfp_temp);
    gfp_n_dsi{ll} = dsi_dg{ll}{1}(~gfp_temp);
end
figure
errorbar(xx, mean(gfp_rho), std(gfp_rho)/sqrt(size(gfp_rho, 1)), 'b')
hold on
errorbar(xx, mean(gfp_n_rho), std(gfp_n_rho)/sqrt(size(gfp_n_rho, 1)), 'r')

rho_dg_bgnd = rho_dg;
dsi_dg_bgnd = dsi_dg;
gfp_rho_bgnd = gfp_rho;
gfp_n_rho_bgnd = gfp_n_rho;
gfp_dsi_bgnd = gfp_dsi;
gfp_n_dsi_bgnd = gfp_n_dsi;