cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);
dirpath = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/';

datadg{1} = load_data(strcat(dirpath, 'data000-map/data000-map'), opt);
datadg{1}.names.stimulus_path = strcat(dirpath, 'stimuli/s00.mat');
datadg{1} = load_stim_matlab(datadg{1}, 'user_defined_trigger_interval', 10);
datadg{2} = load_data(strcat(dirpath, 'data003-map/data003-map'), opt);
datadg{2}.names.stimulus_path = strcat(dirpath, 'stimuli/s03.mat');
datadg{2} = load_stim_matlab(datadg{2}, 'user_defined_trigger_interval', 10);
datadg{3} = load_data(strcat(dirpath, 'data006-map/data006-map'), opt);
datadg{3}.names.stimulus_path = strcat(dirpath, 'stimuli/s06.mat');
datadg{3} = load_stim_matlab(datadg{3}, 'user_defined_trigger_interval', 10);
datadg{4} = load_data(strcat(dirpath, 'data009-map/data009-map'), opt);
datadg{4}.names.stimulus_path = strcat(dirpath, 'stimuli/s09.mat');
datadg{4} = load_stim_matlab(datadg{4}, 'user_defined_trigger_interval', 10);
datadg{5} = load_data(strcat(dirpath, 'data012-map/data012-map'), opt);
datadg{5}.names.stimulus_path = strcat(dirpath, 'stimuli/s12.mat');
datadg{5} = load_stim_matlab(datadg{5}, 'user_defined_trigger_interval', 10);
datadg{5} = load_ei(datadg{5}, 'all', 'array_id', 1551);
%% classification

i = 5; % high light level, NDF 0
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg{i},datadg{i}.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg{i});
close all
% figure 1C
params_idx = [4 5]; % which parameters to use for classification
[ds_id, nonds_id] = classify_ds(datadg{i}, ds_struct, params_idx);
title('Normalized Vector Sum Amplitude')
xlabel('High speed - 480 micron/s')
ylabel('Low speed - 240 micron/s')
% get rasters
load('DS150703-1.mat', 'dg_idx')
n = 5;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg{i},ds_id,0,1);
    DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg{i}));
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
close(1);close(3:100);


%% fig 1A
ll = 5; t = 4; cc = 23;
[idx, xx, yy] = subplot_idx(1, 1);
tt = DG{1}.theta{1}(1, :);
FigHandle = figure;
set(FigHandle, 'Position', [1 1 400 400])
h = subplot(xx, yy, idx(1)); 
u_temp = DG{ll}.U{t}(cc);
v_temp = DG{ll}.V{t}(cc);
P = polar(0, 3);
set(P, 'Visible', 'off')
hold on
compass(DG{ll}.U{t}(cc), DG{ll}.V{t}(cc), 'r');
polar(tt, DG{ll}.rho{t}(cc, :), 'b');
polar_theta_off(h)
for i = 2:9
    subplot(xx, yy, idx(i)); plot_raster(squeeze(raster_dg{ll}{cc}(1, t, i-1, :)), 0, 8);
    if idx(i) == 7
        ylabel('trials')
        xlabel('time (s)')
    end
end
subplot(3, 3, 5)
tt = linspace(0, 2*pi, 100);
polar(tt, ones(1, 100), 'k');
polar(tt, 2*ones(1, 100), 'k');
polar(tt, 3*ones(1, 100), 'k');

%% fig 1B
ll = 5; t = 4; cc = 60;
[idx, xx, yy] = subplot_idx(1, 1);
tt = DG_nds{1}.theta{1}(1, :);
FigHandle = figure;
set(FigHandle, 'Position', [1 1 400 400])
h = subplot(xx, yy, idx(1)); 
u_temp = DG_nds{ll}.U{t}(cc);
v_temp = DG_nds{ll}.V{t}(cc);
P = polar(0, 3);
set(P, 'Visible', 'off')
hold on
compass(DG_nds{ll}.U{t}(cc), DG_nds{ll}.V{t}(cc), 'r');
polar(tt, DG_nds{ll}.rho{t}(cc, :), 'b');
polar_theta_off(h)
for i = 2:9
    subplot(xx, yy, idx(i)); plot_raster(squeeze(raster_dg_nds{ll}{cc}(1, t, i-1, :)), 0, 8);
    if idx(i) == 7
        ylabel('trials')
        xlabel('time (s)')
    end
end
subplot(3, 3, 5)
tt = linspace(0, 2*pi, 100);
polar(tt, ones(1, 100), 'k');
polar(tt, 2*ones(1, 100), 'k');
polar(tt, 3*ones(1, 100), 'k');

%% fig 1D
ll = 5; t = 5;
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

%% fig 1G
clear X


xlabel('1st Principal Component')
ylabel('3rd Principal Component')
title('NDF 0')
close(5)
%% fig 1E-1F
figure
v = 4*datadg{2}.stimulus.params.SPATIAL_PERIOD./datadg{2}.stimulus.params.TEMPORAL_PERIOD;
subplot(1, 2, 1)
semilogx(v, exciseColumn(MAG_all_norm_dg{5}(:, idx_sub{1})), 'b')
xlabel('micron/second')
ylabel('Normalized Response')
xlim([v(end) v(1)])
subplot(1, 2, 2)
semilogx(v, exciseColumn(MAG_all_norm_dg{5}(:, idx_sub{2})), 'r')
xlabel('micron/second')
ylabel('Normalized Response')
xlim([v(end) v(1)])

%% fig 1G inset
t = 4;
figure
subplot(1, 2, 1)
compass(DG{5}.U{t}(idx_sub{1}), DG{5}.V{t}(idx_sub{1}), 'b')
subplot(1, 2, 2)
compass(DG{5}.U{t}(idx_sub{2}), DG{5}.V{t}(idx_sub{2}), 'r')

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
close
%%
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
close
%% dsmap
% see 'dsmap_150703.m' to find out how 'Tform' and 'center_corrected' is generated. 
load('DS150703-1.mat', 'Tform', 'center_corrected')

%% fig 1H
stixel_size_ = 10;
corner_i = [4 126 195 264 386 455 4];
corner_position = datadg{5}.ei.position(corner_i, :);
center_ds = tforminv(Tform, center_corrected * stixel_size_);
radius = 75;

% on-off DSGC
figure
for ct = 1:4
    subplot(2, 2, ct)
    plot(corner_position(:, 1), corner_position(:, 2))
    hold on
    for i = 1:length(idx_dir{ct})
        [x, y] = circle(center_ds(idx_dir{ct}(i), 1), center_ds(idx_dir{ct}(i), 2), radius);
        plot(x, y, 'k')
        hold on
    end
    xlim([-700 700])
    ylim([-500 500])
    axis off
%     title(oo_ds{ct})
end

% on DSGC
figure
for ct = 1:3
    subplot(2, 2, ct)
    plot(corner_position(:, 1), corner_position(:, 2))
    hold on
    for i = 1:length(idx_dir_on{ct})
        [x, y] = circle(center_ds(idx_dir_on{ct}(i), 1), center_ds(idx_dir_on{ct}(i), 2), radius);
        plot(x, y, 'k')
        hold on
    end
    xlim([-700 700])
    ylim([-500 500])
    axis off
%     title(on_ds{ct})
end



%% figure 1I 

datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-03-04-0/data003-sorted/data003-sorted', opt);
datarun.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-03-04-0/stimuli/s03.mat';
datarun = load_stim_matlab(datarun, 'user_defined_trigger_interval', 10);
datarun = load_ei(datarun, 'all', 'array_id', 1551);

pos = datarun.ei.position;
im_hb9 = imread('Hb9_.jpg');
figure
imshow(im_hb9);
array_location_image = ginput; % click 4 corners of array in order 1-->4 (see below)
elec_corner = [195 126 4 455];
array_location_ei = pos(elec_corner,:);
Tform = maketform('projective', array_location_image, array_location_ei);
test = tformfwd(Tform, array_location_image)-array_location_ei; % should be equal or close to zeros
close
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % get coordinates of GFP cell bodies (skip if just want to plot figure)
% soma_location_image = ginput_label('r');
% soma_location_ei = tformfwd(Tform, soma_location_image);
% % get array location in display coordinates
% 
% %                 EI                               DISPLAY
% %
% %               195 (1)                         386(5)  264(6)
% %                 / \                               ______
% %               /     \                            /      \
% %   264 (6)    |       |    126 (2)               /        \
% %   386 (5)    |       |    4   (3)       455(4)  \        / 195(1)
% %               \     /                            \      /
% %                 \ /                               ------
% %                455 (4)                          4(3)  126(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% see 'DS160304_classification.m' to find out id_dir_dg
load('DS160304.mat', 'soma_location_ei', 'id_dir_dg')  
id_dir = id_dir_dg;
ct = 1;
mode = 'neg';

for cc = 1:length(id_dir{ct})
    id = id_dir{ct}(cc);
    idx = get_cell_indices(datarun, id);
    ei = datarun.ei.eis{idx};
    com = ei_com_xy(ei, pos, 30*3, mode);
    com_oo{ct}(cc, :) = com;
    com_oo_image{ct}(cc, :) = tforminv(Tform, com);
    dis = sqrt(sum((soma_location_ei - repmat(com, size(soma_location_ei, 1), 1)).^2, 2));
    [disMin_oo{ct}(cc), matchI{ct}(cc)] = min(dis);
end

figure
imshow(im_hb9);
hold on
plot(com_oo_image{1}(:, 1), com_oo_image{1}(:, 2), 'o', 'color', 'r', 'MarkerSize', 7)
