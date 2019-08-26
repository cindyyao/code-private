% wn_map_path = '/Analysis/xyao/2015-07-03-0/data016-map/data016-map';
wn_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data016/data016';
cd /Users/xyao/matlab/code-private/DS_new/
%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);
% datawn_map = load_data(wn_map_path, opt);
% datawn_map = load_sta(datawn_map);
% datawn_map = get_rf_coms(datawn_map, 'all');

datawn = load_data(wn_path, opt);
datawn = load_sta(datawn);
datawn = get_rf_coms(datawn, 'all');


%% map camera and display coordinates

im_s = imread('WN.jpg'); % load stimulus picture taken by camera
im_array = imread('array_.jpg'); % load array image taken by camera
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
imshow(registered_array);

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
array_location_ei = datawn.ei.position(elec_corner,:);
Tform = maketform('projective', array_location_ei, array_location_display);
test = tformfwd(Tform, array_location_ei)-array_location_display % should be equal or close to zeros

%% get cell location from ei
on1 = get_cell_ids(datawn, 'ON type 1');
off1 = get_cell_ids(datawn, 'OFF type 1');
off2 = get_cell_ids(datawn, 'OFF type 2');
off3 = get_cell_ids(datawn, 'OFF type 3');
off4 = get_cell_ids(datawn, 'OFF type 4');
robust_id = [on1 off1 off2 off3 off4];
robust_idx = get_cell_indices(datawn, robust_id);

distance = 1; % range used to calculate EI com
stixel_size_ = datawn.stimulus.stixel_width;
center_ei = get_ei_com(datawn, robust_id, distance);
center_ei_display = tformfwd(Tform, center_ei);
center_sta = cell2mat(datawn.stas.rf_coms(robust_idx)) * stixel_size_;
center_offset = center_sta - center_ei_display;
[center_offset_pol(:,1), center_offset_pol(:,2)] = cart2pol(center_offset(:,1),center_offset(:,2));

% plot ei-sta offset
figure
quiver(center_ei_display(:,1),center_ei_display(:,2),center_offset(:,1),center_offset(:,2),0)
id_text = cellfun(@num2str, num2cell(robust_id), 'Uni', false);
text(center_ei_display(:,1), center_ei_display(:,2), id_text);

% exclude outliers
figure
hist(center_offset_pol(:,2), 30)
xlabel('offset')
[x,~] = ginput;
outlier_i = find(center_offset_pol(:,2) > x);
center_offset_pol(outlier_i,:) = [];
center_offset(outlier_i,:) = [];
robust_id(outlier_i) = [];
robust_idx(outlier_i) = [];
center_ei_display(outlier_i,:) = [];
center_sta(outlier_i,:) = [];

% replot offset map
figure
quiver(center_ei_display(:,1),center_ei_display(:,2),center_offset(:,1),center_offset(:,2),0)
id_text = cellfun(@num2str, num2cell(robust_id), 'Uni', false);
text(center_ei_display(:,1), center_ei_display(:,2), id_text);

% fit plane
dataX(:,1:2) = center_ei_display;
dataX(:,3) = center_offset(:,1);
dataY(:,1:2) = center_ei_display;
dataY(:,3) = center_offset(:,2);

[normal_x, meanX, pcaExplainedX, sseX] = fit_plane(dataX);
[normal_y, meanY, pcaExplainedY, sseY] = fit_plane(dataY);

%% get DS cell id
% dg_path = '/Analysis/xyao/2015-07-03-0/data012-map/data012-map';
% dg_stimulus_path = '/Analysis/xyao/2015-07-03-0/stimuli/s12.mat';
dg_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data012/data012';
dg_stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s12.mat';
datadg = load_data(dg_path, opt);
datadg.names.stimulus_path = dg_stimulus_path;
datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);

[NumSpikesCell, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0);
ds_struct = dscellanalysis(NumSpikesCell, StimComb, datadg);
params_idx = input('indices of parameters for classification: (eg. [1 2])\n'); 
[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);
ds_idx = get_cell_indices(datadg, ds_id);

%% classify on vs on-off DSGC
[NumSpikesCell, StimComb] = get_spikescellstim(datadg,ds_id,0);
DG = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
MAG_all_norm_dg = normalize_MAG(DG);
mag_pca = MAG_all_norm_dg;
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

figure
plot(scores(idx_sub{1}, pc1), scores(idx_sub{1}, pc2), 'ro', scores(idx_sub{2}, pc1), scores(idx_sub{2}, pc2), 'bo')

figure
v = 4*datadg.stimulus.params.SPATIAL_PERIOD./datadg.stimulus.params.TEMPORAL_PERIOD;
subplot(1, 2, 1)
semilogx(v, exciseColumn(MAG_all_norm_dg(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
xlim([v(end) v(1)])

subplot(1, 2, 2)
semilogx(v, exciseColumn(MAG_all_norm_dg(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
xlim([v(end) v(1)])

t = 4;
figure
compass(DG.U{t}(idx_sub{1}), DG.V{t}(idx_sub{1}), 'r')
hold on
compass(DG.U{t}(idx_sub{2}), DG.V{t}(idx_sub{2}), 'b')

%% classify DSGC into subtypes (directions)
t = 3;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 500 500])
compass(DG.U{t}(idx_sub{1}), DG.V{t}(idx_sub{1}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG.U{t}(idx_sub{1}), DG.V{t}(idx_sub{1}), x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = idx_sub{1}(I);
    id_dir{i} = ds_id(idx_dir{i});
end

%
t = 3;
h = figure;
dirn = 3;
set(h, 'Position', [1 1 500 500])
compass(DG.U{t}(idx_sub{2}), DG.V{t}(idx_sub{2}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG.U{t}(idx_sub{2}), DG.V{t}(idx_sub{2}), x, y);
    [~, I] = find(IN == 1);
    idx_dir_on{i} = idx_sub{2}(I);
    id_dir_on{i} = ds_id(idx_dir_on{i});
end



%% map two dg to wn
[cell_list_map, ~] = map_ei(datadg, datawn);
ds_list_map = cell_list_map(ds_idx);
success_idx = ~cellfun(@isempty, ds_list_map);
ds_id_dgmap = ds_id(success_idx);
ds_idx_dgmap = ds_idx(success_idx);
ds_id_wnmap = cell2mat(ds_list_map);

t = 4;
figure
compass(DG.U{t}(success_idx), DG.V{t}(success_idx), 'b')

for i = 1:4
    id_dir_map{i} = intersect(id_dir{i}, ds_id_dgmap);
    id_dir_nomap{i} = setdiff(id_dir{i}, id_dir_map{i});
end
for i = 1:3
    id_dir_on_map{i} = intersect(id_dir_on{i}, ds_id_dgmap);
    id_dir_on_nomap{i} = setdiff(id_dir_on{i}, id_dir_on_map{i});
end

n = 8;
for i = 1:4
    id_final_oo{i} = id_dir_map{i}(randsample(length(id_dir_map{i}), n));
end
id_final_oo = cell2mat(id_final_oo);

n = 1;
for i = 1:3
    id_final_on{i} = id_dir_on_nomap{i}(randsample(length(id_dir_on_nomap{i}), n));
end
id_final_on = cell2mat(id_final_on);
id_final_on_ = cell2mat(id_dir_on_map);


id_final = sort([id_final_oo id_final_on id_final_on_]);

idx_final = get_cell_indices(datadg, id_final);
id_final_wn = cell_list_map(idx_final);

%% get DS cell ei location (display coordinate)
center_ei_ds = get_ei_com(datadg, id_final, 1);
center_ei_ds = tformfwd(Tform, center_ei_ds);

offset_x_ds = meanX*normal_x*ones(length(id_final),1) - center_ei_ds*normal_x(1:2);
offset_y_ds = meanY*normal_y*ones(length(id_final),1) - center_ei_ds*normal_y(1:2);
offset_ds = [offset_x_ds offset_y_ds];
% center_sta_ds_est = (center_ei_ds + offset_ds)/stixel_size_;
center_sta_ds_est = center_ei_ds + offset_ds;

%% manually correct center location if needed
center_corrected = correct_ei_center(datawn, id_final_wn, center_sta_ds_est, stixel_size_);

%% generate circular masks
radius = 6;
display_width = 800;
display_height = 600;

n = size(center_corrected,1);
radius_all = ones(n,1)*radius;
center_corrected(:,1) = center_corrected(:,1)+abs(display_width-display_height)/2;
masks = make_circular_masks(center_corrected, radius_all, display_width, display_height);
% maploc = '/Volumes/lab/analysis/2015-07-03-0/stimuli/targeted_flash/';
maploc = '/Analysis/xyao/test/stimuli/targeted_flash/';
for i = 1:length(id_final)
    name = ['map-' num2str(i-1, '%04.0f') '.txt'];
    dlmwrite([maploc name], masks{i},  'delimiter', '\t', 'newline', 'pc')
end

%%
% center_sta_ds = cell2mat(datawn_map.stas.rf_coms(ds_idx)) * stixel_size_;
% 
% i = 4;
% figure
% imagesc(masks{i})
% hold on
% plot(center_sta_ds(i, 1)+100, center_sta_ds(i,2), 'ro')
% 
% 
% figure
% for i = 1:length(ds_id)
%     imagesc(masks{i})
%     hold on
% %     plot(center_sta_ds_est(i, 1), center_sta_ds_est(i,2), 'ro')
%     plot(367.05, 271.5, 'ro')
%     pause
% end