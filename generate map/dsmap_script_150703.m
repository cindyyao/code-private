dg_path = '/Volumes/lab/analysis/2015-07-03-0/data012-map/data012-map';
dg_stimulus_path = '/Volumes/lab/analysis/2015-07-03-0/stimuli/s12.mat';
wn_path = '/Volumes/lab/analysis/2015-07-03-0/data016-map/data016-map';
wn_nonmap_path = '/Volumes/lab/analysis/2015-07-03-0/data016/data016';
cd /Users/xyao/matlab/code-private/DS_new/

%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);
datadg = load_data(dg_path, opt);
datadg.names.stimulus_path = dg_stimulus_path;
datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);

datawn = load_data(wn_path, opt);
datawn = load_sta(datawn);
datawn_nonmap = load_data(wn_nonmap_path, opt);
datawn_nonmap = load_sta(datawn_nonmap);
datawn_nonmap = get_rf_coms(datawn_nonmap, 'all');


%% get array location in display coordinates

% load stimulus picture taken by camera
im_s = imread('WN.jpg');

% get stimulus frame in display coordinates
stixel_size = 30;
display_width = 800; display_height = 600;
x_start = 100; x_end = 700; y_start = 0; y_end = 600;
movie_path = '/Volumes/lab/acquisition/movie-xml/BW-30-6-0.48-11111-20x20-60.35.xml';
mov = get_movie(movie_path, datawn_nonmap.triggers, 1);
mov_frame = matrix_scaled_up(squeeze(mov(:,:,1)), stixel_size);

% select control points
clear movingPoints fixedPoints
cpselect(im_s, mov_frame)

%% register two images
tform = fitgeotrans(movingPoints, fixedPoints, 'projective');
registered = imwarp(im_s, tform,'OutputView',imref2d(size(mov_frame)));
figure 
imshow(registered);
figure
imshowpair(mov_frame,registered,'blend');

% load array image taken by camera
im_array = imread('array_.jpg');

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
array_location_ei = datawn_nonmap.ei.position(elec_corner,:);
Tform = maketform('projective', array_location_ei, array_location_display);

test = tformfwd(Tform, array_location_ei)-array_location_display % should be equal or close to zeros

% get cell location from ei
on1 = get_cell_ids(datawn_nonmap, 'ON type 1');
off1 = get_cell_ids(datawn_nonmap, 'OFF type 1');
off2 = get_cell_ids(datawn_nonmap, 'OFF type 2');
off3 = get_cell_ids(datawn_nonmap, 'OFF type 3');
off4 = get_cell_ids(datawn_nonmap, 'OFF type 4');

robust_id = [on1 off1 off2 off3 off4];
robust_idx = get_cell_indices(datawn_nonmap, robust_id);

frame = datawn_nonmap.ei.nrPoints + datawn_nonmap.ei.nlPoints + 1;
elec = size(datawn_nonmap.ei.position, 1);
distance = 1;
center_ei = zeros(length(robust_id),2);
stixel_size_ = datawn_nonmap.stimulus.stixel_width;
for i = 1:length(robust_id)
    distance_temp = distance;
    full_elec_n = sum(1:distance_temp)*6+1; % Assume all arrays have hexagonal-arranged electrodes
    ei = datawn_nonmap.ei.eis{robust_idx(i)};
    ei = ei';
    [~,I] = max(abs(ei(:)));
    elec_n = ceil(I/frame);
    elecs_n = get_ei_neighbors(elec_n, elec, distance_temp);
    while(length(elecs_n) < full_elec_n)
        distance_temp = distance_temp - 1;
        elecs_n = get_ei_neighbors(elec_n, elec, distance_temp);
        full_elec_n = sum(1:distance_temp)*6+1;
    end
    points = datawn_nonmap.ei.position(elecs_n,:);
    mass = max(abs(ei(:,elecs_n)));
    center_ei(i,:) = centroid(points, mass);
end

center_ei_display = tformfwd(Tform, center_ei);
center_sta = cell2mat(datawn_nonmap.stas.rf_coms(robust_idx)) * stixel_size_;
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

% 
figure
scatter3(center_ei_display(:,1),center_ei_display(:,2),center_offset(:,2))

% fit plane
dataX(:,1:2) = center_ei_display;
dataX(:,3) = center_offset(:,1);
dataY(:,1:2) = center_ei_display;
dataY(:,3) = center_offset(:,2);

[normal_x, meanX, pcaExplainedX, sseX] = fit_plane(dataX);
[normal_y, meanY, pcaExplainedY, sseY] = fit_plane(dataY);


%% get DS cell id
% [NumSpikesCell, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0);
% ds_struct = dscellanalysis(NumSpikesCell, StimComb, datadg);
% 
% params_idx = input('indices of parameters for classification: (eg. [1 2])\n'); 
% [ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);
load('DS150703-1.mat')
ds_idx = get_cell_indices(datadg, ds_id);

%% get DS cell location

% cell_location = get_cell_location(datawn, ds_id, 'ei_mode', true);
% [ds_id_wn, cell_location] = Test_Cell_Location(ds_id_test, cell_location);

% get DS cell ei location (display coordinate)
center_ei_ds = zeros(length(ds_id),2);
for i = 1:length(ds_id)
    distance_temp = distance;
    full_elec_n = sum(1:distance_temp)*6+1; % Assume all arrays have hexagonal-arranged electrodes
    ei = datadg.ei.eis{ds_idx(i)};
    ei = ei';
    [~,I] = max(abs(ei(:)));
    elec_n = ceil(I/frame);
    elecs_n = get_ei_neighbors(elec_n, elec, distance_temp);
    while(length(elecs_n) < full_elec_n)
        distance_temp = distance_temp - 1;
        elecs_n = get_ei_neighbors(elec_n, elec, distance_temp);
        full_elec_n = sum(1:distance_temp)*6+1;
    end
    points = datadg.ei.position(elecs_n,:);
    mass = max(abs(ei(:,elecs_n)));
    center_ei_ds(i,:) = centroid(points, mass);
end
center_ei_ds = tformfwd(Tform, center_ei_ds);

offset_x_ds = meanX*normal_x*ones(length(ds_id),1) - center_ei_ds*normal_x(1:2);
offset_y_ds = meanX*normal_y*ones(length(ds_id),1) - center_ei_ds*normal_y(1:2);
offset_ds = [offset_x_ds offset_y_ds];
center_sta_ds_est = (center_ei_ds + offset_ds)/stixel_size_;

%%
% on-off dsgc

oo_ds = {'posterior', 'inferior', 'anterior', 'superior'};
x = [3 3 4 4]; y = [5 5 5 5];
for i = 4:4
    figure
    for cc = 1:length(id_dir{i})
        subplot(x(i), y(i), cc)
        plot_rf(datawn, id_dir{i}(cc), 'fit', false, 'title', false);
        hold on
        plot(center_sta_ds_est(idx_dir{i}(cc), 1), center_sta_ds_est(idx_dir{i}(cc), 2), 'ro')
        plot(center_ei_ds(idx_dir{i}(cc), 1)/stixel_size_, center_ei_ds(idx_dir{i}(cc), 2)/stixel_size_, 'go')
    end
    subplot(x(i), y(i), 1)
    title(oo_ds{i})
    legend('w correction', 'w/o correction')
end

% on dsgc
on_ds = {'inferior', 'anterior', 'superior'};
x = [2 3 2]; y = [2 3 3];
for i = 1:3
    figure
    for cc = 1:length(id_dir_on{i})
        if ismember(id_dir_on{i}(cc), datawn.cell_ids)
            subplot(x(i), y(i), cc)
            plot_rf(datawn, id_dir_on{i}(cc), 'fit', false, 'title', false);
            hold on
            plot(center_sta_ds_est(idx_dir_on{i}(cc), 1), center_sta_ds_est(idx_dir_on{i}(cc), 2), 'ro')
            plot(center_ei_ds(idx_dir_on{i}(cc), 1)/stixel_size_, center_ei_ds(idx_dir_on{i}(cc), 2)/stixel_size_, 'go')
        end
    end
    subplot(x(i), y(i), 2)
    title(on_ds{i})
    legend('w correction', 'w/o correction')
end

%% manually correct center location if needed
center_corrected = correct_ei_center_(datawn, ds_id, center_sta_ds_est*stixel_size_, stixel_size_) / stixel_size_;

center_corrected = center_sta_ds_est;
%% generate circular masks
n = size(cell_location,1);
radius = 3;
masks = make_circular_masks(cell_location, ones(n,1)*radius, 60, 60, 10);
dlmwrite('maptest.txt', masks{1},  'delimiter', '\t', 'newline', 'pc')

%% plot mosaic

corner_i = [4 126 195 264 386 455 4];
corner_position = datadg.ei.position(corner_i, :);
center_ds = tforminv(Tform, center_corrected * stixel_size_);
radius = 75; % micron



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
        [x, y] = circle(center_ds(idx_dir_on{ct}(i), 1), center_ds(idx_dir{ct}(i), 2), radius);
        plot(x, y, 'k')
        hold on
    end
    xlim([-700 700])
    ylim([-500 500])
    axis off
%     title(on_ds{ct})
end



