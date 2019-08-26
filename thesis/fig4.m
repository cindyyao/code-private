dg_path = '/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/data013-map/data013-map';
dg_stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-01-30-0/stimuli/s13.mat';
cd /Users/xyao/matlab/code-private/DS_new/

%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);
datadg = load_data(dg_path, opt);
datadg.names.stimulus_path = dg_stimulus_path;
datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);

%% get array location in display coordinates

% load stimulus picture taken by camera
im_s = imread('/Volumes/lab/Experiments/Array/Images/2016-01-30-0/WN.jpg');

% get stimulus frame in display coordinates
stixel_size = 30;
display_width = 800; display_height = 600;
x_start = 100; x_end = 700; y_start = 0; y_end = 600;
movie_path = '/Volumes/lab/acquisition/movie-xml/BW-30-6-0.48-11111-20x20-60.35.xml';
mov = get_movie(movie_path, 1, 1);
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
im_array = imread('/Volumes/lab/Experiments/Array/Images/2016-01-30-0/array.jpg');

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
array_location_ei = datadg.ei.position(elec_corner,:);
Tform = maketform('projective', array_location_ei, array_location_display);

test = tformfwd(Tform, array_location_ei)-array_location_display % should be equal or close to zeros


%% get DS cell id
% [NumSpikesCell, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0);
% ds_struct = dscellanalysis(NumSpikesCell, StimComb, datadg);
% 
% params_idx = input('indices of parameters for classification: (eg. [1 2])\n'); 
% [ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);
load('DS160130.mat')
% ds_id = intersect(datawn.cell_ids, intersect(datadg.cell_ids, ds_id));
% for i = 1:4
%     [id_dir{i}, idx_dir{i}] = intersect(ds_id, id_dir{i});
% end
% for i = 1:3
%     [id_dir_on{i}, idx_dir_on{i}] = intersect(ds_id, id_dir_on{i});
% end
ds_idx = get_cell_indices(datadg, ds_id);

%% get DS cell location

% cell_location = get_cell_location(datawn, ds_id, 'ei_mode', true);
% [ds_id_wn, cell_location] = Test_Cell_Location(ds_id_test, cell_location);

% get DS cell ei location (display coordinate)
frame = datadg.ei.nrPoints + datadg.ei.nlPoints + 1;
elec = size(datadg.ei.position, 1);
distance = 3;
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



%% plot mosaic

corner_i = [4 126 195 264 386 455 4];
corner_position = datadg.ei.position(corner_i, :);
center_ds = center_ei_ds;
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
        [x, y] = circle(center_ds(idx_dir_on{ct}(i), 1), center_ds(idx_dir_on{ct}(i), 2), radius);
        plot(x, y, 'k')
        hold on
    end
    xlim([-700 700])
    ylim([-500 500])
    axis off
%     title(on_ds{ct})
end



%% NNND
radius = 70;
for dir = 1:4
    com = center_ds(idx_dir{dir}, :);
    nnnd_oo{dir} = get_nnnds_xy(com, repmat(radius,1,2));
end
for dir = 1:3
    com = center_ds(idx_dir_on{dir}, :);
    nnnd_on{dir} = get_nnnds_xy(com, repmat(radius,1,2));
end

figure
for dir = 1:4
    subplot(2,4,dir)
    hist(nnnd_oo{dir}, [0.05:0.1:2.95])
    xlim([0 inf])
end
for dir = 1:3
    subplot(2,4,dir+4)
    hist(nnnd_on{dir}, [0.05:0.1:2.95])
    xlim([0 inf])
end

nnnd_all = [cell2mat(nnnd_oo'); cell2mat(nnnd_on')];
figure
a = hist(nnnd_all, [0.05:0.1:3.45]);
bar([0.05:0.1:3.45], a, 1)
xlabel('nnnd')
ylabel('cell number')