wn_path = '/Volumes/lab/analysis/2015-07-03-0/data016-map/data016-map';
cd /Users/xyao/matlab/code-private/DS_new/

%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);
datawn = load_data(wn_path, opt);
% datawn = load_sta(datawn);


%% get array location in display coordinates

% load stimulus picture taken by camera
im_s = imread('WN.jpg');

% get stimulus frame in display coordinates
stixel_size = 30;
display_width = 800; display_height = 600;
x_start = 100; x_end = 700; y_start = 0; y_end = 600;
movie_path = '/Volumes/lab/acquisition/movie-xml/BW-30-6-0.48-11111-20x20-60.35.xml';
mov = get_movie(movie_path, datawn.triggers, 1);
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
array_location_ei = datawn.ei.position(elec_corner,:);
Tform = maketform('projective', array_location_ei, array_location_display);

test = tformfwd(Tform, array_location_ei)-array_location_display % should be equal or close to zeros

