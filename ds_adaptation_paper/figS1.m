opt = struct('load_params', 1,'load_neurons', 1);
datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-03-04-0/data003-sorted/data003-sorted', opt);
datarun.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-03-04-0/stimuli/s03.mat';
datarun = load_stim_matlab(datarun, 'user_defined_trigger_interval', 10);
datarun = load_ei(datarun, 'all', 'array_id', 1551);

pos = datarun.ei.position;
im_hb9 = imread('Hb9_.jpg');
mode = 'neg';


%% map fluorescent image and ei

figure
imshow(im_hb9);
array_location_image = ginput; % click 4 corners of array in order 1-->4 (see below)
elec_corner = [195 126 4 455];
array_location_ei = pos(elec_corner,:);
Tform = maketform('projective', array_location_image, array_location_ei);
test = tformfwd(Tform, array_location_image)-array_location_ei % should be equal or close to zeros
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
load('DS160304.mat', 'soma_location_ei', 'id_dir_dg', 'id_dir_on_dg')  
id_dir = id_dir_dg;
id_dir_on = id_dir_on_dg;
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

for cc = 1:length(id_dir_on{ct})
    id = id_dir_on{ct}(cc);
    idx = get_cell_indices(datarun, id);
    ei = datarun.ei.eis{idx};
    com = ei_com_xy(ei, pos, 30*3, mode);
    com_on{ct}(cc, :) = com;
    com_on_image{ct}(cc, :) = tforminv(Tform, com);
    dis = sqrt(sum((soma_location_ei - repmat(com, size(soma_location_ei, 1), 1)).^2, 2));
    [disMin_on{ct}(cc), matchI{ct}(cc)] = min(dis);
end


%% figS1A
figure
plot_electrodes(pos, 'scale', 0.05, 'alpha', 0);
axis off

cell_id = [1922 5719 7757];
for i = 1:3
    ds_idx = find(id_dir{1} == cell_id(i));
    cell_idx = get_cell_indices(datarun, cell_id(i));
    ei = datarun.ei.eis{cell_idx};
    figure
    plot_ei_(ei, pos, 0, 'scale', 3.5, 'Alpha', 0);
    hold on
    plot(com_oo{1}(ds_idx, 1), com_oo{1}(ds_idx, 2), 'ro')
    plot(pos(455, 1), pos(455, 2), 'bo')
    axis off
end
%% figS1B
figure
imshow(im_hb9);
hold on
plot(com_oo_image{1}(:, 1), com_oo_image{1}(:, 2), 'o', 'color', 'r', 'MarkerSize', 7)
hold on
plot(com_on_image{1}(:, 1), com_on_image{1}(:, 2), 'o', 'color', 'y', 'MarkerSize', 7)


%% figS1C-D
load('DS160115.mat')
dis_oo{1} = disMin_oo;
dis_on{1} = disMin_on;
dis_oo_permute{1} = disMin_oo_permute;
dis_on_permute{1} = disMin_on_permute;
load('DS160130.mat')
dis_oo{2} = disMin_oo;
dis_on{2} = disMin_on;
dis_oo_permute{2} = disMin_oo_permute;
dis_on_permute{2} = disMin_on_permute;
load('DS160304.mat')
dis_oo{3} = disMin_oo;
dis_on{3} = disMin_on;
dis_oo_permute{3} = disMin_oo_permute;
dis_on_permute{3} = disMin_on_permute;


for i = 1:length(dis_oo{1})
    dis_all_oo{i} = [dis_oo{1}{i} dis_oo{2}{i} dis_oo{3}{i}];
    dis_all_oo_permute{i} = [dis_oo_permute{1}{i} dis_oo_permute{2}{i} dis_oo_permute{3}{i}];
end
for i = 1:length(dis_on{1})
    dis_all_on{i} = [dis_on{1}{i} dis_on{2}{i} dis_on{3}{i}];
    dis_all_on_permute{i} = [dis_on_permute{1}{i} dis_on_permute{2}{i} dis_on_permute{3}{i}];
end

XX = 0:5:120;
figure
subplot(3, 1, 1)
a = hist(dis_all_oo{1}, XX);
bar(XX, a, 1)
xlabel('um')
ylabel('on-off superior')
xlim([0 150])
subplot(3, 1, 2)
a = hist(dis_all_on{1}, XX);
bar(XX, a, 1)
xlabel('um')
ylabel('on superior')
xlim([0 150])
subplot(3, 1, 3)
a = hist(dis_all_oo{2}, XX);
bar(XX, a, 1)
xlabel('um')
ylabel('on-off anterior')
xlim([0 150])

XX = 0:5:120;
figure
subplot(3, 1, 1)
a = hist(dis_all_oo_permute{1}, XX);
bar(XX, a, 1)
xlabel('um')
ylabel('on-off superior')
xlim([0 150])
subplot(3, 1, 2)
a = hist(dis_all_on_permute{1}, XX);
bar(XX, a, 1)
xlabel('um')
ylabel('on superior')
xlim([0 150])
subplot(3, 1, 3)
a = hist(dis_all_oo_permute{2}, XX);
bar(XX, a, 1)
xlabel('um')
ylabel('on-off anterior')
xlim([0 150])

%% figS1E

load('DS150703-1.mat', 'Tform', 'center_corrected', 'idx_dir', 'idx_dir_on')
stixel_size_ = 10;
corner_i = [4 126 195 264 386 455 4];
corner_position = datarun.ei.position(corner_i, :);
center_ds = tforminv(Tform, center_corrected * stixel_size_);
radius = 75;

for dir = 1:4
    com = center_ds(idx_dir{dir}, :);
    nnnd_oo{dir} = get_nnnds_xy(com, repmat(radius,1,2));
end
for dir = 1:3
    com = center_ds(idx_dir_on{dir}, :);
    nnnd_on{dir} = get_nnnds_xy(com, repmat(radius,1,2));
end

nnnd_all = [cell2mat(nnnd_oo'); cell2mat(nnnd_on')];
figure
a = hist(nnnd_all, [0.05:0.1:3.45]);
bar([0.05:0.1:3.45], a, 1)
xlabel('nnnd')
ylabel('cell number')