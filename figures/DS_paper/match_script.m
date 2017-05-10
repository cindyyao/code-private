opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);
datarun{1} = load_data('/Volumes/lab/analysis/2016-01-15-0/data012-sorted/data012-sorted', opt);
datarun{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-01-15-0/stimuli/s12.mat';
datarun{1} = load_stim_matlab(datarun{1}, 'user_defined_trigger_interval', 10);
datarun{2} = load_data('/Volumes/lab/analysis/2016-01-30-0/data013-map/data013-map', opt);
datarun{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-01-30-0/stimuli/s13.mat';
datarun{2} = load_stim_matlab(datarun{2}, 'user_defined_trigger_interval', 10);
datarun{3} = load_data('/Volumes/lab/analysis/2016-03-04-0/data003/data003', opt);
datarun{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-03-04-0/stimuli/s03.mat';
datarun{3} = load_stim_matlab(datarun{3}, 'user_defined_trigger_interval', 10);

pos = datarun{1}.ei.position;
im_hb9 = imread('/Users/xyao/Documents/data/image/2016-03-04/Hb9_.jpg');
% im_hb9 = imread('/Users/xyao/Documents/data/image/2016-01-15/Hb9.jpg');
% im_hb9 = imread('/Users/xyao/Documents/data/image/2016-01-30/Hb9.jpg');
DB = 3;
mode = 'neg';


%% map fluorescent image and ei

figure
imshow(im_hb9);
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
if DB == 1
    load('DS160115.mat')
elseif DB == 2
    load('DS160130.mat')
else
    load('DS160304.mat')
end

for ct = 1:length(id_dir)
    for cc = 1:length(id_dir{ct})
        id = id_dir{ct}(cc);
        idx = get_cell_indices(datarun{DB}, id);
        ei = datarun{DB}.ei.eis{idx};
        com = ei_com_xy(ei, pos, 30*3, mode);
        com_oo{ct}(cc, :) = com;
        com_oo_image{ct}(cc, :) = tforminv(Tform, com);
        dis = sqrt(sum((soma_location_ei - repmat(com, size(soma_location_ei, 1), 1)).^2, 2));
        [disMin_oo{ct}(cc), matchI{ct}(cc)] = min(dis);
    end
end

for ct = 1:length(id_dir_on)
    for cc = 1:length(id_dir_on{ct})
        id = id_dir_on{ct}(cc);
        idx = get_cell_indices(datarun{DB}, id);
        ei = datarun{DB}.ei.eis{idx};
        com = ei_com_xy(ei, pos, 30*3, mode);
        com_on{ct}(cc, :) = com;
        com_on_image{ct}(cc, :) = tforminv(Tform, com);
        dis = sqrt(sum((soma_location_ei - repmat(com, size(soma_location_ei, 1), 1)).^2, 2));
        [disMin_on{ct}(cc), matchI{ct}(cc)] = min(dis);
    end
end

%% plot com of ei on microscopy image
figure
imshow(im_hb9);
hold on
plot(com_oo_image{1}(:, 1), com_oo_image{1}(:, 2), 'o', 'color', 'r', 'MarkerSize', 10)
% plot(com_oo_image{2}(:, 1), com_oo_image{2}(:, 2), 'o', 'color', 'c')
plot(com_on_image{1}(:, 1), com_on_image{1}(:, 2), 'o', 'color', 'y', 'MarkerSize', 10)

%% shortest distance distribution
XX = 0:5:120;
figure
subplot(3, 1, 1)
a = hist(disMin_oo{1}, XX);
bar(XX, a, 1)
xlabel('um')
ylabel('on-off superior')
xlim([0 150])
subplot(3, 1, 2)
a = hist(disMin_on{1}, XX);
bar(XX, a, 1)
xlabel('um')
ylabel('on superior')
xlim([0 150])
subplot(3, 1, 3)
a = hist(disMin_oo{2}, XX);
bar(XX, a, 1)
xlabel('um')
ylabel('on-off anterior')
xlim([0 150])

%%
XX = 0:5:120;
temp = disMin_oo{1};
temp(4) = [];
figure
subplot(3, 1, 1)
a = hist(temp, XX);
bar(XX, a, 1)
xlabel('um')
ylabel('on-off superior')
xlim([0 150])

subplot(3, 1, 2)
a = hist(disMin_on{1}, XX);
bar(XX, a, 1)
xlabel('um')
ylabel('on superior')
xlim([0 150])

subplot(3, 1, 3)
a = hist(disMin_oo{2}, XX);
bar(XX, a, 1)
xlabel('um')
ylabel('on-off anterior')
xlim([0 150])

%% combine 3 datasets
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

%% plot example ei
figure
plot_electrodes(pos, 'scale', 0.05, 'alpha', 0);
axis off

cell_id = [1922 5719 7757];
for i = 1:3
    ds_idx = find(id_dir{1} == cell_id(i));
    cell_idx = get_cell_indices(datarun{3}, cell_id(i));
    ei = datarun{3}.ei.eis{cell_idx};
    figure
    plot_ei_(ei, pos, 0, 'scale', 3.5, 'Alpha', 0);
    hold on
    plot(com_oo{1}(ds_idx, 1), com_oo{1}(ds_idx, 2), 'ro')
    plot(pos(455, 1), pos(455, 2), 'bo')
    axis off
end