cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);


datadg = load_data('/Volumes/lab/analysis/2016-05-06-0/data007-sorted/data007-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-05-06-0/stimuli/s07.mat';
datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);

[NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [3 4]; % which parameters to use for classification

[ds_id, nonds_id, id_init] = classify_ds(datadg, ds_struct, params_idx);

datafs = load_data('/Volumes/lab/analysis/2016-05-06-0/data008-map/data008-map', opt);
datafs.names.stimulus_path = '/Volumes/lab/analysis/2016-05-06-0/stimuli/s08.mat';
datafs = load_stim_matlab(datafs);

datawn = load_data('/Volumes/lab/analysis/2016-05-06-0/data004-map/data004-map', opt);
datawn = load_sta(datawn);

load('/Users/xyao/matlab/data/DS160506.mat')
%% 
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg,ds_id,0,1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);

delta_p = 4; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(DG{i});
rep = datadg.stimulus.repetitions;

%% plot cell summary
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [4 5]);
for cc = 1:length(ds_id);
    id = ds_id(cc);
    plot_ds_raster(DG_cut, raster_dg_cut, cc, id, '', 1, 1, 1)
end

%% speed tuning curves
figure
v = datadg.stimulus.params.SPATIAL_PERIOD./datadg.stimulus.params.TEMPORAL_PERIOD*4;
semilogx(v, exciseColumn(MAG_all_norm_dg{1}), 'b')
xlabel('micron/second')
ylabel('Response')
% title(ll{2})
xlim([v(end) v(1)])

%% classification based on speed tunning
L = 1;
mag_pca = MAG_all_norm_dg{L}(2:end,:);
mag_pca = mag_pca';
[id_sub, idx_sub] = deal(cell(2, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(mag_pca);
pc1 = 2; pc2 = 3;
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
ylabel('3rd Principal Component')
title('NDF 0')

figure
plot(scores(idx_sub{1}, pc1), scores(idx_sub{1}, pc2), 'ro', scores(idx_sub{2}, pc1), scores(idx_sub{2}, pc2), 'bo')
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
title('NDF 0')

v = 4*datadg.stimulus.params.SPATIAL_PERIOD./datadg.stimulus.params.TEMPORAL_PERIOD;
% subplot(1, 2, 2)
figure
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
% title(ll{5})
xlim([v(end) v(1)])

t = 4;
figure
compass(DG{1}.U{t}(idx_sub{1}), DG{1}.V{t}(idx_sub{1}), 'r')
hold on
compass(DG{1}.U{t}(idx_sub{2}), DG{1}.V{t}(idx_sub{2}), 'b')

%% classify ds direction
d = 1;
t = 5;
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

d = 1;
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

%% plot individual cell summary
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [4 5]);
id = 860;
cc = get_cell_indices(datadg, id);
plot_ds_raster(DG_cut, raster_dg_cut, cc, id, '', 1, 1, 0)


%% Moving flashing square ds id
ds_id_fs = ds_id(~idx_fs);

for ct = 1:4
    [id_dir_fs{ct}, idx_dir_fs{ct}] = intersect(ds_id_fs, id_dir{ct});
end
for ct = 1:3
    [id_dir_on_fs{ct}, idx_dir_on_fs{ct}] = intersect(ds_id_fs, id_dir_on{ct});
end

%% wn ds id
ds_id_wn = ds_id(~idx_wn);
%% map camera and display coordinates
cd /Analysis/xyao/2016-05-06-0/images/
im_s = imread('WN.jpg'); % load stimulus picture taken by camera
im_array = imread('array.jpg'); % load array image taken by camera
cd /Users/xyao/matlab/code-private/DS_new/
stixel_size = 30; % frame shown in WN.jpg
movie_path = '/Volumes/lab/acquisition/movie-xml/BW-30-6-0.48-11111-20x20-60.35.xml';

mov = get_movie(movie_path, 0, 1);
mov_frame = matrix_scaled_up(squeeze(mov(:,:,1)), stixel_size);
clear movingPoints fixedPoints
cpselect(im_s, mov_frame) % select 4 control points

%% register two images
x_start = 105;
stixel_size_ = 15;
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

%              EI & Display                              
%
%               195 (1)                         
%                 / \                              
%               /     \                            
%   264 (6)    |       |    126 (2)               
%   386 (5)    |       |    4   (3)       
%               \     /                           
%                 \ /                              
%                455 (4)                          
array_location_display = ginput;

% get array location in ei coordinates
elec_corner = [195 126 4 455];
array_location_ei = datafs.ei.position(elec_corner,:);
Tform = maketform('projective', array_location_ei, array_location_display);
test = tformfwd(Tform, array_location_ei)-array_location_display % should be equal or close to zeros

distance = 1; % range used to calculate EI com
center_ei = get_ei_com(datafs, ds_id_fs, distance);
center_ei_display = (tformfwd(Tform, center_ei) - x_start)/stixel_size_;

%% plot rfs
cell_type = {'superior', 'anterior', 'inferior', 'posterior'};
field_width = 13; field_height = 13;
field_width_sta = 40; field_height_sta = 40;
subregion = 1;

fs_raster = get_fs_raster(datafs, ds_id_fs);
fs_spike = fs_get_spike(fs_raster);
rf_all = get_fs_rf(fs_spike, field_width, field_height,subregion);

for cc = 1:length(ds_id_fs)
    id = ds_id_fs(cc);
    rf = padarray(rf_all{cc},[7,7]);
    ei = center_ei_display(cc,:);
    ei = ei + 7;
    
    figure;
    set(gcf,'Position',[1,1,800,800])
    subplot(2,2,1)
    plot_rf(datawn,id,'foa',-1,'title',false,'ticks',true);
    title(num2str(id))
    
    subplot(2,2,2)
    imagesc(sum(rf,3))
    colormap gray
    hold on 
    plot(ei(1),ei(2),'ro')
    axis image
    
    subplot(2,2,3)
    imagesc(rf(:,:,1))
    colormap gray
    hold on
    plot(ei(1),ei(2),'ro')
    title('on')
    axis image
    
    subplot(2,2,4)
    imagesc(rf(:,:,2))
    colormap gray
    hold on
    plot(ei(1),ei(2),'ro')
    title('off')
    axis image
    print_close(1,[12 12],num2str(id))
end

max_spike = cell(4,1);
figure
for ct = 1:4
    subplot(2,2,ct)
    for cc = 1:length(id_dir_fs{ct})
        max_spike{ct}(cc) = max(rf_all{idx_dir_fs{ct}(cc)}(:));
    end
    hist(max_spike{ct})
    xlabel('spike #')
    ylabel('cell #')
    title(cell_type{ct})
end


%%
fs_spike_temp = cellfun(@(fs_spike) sum(fs_spike,2), fs_spike, 'UniformOutput', false);
fs_spike_temp = cellfun(@(fs_spike_temp) sum(fs_spike_temp,3), fs_spike_temp, 'UniformOutput', false);
[~, I] = cellfun(@(fs_spike_temp) sort(fs_spike_temp, 'descend'), fs_spike_temp, 'UniformOutput', false);trigger = datafs.triggers(2:2:end);
list = datafs.stimulus.trial_list;
repeat = datafs.stimulus.repetitions;
for i = 1:max(list)
    index(i, :) = find(list == i);
end
raster_onoff = cell(length(ds_id_fs), 1);
for i = 1:length(ds_id_fs)
    if ~isempty(intersect(datafs.cell_ids, ds_id_fs(i)))
        idx = get_cell_indices(datafs, ds_id_fs(i));
        raster = get_raster(datafs.spikes{idx}, trigger, 'plot', false);
        raster_onoff{i} = raster(index);
    end
end
max_i = 1;
for dir = 1:4
    figure
    for cc = 1:length(id_dir_fs{dir})
        raster = raster_onoff{idx_dir_fs{dir}(cc)};
        raster_c = raster(I{idx_dir_fs{dir}(cc)}(max_i), :);
        subplot(10, 2, cc)
        plot_raster(raster_c, 0, 2)
    end
end

for dir = 1:3
    figure
    for cc = 1:length(id_dir_on_fs{dir})
        raster = raster_onoff{idx_dir_on_fs{dir}(cc)};
        raster_c = raster(I{idx_dir_on_fs{dir}(cc)}(max_i), :);
        subplot(4, 1, cc)
        plot_raster(raster_c, 0, 2)
    end
end

%%
xx = [0.025:0.05:1.975];
for cc = 1:length(raster_onoff)
    for p = 1:size(raster_onoff{1}, 1)
        raster_all{cc}{p} = sort(cell2mat(raster_onoff{cc}(p,:)'));
        histg{cc}{p} = hist(raster_all{cc}{p}, xx);
    end
end
        
for dir = 1:1
    for cc = 1:1%length(id_dir_fs{dir})
        h = figure;
        set(h, 'Position', [1 1 1080 1080])
        for p = 1:169
            subplot(13, 13, p)
            plot(xx, histg{idx_dir_fs{dir}(cc)}{p})
            ylim([0 max(max(cell2mat(histg{idx_dir_fs{dir}(cc)}')))])
            axis off
        end
%         print_close(1, [14 14], num2str(id_dir_fs{dir}(cc)));
    end
end    

for dir = 1:3
    for cc = 1:length(id_dir_on_fs{dir})
        h = figure;
        set(h, 'Position', [1 1 1080 1080])
        for p = 1:169
            subplot(13, 13, p)
            plot(xx, histg{idx_dir_on_fs{dir}(cc)}{p})
            ylim([0 max(max(cell2mat(histg{idx_dir_on_fs{dir}(cc)}')))])
            axis off
        end
        print_close(1, [14 14], num2str(id_dir_on_fs{dir}(cc)));
    end
end     