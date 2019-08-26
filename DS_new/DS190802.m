cd /Users/xyao/matlab/code-private/DS_new
path = '/Users/xyao/Field-lab/data/2019-08-02-0/';

datapath = strcat(path, 'data004/data004');
stimpath = strcat(path, 'stimuli/s04.txt');
[ds_id, id_dir, id_dir_on, idx_dir, idx_dir_on] = DSClassification(datapath, stimpath, 'params_idx', [2 3], 'delta_p', 3, 'manual', false);

datafs = load_data(strcat(path, 'data006/data006'), opt);
datafs.names.stimulus_path = strcat(path, 'stimuli/s06.mat');
datafs = load_stim_mfs(datafs);

%% parameters 
field_width = 17; field_height = 17; % number of stixel on each dimension
subregion = 0; % if 1: subdivide the stimulus field into 4 regions, show 4 spatially correlated flash squares 
cell_ids = datafs.cell_ids;
bin_size = 0.02; % for PSTH plot


%% get rfs
fs_raster = get_fs_raster(datafs, cell_ids);
fs_spike = fs_get_spike(fs_raster);
[rf_all, rf_std] = get_fs_rf(fs_spike, field_width, field_height,subregion);

%% get histgram
trigger = datafs.stimulus.triggers';
list = datafs.stimulus.trial_list;
for i = 1:max(list)
    index(i, :) = find(list == i);
end
raster_onoff = cell(length(cell_ids), 1);
for i = 1:length(cell_ids)
    idx = get_cell_indices(datafs, cell_ids(i));
    raster = get_raster(datafs.spikes{idx}, trigger, 'plot', false);
    raster_onoff{i} = raster(index);
end

%
xx = [bin_size/2:bin_size:2-bin_size/2];
XX = [bin_size/2:bin_size:1-bin_size/2];
histg_onoff = cell(length(ds_id), 1);
for cc = 1:length(raster_onoff)
    for p = 1:field_width * field_height
        raster_all_onoff{cc}{p} = sort(cell2mat(raster_onoff{cc}(p,:)'));
        histg_onoff{cc}{p} = hist(raster_all_onoff{cc}{p}, xx);
    end
end

%% plot histgram

fig_n = 100;
for cc = 322:322%length(cell_ids)
    plot_mfs_psth(field_width, field_height, bin_size, histg_onoff, cc, fig_n)
    print_close(fig_n, [14 14], [num2str(cell_ids(cc)) '_psth'])
end

%% plot RF

for cc = 1:length(cell_ids)
    figure(fig_n)
    set(gcf, 'Position', [1 1 400 200])
    rf = rf_all{cc};
    subplot(1,2,1)
    imagesc(rf(:,:,1))
    colormap gray
    axis image
    axis off
    title('first second')

    subplot(1,2,2)
    imagesc(rf(:,:,2))
    colormap gray
    axis image
    axis off
    title('last second')
    
    print_close(fig_n, [5 3], [num2str(cell_ids(cc)) '_RF'])
end
