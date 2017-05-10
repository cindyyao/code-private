cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);


datarun = load_data('/Volumes/lab/analysis/2016-10-10-0/data000/data000', opt);
datarun.names.stimulus_path = '/Volumes/lab/analysis/2016-10-10-0/stimuli/s00.mat';
datarun = load_stim_mfs(datarun);


%% plot rfs
field_width = 13; field_height = 13;
field_width_sta = 40; field_height_sta = 40;
subregion = 1;
stop = 0.5; %second
fs_raster = get_fs_raster(datarun, datarun.cell_ids);
fs_spike = fs_get_spike(fs_raster);
[rf_all, rf_std] = get_fs_rf(fs_spike, field_width, field_height,subregion);



for cc = 1:length(datarun.cell_ids)
    figure(1)
    set(gcf, 'Position', [1 1 1000 500])
    id = datarun.cell_ids(cc);
    rf = padarray(rf_all{cc},[7,7]);

    subplot(1,3,1)
    imagesc(sum(rf,3))
    colormap gray
    axis image

    subplot(1,3,2)
    imagesc(rf(:,:,1))
    colormap gray
    title('on')
    axis image

    subplot(1,3,3)
    imagesc(rf(:,:,2))
    colormap gray
    title('off')
    axis image
    pause
    close(1)
%     print_close(1,[12 12],num2str(id))
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
