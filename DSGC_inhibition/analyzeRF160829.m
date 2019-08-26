function [rf_all] = analyzeRF160829

cd /Volumes/dusom_fieldlab/All_Staff/lab/Development/matlab/private/xyao/matlab/code/DS_new
path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2016-08-29-0/';
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

datafs{1} = load_data(strcat(path, 'data000-map/data000-map'), opt);
datafs{1}.names.stimulus_path = strcat(path, 'stimuli/s00.mat');
datafs{1} = load_stim_matlab(datafs{1});
datafs{1}.triggers = datafs{1}.triggers(2:end);

datafs{2} = load_data(strcat(path, 'data001-map/data001-map'), opt);
datafs{2}.names.stimulus_path = strcat(path, 'stimuli/s01.mat');
datafs{2} = load_stim_matlab(datafs{2});
datafs{2}.triggers = datafs{2}.triggers(2:end);

datafs{3} = load_data(strcat(path, 'data004-map/data004-map'), opt);
datafs{3}.names.stimulus_path = strcat(path, 'stimuli/s04.mat');
datafs{3} = load_stim_matlab(datafs{3});
datafs{3}.triggers = datafs{3}.triggers(2:end);

% classification code below. Load pre-saved data to keep consistancy

% datapath = strcat(path, 'data003-sorted/data003-sorted');
% stimpath = strcat(path, 'stimuli/s03.mat');
% [ds_id, id_dir, id_dir_on, idx_dir, idx_dir_on] = DSClassification(datapath, stimpath, 'params_idx', [2 5], 'delta_p', 4, 'manual', true);
load('DS160829.mat', 'ds_id', 'id_dir', 'idx_dir', 'id_dir_on', 'idx_dir_on', 'fs_idx')


%% get rfs
field_width = 20; field_height = 20;
subregion = 0;
for i = 1:3
    fs_raster{i} = get_fs_raster(datafs{i}, ds_id);
    for cc = 1:length(ds_id)
        if fs_idx(cc,i)
            fs_raster{i}{cc} = [];
        end
    end
    fs_spike{i} = fs_get_spike(fs_raster{i});
    [rf_all{i}, rf_std{i}] = get_fs_rf(fs_spike{i}, field_width, field_height,subregion);
end

end
