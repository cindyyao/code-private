function [rf_all, histg_onoff, rf_area_clean_mean, rf_area_clean_ste] = analyzeRF161122

cd /Volumes/dusom_fieldlab/All_Staff/lab/Development/matlab/private/xyao/matlab/code/DS_new
path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2016-11-22-0/';
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

datafs{1} = load_data(strcat(path, 'data000-map/data000-map'), opt);
datafs{1}.names.stimulus_path = strcat(path, 'stimuli/s00.mat');
datafs{1} = load_stim_mfs(datafs{1});
datafs{2} = load_data(strcat(path, 'data001-map/data001-map'), opt);
datafs{2}.names.stimulus_path = strcat(path, 'stimuli/s01.mat');
datafs{2} = load_stim_mfs(datafs{2});
datafs{3} = load_data(strcat(path, 'data002-map/data002-map'), opt);
datafs{3}.names.stimulus_path = strcat(path, 'stimuli/s02.mat');
datafs{3} = load_stim_mfs(datafs{3});
datafs{4} = load_data(strcat(path, 'data003-map/data003-map'), opt);
datafs{4}.names.stimulus_path = strcat(path, 'stimuli/s03.mat');
datafs{4} = load_stim_mfs(datafs{4});
datafs{5} = load_data(strcat(path, 'data006-map/data006-map'), opt);
datafs{5}.names.stimulus_path = strcat(path, 'stimuli/s06.mat');
datafs{5} = load_stim_mfs(datafs{5});

% classification code below. Load pre-saved data to keep consistancy

% datapath = strcat(path, 'data007-sorted/data007-sorted');
% stimpath = strcat(path, 'stimuli/s07.txt');
% [ds_id, id_dir, id_dir_on, idx_dir, idx_dir_on] = DSClassification(datapath, stimpath, 'params_idx', [2 3], 'delta_p', 3);

load('DS161122.mat', 'ds_id', 'id_dir', 'idx_dir', 'id_dir_on', 'idx_dir_on', 'fs_idx')

%% get rfs
field_width = 17; field_height = 17;
subregion = 0; % 0: one square each time; 1: 4 synchronized squares each time
% stop = 0.5; %second
for i = 1:5
%     fs_raster{i} = get_fs_raster(datafs{i}, ds_id, 'stop', 0.5);
    fs_raster{i} = get_fs_raster(datafs{i}, ds_id);
    for cc = 1:length(ds_id)
        if fs_idx(cc,i)
            fs_raster{i}{cc} = [];
        end
    end
    fs_spike{i} = fs_get_spike(fs_raster{i});
    [rf_all{i}, rf_std{i}] = get_fs_rf(fs_spike{i}, field_width, field_height,subregion);
end

%% get histgram
for ll = 1:5 % conditions
    trigger = datafs{ll}.triggers(1:2:end);
    list = datafs{ll}.stimulus.trial_list;
    for i = 1:max(list)
        index(i, :) = find(list == i);
    end
    raster_onoff{ll} = cell(length(ds_id), 1);
    for i = 1:length(ds_id)
        if ~fs_idx(i,ll) % a boolean matrix of cell# X condition# 1: delete 0: keep; I manually fill fs_idx
            idx = get_cell_indices(datafs{ll}, ds_id(i));
            raster = get_raster(datafs{ll}.spikes{idx}, trigger, 'plot', false);
            raster_onoff{ll}{i} = raster(index);
        end
    end
end

%
bin_size = 0.02;
xx = [bin_size/2:bin_size:2-bin_size/2];
XX = [bin_size/2:bin_size:1-bin_size/2];
for ll = 1:5
    histg_onoff{ll} = cell(length(ds_id), 1);
    for cc = 1:length(raster_onoff{ll})
        if ~isempty(raster_onoff{ll}{cc})
            for p = 1:field_width * field_height
                raster_all_onoff{ll}{cc}{p} = sort(cell2mat(raster_onoff{ll}{cc}(p,:)'));
                histg_onoff{ll}{cc}{p} = hist(raster_all_onoff{ll}{cc}{p}, xx);
                for onoff = 1:2
                    raster_all{ll}{cc}{onoff}{p} = sort(cell2mat(fs_raster{ll}{cc}(p,:,onoff)'));
                    histg_all{ll}{cc}{onoff}{p} = hist(raster_all{ll}{cc}{onoff}{p}, XX);
                end
            end
        end
    end
end

%% gaussian filter
stixel_size = 30;
PixelArea = (stixel_size*4)^2/10^6; % 4µm/pix
threshold = 0.3;

clear rf_wt rf_wt_area rf_wt_area_mean rf_wt_area_ste rf_area rf_area_clean
tau = 0.05;
tt = -3*tau:bin_size:3*tau;
filter = exp(-tt.^2/(2*tau^2));

filter = filter/norm(filter);
npixel = 5;
for dir = 1:4
    for ll = 1:5
        for cc = 1:length(id_dir{dir})
            if ~fs_idx(idx_dir{dir}(cc), ll)
                for onoff = 1:2
                    % use the npixel brightest pixels to calculate a psth template
                    [~, idx] = sort(cellfun(@max,histg_all{ll}{idx_dir{dir}(cc)}{onoff}), 'descend');
                    idx = idx(1:npixel);
                    temp = histg_all{ll}{idx_dir{dir}(cc)}{onoff}(idx);
                    for i = 1:npixel
                        temp{i} = temp{i}/norm(temp{i});
                    end
                    temp_mean = mean(cell2mat(temp'));
                    temp_mean_norm = temp_mean/norm(temp_mean);

    %                 figure
    %                 subplot(1, 2, 1); plot(temp_mean_norm)

                    temp_mean_norm = conv(temp_mean_norm, filter, 'same');

    %                 subplot(1, 2, 2); plot(temp_mean_norm)
    %                 pause

                    for p = 1:field_width * field_height
                        rf_wt{ll}{dir}{cc}{onoff}(p) = histg_all{ll}{idx_dir{dir}(cc)}{onoff}{p}*temp_mean_norm';
                    end
                    rf_wt{ll}{dir}{cc}{onoff} = reshape(rf_wt{ll}{dir}{cc}{onoff}, field_width, field_height);
                    rf_wt_area{dir}{ll}{onoff}{cc} = sum(sum(rf_wt{ll}{dir}{cc}{onoff} > max(rf_wt{ll}{dir}{cc}{onoff}(:))*threshold))*PixelArea;
                end
            end
        end
        for onoff = 1:2
            rf_wt_area{dir}{ll}{onoff} = cell2mat(rf_wt_area{dir}{ll}{onoff});
            rf_wt_area_mean(dir, ll, onoff) = mean(rf_wt_area{dir}{ll}{onoff});
            rf_wt_area_ste(dir, ll, onoff) = std(rf_wt_area{dir}{ll}{onoff})/sqrt(length(rf_wt_area{dir}{ll}{onoff}));
        end
    end
end

% fit and compute rf area
for ll = 1:5
    for dir = 1:4
        clear rf_area_temp
        for onoff = 1:3
            rf_area_temp = [];
%             rf_center_temp = [];
            for cc = 1:length(id_dir{dir})
                if ~fs_idx(idx_dir{dir}(cc), ll)
                    if onoff < 3
                        data = rf_all{ll}{idx_dir{dir}(cc)}(:, :, onoff);
                    else
                        data = sum(rf_all{ll}{idx_dir{dir}(cc)}, 3);
                    end
                    if sum(sum(data > mean(data(:))+5*std(data(:))))>0
%                         figure(100)
%                         imagesc(data)
%                         colormap gray
%                         pause
                        
                        params = fit_2d_gaussian(data);
    %                     Gaussian_params{ll}{dir}{cc}{onoff} = params;
                        rf_area_temp = [rf_area_temp params.xd * params.yd * pi * PixelArea];
                        params_all{ll}{dir}{onoff}{cc} = params;
                    end
                end
            end
            
            rf_area{ll}{dir}{onoff} = rf_area_temp;
        end
    end
end


% exclude outliers
stdn = 3;
for ll = 1:5
    for dir = 1:4
        for onoff = 1:2
            notdone = 1;
            rf_area_temp = rf_area{ll}{dir}{onoff};
            while notdone
                a = length(rf_area_temp);
                rf_area_temp(rf_area_temp > std(rf_area_temp)*stdn + mean(rf_area_temp)) = [];
                b = length(rf_area_temp);
                if a == b
                    notdone = 0;
                    rf_area_clean{ll}{dir}{onoff} = rf_area_temp;
                end
            end
            rf_area_clean_mean{onoff}(ll, dir) = mean(rf_area_clean{ll}{dir}{onoff});
            rf_area_clean_ste{onoff}(ll, dir) = std(rf_area_clean{ll}{dir}{onoff})/sqrt(length(rf_area_clean{ll}{dir}{onoff}));
        end
    end
end

