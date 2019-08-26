function [rf_all] = analyzeRF161115

cd /Volumes/dusom_fieldlab/All_Staff/lab/Development/matlab/private/xyao/matlab/code/DS_new
path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2016-11-15-0/';
opt = struct('load_params', 1,'load_neurons', 1);

datafs{1} = load_data(strcat(path, 'data001-map/data001-map'), opt);
datafs{1}.names.stimulus_path = strcat(path, 'stimuli/s01.mat');
datafs{1} = load_stim_mfs(datafs{1});

datafs{2} = load_data(strcat(path, 'data002-map/data002-map'), opt);
datafs{2}.names.stimulus_path = strcat(path, 'stimuli/s02.mat');
datafs{2} = load_stim_mfs(datafs{2});

datafs{3} = load_data(strcat(path, 'data004-map/data004-map'), opt);
datafs{3}.names.stimulus_path = strcat(path, 'stimuli/s04.mat');
datafs{3} = load_stim_mfs(datafs{3});

datafs{4} = load_data(strcat(path, 'data005-map/data005-map'), opt);
datafs{4}.names.stimulus_path = strcat(path, 'stimuli/s05.mat');
datafs{4} = load_stim_mfs(datafs{4});

datafs{5} = load_data(strcat(path, 'data008-map/data008-map'), opt);
datafs{5}.names.stimulus_path = strcat(path, 'stimuli/s08.mat');
datafs{5} = load_stim_mfs(datafs{5});

% classification code below. Load pre-saved data to keep consistancy

% datapath = strcat(path, 'data007-sorted/data007-sorted');
% stimpath = strcat(path, 'stimuli/s07.txt');
% [ds_id, id_dir, id_dir_on, idx_dir, idx_dir_on] = DSClassification(datapath, stimpath, 'params_idx', [4 5], 'delta_p', 3);

load('DS161115.mat', 'ds_id', 'id_dir', 'idx_dir', 'id_dir_on', 'idx_dir_on', 'fs_idx', 'fs_idx_onoff', 'center')

%% get rfs
field_width = 13; field_height = 13;
subregion = 1;
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

%% the code below get center of rf (any one) since there are 4 identical subregion
% skip clicking by loading pre-saved data


% for cc = 1:length(ds_id)
%     figure(1)
%     id = ds_id(cc);
%     i = 5;
%     while isempty(rf_all{i}{cc})
%         i = i - 1;
%         if i == 0
%             break
%         end
%     end
%     if i>0
%         rf = rf_all{i}{cc};
% 
%         imagesc(sum(rf,3))
%         colormap gray
%         axis image
%         axis off
%         center(cc, :) = ginput;
%         close(1)
%     end
% end
% center = floor(center);
% center = max(center, ones(length(ds_id),2)*7);
% center = min(center, ones(length(ds_id),2)*20);


for ll = 1:5
    trigger = datafs{ll}.triggers(1:2:end);
    list = datafs{ll}.stimulus.trial_list;
    repeat = datafs{ll}.stimulus.repetitions;
    for i = 1:max(list)
        index(i, :) = find(list == i);
    end
    raster_onoff{ll} = cell(length(ds_id), 1);
    for i = 1:length(ds_id)
        if ~fs_idx(i,ll)
            idx = get_cell_indices(datafs{ll}, ds_id(i));
            raster = get_raster(datafs{ll}.spikes{idx}, trigger, 'plot', false);
            raster_onoff{ll}{i} = raster(index);
        end
    end
end

LL = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};
bin_size = 0.02;
xx = [bin_size/2:bin_size:2-bin_size/2];
XX = [bin_size/2:bin_size:1-bin_size/2];
for ll = 1:5
    for cc = 1:length(raster_onoff{ll})
        if ~isempty(raster_onoff{ll}{cc})
            for p = 1:169%size(raster_onoff{ll}{1}, 1)
                raster_all_onoff{ll}{cc}{p} = sort(cell2mat(raster_onoff{ll}{cc}(p,:)'));
                histg_onoff{ll}{cc}{p} = hist(raster_all_onoff{ll}{cc}{p}, xx);
                for onoff = 1:2
                    raster_all{ll}{cc}{onoff}{p} = sort(cell2mat(fs_raster{ll}{cc}(p,:,onoff)'));
                    histg_all{ll}{cc}{onoff}{p} = hist(raster_all{ll}{cc}{onoff}{p}, XX);
                end
            end
            histg_temp = reshape(histg_onoff{ll}{cc}, field_height, field_width);
            histg_temp = repmat(histg_temp, 2, 2);
            histg_temp = histg_temp(center(cc, 1)-6:center(cc, 1)+6, center(cc, 2)-6:center(cc, 2)+6);
            histg_temp = reshape(histg_temp, 1, 169);
            histg_onoff_center{ll}{cc} = histg_temp;
            for p = 1:169
                histg_center{ll}{cc}{1}{p} = histg_temp{p}(1:1/bin_size);
                histg_center{ll}{cc}{2}{p} = histg_temp{p}(1/bin_size+1:end);
            end
        end
    end
end
%% gaussian filter
PixelArea = (15*4)^2/10^6;
threshold = 0.3;
bin_size = 0.02;

clear rf_wt rf_wt_area rf_wt_area_mean rf_wt_area_ste rf_area rf_area_clean
tau = 0.5;
tt = -3*tau:bin_size:3*tau;
filter = exp(-tt.^2/(2*tau^2));

filter = filter/norm(filter);
npixel = 3;
for dir = 1:4

for cc = 1:length(id_dir{dir})
    for onoff = 1:2
        for ll = 5:-1:1
            if ~fs_idx(idx_dir{dir}(cc), ll)
                if ll == 5
                    % use the npixel brightest pixels to calculate a psth template
                    [~, idx] = sort(cellfun(@max,histg_center{ll}{idx_dir{dir}(cc)}{onoff}), 'descend');
                    idx = idx(1:npixel);
                    temp = histg_center{ll}{idx_dir{dir}(cc)}{onoff}(idx);
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
                end
                
                for p = 1:169
                    rf_wt{ll}{dir}{cc}{onoff}(p) = histg_center{ll}{idx_dir{dir}(cc)}{onoff}{p}*temp_mean_norm';
                end
                rf_wt{ll}{dir}{cc}{onoff} = reshape(rf_wt{ll}{dir}{cc}{onoff}, field_width, field_height);
%                 rf_wt_area{dir}{ll}{onoff}{cc} = sum(sum(rf_wt{ll}{dir}{cc}{onoff} > max(rf_wt{ll}{dir}{cc}{onoff}(:))*threshold))*PixelArea;
            end
        end
    end
%     for onoff = 1:2
%         rf_wt_area{dir}{ll}{onoff} = cell2mat(rf_wt_area{dir}{ll}{onoff});
%         rf_wt_area_mean(dir, ll, onoff) = mean(rf_wt_area{dir}{ll}{onoff});
%         rf_wt_area_ste(dir, ll, onoff) = std(rf_wt_area{dir}{ll}{onoff})/sqrt(length(rf_wt_area{dir}{ll}{onoff}));
%     end
end
end



%% fit RF with Gaussian

for ll = 1:5
    for cc = 1:length(ds_id)
        if ~isempty(rf_all{ll}{cc})
            rf_all_center{ll}{cc} = rf_all{ll}{cc}(center(cc, 2)-6:center(cc, 2)+6, center(cc, 1)-6:center(cc, 1)+6, :);
        end
    end
end
                
PixelArea = (15*4)^2/10^6;
% fit and compute rf area
for ll = 1:5
    for dir = 1:3
        clear rf_area_temp
        for onoff = 1:2
            rf_area_temp = [];
            for cc = 1:length(id_dir{dir})
%                 if ~isempty(rf_all_center{ll}{idx_dir{dir}(cc)})
                if ~fs_idx_onoff{onoff}(idx_dir{dir}(cc), ll)
                    data = rf_all_center{ll}{idx_dir{dir}(cc)}(:, :, onoff);
%                     data = rf_wt{ll}{dir}{cc}{onoff};
                    if sum(sum(data > mean(data(:))+2*std(data(:))))>0
%                         figure(100)
%                         imagesc(data)
%                         colormap gray

                        params = fit_2d_gaussian(data);
    %                     Gaussian_params{ll}{dir}{cc}{onoff} = params;
                        rf_area_temp = [rf_area_temp params.xd * params.yd * pi * PixelArea];
                        
%                         params.xd * params.yd * pi * PixelArea
%                         id_dir{dir}(cc)
%                         pause

                    end
                end
            end
            rf_area{ll}{dir}{onoff} = rf_area_temp;
        end
    end
end


% exclude outliers

stdn = 2;
for ll = 1:5
    for onoff = 1:2
        rf_area_temp = [];
        for dir = 1:3
            rf_area_temp = [rf_area_temp rf_area{ll}{dir}{onoff}];
        end
        notdone = 1;
        while notdone
            a = length(rf_area_temp);
            rf_area_temp(rf_area_temp > mean(rf_area_temp) + std(rf_area_temp)*stdn) = [];
            rf_area_temp(rf_area_temp < mean(rf_area_temp) - std(rf_area_temp)*stdn) = [];
            b = length(rf_area_temp);
            if a == b
                notdone = 0;
                rf_area_clean{ll}{onoff} = rf_area_temp;
            end
        end
        rf_area_clean_mean{onoff}(ll) = mean(rf_area_clean{ll}{onoff});
        rf_area_clean_ste{onoff}(ll) = std(rf_area_clean{ll}{onoff})/sqrt(length(rf_area_clean{ll}{onoff}));
    end
end


