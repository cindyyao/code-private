%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun{1} = load_data('/Volumes/lab/analysis/2013-02-14-0/data006/data006', opt);
datarun{2} = load_data('/Volumes/lab/analysis/2013-02-21-0/data004/data004', opt);
datarun{1} = get_rfs(datarun{1}, 'all');
datarun{2} = get_rfs(datarun{2}, 'all');

%% get cell ids
cell_type = {'ON transient', 'ON brisk transient', 'OFF brisk transient', 'OFF transient', 'OFF sustained'};
n = length(datarun);
m = length(cell_type);
cell_ids = cell(n, 1);
for i = 1:n
    cell_id_temp = cell(m, 1); 
    for j = 1:m
        id_temp = get_cell_ids(datarun{i}, cell_type{j});
        idx_temp = get_cell_indices(datarun{i}, cell_type{j});
        cell_id_temp{j} = struct('name', cell_type{j}, 'cell_ids', id_temp, 'cell_indices', idx_temp);
    end
    cell_ids{i} = cell_id_temp;
end

%% fit individual cell
id_all = cell(n, 1);
idx_all = cell(n, 1);
for i = 1:n
    id_all_temp = [];
    idx_all_temp = [];
    for j = 1:m
        id_all_temp = [id_all_temp cell_ids{i}{j}.cell_ids];
        idx_all_temp = [idx_all_temp cell_ids{i}{j}.cell_indices];
    end
    id_all{i} = sort(id_all_temp);
    idx_all{i} = sort(idx_all_temp);
end

% % fit
% fit_instruction = struct('verbose', true);
% for i = 1:n
%     datarun{i} = compute_sta_fits_sequence(datarun{i}, id_all{i}, 'verbose', true);
% end

%% Gaussian Volume Ratio
load('vgat_130214.mat')
datarun{1}.matlab.sta_fits = datarun1_sta_fits;
datarun{2}.matlab.sta_fits = datarun2_sta_fits;
sta_fits{1} = datarun1_sta_fits;
sta_fits{2} = datarun2_sta_fits;

% surround gaussian / center gaussian
VolumeRatio = cell(n, 1);
for i = 1:n
    volumeratio_temp = cell(m, 1);
    for j = 1:m
        id_temp = cell_ids{i}{j}.cell_ids;
        idx_temp = cell_ids{i}{j}.cell_indices;
        ratio_temp = [];
        for cc = 1:length(id_temp)
            fit_temp = sta_fits{i}{idx_temp(cc)};
%             if ~isempty(fit_temp)
                ratio_temp = [ratio_temp fit_temp.surround_sd_scale^2*fit_temp.surround_amp_scale];
%             end
        end
        volumeratio_temp{j} = ratio_temp;
    end
    VolumeRatio{i} = volumeratio_temp;
end

% surround region volume / center region volume
VolumeRatio = cell(n, 1);
for i = 1:n
    volumeratio_temp = cell(m, 1);
    for j = 1:m
        id_temp = cell_ids{i}{j}.cell_ids;
        idx_temp = cell_ids{i}{j}.cell_indices;
        ratio_temp = [];
        for cc = 1:length(id_temp)
            params = get_params(datarun{i}, id_temp(cc));
            fit_para_center = struct('center_point_x', params(1), 'center_point_y', ...
                params(2), 'sd_x', params(3), 'sd_y', params(4), ...
                'rotation_angle', params(5), 'x_dim', params(9), 'y_dim', params(10));
            fit_para_surround = struct('center_point_x', params(1), 'center_point_y', ...
                params(2), 'sd_x', params(3)*params(13), 'sd_y', params(4)*params(13),...
                'amp_scale', params(14), 'rotation_angle', params(5), 'x_dim', params(9),...
                'y_dim', params(10));
            fit1 = make_Gaussian_two_d(fit_para_center);
            fit2 = make_Gaussian_two_d(fit_para_surround);
            fit = fit1 - fit2;
            center_v = sum(fit(fit > 0));
            surround_v = sum(fit(fit < 0));
            ratio = -surround_v/center_v;
            ratio_temp = [ratio_temp ratio];
        end
        volumeratio_temp{j} = ratio_temp;
    end
    VolumeRatio{i} = volumeratio_temp;
end

binN = 10;
figure
for j = 1:m
    ratio_temp = [VolumeRatio{1}{j} VolumeRatio{2}{j}];
    XX = linspace(min(ratio_temp), max(ratio_temp), binN);
    h_ko = hist(VolumeRatio{1}{j}, XX);
    h_wt = hist(VolumeRatio{2}{j}, XX);
    h = [h_wt' h_ko'];
    subplot(2, 3, j)
    bar(XX, h)
    title(cell_type{j})
    legend('WT', 'KO')
    xlabel('Volume Ratio')
    ylabel('# of cells')
    for i = 1:n
        VolumeRatioMean{i}(j) = mean(VolumeRatio{i}{j});
        VolumeRatioSte{i}(j) = std(VolumeRatio{i}{j})/sqrt(length(VolumeRatio{i}{j}));
    end
    [~,p] = ttest2(VolumeRatio{1}{j}, VolumeRatio{2}{j});
    PValue(j) = p;
end

% for i = 1:n
%     for cc = 1:length(id_all{i})
%         plot_fit(datarun{i}, id_all{i}(cc));
%         pause
%     end
% end
%% Bootstrap sta
% get centered (and cropped) sta
cell_type_ = {'ON_transient', 'ON_brisk_transient', 'OFF_brisk_transient', 'OFF_transient', 'OFF_sustained'};
for i = 1:n
    centered_sta{i} = struct();
    for j = 1:m
        id_temp = cell_ids{i}{j}.cell_indices;
        sta = datarun{i}.stas.stas(id_temp);
        new_sta = [];
        for k = 1:length(sta)
            new_sta_temp = {crop_sta(sta{k})};
            if ~isempty(new_sta_temp{1})
                new_sta = [new_sta new_sta_temp];
            end
        end
        centered_sta{i} = setfield(centered_sta{i}, cell_type_{j}, new_sta);
    end
end
        
% calculate bootstrap mean

clear bootsam mean_sta_all
scale = 10;
for i = 1:n
    for j = 1:m
        sta = getfield(centered_sta{i}, cell_type_{j});
        celln = length(sta);
        sta_all = cat(5, sta{:});
        for k = 1:celln
            bootsam{i}{j}(:, k) = randsample(celln, celln, true);
            mean_sta = mean(sta_all(:, :, :, :, bootsam{i}{j}(:, k)), 5);
            mean_sta_all{i}{j}{k} = downscale_sta(mean_sta, scale);
        end
    end
end

% fit mean sta
center = size(mean_sta_all{1}{1}{1}, 1)/2+0.5;
for i = 1:n
    for j = 1:m
        celln = length(mean_sta_all{i}{j});
        for cc = 1:celln
            sig_sti = significant_stixels(mean_sta_all{i}{j}{cc}, 'thresh', 3.5);
            fit_instructions = struct('initial_center_point_x', center, 'initial_center_point_y', center, 'sig_stixels', sig_sti);
            mean_fit{i}{j}{cc} = fit_sta_sequence(mean_sta_all{i}{j}{cc}, 'verbose', true, 'fit_instructions', fit_instructions);
        end
    end
end


% surround gaussian / center gaussian
VolumeRatio = cell(n, 1);
for i = 1:n
    volumeratio_temp = cell(m, 1);
    for j = 1:m
        ratio_temp = [];
        celln = length(mean_sta_all{i}{j});
        for cc = 1:celln
            fit_temp = mean_fit{i}{j}{cc};
            ratio_temp = [ratio_temp fit_temp.surround_sd_scale^2*fit_temp.surround_amp_scale];
        end
        volumeratio_temp{j} = ratio_temp;
    end
    VolumeRatio{i} = volumeratio_temp;
end

% surround region volume / center region volume
VolumeRatio = cell(n, 1);
for i = 1:n
    volumeratio_temp = cell(m, 1);
    for j = 1:m
        ratio_temp = [];
        celln = length(mean_sta_all{i}{j});
        for cc = 1:celln
            params = get_params_(mean_fit{i}{j}{cc});
            fit_para_center = struct('center_point_x', params(1), 'center_point_y', ...
                params(2), 'sd_x', params(3), 'sd_y', params(4), ...
                'rotation_angle', params(5), 'x_dim', params(9), 'y_dim', params(10));
            fit_para_surround = struct('center_point_x', params(1), 'center_point_y', ...
                params(2), 'sd_x', params(3)*params(13), 'sd_y', params(4)*params(13),...
                'amp_scale', params(14), 'rotation_angle', params(5), 'x_dim', params(9),...
                'y_dim', params(10));
            fit1 = make_Gaussian_two_d(fit_para_center);
            fit2 = make_Gaussian_two_d(fit_para_surround);
            fit = fit1 - fit2;
            center_v = sum(fit(fit > 0));
            surround_v = sum(fit(fit < 0));
            ratio = -surround_v/center_v;
            ratio_temp = [ratio_temp ratio];
        end
        volumeratio_temp{j} = ratio_temp;
    end
    VolumeRatio{i} = volumeratio_temp;
end

binN = 10;
figure
for j = 1:m
    ratio_temp = [VolumeRatio{1}{j} VolumeRatio{2}{j}];
    XX = linspace(min(ratio_temp), max(ratio_temp), binN);
    h_ko = hist(VolumeRatio{1}{j}, XX);
    h_wt = hist(VolumeRatio{2}{j}, XX);
    h = [h_wt' h_ko'];
    subplot(2, 3, j)
    bar(XX, h)
    title(cell_type{j})
    legend('WT', 'KO')
    xlabel('Volume Ratio')
    ylabel('# of cells')
    for i = 1:n
        VolumeRatioMean{i}(j) = mean(VolumeRatio{i}{j});
        VolumeRatioSte{i}(j) = std(VolumeRatio{i}{j})/sqrt(length(VolumeRatio{i}{j}));
    end
    [~,p] = ttest2(VolumeRatio{1}{j}, VolumeRatio{2}{j});
    PValue(j) = p;
end

%% SNR/volume ratio correlation

% individual cell
clear snr_i
for i = 1:n
    for j = 1:m
        id_temp = cell_ids{i}{j}.cell_ids;
        idx_temp = cell_ids{i}{j}.cell_indices;
        for cc = 1:length(id_temp)
            rf = datarun{i}.stas.rfs{idx_temp(cc)};
            if j < 3
                rf_sort = sort(rf(:), 'descend');
                snr_i{i}{j}(cc) = mean(rf_sort(1:4))/std(rf_sort);
            else
                rf_sort = sort(rf(:));
                snr_i{i}{j}(cc) = -mean(rf_sort(1:4))/std(rf_sort);
            end
        end
    end
end

for i = 1:n
    figure
    for j = 1:m
        subplot(2, 3, j)
        plot(VolumeRatio{i}{j}, snr_i{i}{j}, 'o')
        coef = corrcoef(VolumeRatio{i}{j}, snr_i{i}{j});
        corr_coef(i, j) = coef(1, 2);
        xlabel('VolumeRatio')
        ylabel('SNR')
        title([cell_type{j} '  r^2 = ' num2str(coef(1,2))])
    end
end

% bootstrapped mean
clear snr_mean
for i = 1:n
    for j = 1:m
        for cc = 1:length(mean_sta_all{i}{j})
            rf = rf_from_sta(mean_sta_all{i}{j}{cc});
            if j < 3
                rf_sort = sort(rf(:), 'descend');
                snr_mean{i}{j}(cc) = mean(rf_sort(1:4))/std(rf_sort);
            else
                rf_sort = sort(rf(:));
                snr_mean{i}{j}(cc) = -mean(rf_sort(1:4))/std(rf_sort);
            end
        end
    end
end

for i = 1:n
    figure
    for j = 1:m
        subplot(2, 3, j)
        plot(VolumeRatio{i}{j}, snr_mean{i}{j}, 'o')
        coef = corrcoef(VolumeRatio{i}{j}, snr_mean{i}{j});
        corr_coef(i, j) = coef(1, 2);
        xlabel('VolumeRatio')
        ylabel('SNR')
        title([cell_type{j} '  r^2 = ' num2str(coef(1,2))])
    end
end

%% rmse volume ratio correlation
% individual cell
clear rmse_i
for i = 1:n
    for j = 1:m
        id_temp = cell_ids{i}{j}.cell_ids;
        idx_temp = cell_ids{i}{j}.cell_indices;
        for cc = 1:length(id_temp)
            rmse_i{i}{j}(cc) = datarun{i}.matlab.sta_fits{idx_temp(cc)}.rmse;
        end
    end
end
for i = 1:n
    figure
    for j = 1:m
        subplot(2, 3, j)
        plot(VolumeRatio{i}{j}, rmse_i{i}{j}, 'o')
        coef = corrcoef(VolumeRatio{i}{j}, rmse_i{i}{j});
        corr_coef(i, j) = coef(1, 2);
        xlabel('VolumeRatio')
        ylabel('rmse')
        title([cell_type{j} '  r^2 = ' num2str(coef(1,2))])
    end
end

% bootstrapped mean sta
clear rmse_mean
for i = 1:n
    for j = 1:m
        for cc = 1:length(mean_fit{i}{j})
            rmse_mean{i}{j}(cc) = mean_fit{i}{j}{cc}.rmse;
        end
    end
end
for i = 1:n
    figure
    for j = 1:m
        subplot(2, 3, j)
        plot(VolumeRatio{i}{j}, rmse_mean{i}{j}, 'o')
        coef = corrcoef(VolumeRatio{i}{j}, rmse_mean{i}{j});
        corr_coef(i, j) = coef(1, 2);
        xlabel('VolumeRatio')
        ylabel('rmse')
        title([cell_type{j} '  r^2 = ' num2str(coef(1,2))])
    end
end