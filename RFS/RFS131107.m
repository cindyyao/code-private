%% VGAT-/- rod cone light level comparison

opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);

% % fit individual neurons

% ndf 3
datarun{1} = load_data('/Analysis/xyao/2013-11-07-0/data000/data000', opt);

cell_type = {'ON brisk transient', 'ON transient', 'ON transient2', 'OFF brisk transient', ...
    'OFF transient', 'OFF transient2', 'OFF sustained', 'OFF sustained2', 'OFF slow'};
n = length(cell_type);

cell_id = [];
cell_idx = [];

for i =  1:n
    cell_id_temp = get_cell_ids(datarun{1}, cell_type{i});
    cell_idx_temp = get_cell_indices(datarun{1}, cell_type{i});
    
    cell_id = [cell_id cell_id_temp];
    cell_idx = [cell_idx cell_idx_temp];
end

datarun{1} = compute_sta_fits_sequence(datarun{1}, cell_id, 'verbose', true);

% ndf 0

cell_type = {'ON brisk transient', 'ON transient', 'ON sustained', 'OFF brisk transient', ...
    'OFF transient', 'OFF sustained', 'OFF slow1', 'OFF slow2'};

datarun{2} = load_data('/Analysis/xyao/2013-11-07-0/data002/data002', opt);
n = length(cell_type);

cell_id = [];
cell_idx = [];

for i =  1:n
    cell_id_temp = get_cell_ids(datarun{2}, cell_type{i});
    cell_idx_temp = get_cell_indices(datarun{2}, cell_type{i});
    
    cell_id = [cell_id cell_id_temp];
    cell_idx = [cell_idx cell_idx_temp];
end

datarun{2} = compute_sta_fits_sequence(datarun{2}, cell_id, 'verbose', true);

% % get snls

% ndf 3
datarun{1} = load_java_movie(datarun{1}, '/Volumes/lab/acquisition/movie-xml/BW-20-4-0.48-11111-30x30-60.35.xml');
datarun{1} = get_sta_summaries(datarun{1}, 'all');

datarun{1} = get_snls(datarun{1}, 'all');

% ndf 0
datarun{2} = load_java_movie(datarun{2}, '/Volumes/lab/acquisition/movie-xml/BW-10-2-0.48-11111-60x60-60.35.xml');
datarun{2} = get_sta_summaries(datarun{2}, 'all');

datarun{2} = get_snls(datarun{2}, 'all');


save('RFS131107.mat', 'datarun', '-append')


%% ei mapping

cell_type = {'ON brisk transient', 'ON transient', 'ON sustained', 'OFF brisk transient', ...
    'OFF transient', 'OFF sustained', 'OFF slow1'};

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun{1}, ...
    datarun{2}, cell_type);

[sta_f_weighted, sta_f_n, sta_c] = mean_rf2(datarun{1}, datarun{2}, cell_type, ...
    cell_id, cell_idx);

%% RF fitting

load('RFS131107.mat')
% fit mean RF

n = length(cell_type);
fit_sta_c = cell(n, length(datarun));
params_c = cell(n, length(datarun));

for j = 1:length(datarun)
    d = size(datarun{j}.stas.stas{1});
    for i = 1:n
        f = 40/d(1);
%         temp_marks_sta = significant_stixels(sta_c{i, j}, 'thresh', 5, 'time', 'max');
%         fit_ins = struct('sig_stixels', temp_marks_sta);
        fit_temp = fit_sta_sequence(sta_c{i, j}); %, 'fit_instructions', fit_ins);
        fit_sta_c{i, j} = fit_temp;
        params_temp = zeros(1, 8);
        params_temp(1) = fit_temp.center_point_x*f;
        params_temp(2) = fit_temp.center_point_y*f;
        params_temp(3) = fit_temp.center_sd_x*f;
        params_temp(4) = fit_temp.center_sd_y*f;
        params_temp(5) = fit_temp.center_rotation_angle;
        params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale*f;
        params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale*f;
        params_temp(8) = fit_temp.surround_amp_scale;
        params_c{i, j} = params_temp;
        
        fprintf(['i = ' num2str(i) '  j = ' num2str(j) '\n'])
    end
   

end


save('RFS131107.mat', 'sta_f_weighted', 'sta_f_n', 'sta_c', 'fit_sta_c', 'params_c', '-append')


%% individual cell summaries

display_rate = 60.35;
refresh_rate = [4 2];

ct = 7;
cell_num = size(cell_id{ct}, 1);

    for cn = 1:cell_num
          
        FigHandle = figure;
        set(FigHandle,  'Position', [100, 100, 1200, 730]);
        
        for i = 1:length(datarun)
            idx = cell_idx{ct}(cn, i);
            com = datarun{i}.stas.rf_coms{idx};
            gen_signals{i} = datarun{i}.stas.snls{idx}.gen_signal;
            spikes{i} = datarun{i}.stas.snls{idx}.spikes;
            tc(:, i) = datarun{i}.stas.time_courses{idx};
            r = datarun{i}.stimulus.field_height;
            
            rf{i} = datarun{i}.stas.rfs{idx};
            pix = [1:r] + 0.5;
            [XX YY] = meshgrid(pix);
            xdis = com(1) - XX;
            ydis = com(2) - YY; 

            dis = zeros(r, r);
            for k = 1:r
                for j = 1:r
                    dis(k, j) = norm([xdis(k, j) ydis(k, j)]);
                end
            end
            
            dis_1 = reshape(dis, 1, r^2);
            dis_1 = dis_1*40/r;
            rf_1 = reshape(rf{i}, 1, r^2);
            [Dis1{i} RF1{i}] = curve_from_binning(dis_1, rf_1, 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:20);
            Dis1{i} = [-Dis1{i}(end:-1:1); Dis1{i}];
            RF1{i} = [RF1{i}(end:-1:1); RF1{i}];

            [X(:, i), Y(:, i)] = curve_from_binning(gen_signals{i}, spikes{i}, 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
            Y(:, i) = Y(:, i)*display_rate/refresh_rate(i);
        

       
        end
        
        % unscaled
        subplot(2, 3, 1)
        plot_rf(datarun{1}, cell_id{ct}(cn, 1), 'scale', 8);
        title([num2str(cell_id{ct}(cn, 1)) '  NDF 3'])

 
        subplot(2, 3, 2)
        plot_rf(datarun{2}, cell_id{ct}(cn, 2), 'scale', 8);
        title([num2str(cell_id{ct}(cn, 2)) '  NDF 0'])
        

        % calculate the scaling factor
        sigma_1 = std(gen_signals{1});
        sigma_2 = std(gen_signals{2});
        f_gs = sigma_1/sigma_2;
        
        % apply this factor to strf and recalculate the sf and tf
        datarun_2 = datarun{2};
        datarun_2.stas.stas{idx} = datarun{2}.stas.stas{idx}*f_gs;
        datarun_2 = get_sta_summaries(datarun_2, datarun_2.cell_ids(idx));
        
        
        gen_signals{2} = gen_signals{2}*f_gs;
        tc(:, 2) = datarun_2.stas.time_courses{idx};
        rf{2} = datarun_2.stas.rfs{idx};
        pix = [1:r] + 0.5;
        [XX YY] = meshgrid(pix);
        xdis = com(1) - XX;
        ydis = com(2) - YY; 

        dis = zeros(r, r);
        for k = 1:r
            for j = 1:r
                dis(k, j) = norm([xdis(k, j) ydis(k, j)]);
            end
        end
            
        dis_1 = reshape(dis, 1, r^2);
        dis_1 = dis_1*40/r;
        rf_1 = reshape(rf{2}, 1, r^2);
        [Dis1{2} RF1{2}] = curve_from_binning(dis_1, rf_1, 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:20);
        Dis1{2} = [-Dis1{2}(end:-1:1); Dis1{2}];
        RF1{2} = [RF1{2}(end:-1:1); RF1{2}];

        RF1_std{1} = RF1{1}/std(RF1{1});
        RF1_std{2} = RF1{2}/std(RF1{2});
        
        t1 = -1/display_rate*refresh_rate(1)*[14:-1:0];
        t2 = -1/display_rate*refresh_rate(2)*[14:-1:0];


        tc_std(:, 1) = tc(:, 1)/std(tc(:, 1));
        tc_std(:, 2) = tc(:, 2)/std(tc(:, 2));

        [X(:, 2), Y(:, 2)] = curve_from_binning(gen_signals{2}, spikes{2}, 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
        Y(:, 2) = Y(:, 2)*display_rate/refresh_rate(2);
        
        subplot(2, 3, 3);
        plot(Dis1{1}, RF1_std{1}, 'b');
        hold on
        plot(Dis1{2}, RF1_std{2}, 'r');
        xlim([-10 10])
        title('Spatial Filter')
        h_legend = legend('NDF 3', 'NDF 0', 'location', 'northwest');
%         set(h_legend,'FontSize',10);

        subplot(2, 3, [4 5]);
        plot(t1, tc_std(:, 1), 'b');
        hold on
        plot(t2, tc_std(:, 2), 'r');
        xlabel('time to spike (s)')
        ylabel('STA contrast')
        xlim([-0.75 0])
        ylim([min(tc_std(:))*1.2 max(tc_std(:))*1.2])
        title('Temporal Filter')

        subplot(2, 3, 6);
        plot(X(:, 1), Y(:, 1), 'b');
        hold on
        plot(X(:, 2), Y(:, 2), 'r');
        xlabel('generator signal')
        ylabel('spike rate(HZ)')
        title('Nonlinearity')
 
%         subplot(3, 3, 7)
%         plot_ei(datarun{1}, cell_id{ct}(cn, 1))
%         subplot(3, 3, 8)
%         plot_ei(datarun{2}, cell_id{ct}(cn, 2))
       

    end


%% cell type summaries


% get mean spatial profile
cell_type_n = length(cell_type);

for ct = 1:7;
cell_num = size(cell_id{ct}, 1);

temp_sf = zeros(2, 40, cell_num, 2);
temp_tf = zeros(2, 15, cell_num);



    for cn = 1:cell_num
          
        
        for i = 1:2
            idx = cell_idx{ct}(cn, i);
            com = datarun{i}.stas.rf_coms{idx};
            gen_signals{i} = datarun{i}.stas.snls{idx}.gen_signal;
            spikes{i} = datarun{i}.stas.snls{idx}.spikes;
            tc(:, i) = datarun{i}.stas.time_courses{idx};
            r = datarun{i}.stimulus.field_height;
            
            rf{i} = datarun{i}.stas.rfs{idx};
            pix = [1:r] + 0.5;
            [XX YY] = meshgrid(pix);
            xdis = com(1) - XX;
            ydis = com(2) - YY; 

            dis = zeros(r, r);
            for k = 1:r
                for j = 1:r
                    dis(k, j) = norm([xdis(k, j) ydis(k, j)]);
                end
            end
            
            dis_1 = reshape(dis, 1, r^2);
            dis_1 = dis_1*40/r;
            rf_1 = reshape(rf{i}, 1, r^2);
            [Dis1{i} RF1{i}] = curve_from_binning(dis_1, rf_1, 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:20);
            Dis1{i} = [-Dis1{i}(end:-1:1); Dis1{i}];
            RF1{i} = [RF1{i}(end:-1:1); RF1{i}];

            [X(:, i), Y(:, i)] = curve_from_binning(gen_signals{i}, spikes{i}, 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
            Y(:, i) = Y(:, i)*display_rate/refresh_rate(i);
        

       
        end
        

        % calculate the scaling factor
        sigma_1 = std(gen_signals{1});
        sigma_2 = std(gen_signals{2});
        f_gs = sigma_1/sigma_2;
        
        % apply this factor to strf and recalculate the sf and tf
        datarun_2 = datarun{2};
        datarun_2.stas.stas{idx} = datarun{2}.stas.stas{idx}*f_gs;
        datarun_2 = get_sta_summaries(datarun_2, datarun_2.cell_ids(idx));
        
        
        gen_signals{2} = gen_signals{2}*f_gs;
        tc(:, 2) = datarun_2.stas.time_courses{idx};
        rf{2} = datarun_2.stas.rfs{idx};
        pix = [1:r] + 0.5;
        [XX YY] = meshgrid(pix);
        xdis = com(1) - XX;
        ydis = com(2) - YY; 

        dis = zeros(r, r);
        for k = 1:r
            for j = 1:r
                dis(k, j) = norm([xdis(k, j) ydis(k, j)]);
            end
        end
            
        dis_1 = reshape(dis, 1, r^2);
        dis_1 = dis_1*40/r;
        rf_1 = reshape(rf{2}, 1, r^2);
        [Dis1{2} RF1{2}] = curve_from_binning(dis_1, rf_1, 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:20);
        Dis1{2} = [-Dis1{2}(end:-1:1); Dis1{2}];
        RF1{2} = [RF1{2}(end:-1:1); RF1{2}];

        RF1_std{1} = RF1{1}/std(RF1{1});
        RF1_std{2} = RF1{2}/std(RF1{2});

        tc_std(:, 1) = tc(:, 1)/std(tc(:, 1));
        tc_std(:, 2) = tc(:, 2)/std(tc(:, 2));

        [X(:, 2), Y(:, 2)] = curve_from_binning(gen_signals{2}, spikes{2}, 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
        Y(:, 2) = Y(:, 2)*display_rate/refresh_rate(2);
        
        
        temp_sf(1, :, cn, 1) = Dis1{1}; 
        temp_sf(1, :, cn, 2) = Dis1{2}; 
        temp_sf(2, :, cn, 1) = RF1_std{1};
        temp_sf(2, :, cn, 2) = RF1_std{2};
        
        temp_tf(1, :, cn) = tc_std(:, 1);
        temp_tf(2, :, cn) = tc_std(:, 2);
        
        
    end
    
    spatial_profile{ct} = temp_sf;
    temporal_filter{ct} = temp_tf;
end
    
for ct = 1:7
    Temporal_Filter{ct} = mean(temporal_filter{ct}, 3);
    cell_num = size(cell_id{ct}, 1);
    Spatial_Profile_temp = zeros(21, cell_num, 2);

    for cn = 1:cell_num
        for cd = 1:2
            Spatial_Profile_temp(:, cn, cd) = interp1(spatial_profile{ct}(1, :, cn, cd), spatial_profile{ct}(2, :, cn, cd), -10:10);
        end
    end
    Spatial_Profile{ct} = squeeze(mean(Spatial_Profile_temp, 2));
end

% save('snl130221_avg.mat', 'Spatial_Profile', 'Temporal_Filter', 'sta_c')

%% get average snls

display_rate = 60.35;
refresh_rate = [4 2];

NL = cell(cell_type_n, 1);

for ct = 1:cell_type_n
    n = size(cell_id{ct}, 1);
    NL_temp = zeros(2, 20, n, 2);
    for cc = 1:n
        for cd = 1:2
            idx = cell_idx{ct}(cc, cd);
            gen_signals = datarun{cd}.stas.snls{idx}.gen_signal;
            gen_signals = gen_signals/std(gen_signals);
            spikes = datarun{cd}.stas.snls{idx}.spikes;

            [X, Y] = curve_from_binning(gen_signals, spikes, ...
                'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
            Y = Y*display_rate/refresh_rate(cd);
            NL_temp(1, :, cc, cd) = X;
            NL_temp(2, :, cc, cd) = Y;
        end
        temp = NL_temp(2, :, cc, :);
        NL_temp(2, :, cc, :) = temp/max(temp(:));
    end
    NL{ct} = NL_temp;
end

XX = zeros(20, 2, cell_type_n);
for ct = 1:cell_type_n
    for cd = 1:2
        xx = squeeze(NL{ct}(1, :, :, cd));
        xx = xx';
        XX_temp = mean(xx);
        XX_temp(1) = max(xx(:, 1));
        XX_temp(20) = min(xx(:, 20));
        XX(:, cd, ct) = XX_temp;
    end
end

NL_Y_mean = zeros(20, 2, cell_type_n);
Nonlinearity_stev = zeros(20, 2, cell_type_n);

for ct = 1:cell_type_n
    n = size(cell_id{ct}, 1);
    NL_Y_temp = zeros(20, n, 2);
    for cd = 1:2
        for cc = 1:n
            YY_temp = interp1(NL{ct}(1, :, cc, cd), NL{ct}(2, :, cc, cd), XX(:, cd, ct));
            NL_Y_temp(:, cc, cd) = YY_temp;
        end
    end
    NL_Y_mean_temp = squeeze(mean(NL_Y_temp, 2));
    Nonlinearity_stev_temp = squeeze(std(NL_Y_temp, 0, 2))/sqrt(n);
    NL_Y_mean(:, :, ct) = NL_Y_mean_temp;
    Nonlinearity_stev(:, :, ct) = Nonlinearity_stev_temp;
end

Nonlinearity_mean(:, :, :, 1) = XX;
Nonlinearity_mean(:, :, :, 2) = NL_Y_mean;

%% surround strengths & center size

% datarun{1} = load_sta_fits(datarun{1});
% datarun{2} = load_sta_fits(datarun{2});

n = length(cell_type);
surround_strengths_mean = zeros(n, 2);
surround_strengths_var = zeros(n, 2);
surround_strengths_stev = zeros(n, 2);
surround_strengths_p = zeros(n, 1);

center_size_mean = zeros(n, 2);
center_size_var = zeros(n, 2);
center_size_stev = zeros(n, 2);
center_size_p = zeros(n, 1);

for ct = 1:n
    idx = cell_idx{ct};
    cell_num = size(idx, 1);
    for cc = 1:cell_num
        for i = 1:2
        r = datarun{i}.stimulus.field_height;
        params_temp = zeros(7, 1);
        params_temp(1) = datarun{i}.matlab.sta_fits{cell_idx{ct}(cc, i)}.center_point_x*40/r;
        params_temp(2) = datarun{i}.matlab.sta_fits{cell_idx{ct}(cc, i)}.center_point_y*40/r;
        params_temp(3) = datarun{i}.matlab.sta_fits{cell_idx{ct}(cc, i)}.center_sd_x*40/r;
        params_temp(4) = datarun{i}.matlab.sta_fits{cell_idx{ct}(cc, i)}.center_sd_y*40/r;
        params_temp(5) = datarun{i}.matlab.sta_fits{cell_idx{ct}(cc, i)}.center_rotation_angle;
        params_temp(6) = datarun{i}.matlab.sta_fits{cell_idx{ct}(cc, i)}.surround_sd_scale;
        params_temp(7) = datarun{i}.matlab.sta_fits{cell_idx{ct}(cc, i)}.surround_amp_scale;
        params{ct}{i}{cc} = params_temp;
        surround_strengths{ct}{i}(cc) = params_temp(7) * params_temp(6)^2;
        center_size_temp{ct}{i}(cc) = params_temp(3) * params_temp(4);
        end
    end
    
    for i = 1:2
        surround_strengths_mean(ct, i) = mean(surround_strengths{ct}{i});
        surround_strengths_var(ct, i) = var(surround_strengths{ct}{i});
        surround_strengths_stev(ct, i) = std(surround_strengths{ct}{i})/sqrt(cell_num);
        
        center_size_mean(ct, i) = mean(center_size_temp{ct}{i});
        center_size_var(ct, i) = var(center_size_temp{ct}{i});
        center_size_stev(ct, i) = std(center_size_temp{ct}{i})/sqrt(cell_num);   

    end
    
    [~, surround_strengths_p(ct)] = ttest(surround_strengths{ct}{1} - surround_strengths{ct}{2});
    [~, center_size_p(ct)] = ttest(center_size_temp{ct}{1} - center_size_temp{ct}{2});

end


surround = struct('mean', surround_strengths_mean, 'variance', surround_strengths_var, ...
    'stev', surround_strengths_stev, 'p_value', surround_strengths_p);
center = struct('mean', center_size_mean, 'variance', center_size_var, 'stev', ...
    center_size_stev, 'p_value', center_size_p);

%% degree of transience(DOT)

smooth_f = 0.01;
DOT = cell(cell_type_n, 1);
t = -14:0;
for ct = 1:cell_type_n
    n = size(cell_id{ct}, 1);
    DOT_temp = zeros(n, 2);
    for cn = 1:n
        for cd = 1:2
            tf_smooth = interp1(-14:0, temporal_filter{ct}(cd, :, cn), -14:smooth_f:0);
            area = smooth_f * trapz(tf_smooth);
            tf_smooth_p = tf_smooth;
            tf_smooth_p(tf_smooth_p<0) = 0;
            area_p = smooth_f * trapz(tf_smooth_p);
            area_n = area_p - area;
            DOT_temp(cn, cd) = 1 - abs(area/(area_p + area_n));
        end
    end
    DOT{ct} = DOT_temp;
end

DOT_mean = zeros(4, 2);
DOT_stev = zeros(4, 2);
DOT_pvalue = zeros(4, 1);

for ct = 1:cell_type_n
    n = size(cell_id{ct}, 1);
    DOT_mean(ct, :) = mean(DOT{ct});
    DOT_stev(ct, :) = std(DOT{ct})/sqrt(n);
    [~, DOT_pvalue(ct)] = ttest(DOT{ct}(:, 1), DOT{ct}(:, 2));
end

dot = struct('mean', DOT_mean, 'stev', DOT_stev, 'p_value', DOT_pvalue);


%% plot cell type summaries

% load('snls130221_avg.mat')

for ct = 1:cell_type_n
    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 1550, 700]);
    
    % NDF 4 sta image
    r = datarun{1}.stimulus.field_height;
    sta = sta_c{ct, 1};
    if cell_type{ct}(2) == 'N' 
        [~, b] = max(sta(:));
    else
        [~, b] = min(sta(:));  
    end

    z = ceil(b/r^2);
    sta_1 = sta(:, :, 1, z);
    sta_1 = matrix_scaled_up(sta_1, 8);
    
    subplot(2, 4, 1)
    im = norm_image(sta_1);
    image(im)
    title('NDF 4')
    axis off
    
    % NDF 0 sta image

    r = datarun{2}.stimulus.field_height;
    sta = sta_c{ct, 2};
    if cell_type{ct}(2) == 'N' 
        [~, b] = max(sta(:));
    else
        [~, b] = min(sta(:));  
    end

    z = ceil(b/r^2);
    sta_1 = sta(:, :, 1, z);
    sta_1 = matrix_scaled_up(sta_1, 8);
    
    subplot(2, 4, 2)
    im = norm_image(sta_1);
    image(im)
    title('NDF 0')
    axis off
    
    subplot(2, 4, 3)
    plot_rf_summaries(datarun{1}, cell_id{ct}(:, 1), 'fit_color', 'b', 'scale', 2)
    hold on
    plot_rf_summaries(datarun{2}, cell_id{ct}(:, 2), 'fit_color', 'r')
    title('Receptive Fields Mosaic')
    
%     % bar
%     subplot(2, 3, [7:9])
%     
%     xtick = {'Surround Strengths', 'Center Size', 'Degree of Transience'};
%     model_series = [surround.mean(ct, 1) surround.mean(ct, 2); center.mean(ct, 1)/2 center.mean(ct, 2)/2; dot.mean(ct, 1) dot.mean(ct, 2)];   
%     model_error = [surround.stev(ct, 1) surround.stev(ct, 2); center.stev(ct, 1)/2 center.stev(ct, 2)/2; dot.stev(ct, 1) dot.stev(ct, 2)];
%     h = bar(model_series);
%     set(h,'BarWidth',1); % The bars will now touch each other
% 
%     set(gca,'XTicklabel',xtick)
%     % ylabel('Surround Strengths')
%     legend('NDF 4','NDF 0', 'location', 'northwest');
%     hold on;
%  
%     numgroups = size(model_series, 1); 
%     numbars = size(model_series, 2); 
% 
%     groupwidth = min(0.8, numbars/(numbars+1.5));
%     
%     for i = 1:numbars
%     % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
%     x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
%     errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
%     end
%     
%     title(cell_type{ct})

    % spatial profile
    subplot(2, 4, 4)
    dis = -10:10;
    plot(dis, Spatial_Profile{ct}(:, 1), 'b', dis, Spatial_Profile{ct}(:, 2), 'r');
    xlim([-10 10])
    title('Spatial Filter')
    h_legend = legend('NDF 4', 'NDF 0', 'location', 'northwest');
%     set(h_legend,'FontSize',14);


    % temporal profile
    subplot(2, 4, [5 6])
    t1 = -0.066*[14:-1:0];
    t2 = -0.033*[14:-1:0];
    plot(t1, Temporal_Filter{ct}(1, :), 'b', t2, Temporal_Filter{ct}(2, :), 'r');
    xlim([-0.75 0])
    title('Temporal Filter')
    
    % Nonlinearity
    subplot(2, 4, 7)
    errorbar(Nonlinearity_mean(:, 1, ct, 1), Nonlinearity_mean(:, 1, ct, 2), Nonlinearity_stev(:, 1, ct), 'b');
    hold on
    errorbar(Nonlinearity_mean(:, 2, ct, 1), Nonlinearity_mean(:, 2, ct, 2), Nonlinearity_stev(:, 2, ct), 'r');
    title('Average Nonlinearity')
    xlim([-2.5 2.5])
    
    % individual nonlinearity
    subplot(2, 4, 8)
    n = size(cell_id{ct}, 1);
    for j = 1:n
        plot(NL{ct}(1, :, j, 1), NL{ct}(2, :, j, 1), 'b', NL{ct}(1, :, j, 2), NL{ct}(2, :, j, 2), 'r');
        xlim([-2.5 2.5]) 
        hold on
    end
    xlim([-2.5 2.5])    
    title('Nonlinearity')
end



%% bar graph
% surround strengths
    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 1400, 750]);

    xtick = cell_type;
    model_series = [surround.mean(1, 1) surround.mean(1, 2); surround.mean(2, 1) surround.mean(2, 2); surround.mean(3, 1) surround.mean(3, 2); surround.mean(4, 1) surround.mean(4, 2); surround.mean(5, 1) surround.mean(5, 2); surround.mean(6, 1) surround.mean(6, 2); surround.mean(7, 1) surround.mean(7, 2)];   
    model_error = [surround.stev(1, 1) surround.stev(1, 2); surround.stev(2, 1) surround.stev(2, 2);surround.stev(3, 1) surround.stev(3, 2);surround.stev(4, 1) surround.stev(4, 2); surround.stev(5, 1) surround.stev(5, 2);surround.stev(6, 1) surround.stev(6, 2);surround.stev(7, 1) surround.stev(7, 2)];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('Surround Center Volume Ratio')
    legend('NDF 4','NDF 0', 'location', 'northwest');
    hold on;
 
    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 

    groupwidth = min(0.8, numbars/(numbars+1.5));
    
    for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end
    
    title('Surround Strengths')

% center size

    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 1400, 750]);

    xtick = cell_type;
    model_series = [center.mean(1, 1) center.mean(1, 2); center.mean(2, 1) center.mean(2, 2); center.mean(3, 1) center.mean(3, 2); center.mean(4, 1) center.mean(4, 2); center.mean(5, 1) center.mean(5, 2); center.mean(6, 1) center.mean(6, 2); center.mean(7, 1) center.mean(7, 2)];   
    model_error = [center.stev(1, 1) center.stev(1, 2); center.stev(2, 1) center.stev(2, 2); center.stev(3, 1) center.stev(3, 2);center.stev(4, 1) center.stev(4, 2);center.stev(5, 1) center.stev(5, 2);center.stev(6, 1) center.stev(6, 2);center.stev(7, 1) center.stev(7, 2)];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('Center Size(a.u.)')
    legend('NDF 4','NDF 0', 'location', 'northwest');
    hold on;
 
    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 

    groupwidth = min(0.8, numbars/(numbars+1.5));
    
    for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end
    
    title('Center Size')

    % DOT
    
    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 1400, 750]);

    xtick = cell_type;
    model_series = [dot.mean(1, 1) dot.mean(1, 2); dot.mean(2, 1) dot.mean(2, 2); dot.mean(3, 1) dot.mean(3, 2); dot.mean(4, 1) dot.mean(4, 2); dot.mean(5, 1) dot.mean(5, 2); dot.mean(6, 1) dot.mean(6, 2); dot.mean(7, 1) dot.mean(7, 2)];   
    model_error = [dot.stev(1, 1) dot.stev(1, 2); dot.stev(2, 1) dot.stev(2, 2); dot.stev(3, 1) dot.stev(3, 2);dot.stev(4, 1) dot.stev(4, 2);dot.stev(5, 1) dot.stev(5, 2);dot.stev(6, 1) dot.stev(6, 2);dot.stev(7, 1) dot.stev(7, 2);];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('DOT')
    legend('NDF 4','NDF 0', 'location', 'northwest');
    hold on;
 
    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 

    groupwidth = min(0.8, numbars/(numbars+1.5));
    
    for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end
    
    ylim([0 1.2])
    title('Degree of Transience')


