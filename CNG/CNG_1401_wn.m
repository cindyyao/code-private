%% load data, calculate fitting and snls
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);

cell_type = {'ON 1', 'ON 2', 'OFF 1', 'OFF 2', 'OFF 3', 'OFF 4'};
datarun{1} = load_data('/Analysis/xyao/2014-01-14-0/data008-map/data008-map', opt);

datarun{1} = compute_sta_fits_sequence(datarun{1}, 'OFF 4');
datarun{1} = load_java_movie(datarun{1}, '/Volumes/lab/acquisition/movie-xml/BW-15-1-0.48-11111-40x40-60.35.xml');
datarun{1} = get_sta_summaries(datarun{1}, 'all');

datarun{1} = get_snls(datarun{1}, 'OFF 4');



datarun{2} = load_data('/Analysis/xyao/2014-01-16-0/data008-map/data008-map', opt);
datarun{2} = compute_sta_fits_sequence(datarun{2}, 'OFF 4');
datarun{2} = load_java_movie(datarun{2}, '/Volumes/lab/acquisition/movie-xml/BW-15-1-0.48-11111-40x40-60.35.xml');
datarun{2} = get_sta_summaries(datarun{2}, 'all');

datarun{2} = get_snls(datarun{2}, 'OFF 4');

%% get cell ids
n = length(cell_type);

cell_id = cell(n, 1);
cell_idx = cell(n, 1);

for i =  1:n
    for j = 1:2
        cell_id{i}{j} = get_cell_ids(datarun{j}, cell_type{i});
        cell_idx{i}{j} = get_cell_indices(datarun{j}, cell_type{i});
    end
end

% get mean sta image

d(1, :) = size(datarun{1}.stas.stas{1});
d(2, :) = size(datarun{2}.stas.stas{1});
f(1) = 1200/d(1, 1);
f(2) = 1200/d(2, 1);

sta_f_weighted = cell(n, 2);
sta_f_n = zeros(n, 2);

for cd = 1:2
for k = 1:n
    index = cell_idx{k}{cd};
    cell_numb = length(index);
    sta_f{cd} = zeros(f(cd)*d(cd, 1)/2, f(cd)*d(cd, 2)/2, d(cd, 3), d(cd, 4));
    if isempty(index) == 0
       for i = 1:cell_numb
           sta = datarun{cd}.stas.stas{index(i)};
           if cell_type{k}(2) == 'N' 
              [~, b] = max(sta(:));
           else
              [~, b] = min(sta(:));
           end
           z = ceil(b/(d(cd, 1)*d(cd, 2)));
           sta_1 = sta(:, :, 1, z);
           noise = robust_std(sta_1(:), 3);
           sta_temp = zeros(f(cd)*d(cd, 1), f(cd)*d(cd, 2), d(cd, 3), d(cd, 4));
           for j = 1:d(cd, 4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f(cd), f(cd)));
           end
           xy = floor(f(cd)*rf_com(sta));
           if sum(xy<=f(cd)*d(cd, 1)/4) == 0 && sum(xy>=3*f(cd)*d(cd, 1)/4) == 0 && isempty(xy) == 0
           sta_f{cd} = sta_f{cd} + sta_temp(xy(2)-f(cd)*d(cd, 1)/4+1:xy(2)+f(cd)*d(cd, 1)/4, xy(1)-f(cd)*d(cd, 2)/4+1:xy(1)+f(cd)*d(cd, 2)/4, :, :)/noise;
           sta_f_n(k, cd) = sta_f_n(k, cd) + 1;
           end
           
       end       
    end
    sta_f_weighted{k, cd} = sta_f{cd};
end
end

sta_c = cell(n, 2);
for i = 1:n
    for cd = 1:2
        sta = sta_f_weighted{i, cd};
        sta_c_temp = zeros(d(cd, 1), d(cd, 2), d(cd, 3), d(cd, 4));
        for j = 1:d(cd, 1)
            for k = 1:d(cd, 2)
                for m = 1:d(cd, 4)
                    sta_temp = sta(f(cd)/2*(j-1)+1:f(cd)/2*j, f(cd)/2*(k-1)+1:f(cd)/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c{i, cd} = sta_c_temp;
    end
end


%% get mean spatial profile


for ct = 1:n;
    for i = 1:2
        cell_num = length(cell_id{ct}{i});
        temp_sf{i} = zeros(2, 40, cell_num);
        temp_tf{i} = zeros(30, cell_num);
        ns = [];

        for cn = 1:cell_num
            idx = cell_idx{ct}{i}(cn);
            com = datarun{i}.stas.rf_coms{idx};
            if isempty(com) == 0
       
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
            RF1_std{i} = RF1{i}/std(RF1{i});
            
            tc_std(:, i) = tc(:, i)/std(tc(:, i));

            temp_sf{i}(1, :, cn) = Dis1{i}; 
            temp_sf{i}(2, :, cn) = RF1_std{i};
        
            temp_tf{i}(:, cn) = tc_std(:, i);
            else
                ns = [ns cn];
            end
        end
        
        temp_sf{i}(:, :, ns) = [];
        temp_tf{i}(:, ns) = [];
        cell_id{ct}{i}(ns) = [];
        cell_idx{ct}{i}(ns) = [];

        spatial_profile{ct}{i} = temp_sf{i};
        temporal_filter{ct}{i} = temp_tf{i};

        
    end
    
end
    



for ct = 1:n    
    for cd = 1:2
    Temporal_Filter{ct}(:, cd) = mean(temporal_filter{ct}{cd}, 2);

    cell_num = length(cell_id{ct}{cd});
    Spatial_Profile_temp{cd} = zeros(21, cell_num);
    for cn = 1:cell_num
        Spatial_Profile_temp{cd}(:, cn) = interp1(spatial_profile{ct}{cd}(1, :, cn), spatial_profile{ct}{cd}(2, :, cn), -10:10);
    end
    Spatial_Profile{ct}(:, cd) = squeeze(mean(Spatial_Profile_temp{cd}, 2));

    end
end

% save('snl130221_avg.mat', 'Spatial_Profile', 'Temporal_Filter', 'sta_c')

%% get average snls

display_rate = 60.35;
refresh_rate = [1 1];

NL = cell(n, 1);

for ct = 1:n
    for cd = 1:2
    cell_num = length(cell_id{ct}{cd});
    NL_temp{cd} = zeros(2, 20, cell_num);    
        for cc = 1:cell_num
            idx = cell_idx{ct}{cd}(cc);
            gen_signals = datarun{cd}.stas.snls{idx}.gen_signal;
            gen_signals = gen_signals/std(gen_signals);
            spikes = datarun{cd}.stas.snls{idx}.spikes;

            [X, Y] = curve_from_binning(gen_signals, spikes, ...
                'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
            Y = Y*display_rate/refresh_rate(cd);
            NL_temp{cd}(1, :, cc) = X;
            NL_temp{cd}(2, :, cc) = Y/max(Y);
        end
        
    end
    NL{ct} = NL_temp;
end

XX = zeros(20, 2, n);
for ct = 1:n
    for cd = 1:2
        xx = squeeze(NL{ct}{cd}(1, :, :));
        XX_temp = mean(xx, 2);
        XX_temp(1) = max(xx(1, :));
        XX_temp(20) = min(xx(20, :));
        XX(:, cd, ct) = XX_temp;
    end
end

NL_Y_mean = zeros(20, 2, n);
Nonlinearity_stev = zeros(20, 2, n);

for ct = 1:n
    for cd = 1:2
    cell_num = length(cell_id{ct}{cd});
    NL_Y_temp{cd} = zeros(20, cell_num);
        for cc = 1:cell_num
            YY_temp = interp1(NL{ct}{cd}(1, :, cc), NL{ct}{cd}(2, :, cc), XX(:, cd, ct));
            NL_Y_temp{cd}(:, cc) = YY_temp;
        end
    NL_Y_mean_temp(cd, :) = squeeze(mean(NL_Y_temp{cd}, 2));
    Nonlinearity_stev_temp(cd, :) = squeeze(std(NL_Y_temp{cd}, 0, 2))/sqrt(cell_num);

    end
    NL_Y_mean(:, :, ct) = NL_Y_mean_temp';
    Nonlinearity_stev(:, :, ct) = Nonlinearity_stev_temp';
end

Nonlinearity_mean(:, :, :, 1) = XX;
Nonlinearity_mean(:, :, :, 2) = NL_Y_mean;

%% surround strengths & center size


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
    for i = 1:2
        idx = cell_idx{ct}{i};
        cell_num = length(idx);
        for cc = 1:cell_num
        r = datarun{i}.stimulus.field_height;
        params_temp = zeros(7, 1);
        params_temp(1) = datarun{i}.matlab.sta_fits{cell_idx{ct}{i}(cc)}.center_point_x*40/r;
        params_temp(2) = datarun{i}.matlab.sta_fits{cell_idx{ct}{i}(cc)}.center_point_y*40/r;
        params_temp(3) = datarun{i}.matlab.sta_fits{cell_idx{ct}{i}(cc)}.center_sd_x*40/r;
        params_temp(4) = datarun{i}.matlab.sta_fits{cell_idx{ct}{i}(cc)}.center_sd_y*40/r;
        params_temp(5) = datarun{i}.matlab.sta_fits{cell_idx{ct}{i}(cc)}.center_rotation_angle;
        params_temp(6) = datarun{i}.matlab.sta_fits{cell_idx{ct}{i}(cc)}.surround_sd_scale;
        params_temp(7) = datarun{i}.matlab.sta_fits{cell_idx{ct}{i}(cc)}.surround_amp_scale;
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
    
    [~, surround_strengths_p(ct)] = ttest2(surround_strengths{ct}{1}, surround_strengths{ct}{2});
    [~, center_size_p(ct)] = ttest2(center_size_temp{ct}{1}, center_size_temp{ct}{2});

end


surround = struct('mean', surround_strengths_mean, 'variance', surround_strengths_var, ...
    'stev', surround_strengths_stev, 'p_value', surround_strengths_p);
center = struct('mean', center_size_mean, 'variance', center_size_var, 'stev', ...
    center_size_stev, 'p_value', center_size_p);

%% degree of transience(DOT)

smooth_f = 0.01;
DOT = cell(n, 1);
t = -29:0;
for ct = 1:n
    for cd = 1:2
    cell_num = length(cell_id{ct}{cd});
    DOT_temp{cd} = zeros(cell_num, 1);    
        for cn = 1:cell_num
            tf_smooth = interp1(-29:0, temporal_filter{ct}{cd}(:, cn), -14:smooth_f:0);
            area = smooth_f * trapz(tf_smooth);
            tf_smooth_p = tf_smooth;
            tf_smooth_p(tf_smooth_p<0) = 0;
            area_p = smooth_f * trapz(tf_smooth_p);
            area_n = area_p - area;
            DOT_temp{cd}(cn) = 1 - abs(area/(area_p + area_n));
        end
    end
    DOT{ct} = DOT_temp;
end

DOT_mean = zeros(n, 2);
DOT_stev = zeros(n, 2);
DOT_pvalue = zeros(n, 1);

for ct = 1:n
    for cd = 1:2
    cell_num = length(cell_id{ct}{cd});
    DOT_mean(ct, cd) = mean(DOT{ct}{cd});
    DOT_stev(ct, cd) = std(DOT{ct}{cd})/sqrt(cell_num);
    end
    
    [~, DOT_pvalue(ct)] = ttest2(DOT{ct}{1}, DOT{ct}{2});
end

dot = struct('mean', DOT_mean, 'stev', DOT_stev, 'p_value', DOT_pvalue);


%% plot cell type summaries

% load('snls130221_avg.mat')

for ct = 1:n
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
    title('control')
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
    title('rescue')
    axis off
    
    subplot(2, 4, 3)
    plot_rf_summaries(datarun{1}, cell_id{ct}{1})
    title('control')
    
    subplot(2, 4, 4)
    plot_rf_summaries(datarun{2}, cell_id{ct}{2})
    title('rescue')
    

    % spatial profile
    subplot(2, 4, 5)
    dis = -10:10;
    plot(dis, Spatial_Profile{ct}(:, 1), 'b', dis, Spatial_Profile{ct}(:, 2), 'r');
    xlim([-10 10])
    title('Spatial Filter')
%     set(h_legend,'FontSize',14);


    % temporal profile
    subplot(2, 4, 6)
    t1 = -0.0166*[29:-1:0];
    t2 = -0.0166*[29:-1:0];
    plot(t1, Temporal_Filter{ct}(:, 1), 'b', t2, Temporal_Filter{ct}(:, 2), 'r');
    xlim([-0.5 0])
    title('Temporal Filter')
    
    % Nonlinearity
    subplot(2, 4, 8)
    errorbar(Nonlinearity_mean(:, 1, ct, 1), Nonlinearity_mean(:, 1, ct, 2), Nonlinearity_stev(:, 1, ct), 'b');
    hold on
    errorbar(Nonlinearity_mean(:, 2, ct, 1), Nonlinearity_mean(:, 2, ct, 2), Nonlinearity_stev(:, 2, ct), 'r');
    title('Average Nonlinearity')
    xlim([-2.5 2.5])
    h_legend = legend('control', 'rescue', 'location', 'northwest');

    % individual nonlinearity
    subplot(2, 4, 7)
    cl = 'br';
    for cd = 1:2
    cell_num = length(cell_id{ct}{cd});
    for j = 1:cell_num
        plot(NL{ct}{cd}(1, :, j), NL{ct}{cd}(2, :, j), cl(cd));
        xlim([-2.5 2.5])
        hold on
    end
    end
    xlim([-2.5 2.5])
    title('Nonlinearity')
end



%% bar graph
% surround strengths
    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 1000, 500]);

    xtick = cell_type;
    model_series = [surround.mean(1, 1) surround.mean(1, 2); surround.mean(2, 1) surround.mean(2, 2); surround.mean(3, 1) surround.mean(3, 2); surround.mean(4, 1) surround.mean(4, 2); surround.mean(5, 1) surround.mean(5, 2); surround.mean(6, 1) surround.mean(6, 2)];   
    model_error = [surround.stev(1, 1) surround.stev(1, 2); surround.stev(2, 1) surround.stev(2, 2);surround.stev(3, 1) surround.stev(3, 2);surround.stev(4, 1) surround.stev(4, 2); surround.stev(5, 1) surround.stev(5, 2); surround.stev(6, 1) surround.stev(6, 2)];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('Surround Center Volume Ratio')
    legend('control','rescue', 'location', 'northwest');
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
    set(FigHandle, 'Position', [100, 100, 1000, 500]);

    xtick = cell_type;
    model_series = [center.mean(1, 1) center.mean(1, 2); center.mean(2, 1) center.mean(2, 2); center.mean(3, 1) center.mean(3, 2); center.mean(4, 1) center.mean(4, 2); center.mean(5, 1) center.mean(5, 2); center.mean(6, 1) center.mean(6, 2)];   
    model_error = [center.stev(1, 1) center.stev(1, 2); center.stev(2, 1) center.stev(2, 2); center.stev(3, 1) center.stev(3, 2);center.stev(4, 1) center.stev(4, 2); center.stev(5, 1) center.stev(5, 2); center.stev(6, 1) center.stev(6, 2)];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('Center Size(a.u.)')
    legend('control','rescue', 'location', 'northwest');
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
    set(FigHandle, 'Position', [100, 100, 1000, 500]);

    xtick = cell_type;
    model_series = [dot.mean(1, 1) dot.mean(1, 2); dot.mean(2, 1) dot.mean(2, 2); dot.mean(3, 1) dot.mean(3, 2); dot.mean(4, 1) dot.mean(4, 2); dot.mean(5, 1) dot.mean(5, 2); dot.mean(6, 1) dot.mean(6, 2)];   
    model_error = [dot.stev(1, 1) dot.stev(1, 2); dot.stev(2, 1) dot.stev(2, 2); dot.stev(3, 1) dot.stev(3, 2);dot.stev(4, 1) dot.stev(4, 2); dot.stev(5, 1) dot.stev(5, 2); dot.stev(6, 1) dot.stev(6, 2)];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('DOT')
    legend('control','rescue', 'location', 'northwest');
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



%% Plot Spatial

figure
for i = 1:n
    for j = 1:2
        subplot(6, 6, 6*(i-1)+j)
        plot_rf_summaries(datarun{j}, cell_id{i}{j})
        axis([10 30 10 30])
        set(gca, 'XTick', [], 'YTick', [])
    end
end
subplot(6, 6, [3:6 9:12 15:18])
xtick = cell_type;
    model_series = [center.mean(1, 1) center.mean(1, 2); center.mean(2, 1) center.mean(2, 2); center.mean(3, 1) center.mean(3, 2); center.mean(4, 1) center.mean(4, 2); center.mean(5, 1) center.mean(5, 2); center.mean(6, 1) center.mean(6, 2)];   
    model_error = [center.stev(1, 1) center.stev(1, 2); center.stev(2, 1) center.stev(2, 2); center.stev(3, 1) center.stev(3, 2);center.stev(4, 1) center.stev(4, 2); center.stev(5, 1) center.stev(5, 2); center.stev(6, 1) center.stev(6, 2)];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',[])
    ylabel('Center Size(a.u.)')
    legend('control','rescue');
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
    
    subplot(6, 6, [21:24 27:30 33:36])
    xtick = cell_type;
    model_series = [surround.mean(1, 1) surround.mean(1, 2); surround.mean(2, 1) surround.mean(2, 2); surround.mean(3, 1) surround.mean(3, 2); surround.mean(4, 1) surround.mean(4, 2); surround.mean(5, 1) surround.mean(5, 2); surround.mean(6, 1) surround.mean(6, 2)];   
    model_error = [surround.stev(1, 1) surround.stev(1, 2); surround.stev(2, 1) surround.stev(2, 2);surround.stev(3, 1) surround.stev(3, 2);surround.stev(4, 1) surround.stev(4, 2); surround.stev(5, 1) surround.stev(5, 2); surround.stev(6, 1) surround.stev(6, 2)];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('Surround Center Volume Ratio')
    legend('control','rescue');
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

    
%% time to peak
for ct = 1:n
    for cd = 1:2
        if cell_type{ct}(2) == 'N'
            [~, I] = max(temporal_filter{ct}{cd});
        else
            [~, I] = min(temporal_filter{ct}{cd});
        end
        TTP{ct}{cd} = (30 - I)*0.0166;
        TTP_mean(ct, cd) = mean(TTP{ct}{cd});
        TTP_ste(ct, cd) = std(TTP{ct}{cd})/sqrt(length(TTP{ct}{cd}));

    end
    [~, TTP_p(ct)] = ttest2(TTP{ct}{1}, TTP{ct}{2});

end

time_to_peak = struct('mean', TTP_mean, 'stev', TTP_ste, 'p_value', TTP_p);

%% Plot Temporal

figure
for i = 1:n
        subplot(6, 6, 6*(i-1)+1)
        plot(t1, Temporal_Filter{i}(:, 1), 'b', t2, Temporal_Filter{i}(:, 2), 'r');
        xlim([-0.5 0])
        set(gca, 'XTick', [], 'YTick', [])
end
subplot(6, 6, [3:6 9:12 15:18])
xtick = cell_type;
    model_series = [time_to_peak.mean(1, 1) time_to_peak.mean(1, 2); time_to_peak.mean(2, 1) ...
        time_to_peak.mean(2, 2); time_to_peak.mean(3, 1) time_to_peak.mean(3, 2); ... 
        time_to_peak.mean(4, 1) time_to_peak.mean(4, 2); time_to_peak.mean(5, 1) ... 
        time_to_peak.mean(5, 2); time_to_peak.mean(6, 1) time_to_peak.mean(6, 2)];   
    model_error = [time_to_peak.stev(1, 1) time_to_peak.stev(1, 2); time_to_peak.stev(2, 1) ... 
        time_to_peak.stev(2, 2); time_to_peak.stev(3, 1) time_to_peak.stev(3, 2); ...
        time_to_peak.stev(4, 1) time_to_peak.stev(4, 2); time_to_peak.stev(5, 1) ... 
        time_to_peak.stev(5, 2); time_to_peak.stev(6, 1) time_to_peak.stev(6, 2)];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',[])
    ylabel('Time to Peak(Sec)')
    legend('control','rescue');
    hold on;
 
    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 

    groupwidth = min(0.8, numbars/(numbars+1.5));
    
    for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end
    
    title('Time to Peak')
    
    subplot(6, 6, [21:24 27:30 33:36])
    xtick = cell_type;
    model_series = [dot.mean(1, 1) dot.mean(1, 2); dot.mean(2, 1) dot.mean(2, 2); dot.mean(3, 1) dot.mean(3, 2); dot.mean(4, 1) dot.mean(4, 2); dot.mean(5, 1) dot.mean(5, 2); dot.mean(6, 1) dot.mean(6, 2)];   
    model_error = [dot.stev(1, 1) dot.stev(1, 2); dot.stev(2, 1) dot.stev(2, 2); dot.stev(3, 1) dot.stev(3, 2);dot.stev(4, 1) dot.stev(4, 2); dot.stev(5, 1) dot.stev(5, 2); dot.stev(6, 1) dot.stev(6, 2)];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('DOT')
    legend('control','rescue', 'location', 'northwest');
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

    %% plot nonlinearity
    figure
    for ct = 1:n
        subplot(2, 3, ct)
        errorbar(Nonlinearity_mean(:, 1, ct, 1), Nonlinearity_mean(:, 1, ct, 2), Nonlinearity_stev(:, 1, ct), 'b');
        hold on
        errorbar(Nonlinearity_mean(:, 2, ct, 1), Nonlinearity_mean(:, 2, ct, 2), Nonlinearity_stev(:, 2, ct), 'r');
        title(cell_type{ct})
        xlim([-2.5 2.5])
        if ct == 1
            h_legend = legend('control', 'rescue', 'location', 'northwest');
        end
    end