% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun2 = load_data('/Analysis/xyao/2013-04-24-0/data002/data002', opt);
datarun3 = load_data('/Analysis/xyao/2013-04-24-0/data003/data003', opt);

% map ei
[cell_list_map32, ~] = map_ei(datarun3, datarun2);


% map data002 and data003
cell_type = {'ON transient', 'OFF brisk transient', 'OFF transient', 'OFF sustained', 'OFF unknown'};
n = length(cell_type);
cell_id23 = cell(n, 1);
cell_idx23 = cell(n, 1);
cell_id_matched23 = [];
cell_idx_matched23 = [];

for i = 1:n
    [cell_idx3, ~, ~] = get_cell_indices(datarun3, cell_type{i});
    cell_id2 = cell_list_map32(cell_idx3);
    ept = 1-cellfun(@isempty, cell_id2);
    cell_idx3 = cell_idx3.*ept;
    cell_idx3(cell_idx3 == 0) = [];
    cell_id2 = cell2mat(cell_id2);
    cell_idx2 = get_cell_indices(datarun2, cell_id2);
    cell_id3 = datarun3.cell_ids(cell_idx3);
    id = [cell_id2' cell_id3'];
    idx = [cell_idx2' cell_idx3'];
    cell_id23{i} = id;
    cell_idx23{i} = idx;
    cell_id_matched23 = [cell_id_matched23; id];
    cell_idx_matched23 = [cell_idx_matched23; idx];

end



% get parameters
datarun2 = compute_sta_fits_sequence(datarun2, cell_id_matched23(:, 1), 'verbose', true);
datarun3 = compute_sta_fits_sequence(datarun3, cell_id_matched23(:, 2), 'verbose', true);

save_sta_fits(datarun2);
save_sta_fits(datarun3);

datarun2 = load_sta_fits(datarun2);
datarun3 = load_sta_fits(datarun3);

params23 = cell(n, 1);
for i = 1:n
    index = cell_idx23{i};
    cell_numb = size(index, 1);
    params_temp = zeros(cell_numb, 8, 2);
    for j = 1:cell_numb
        params2 = datarun2.matlab.sta_fits{index(j, 1)};
        params3 = datarun3.matlab.sta_fits{index(j, 2)};
        if isempty(params2) == 1 || isempty(params3) == 1;
           params_temp(j, :, 1) = zeros(1, 8);
           params_temp(j, :, 2) = zeros(1, 8);

        else  
        params_temp(j, 1, 1) = params2.center_point_x;
        params_temp(j, 2, 1) = params2.center_point_y;
        params_temp(j, 3, 1) = params2.center_sd_x;
        params_temp(j, 4, 1) = params2.center_sd_y;
        params_temp(j, 5, 1) = params2.center_rotation_angle;
        params_temp(j, 6, 1) = params2.surround_sd_scale*params2.center_sd_x;
        params_temp(j, 7, 1) = params2.surround_sd_scale*params2.center_sd_y;
        params_temp(j, 8, 1) = params2.surround_amp_scale;
        
        params_temp(j, 1, 2) = params3.center_point_x;
        params_temp(j, 2, 2) = params3.center_point_y;
        params_temp(j, 3, 2) = params3.center_sd_x;
        params_temp(j, 4, 2) = params3.center_sd_y;
        params_temp(j, 5, 2) = params3.center_rotation_angle;
        params_temp(j, 6, 2) = params3.surround_sd_scale*params3.center_sd_x;
        params_temp(j, 7, 2) = params3.surround_sd_scale*params3.center_sd_y;
        params_temp(j, 8, 2) = params3.surround_amp_scale;
        end
        
    end
    params23{i} = params_temp;
end

%% plot mosaic

cell_type_n = 5;
params_temp = params23{cell_type_n};
cell_numb = size(params_temp, 1);
figure;
subplot(1, 2, 1);
for i = 1:cell_numb
    [X1, Y1] = drawEllipse(params_temp(i, 1:5, 1));
    [X2, Y2] = drawEllipse(params_temp(i, 1:5, 2));
    if i == 1
        h1 = plot(X1, Y1, 'r');
        hold on
        h2 = plot(X2, Y2, 'b');
    else
        plot(X1, Y1, 'r');
        plot(X2, Y2, 'b');
    end
end

axis([0 40 0 40]);
legend([h1 h2], 'NDF 4', 'NDF 0');
title([cell_type{cell_type_n} ' center']);

subplot(1, 2, 2);
for i = 1:cell_numb
    [X1, Y1] = drawEllipse(params_temp(i, [1, 2, 6, 7, 5], 1));
    [X2, Y2] = drawEllipse(params_temp(i, [1, 2, 6, 7, 5], 2));
    if i == 1
        h1 = plot(X1, Y1, 'r');
        hold on
        h2 = plot(X2, Y2, 'b');
    else
        plot(X1, Y1, 'r');
        plot(X2, Y2, 'b');
    end
end

axis([0 40 0 40]);
legend([h1 h2], 'NDF 4', 'NDF 0');
title([cell_type{cell_type_n} ' surround']);


%% weighted addition of RF by SNR

d2 = size(datarun2.stas.stas{1});
d3 = size(datarun3.stas.stas{1});
f2 = 30;
f3 = 30;


sta_f_weighted = cell(n, 2);
sta_f_n = zeros(n, 2);
for k = 1:n;
    index = cell_idx23{k};
    cell_numb = size(cell_id23{k}, 1);
    sta_f2 = zeros(f2*d2(1)/2, f2*d2(1)/2, 1, d2(4));
    sta_f3 = zeros(f3*d3(1)/2, f3*d3(1)/2, 1, d3(4));
    sta_f_n2 = 0;
    sta_f_n3 = 0;
    if isempty(index) == 0
       for i = 1:cell_numb
           sta = datarun2.stas.stas{index(i, 1)};
           if k == 1 
              [~, b] = max(sta(:));
           else
              [~, b] = min(sta(:));
           end
           z = ceil(b/(d2(1)*d2(2)));
           sta_1 = sta(:, :, 1, z);
           noise = robust_std(sta_1(:), 3);
           sta_temp = zeros(f2*d2(1), f2*d2(2), d2(3), d2(4));
           for j = 1:d2(4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f2, f2));
           end
           xy = floor(f2*rf_com(sta));
           if sum(xy<=f2*d2(1)/8) == 0 && sum(xy>=7*f2*d2(1)/8) == 0 && isempty(xy) == 0
               x_r = min(f2*d2(2), xy(1)+f2*d2(2)/4);
               x_l = max(xy(1)-f2*d2(2)/4+1, 0);
               y_d = min(f2*d2(1), xy(2)+f2*d2(1)/4);
               y_u = max(xy(2)-f2*d2(1)/4+1, 0);
               sta_f2_temp = zeros(f2*d2(1)/2, f2*d2(2)/2, d2(3), d2(4));
               x_l_temp = min(xy(1)-f2*d2(2)/4+1, 1);
               x_r_temp = x_l_temp+x_r-x_l;
               y_u_temp = min(xy(2)-f2*d2(1)/4+1, 1);
               y_d_temp = y_u_temp+y_d-y_u;
               sta_f2_temp(y_u_temp:y_d_temp, x_l_temp:x_r_temp, :, :) = sta_temp(y_u:y_d, x_l:x_r, :, :); 
               sta_f2 = sta_f2 + sta_f2_temp/noise;
               sta_f_n2 = sta_f_n2 + 1;
           end
           
    
           sta = datarun3.stas.stas{index(i, 2)};
           if k == 1 
              [~, b] = max(sta(:));
           else
              [~, b] = min(sta(:));
           end
           z = ceil(b/(d3(1)*d3(2)));
           sta_1 = sta(:, :, 1, z);
           noise = robust_std(sta_1(:), 3);
           sta_temp = zeros(f3*d3(1), f3*d3(2), d3(3), d3(4));
           for j = 1:d3(4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f3, f3));
           end
           xy = floor(f3*rf_com(sta));
           if sum(xy<=f3*d3(1)/8) == 0 && sum(xy>=7*f3*d3(1)/8) == 0 && isempty(xy) == 0
               x_r = min(f3*d3(2), xy(1)+f2*d3(2)/4);
               x_l = max(xy(1)-f3*d3(2)/4+1, 0);
               y_d = min(f3*d3(1), xy(2)+f3*d3(1)/4);
               y_u = max(xy(2)-f3*d3(1)/4+1, 0);
               sta_f3_temp = zeros(f3*d3(1)/2, f3*d3(2)/2, d3(3), d3(4));
               x_l_temp = min(xy(1)-f3*d3(2)/4+1, 1);
               x_r_temp = x_l_temp+x_r-x_l;
               y_u_temp = min(xy(2)-f3*d3(1)/4+1, 1);
               y_d_temp = y_u_temp+y_d-y_u;
               sta_f3_temp(y_u_temp:y_d_temp, x_l_temp:x_r_temp, :, :) = sta_temp(y_u:y_d, x_l:x_r, :, :); 
               sta_f3 = sta_f3 + sta_f3_temp/noise;
               sta_f_n3 = sta_f_n3 + 1;
           end
           
       end       
    end
    sta_f_weighted{k, 1} = sta_f2;
    sta_f_weighted{k, 2} = sta_f3;
    sta_f_n(k, 1) = sta_f_n2;
    sta_f_n(k, 2) = sta_f_n3;
    
end


%% mean RF for different sampling


sta_c = cell(n, 2);
for i = 1:n
        sta = sta_f_weighted{i, 1};
        sta_c_temp = zeros(d2(1), d2(2), d2(3), d2(4));
        for j = 1:d2(1)
            for k = 1:d2(2)
                for m = 1:d2(4)
                    sta_temp = sta(f2/2*(j-1)+1:f2/2*j, f2/2*(k-1)+1:f2/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c{i, 1} = sta_c_temp;
    
        sta = sta_f_weighted{i, 2};
        sta_c_temp = zeros(d3(1), d3(2), d3(3), d3(4));
        for j = 1:d3(1)
            for k = 1:d3(2)
                for m = 1:d3(4)
                    sta_temp = sta(f3/2*(j-1)+1:f3/2*j, f3/2*(k-1)+1:f3/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c{i, 2} = sta_c_temp;
end

figure
sta = sta_c{3, 1};
for i = 1:15
    im = norm_image(sta(:, :, 1, i));
    image(im);
    pause
end

% sta15 = sta_c_15{4, 1};
% [~, b] = max(sta15(:));
% 
% z = ceil(b/225);
% sta15 = sta15(:, :, 1, z);
% 
% sta30 = sta_c_30{4, 1};
% [~, b] = max(sta30(:));
% 
% z = ceil(b/900);
% sta30 = sta30(:, :, 1, z);
% 
% figure
% hist(sta15(:), 20)
% xlabel('stixel value')
% ylabel('number of stixels')
% title('15*15')
% figure
% hist(sta30(:), 20)
% xlabel('stixel value')
% ylabel('number of stixels')
% title('30*30')
% 

%% fit sta

fit_sta_c = cell(n, 2);
params_c = cell(n, 2);

for i = 1:n
        
    fit_temp = fit_sta_sequence(sta_c{i, 1});
    fit_sta_c{i, 1} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x;
    params_temp(2) = fit_temp.center_point_y;
    params_temp(3) = fit_temp.center_sd_x;
    params_temp(4) = fit_temp.center_sd_y;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_c{i, 1} = params_temp;


    fit_temp = fit_sta_sequence(sta_c{i, 2});
    fit_sta_c{i, 2} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x;
    params_temp(2) = fit_temp.center_point_y;
    params_temp(3) = fit_temp.center_sd_x;
    params_temp(4) = fit_temp.center_sd_y;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_c{i, 2} = params_temp;
    
    i

end






%% plot mean RF versus fitting

cell_type_n = 1;
light_level = 2;

sta = sta_c{cell_type_n, light_level};

if light_level == 1
    d = d2(1);
else
    d = d3(1);
end

if cell_type_n == 1 
    [~, b] = max(sta(:));
else
    [~, b] = min(sta(:));
end

z = ceil(b/d^2);
sta_1 = sta(:, :, 1, z);


figure
subplot(1, 2, 1);
im = norm_image(sta_1);
image(im)
hold on
if light_level == 1
    clear params_temp
    params_temp = params_c{cell_type_n, light_level}(1:7);
    params_temp(5) = params_c{cell_type_n, light_level}(5);
else
    params_temp = params_c{cell_type_n, light_level};
end

[X1, Y1] = drawEllipse(params_temp(1:5));
[X2, Y2] = drawEllipse(params_temp([1, 2, 6, 7, 5]));
plot(X1, Y1, 'r')
plot(X2, Y2, 'b')
legend('center', 'surround')
NDF = [4 0.3];
title([cell_type{cell_type_n} ' NDF = ' num2str(NDF(light_level))])




if cell_type_n == 1 
    com = rf_com(sta_1);
else
    com = rf_com(-sta_1);
end


pix = [1:d] + 0.5;
[X Y] = meshgrid(pix);
xdis = com(1) - X;
ydis = com(2) - Y;

dis = zeros(d, d);
for i = 1:d
    for j = 1:d
        dis(i, j) = norm([xdis(i, j) ydis(i, j)]);
    end
end

dis1 = reshape(dis, 1, d^2);
sta_temp = reshape(sta_1, 1, d^2);

para_temp = params_c{cell_type_n, light_level};
if light_level == 1
    fit_para1 = struct('center_point_x', para_temp(1), 'center_point_y', para_temp(2), ...
    'sd_x', para_temp(3), 'sd_y', para_temp(4), 'rotation_angle', para_temp(5), ...
    'x_dim', d, 'y_dim', d);

    fit_para2 = struct('center_point_x', para_temp(1), 'center_point_y', para_temp(2), ...
    'sd_x', para_temp(6), 'sd_y', para_temp(7), 'amp_scale', para_temp(8), ...
    'rotation_angle', para_temp(5), 'x_dim', d, 'y_dim', d);

else
    fit_para1 = struct('center_point_x', para_temp(1), 'center_point_y', para_temp(2), ...
    'sd_x', para_temp(3), 'sd_y', para_temp(4), 'rotation_angle', para_temp(5), ...
    'x_dim', d, 'y_dim', d);

    fit_para2 = struct('center_point_x', para_temp(1), 'center_point_y', para_temp(2), ...
    'sd_x', para_temp(6), 'sd_y', para_temp(7), 'amp_scale', para_temp(8), ...
    'rotation_angle', para_temp(5), 'x_dim', d, 'y_dim', d);
end

fit1 = make_Gaussian_two_d(fit_para1);
fit2 = make_Gaussian_two_d(fit_para2);
if cell_type_n == 1 
    fit = fit1 - fit2;
    comf = rf_com(fit);
    r = max(sta_temp(:))/max(fit(:));
    fit_center = fit;
    fit_center(fit_center < 0) = 0;
    center = sum(fit_center(:));
    fit_surround = fit;
    fit_surround(fit_surround > 0) = 0;
    surround = sum(fit_surround(:));
    ratio = -center/surround
else
    fit = fit2 - fit1;
    comf = rf_com(-fit);
    r = min(sta_temp(:))/min(fit(:));
    fit_center = fit;
    fit_center(fit_center > 0) = 0;
    center = sum(fit_center(:));
    fit_surround = fit;
    fit_surround(fit_surround < 0) = 0;
    surround = sum(fit_surround(:));
    ratio = -center/surround
end



pix = [1:d] + 0.5;
[X Y] = meshgrid(pix);
xdisf = comf(1) - X;
ydisf = comf(2) - Y;

disf = zeros(d, d);
for i = 1:d
    for j = 1:d
        disf(i, j) = norm([xdisf(i, j) ydisf(i, j)]);
    end
end

dis1f = reshape(disf, 1, d^2);
fit1 = reshape(fit, 1, d^2);

subplot(1, 2, 2);
if cell_type_n == 1 
    plot(dis1, sta_temp, '.', 'MarkerSize', 7)
    hold on 
    plot(dis1f, fit1*r, '.', 'MarkerSize', 7, 'Color', 'r')
    plot(dis1, zeros(1, d^2), 'm')
else
    plot(dis1, -sta_temp, '.', 'MarkerSize', 7)
    hold on
    plot(dis1f, -fit1*r, '.', 'MarkerSize', 7, 'Color', 'r')
    plot(dis1, zeros(1, d^2), 'm')
end
xlim([0 25])
legend('data', 'fit')
title([cell_type{cell_type_n} '   NDF ' num2str(NDF(light_level))]);

%% real data of mean RF

cell_type_n = 3;
distance = cell(2, 1);
STA = cell(2, 1);

for k = 1:2
    
sta = sta_c{cell_type_n, k};

if k == 1
    d = d2(1);
else
    d = d3(1);
end

if cell_type_n == 1 
    [~, b] = max(sta(:));
else
    [~, b] = min(sta(:));
end

z = ceil(b/d^2);
sta_1 = sta(:, :, 1, z);


if cell_type_n == 1 
    com = rf_com(sta_1);
else
    com = rf_com(-sta_1);
end


pix = [1:d] + 0.5;
[X Y] = meshgrid(pix);
xdis = com(1) - X;
ydis = com(2) - Y;

dis = zeros(d, d);
for i = 1:d
    for j = 1:d
        dis(i, j) = norm([xdis(i, j) ydis(i, j)]);
    end
end

dis1 = reshape(dis, 1, d^2);
sta_temp = reshape(sta_1, 1, d^2);


if k == 1
    figure
    if cell_type_n == 1 
        plot(dis1, sta_temp, '.', 'MarkerSize', 7)
        hold on 
    else
        plot(dis1, -sta_temp, '.', 'MarkerSize', 7)
        hold on
    end
else
    if cell_type_n == 1 
        plot(dis1, sta_temp, '.', 'MarkerSize', 7, 'color', 'r')
        plot(dis1, zeros(1, d^2), 'm')
    else
        plot(dis1, -sta_temp, '.', 'MarkerSize', 7, 'color', 'r')
        plot(dis1, zeros(1, d^2), 'm')
    end
end


end
xlim([0 25])
legend('NDF 4', 'NDF 0.3')
title([cell_type{cell_type_n} '   NDF ' num2str(NDF(light_level))]);

