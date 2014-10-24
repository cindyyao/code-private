% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun0 = load_data('/Analysis/xyao/2013-02-14-0/data000/data000', opt);
datarun6 = load_data('/Analysis/xyao/2013-02-14-0/data006/data006', opt);

% map ei
[cell_list_map60, ~] = map_ei(datarun6, datarun0);

% map data000 and data006
cell_type = {'ON transient', 'ON brisk transient', 'OFF brisk transient', 'OFF transient', 'OFF sustained'};
n = length(cell_type);
cell_id06 = cell(n, 1);
cell_idx06 = cell(n, 1);
cell_id_matched06 = [];
cell_idx_matched06 = [];

for i = 1:n
    [cell_idx6, ~, ~] = get_cell_indices(datarun6, cell_type{i});
    cell_id0 = cell_list_map60(cell_idx6);
    ept = 1-cellfun(@isempty, cell_id0);
    cell_idx6 = cell_idx6.*ept;
    cell_idx6(cell_idx6 == 0) = [];
    cell_id0 = cell2mat(cell_id0);
    cell_idx0 = get_cell_indices(datarun0, cell_id0);
    cell_id6 = datarun6.cell_ids(cell_idx6);
    id = [cell_id0' cell_id6'];
    idx = [cell_idx0' cell_idx6'];
    cell_id06{i} = id;
    cell_idx06{i} = idx;
    cell_id_matched06 = [cell_id_matched06; id];
    cell_idx_matched06 = [cell_idx_matched06; idx];

end






load('map130214.mat')


% get parameters
datarun0 = compute_sta_fits_sequence(datarun0, cell_id_matched06(:, 1), 'verbose', true);
datarun6 = compute_sta_fits_sequence(datarun6, cell_id_matched06(:, 2), 'verbose', true);

save_sta_fits(datarun0);
save_sta_fits(datarun6);

datarun0 = load_sta_fits(datarun0);
datarun6 = load_sta_fits(datarun6);

params06 = cell(n, 1);
for i = 1:n
    index = cell_idx06{i};
    cell_numb = size(index, 1);
    params_temp = zeros(cell_numb, 8, 2);
    for j = 1:cell_numb
        params0 = datarun0.matlab.sta_fits{index(j, 1)};
        params6 = datarun6.matlab.sta_fits{index(j, 2)};
        if isempty(params0) == 1 || isempty(params6) == 1;
           params_temp(j, :, 1) = zeros(1, 8);
           params_temp(j, :, 2) = zeros(1, 8);

        else  
        params_temp(j, 1, 1) = params0.center_point_x*4/3;
        params_temp(j, 2, 1) = params0.center_point_y*4/3;
        params_temp(j, 3, 1) = params0.center_sd_x*4/3;
        params_temp(j, 4, 1) = params0.center_sd_y*4/3;
        params_temp(j, 5, 1) = params0.center_rotation_angle;
        params_temp(j, 6, 1) = params0.surround_sd_scale*params0.center_sd_x*4/3;
        params_temp(j, 7, 1) = params0.surround_sd_scale*params0.center_sd_y*4/3;
        params_temp(j, 8, 1) = params0.surround_amp_scale;
        
        params_temp(j, 1, 2) = params6.center_point_x - 1;
        params_temp(j, 2, 2) = params6.center_point_y + 1.37;
        params_temp(j, 3, 2) = params6.center_sd_x;
        params_temp(j, 4, 2) = params6.center_sd_y;
        params_temp(j, 5, 2) = params6.center_rotation_angle;
        params_temp(j, 6, 2) = params6.surround_sd_scale*params6.center_sd_x;
        params_temp(j, 7, 2) = params6.surround_sd_scale*params6.center_sd_y;
        params_temp(j, 8, 2) = params6.surround_amp_scale;
        end
        
    end
    params06{i} = params_temp;
end

%% plot mosaic

cell_type_n = 5;
params_temp = params06{cell_type_n};
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

d0 = size(datarun0.stas.stas{1});
d6 = size(datarun6.stas.stas{1});
f0 = 40;
f6 = 30;


sta_f_weighted = cell(n, 2);
for k = 1:n;
    index = cell_idx06{k};
    cell_numb = size(cell_id06{k}, 1);
    sta_f0 = zeros(f0*d0(1)/2, f0*d0(1)/2, 1, d0(4));
    sta_f6 = zeros(f6*d6(1)/2, f6*d6(1)/2, 1, d6(4));
    if isempty(index) == 0
       for i = 1:cell_numb
           sta = datarun0.stas.stas{index(i, 1)};
           if k == 1 || k == 2 
              [~, b] = max(sta(:));
           else
              [~, b] = min(sta(:));
           end
           z = ceil(b/(d0(1)*d0(2)));
           sta_1 = sta(:, :, 1, z);
           noise = robust_std(sta_1(:), 3);
           sta_temp = zeros(f0*d0(1), f0*d0(2), d0(3), d0(4));
           for j = 1:d0(4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f0, f0));
           end
           xy = floor(f0*rf_com(sta));
           if sum(xy<=f0*d0(1)/4) == 0 && sum(xy>=3*f0*d0(1)/4) == 0 && isempty(xy) == 0
           sta_f0 = sta_f0 + sta_temp(xy(2)-f0*d0(1)/4+1:xy(2)+f0*d0(1)/4, xy(1)-f0*d0(2)/4+1:xy(1)+f0*d0(2)/4, :, :)/noise;
           end
           
    
           sta = datarun6.stas.stas{index(i, 2)};
           if k == 1 || k == 2
              [~, b] = max(sta(:));
           else
              [~, b] = min(sta(:));
           end
           z = ceil(b/(d6(1)*d6(2)));
           sta_1 = sta(:, :, 1, z);
           noise = robust_std(sta_1(:), 3);
           sta_temp = zeros(f6*d6(1), f6*d6(2), d6(3), d6(4));
           for j = 1:d6(4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f6, f6));
           end
           xy = floor(f6*rf_com(sta));
           if sum(xy<=f0*d0(1)/4) == 0 && sum(xy>=3*f0*d0(1)/4) == 0 && isempty(xy) == 0
           sta_f6 = sta_f6 + sta_temp(xy(2)-f6*d6(1)/4+1:xy(2)+f6*d6(1)/4, xy(1)-f6*d6(2)/4+1:xy(1)+f6*d6(2)/4, :, :)/noise;
           end
           
       end       
    end
    sta_f_weighted{k, 1} = sta_f0;
    sta_f_weighted{k, 2} = sta_f6;
end


%% mean RF for different sampling


sta_c_15 = cell(n, 2);
for i = 1:n
        sta = sta_f_weighted{i, 1};
        sta_c_temp = zeros(d0(1)/2, d0(2)/2, d0(3), d0(4));
        for j = 1:d0(1)/2
            for k = 1:d0(2)/2
                for m = 1:d0(4)
                    sta_temp = sta(f0*(j-1)+1:f0*j, f0*(k-1)+1:f0*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c_15{i, 1} = sta_c_temp;
    
    sta = sta_f_weighted{i, 2};
    sta_c_temp = zeros(d6(1)/2, d6(2)/2, d6(3), d6(4));
        for j = 1:d6(1)/2
            for k = 1:d6(2)/2
                for m = 1:d6(4)
                    sta_temp = sta(f6*(j-1)+1:f6*j, f6*(k-1)+1:f6*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c_15{i, 2} = sta_c_temp;
end


sta_c_30 = cell(n, 2);
for i = 1:n
        sta = sta_f_weighted{i, 1};
        sta_c_temp = zeros(d0(1), d0(2), d0(3), d0(4));
        for j = 1:d0(1)
            for k = 1:d0(2)
                for m = 1:d0(4)
                    sta_temp = sta(f0/2*(j-1)+1:f0/2*j, f0/2*(k-1)+1:f0/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c_30{i, 1} = sta_c_temp;
    
    sta = sta_f_weighted{i, 2};
    sta_c_temp = zeros(d6(1), d6(2), d6(3), d6(4));
        for j = 1:d6(1)
            for k = 1:d6(2)
                for m = 1:d6(4)
                    sta_temp = sta(f6/2*(j-1)+1:f6/2*j, f6/2*(k-1)+1:f6/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c_30{i, 2} = sta_c_temp;
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

fit_sta_15 = cell(n, 2);
params_15 = cell(n, 2);
for i = 1:n
        
    fit_temp = fit_sta_sequence(sta_c_15{i, 1});
    fit_sta_15{i, 1} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x*4/3;
    params_temp(2) = fit_temp.center_point_y*4/3;
    params_temp(3) = fit_temp.center_sd_x*4/3;
    params_temp(4) = fit_temp.center_sd_y*4/3;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale*4/3;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale*4/3;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_15{i, 1} = params_temp;


    fit_temp = fit_sta_sequence(sta_c_15{i, 2});
    fit_sta_15{i, 2} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x;
    params_temp(2) = fit_temp.center_point_y;
    params_temp(3) = fit_temp.center_sd_x;
    params_temp(4) = fit_temp.center_sd_y;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_15{i, 2} = params_temp;
    
    i
    

end


fit_sta_30 = cell(n, 2);
params_30 = cell(n, 2);

for i = 1:n
        
    fit_temp = fit_sta_sequence(sta_c_30{i, 1});
    fit_sta_30{i, 1} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x*4/3;
    params_temp(2) = fit_temp.center_point_y*4/3;
    params_temp(3) = fit_temp.center_sd_x*4/3;
    params_temp(4) = fit_temp.center_sd_y*4/3;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale*4/3;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale*4/3;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_30{i, 1} = params_temp;


    fit_temp = fit_sta_sequence(sta_c_30{i, 2});
    fit_sta_30{i, 2} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x;
    params_temp(2) = fit_temp.center_point_y;
    params_temp(3) = fit_temp.center_sd_x;
    params_temp(4) = fit_temp.center_sd_y;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_30{i, 2} = params_temp;
    
    i

end






%% plot mean RF versus fitting

cell_type_n = 1;
light_level = 2;

sta = sta_c_15{cell_type_n, light_level};

if light_level == 1
    d = d0(1)/2;
else
    d = d6(1)/2;
end

if cell_type_n == 1 || cell_type_n == 2 
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
    params_temp = params_15{cell_type_n, light_level}(1:7)/4*3;
    params_temp(5) = params_15{cell_type_n, light_level}(5);
else
    params_temp = params_15{cell_type_n, light_level};
end

[X1, Y1] = drawEllipse(params_temp(1:5));
[X2, Y2] = drawEllipse(params_temp([1, 2, 6, 7, 5]));
plot(X1, Y1, 'r')
plot(X2, Y2, 'b')
legend('center', 'surround')
NDF = [4 0];
title([cell_type{cell_type_n} ' NDF = ' num2str(NDF(light_level))])




if cell_type_n == 1 || cell_type_n == 2 
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

para_temp = params_15{cell_type_n, light_level};
if light_level == 1
    fit_para1 = struct('center_point_x', para_temp(1)*3/4, 'center_point_y', para_temp(2)*3/4, ...
    'sd_x', para_temp(3)*3/4, 'sd_y', para_temp(4)*3/4, 'rotation_angle', para_temp(5), ...
    'x_dim', d, 'y_dim', d);

    fit_para2 = struct('center_point_x', para_temp(1)*3/4, 'center_point_y', para_temp(2)*3/4, ...
    'sd_x', para_temp(6)*3/4, 'sd_y', para_temp(7)*3/4, 'amp_scale', para_temp(8), ...
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
if cell_type_n == 1 || cell_type_n == 2 
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
if cell_type_n == 1 || cell_type_n == 2 
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
xlim([0 10])
legend('data', 'fit')
title([cell_type{cell_type_n} '   NDF ' num2str(NDF(light_level))]);

%% real data of mean RF

cell_type_n = 3;
distance = cell(2, 1);
STA = cell(2, 1);

for k = 1:2
    
sta = sta_c_15{cell_type_n, k};

if k == 1
    d = d0(1)/2;
else
    d = d6(1)/2;
end

if cell_type_n == 1 || cell_type_n == 2 
    [~, b] = max(sta(:));
else
    [~, b] = min(sta(:));
end

z = ceil(b/d^2);
sta_1 = sta(:, :, 1, z);


if cell_type_n == 1 || cell_type_n == 2 
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
    if cell_type_n == 1 || cell_type_n == 2 
        plot(dis1, sta_temp, '.', 'MarkerSize', 7)
        hold on 
    else
        plot(dis1, -sta_temp, '.', 'MarkerSize', 7)
        hold on
    end
else
    if cell_type_n == 1 || cell_type_n == 2 
        plot(dis1, sta_temp, '.', 'MarkerSize', 7, 'color', 'r')
        plot(dis1, zeros(1, d^2), 'm')
    else
        plot(dis1, -sta_temp, '.', 'MarkerSize', 7, 'color', 'r')
        plot(dis1, zeros(1, d^2), 'm')
    end
end


end
xlim([0 10])
legend('NDF 4', 'NDF 0')
title([cell_type{cell_type_n} '   NDF ' num2str(NDF(light_level))]);

