% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun1 = load_data('/Analysis/xyao/2013-02-13-0/data001/data001', opt);
datarun3 = load_data('/Analysis/xyao/2013-02-13-0/data003/data003', opt);
datarun5 = load_data('/Analysis/xyao/2013-02-13-0/data005/data005', opt);

% map ei
[cell_list_map51, ~] = map_ei(datarun5, datarun1);
[cell_list_map53, ~] = map_ei(datarun5, datarun3);

% map data000 and data006
cell_type = {'ON transient', 'OFF brisk transient', 'OFF brisk transient large', 'OFF transient'};
n = length(cell_type);
cell_id135 = cell(n, 1);
cell_idx135 = cell(n, 1);
cell_id_matched135 = [];
cell_idx_matched135 = [];

for i = 1:n
    [cell_idx5, ~, ~] = get_cell_indices(datarun5, cell_type{i});
    cell_id3 = cell_list_map53(cell_idx5);
    ept = 1-cellfun(@isempty, cell_id3);
    cell_idx5 = cell_idx5.*ept;
    cell_idx5(cell_idx5 == 0) = [];
    cell_id3 = cell2mat(cell_id3);
    cell_idx3 = get_cell_indices(datarun3, cell_id3);
    cell_id5 = datarun5.cell_ids(cell_idx5);
    
    cell_id1 = cell_list_map51(cell_idx5);
    ept = 1-cellfun(@isempty, cell_id1);
    cell_idx5 = cell_idx5.*ept;
    cell_idx3 = cell_idx3.*ept;
    cell_idx5(cell_idx5 == 0) = [];
    cell_idx3(cell_idx3 == 0) = [];
    cell_id1 = cell2mat(cell_id1);
    cell_idx1 = get_cell_indices(datarun1, cell_id1);
    cell_id5 = datarun5.cell_ids(cell_idx5);
    cell_id3 = datarun3.cell_ids(cell_idx3);
    
    id = [cell_id1' cell_id3' cell_id5'];
    idx = [cell_idx1' cell_idx3' cell_idx5'];
    cell_id135{i} = id;
    cell_idx135{i} = idx;
    cell_id_matched135 = [cell_id_matched135; id];
    cell_idx_matched135 = [cell_idx_matched135; idx];

end





%% weighted addition of RF by SNR

d1 = size(datarun1.stas.stas{1});
d3 = size(datarun3.stas.stas{1});
d5 = size(datarun5.stas.stas{1});
f1 = 40;
f3 = 30;
f5 = 30;


sta_f_weighted = cell(n, 3);
sta_f_n = zeros(n, 3);

for k = 1:n;
    index = cell_idx135{k};
    cell_numb = size(cell_id135{k}, 1);
    sta_f1 = zeros(f1*d1(1)/2, f1*d1(1)/2, 1, d1(4));
    sta_f3 = zeros(f3*d3(1)/2, f3*d3(1)/2, 1, d3(4));
    sta_f5 = zeros(f5*d5(1)/2, f5*d5(1)/2, 1, d5(4));
    
    sta_f_n1 = 0;
    sta_f_n3 = 0;
    sta_f_n5 = 0;
    if isempty(index) == 0
       for i = 1:cell_numb
           
           sta = datarun1.stas.stas{index(i, 1)};
           if k == 1 
              [~, b] = max(sta(:));
           else
              [~, b] = min(sta(:));
           end
           z = ceil(b/(d1(1)*d1(2)));
           sta_1 = sta(:, :, 1, z);
           noise = robust_std(sta_1(:), 3);
           sta_temp = zeros(f1*d1(1), f1*d1(2), d1(3), d1(4));
           for j = 1:d1(4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f1, f1));
           end
           xy = floor(f1*rf_com(sta));
           if sum(xy<=f1*d1(1)/4) == 0 && sum(xy>=3*f1*d1(1)/4) == 0 && isempty(xy) == 0
           sta_f1 = sta_f1 + sta_temp(xy(2)-f1*d1(1)/4+1:xy(2)+f1*d1(1)/4, xy(1)-f1*d1(2)/4+1:xy(1)+f1*d1(2)/4, :, :)/noise;
           sta_f_n1 = sta_f_n1 + 1;
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
           if sum(xy<=f3*d3(1)/4) == 0 && sum(xy>=3*f3*d3(1)/4) == 0 && isempty(xy) == 0
           sta_f3 = sta_f3 + sta_temp(xy(2)-f3*d3(1)/4+1:xy(2)+f3*d3(1)/4, xy(1)-f3*d3(2)/4+1:xy(1)+f3*d3(2)/4, :, :)/noise;
           sta_f_n3 = sta_f_n3 + 1;
           end
           
    
           sta = datarun5.stas.stas{index(i, 3)};
           if k == 1 
              [~, b] = max(sta(:));
           else
              [~, b] = min(sta(:));
           end
           z = ceil(b/(d5(1)*d5(2)));
           sta_1 = sta(:, :, 1, z);
           noise = robust_std(sta_1(:), 3);
           sta_temp = zeros(f5*d5(1), f5*d5(2), d5(3), d5(4));
           for j = 1:d5(4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f5, f5));
           end
           xy = floor(f5*rf_com(sta));
           if sum(xy<=f3*d3(1)/4) == 0 && sum(xy>=3*f3*d3(1)/4) == 0 && isempty(xy) == 0
           sta_f5 = sta_f5 + sta_temp(xy(2)-f5*d5(1)/4+1:xy(2)+f5*d5(1)/4, xy(1)-f5*d5(2)/4+1:xy(1)+f5*d5(2)/4, :, :)/noise;
           sta_f_n5 = sta_f_n5 + 1;
           end
           
       end       
    end
    sta_f_weighted{k, 1} = sta_f1;
    sta_f_weighted{k, 2} = sta_f3;
    sta_f_weighted{k, 3} = sta_f5;
    
    sta_f_n(k, 1) = sta_f_n1;
    sta_f_n(k, 2) = sta_f_n3;
    sta_f_n(k, 3) = sta_f_n5;
end


%% mean RF for different sampling


sta_c = cell(n, 3);
for i = 1:n
    
     sta = sta_f_weighted{i, 1};
        sta_c_temp = zeros(d1(1), d1(2), d1(3), d1(4));
        for j = 1:d1(1)
            for k = 1:d1(2)
                for m = 1:d1(4)
                    sta_temp = sta(f1/2*(j-1)+1:f1/2*j, f1/2*(k-1)+1:f1/2*k, 1, m);
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
    
    sta = sta_f_weighted{i, 3};
    sta_c_temp = zeros(d5(1), d5(2), d5(3), d5(4));
        for j = 1:d5(1)
            for k = 1:d5(2)
                for m = 1:d5(4)
                    sta_temp = sta(f5/2*(j-1)+1:f5/2*j, f5/2*(k-1)+1:f5/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c{i, 3} = sta_c_temp;
end


%% fit sta

fit_sta_c = cell(n, 3);
params_c = cell(n, 3);

for i = 4:n
   
    fit_temp = fit_sta_sequence(sta_c{i, 1}, 'verbose', true);
    fit_sta_c{i, 1} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x*4/3;
    params_temp(2) = fit_temp.center_point_y*4/3;
    params_temp(3) = fit_temp.center_sd_x*4/3;
    params_temp(4) = fit_temp.center_sd_y*4/3;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale*4/3;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale*4/3;
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


    fit_temp = fit_sta_sequence(sta_c{i, 3});
    fit_sta_c{i, 3} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x;
    params_temp(2) = fit_temp.center_point_y;
    params_temp(3) = fit_temp.center_sd_x;
    params_temp(4) = fit_temp.center_sd_y;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_c{i, 3} = params_temp;
    
    i

end


save('RFS130213_map3.mat', 'fit_sta_c', 'params_c', 'sta_f_weighted',...
    'sta_c', 'n', 'd1', 'd3', 'd5', 'f1', 'f3', 'f5', 'cell_type')



%% plot mean RF versus fitting

dim = [d1; d3; d5];
factor = [f1 f3 f5];
cell_type_n = 4;

for k = 1:3
   
sta = sta_c{cell_type_n, k};
d = dim(k, 1);

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
clear params_temp

if k == 1
    params_temp = params_c{cell_type_n, k}(1:7)/4*3;
    params_temp(5) = params_c{cell_type_n, k}(5);
else
    params_temp = params_c{cell_type_n, k};
end

[X1, Y1] = drawEllipse(params_temp(1:5));
[X2, Y2] = drawEllipse(params_temp([1, 2, 6, 7, 5]));
plot(X1, Y1, 'r')
plot(X2, Y2, 'b')
legend('center', 'surround')
NDF = [4 2 0];
title([cell_type{cell_type_n} ' NDF = ' num2str(NDF(k))])




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

para_temp = params_c{cell_type_n, k};
if k == 1
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
    ratio = -surround/center
    
    ratio_all(cell_type_n, k) = ratio;

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
    ratio = -surround/center
    
    ratio_all(cell_type_n, k) = ratio;

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
xlim([0 20])
legend('data', 'fit')
title([cell_type{cell_type_n} '   NDF ' num2str(NDF(k))]);
end

%% real data of mean RF

cell_type_n = 5;
distance = cell(3, 1);
STA = cell(3, 1);

for k = 1:3
    
sta = sta_c{cell_type_n, k};
d = dim(k, 1);

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
if k ~= 1
    if cell_type_n == 1 || cell_type_n == 2 
        r = max(sta_temp1(:))/max(sta_1(:));
    else
        r = min(sta_temp1(:))/min(sta_1(:));
    end
end
    
sta_temp = reshape(sta_1, 1, d^2);

if k == 1
    sta_temp1 = sta_temp;
end

if k == 1
    figure
    if cell_type_n == 1 || cell_type_n == 2 
        plot(dis1, sta_temp, '.', 'MarkerSize', 7)
        hold on 
    else
        plot(dis1, -sta_temp, '.', 'MarkerSize', 7)
        hold on
    end
elseif k == 2
    if cell_type_n == 1 || cell_type_n == 2 
        plot(dis1, sta_temp*r, '.', 'MarkerSize', 7, 'color', 'r')
    else
        plot(dis1, -sta_temp*r, '.', 'MarkerSize', 7, 'color', 'r')
    end
else
    if cell_type_n == 1 || cell_type_n == 2 
        plot(dis1, sta_temp*r, '.', 'MarkerSize', 7, 'color', 'g')
    else
        plot(dis1, -sta_temp*r, '.', 'MarkerSize', 7, 'color', 'g')
    end
end


end
plot(dis1, zeros(1, d^2), 'm')
xlim([0 20])
legend('NDF 4', 'NDF 2', 'NDF 0')
title(cell_type{cell_type_n});

