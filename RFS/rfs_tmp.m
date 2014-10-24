%%
a = clock;
    a(4:5)
    fit_temp = fit_sta_sequence(sta_c_30{4, 2}, 'verbose', true);
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x;
    params_temp(2) = fit_temp.center_point_y;
    params_temp(3) = fit_temp.center_sd_x;
    params_temp(4) = fit_temp.center_sd_y;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale;
    params_temp(8) = fit_temp.surround_amp_scale;
    a = clock;
    a(4:5)
    
%%

cell_type_n = 4;
light_level = 2;

sta = sta_c_30{cell_type_n, light_level};

if light_level == 1
    d = d0(1);
else
    d = d4(1);
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

para_temp = params_30{cell_type_n, light_level};
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
    ratio = -surround/center
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
title([cell_type{cell_type_n} '   NDF ' num2str(NDF(light_level))]);


%% 

rmse = zeros(8, 2);
for i = 1:2
    for j = 1:4
        fit_temp = fit_sta_30{j, i};
        rmse((i-1)*4+j, 1) = fit_temp.rmse;
        fit_temp = fit_sta_30_2{j, i};
        rmse((i-1)*4+j, 2) = fit_temp.rmse;
    end
end

