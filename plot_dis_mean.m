cell_type_n = 1;

sta = sta_c{cell_type_n, 1};

[~, b] = max(sta(:));

% [~, b] = min(sta(:));

z = ceil(b/900);
sta = sta(:, :, 1, z);

com = rf_com(sta);

% com = rf_com(-sta);

pix = [1:30] + 0.5;
[X Y] = meshgrid(pix);
xdis = com(1) - X;
ydis = com(2) - Y;

dis = zeros(30, 30);
for i = 1:30
    for j = 1:30
        dis(i, j) = norm([xdis(i, j) ydis(i, j)]);
    end
end

dis1 = reshape(dis, 1, 900);


sta_temp = reshape(sta, 1, 900);

para_temp = params_c{cell_type_n, 1};
fit_para1 = struct('center_point_x', para_temp(1)*3/4, 'center_point_y', para_temp(2)*3/4, ...
    'sd_x', para_temp(3)*3/4, 'sd_y', para_temp(4)*3/4, 'rotation_angle', para_temp(5), ...
    'x_dim', 30, 'y_dim', 30);

fit_para2 = struct('center_point_x', para_temp(1)*3/4, 'center_point_y', para_temp(2)*3/4, ...
    'sd_x', para_temp(6)*3/4, 'sd_y', para_temp(7)*3/4, 'amp_scale', abs(para_temp(8)), ...
    'rotation_angle', para_temp(5), 'x_dim', 30, 'y_dim', 30);


fit1 = make_Gaussian_two_d(fit_para1);
fit2 = make_Gaussian_two_d(fit_para2);
fit = fit1 - fit2;
% fit = fit2 - fit1;

% figure
% colormap gray
% imagesc(fit)

comf = rf_com(fit);
% comf = rf_com(-fit);

pix = [1:30] + 0.5;
[X Y] = meshgrid(pix);
xdisf = comf(1) - X;
ydisf = comf(2) - Y;

disf = zeros(30, 30);
for i = 1:30
    for j = 1:30
        disf(i, j) = norm([xdisf(i, j) ydisf(i, j)]);
    end
end

dis1f = reshape(disf, 1, 900);
fit1 = reshape(fit, 1, 900);




r = max(sta_temp(:))/max(fit(:));
% r = min(sta_temp(:))/min(fit(:));



figure
plot(dis1, sta_temp, '.', 'MarkerSize', 7)
% plot(dis1, -sta_temp, '.', 'MarkerSize', 7)

hold on 
plot(dis1f, fit1*r, '.', 'MarkerSize', 7, 'Color', 'r')
% plot(dis1f, -fit1*r, '.', 'MarkerSize', 7, 'Color', 'r')

plot(dis1, zeros(1, 900), 'm')
xlim([0 20])
legend('data', 'fit')
title([cell_type{cell_type_n} '   NDF 4   ']);










% NDF 0

% cell_type_n = 1;
% cell_n = 6;

sta = sta_c{cell_type_n, 2};

[~, b] = max(sta(:));
% [~, b] = min(sta(:));

z = ceil(b/1600);
sta = sta(:, :, 1, z);

com = rf_com(sta);
% com = rf_com(-sta);

pix = [1:40] + 0.5;
[X Y] = meshgrid(pix);
xdis = com(1) - X;
ydis = com(2) - Y;

dis = zeros(40, 40);
for i = 1:40
    for j = 1:40
        dis(i, j) = norm([xdis(i, j) ydis(i, j)]);
    end
end
% dis = dis*0.75; % compare NDF 4 with NDF 0

dis1 = reshape(dis, 1, 1600);


sta_temp = reshape(sta, 1, 1600);

para_temp = params_c{cell_type_n, 2};
fit_para1 = struct('center_point_x', para_temp(1), 'center_point_y', para_temp(2), ...
    'sd_x', para_temp(3), 'sd_y', para_temp(4), 'rotation_angle', para_temp(5), ...
    'x_dim', 40, 'y_dim', 40);

fit_para2 = struct('center_point_x', para_temp(1), 'center_point_y', para_temp(2), ...
    'sd_x', para_temp(6), 'sd_y', para_temp(7), 'amp_scale', abs(para_temp(8)), ...
    'rotation_angle', para_temp(5), 'x_dim', 40, 'y_dim', 40);


fit1 = make_Gaussian_two_d(fit_para1);
fit2 = make_Gaussian_two_d(fit_para2);
fit = fit1 - fit2;
% fit = fit2 - fit1;

% figure
% colormap gray
% imagesc(fit)

comf = rf_com(fit);
% comf = rf_com(-fit);

pix = [1:40] + 0.5;
[X Y] = meshgrid(pix);
xdisf = comf(1) - X;
ydisf = comf(2) - Y;

disf = zeros(40, 40);
for i = 1:40
    for j = 1:40
        disf(i, j) = norm([xdisf(i, j) ydisf(i, j)]);
    end
end

dis1f = reshape(disf, 1, 1600);
fit1 = reshape(fit, 1, 1600);




r = max(sta_temp(:))/max(fit(:));
% r = min(sta_temp(:))/min(fit(:));



figure
plot(dis1, sta_temp, '.', 'MarkerSize', 7)
% plot(dis1, -sta_temp, '.', 'MarkerSize', 7)
hold on 
plot(dis1f, fit1*r, '.', 'MarkerSize', 7, 'Color', 'r')
% plot(dis1f, -fit1*r, '.', 'MarkerSize', 7, 'Color', 'r')

plot(dis1, zeros(1, 1600), 'm')
xlim([0 20])
legend('data', 'fit')
title([cell_type{cell_type_n} '  NDF 0   ']);
