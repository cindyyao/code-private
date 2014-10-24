% NDF 4

cell_type_n = 1;
cell_n = 3;

idx_temp = cell_idx04{cell_type_n}(cell_n, 1);
sta = datarun0.stas.stas{idx_temp};
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

para_temp = params04{cell_type_n}(cell_n, :, 1);
fit_para1 = struct('center_point_x', para_temp(1), 'center_point_y', para_temp(2), ...
    'sd_x', para_temp(3), 'sd_y', para_temp(4), 'rotation_angle', para_temp(5), ...
    'x_dim', 30, 'y_dim', 30);

fit_para2 = struct('center_point_x', para_temp(1), 'center_point_y', para_temp(2), ...
    'sd_x', para_temp(6), 'sd_y', para_temp(7), 'amp_scale', para_temp(8), ...
    'rotation_angle', para_temp(5), 'x_dim', 30, 'y_dim', 30);


fit1 = make_Gaussian_two_d(fit_para1);
fit2 = make_Gaussian_two_d(fit_para2);
fit = fit1 - fit2;

% figure
% colormap gray
% imagesc(fit)

comf = rf_com(fit);
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


figure
plot(dis1, sta_temp, '.', 'MarkerSize', 7)
% plot(dis1, -sta_temp, '.', 'MarkerSize', 7)

hold on 
plot(dis1, zeros(1, 900), 'm')
plot(dis1f, fit1*r, '.', 'MarkerSize', 7, 'Color', 'r')
xlim([0 10])
title([cell_type{cell_type_n} '   NDF 4   cell n = ' num2str(cell_n)]);










% NDF 0

% cell_type_n = 1;
% cell_n = 6;

idx_temp = cell_idx04{cell_type_n}(cell_n, 2);
sta = datarun4.stas.stas{idx_temp};

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
dis = dis*0.75;
dis1 = reshape(dis, 1, 1600);


sta_temp = reshape(sta, 1, 1600);

para_temp = params04{cell_type_n}(cell_n, :, 2);
fit_para1 = struct('center_point_x', para_temp(1), 'center_point_y', para_temp(2), ...
    'sd_x', para_temp(3), 'sd_y', para_temp(4), 'rotation_angle', para_temp(5), ...
    'x_dim', 40, 'y_dim', 40);

fit_para2 = struct('center_point_x', para_temp(1), 'center_point_y', para_temp(2), ...
    'sd_x', para_temp(6), 'sd_y', para_temp(7), 'amp_scale', para_temp(8), ...
    'rotation_angle', para_temp(5), 'x_dim', 40, 'y_dim', 40);


fit1 = make_Gaussian_two_d(fit_para1);
fit2 = make_Gaussian_two_d(fit_para2);
fit = fit1 - fit2;

% figure
% colormap gray
% imagesc(fit)

comf = rf_com(fit);
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


figure
plot(dis1, sta_temp, '.', 'MarkerSize', 7)
% plot(dis1, -sta_temp, '.', 'MarkerSize', 7)
hold on 
plot(dis1, zeros(1, 1600), 'm')
plot(dis1f, fit1*r, '.', 'MarkerSize', 7, 'Color', 'r')
xlim([0 10])
title([cell_type{cell_type_n} '   NDF 0   cell n = ' num2str(cell_n)]);











% %% plot mosiac
% 
% cell_type_n = 1;
% params_temp = params04{cell_type_n};
% cell_numb = size(params_temp, 1);
% figure;
% for i = 1:cell_numb
%     [X1, Y1] = drawEllipse(params_temp(i, 1:5, 1));
%     [X2, Y2] = drawEllipse(params_temp(i, [1, 2, 6, 7, 5], 1));
%     if i == 1
%         h1 = plot(X1, Y1, 'r');
%         hold on
%         h2 = plot(X2, Y2, 'b');
%     else
%         plot(X1, Y1, 'r');
%         plot(X2, Y2, 'b');
%     end
% end
% 
% axis([0 40 0 40]);
% legend([h1 h2], 'center', 'surround');
% title([cell_type_n{cell_type_n} ' NDF 4']);
% 
% figure;
% for i = 1:cell_numb
%     [X1, Y1] = drawEllipse(params_temp(i, 1:5, 2));
%     [X2, Y2] = drawEllipse(params_temp(i, [1, 2, 6, 7, 5], 2));
%     if i == 1
%         h1 = plot(X1, Y1, 'r');
%         hold on
%         h2 = plot(X2, Y2, 'b');
%     else
%         plot(X1, Y1, 'r');
%         plot(X2, Y2, 'b');
%     end
% end
% 
% axis([0 40 0 40]);
% legend([h1 h2], 'center', 'surround');
% title([cell_type_n{cell_type_n} ' NDF 0']);
