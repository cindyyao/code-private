function [ratio] = plot_mean_rf3(sta_c, cell_type, cell_type_n, params_c)
% [ratio] = plot_mean_rf3(sta_c, cell_type, cell_type_n, params_c)
ratio = zeros(1, 3);

for k = 1:3
   
sta = sta_c{cell_type_n, k};
d = size(sta, 1);

if cell_type{cell_type_n}(2) == 'N' 
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

params_temp = params_c{cell_type_n, k}(1:7)/40*d;
params_temp(5) = params_c{cell_type_n, k}(5);

[X1, Y1] = drawEllipse(params_temp(1:5));
[X2, Y2] = drawEllipse(params_temp([1, 2, 6, 7, 5]));
plot(X1, Y1, 'r')
plot(X2, Y2, 'b')
legend('center', 'surround')
NDF = [0 2 4];
title([cell_type{cell_type_n} ' NDF = ' num2str(NDF(k))])




if cell_type{cell_type_n}(2) == 'N' 
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
dis = dis/d*40;

dis1 = reshape(dis, 1, d^2);
sta_temp = reshape(sta_1, 1, d^2);

para_temp = params_c{cell_type_n, k};

fit_para1 = struct('center_point_x', para_temp(1)*d/40, 'center_point_y', para_temp(2)*d/40, ...
'sd_x', para_temp(3)*d/40, 'sd_y', para_temp(4)*d/40, 'rotation_angle', para_temp(5), ...
'x_dim', d, 'y_dim', d);

fit_para2 = struct('center_point_x', para_temp(1)*d/40, 'center_point_y', para_temp(2)*d/40, ...
'sd_x', para_temp(6)*d/40, 'sd_y', para_temp(7)*d/40, 'amp_scale', para_temp(8), ...
'rotation_angle', para_temp(5), 'x_dim', d, 'y_dim', d);


fit1 = make_Gaussian_two_d(fit_para1);
fit2 = make_Gaussian_two_d(fit_para2);
if cell_type{cell_type_n}(2) == 'N' 
    fit = fit1 - fit2;
    comf = rf_com(fit);
    r = max(sta_temp(:))/max(fit(:));
    ratio(k) = para_temp(8)*(para_temp(6)/para_temp(3))^2;
    
else
    fit = fit2 - fit1;
    comf = rf_com(-fit);
    r = min(sta_temp(:))/min(fit(:));
    ratio(k) = para_temp(8)*(para_temp(6)/para_temp(3))^2;
    
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
if cell_type{cell_type_n}(2) == 'N' 
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

for k = 1:3
    
sta = sta_c{cell_type_n, k};
d = size(sta, 1);

if cell_type{cell_type_n}(2) == 'N' 
    [~, b] = max(sta(:));
else
    [~, b] = min(sta(:));
end

z = ceil(b/d^2);
sta_1 = sta(:, :, 1, z);


if cell_type{cell_type_n}(2) == 'N' 
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
    if cell_type{cell_type_n}(2) == 'N'  
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
    if cell_type{cell_type_n}(2) == 'N' 
        plot(dis1, sta_temp, '.', 'MarkerSize', 7)
        hold on 
    else
        plot(dis1, -sta_temp, '.', 'MarkerSize', 7)
        hold on
    end
elseif k == 2
    if cell_type{cell_type_n}(2) == 'N' 
        plot(dis1, sta_temp*r, '.', 'MarkerSize', 7, 'color', 'r')
    else
        plot(dis1, -sta_temp*r, '.', 'MarkerSize', 7, 'color', 'r')
    end
else
    if cell_type{cell_type_n}(2) == 'N' 
        plot(dis1, sta_temp*r, '.', 'MarkerSize', 7, 'color', 'g')
    else
        plot(dis1, -sta_temp*r, '.', 'MarkerSize', 7, 'color', 'g')
    end
end


end
plot(dis1, zeros(1, d^2), 'm')
xlim([0 20])
legend('NDF 0', 'NDF 2', 'NDF 4')
title(cell_type{cell_type_n});
