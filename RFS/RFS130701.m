clear all
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun3 = load_data('/Volumes/lab/Analysis/2013-07-01-0/data003/data003', opt);
datarun8 = load_data('/Volumes/lab/Analysis/2013-07-01-0/data008/data008', opt);

cell_type = {'ON transient', 'ON brisk transient', 'OFF brisk transient', 'OFF sustained', 'OFF transient'};

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun3, ...
    datarun8, cell_type);

[sta_f_weighted, sta_f_n, sta_c] = mean_rf2(datarun3, datarun8, cell_type, ...
    cell_id, cell_idx);


d0 = size(datarun3.stas.stas{1});
d1 = size(datarun8.stas.stas{1});
n = size(cell_type, 2);

fit_sta_c = cell(n, 2);
params_c = cell(n, 2);

for i = 1:n
    
    f = 40/d0(1);
%     temp_marks_sta = significant_stixels(sta_c{i, 1}, 'thresh', 3.5, 'time', 'max');
%     fit_ins = struct('sig_stixels', temp_marks_sta);
    fit_temp = fit_sta_sequence(sta_c{i, 1}); %, 'fit_instructions', fit_ins);
    fit_sta_c{i, 1} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x*f;
    params_temp(2) = fit_temp.center_point_y*f;
    params_temp(3) = fit_temp.center_sd_x*f;
    params_temp(4) = fit_temp.center_sd_y*f;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale*f;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale*f;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_c{i, 1} = params_temp;
    
    f = 40/d1(1);
%     temp_marks_sta = significant_stixels(sta_c{i, 2}, 'thresh', 5, 'time', 'max');
%     fit_ins = struct('sig_stixels', temp_marks_sta);
    fit_temp = fit_sta_sequence(sta_c{i, 2}); % 'fit_instructions', fit_ins);
    fit_sta_c{i, 2} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x*f;
    params_temp(2) = fit_temp.center_point_y*f;
    params_temp(3) = fit_temp.center_sd_x*f;
    params_temp(4) = fit_temp.center_sd_y*f;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale*f;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale*f;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_c{i, 2} = params_temp;


    i

end

save('RFS130701.mat', 'sta_f_weighted', 'sta_f_n', 'sta_c', 'fit_sta_c', 'params_c')

%%

% ratio_s = zeros(4, 4);

cell_type_n = 4;
figure
for k = 1:4
   if k < 3
       load('RFS130701.mat')
   else
       load('RFS130221_map3.mat')
       sta_c = sta_c(:, [1 3]);
       fit_sta_c = fit_sta_c(:, [1 3]);
       params_c = params_c(:, [1 3]);
       sta_f_n = sta_f_n(:, [1 3]);
   end
   m = [1 2 1 2];
sta = sta_c{cell_type_n, m(k)};
d = size(sta, 1);

if cell_type{cell_type_n}(2) == 'N' 
    [~, b] = max(sta(:));
else
    [~, b] = min(sta(:));
end

z = ceil(b/d^2);
sta_1 = sta(:, :, 1, z);



subplot(2, 4, 2*k-1);
im = norm_image(sta_1);
image(im)
hold on

params_temp = params_c{cell_type_n, m(k)}(1:7)/40*d;
params_temp(5) = params_c{cell_type_n, m(k)}(5);

[X1, Y1] = drawEllipse(params_temp(1:5));
[X2, Y2] = drawEllipse(params_temp([1, 2, 6, 7, 5]));
plot(X1, Y1, 'r')
plot(X2, Y2, 'b')
legend('center', 'surround')
NDF = [4 0];
condition = {'night', 'day'};
title([cell_type{cell_type_n} ' NDF ' num2str(NDF(m(k))) '  ' condition{ceil(k/2)}])




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
sta_temp = reshape(sta_1, 1, d^2);

para_temp = params_c{cell_type_n, m(k)};
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
    ratio_s(cell_type_n, k) = para_temp(8)*(para_temp(6)/para_temp(3))^2;
    

else
    fit = fit2 - fit1;
    comf = rf_com(-fit);
    r = min(sta_temp(:))/min(fit(:));
    ratio_s(cell_type_n, k) = para_temp(8)*(para_temp(6)/para_temp(3))^2;
    

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

subplot(2, 4, 2*k);
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
xlim([0 d/2])
legend('data', 'fit')
title([cell_type{cell_type_n} '   NDF ' num2str(NDF(m(k))) '  ' condition{ceil(k/2)}]);
end

%%

for k = 1:2
    
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

dis = dis/d*40;
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
else
    if cell_type{cell_type_n}(2) == 'N' 
        plot(dis1, sta_temp*r, '.', 'MarkerSize', 7, 'color', 'r')
    else
        plot(dis1, -sta_temp*r, '.', 'MarkerSize', 7, 'color', 'r')
    end
end


end
plot(dis1, zeros(1, d^2), 'm')
xlim([0 20])
legend(['NDF ' num2str(NDF(1))], ['NDF ' num2str(NDF(2))])
title(cell_type{cell_type_n});

%%



figure
x = [1:4]; y = ratio_s; 
bar(x,y)
set(gca,'xticklabel',cell_type(1:4))
legend('NDF 4 night', 'NDF 0 night', 'NDF 4 day', 'NDF 0 day', 'location', 'northwest')
ylabel('surround center volume ratio')