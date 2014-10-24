%% 2013-07-22-0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changed the classification of data005-map while comparing with
% 2013-11-07-0
% off brisk transient --> off transient
% off transient --> off brisk transient
%
% 2013-11-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun{1} = load_data('/Analysis/xyao/2013-07-22-0/data002/data002', opt);
datarun{2} = load_data('/Analysis/xyao/2013-07-22-0/data005/data005', opt);
datarun{3} = load_data('/Analysis/xyao/2013-07-22-0/data010/data010', opt);

cell_type = {'ON brisk transient', 'ON transient', 'OFF brisk transient', ...
    'OFF transient', 'OFF slow', 'OFF sustained'};

%% compare between ndf4 melatonin and non-melatonin
[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun{1}, ...
    datarun{2}, cell_type);


[sta_f_weighted, sta_f_n, sta_c] = mean_rf2(datarun{1}, datarun{2}, cell_type, ...
    cell_id, cell_idx);

%% compare between non-melatonin ndf0 and ndf4
[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun{3}, ...
    datarun{2}, cell_type);


[sta_f_weighted, sta_f_n, sta_c] = mean_rf2(datarun{3}, datarun{2}, cell_type, ...
    cell_id, cell_idx);

%% fitting

n = size(cell_type, 2);
fit_sta_c = cell(n, 3);
params_c = cell(n, 3);

for j = 1:1
    d = size(datarun{3}.stas.stas{1});
    D{3} = d;
    F(3) = 1200/d(1);
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


save('RFS130722.mat', 'sta_f_weighted', 'sta_f_n', 'sta_c', 'fit_sta_c', 'params_c')


%% compare ndf 4 and ndf 0

% ratio = zeros(n, 2);

cell_type_n = 1;

figure
for k = 2:3
   
sta = sta_c{cell_type_n, k};
d = size(sta, 1);

if cell_type{cell_type_n}(2) == 'N' 
    [~, b] = max(sta(:));
else
    [~, b] = min(sta(:));
end

z = ceil(b/d^2);
sta_1 = sta(:, :, 1, z);



subplot(2, 2, 2*(k-1)-1);
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
condition = {'ndf 4', 'ndf 0'};
title([cell_type{cell_type_n} '   ' condition{k-1}])




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
    ratio(cell_type_n, k) = para_temp(8)*(para_temp(6)/para_temp(3))^2;
    

else
    fit = fit2 - fit1;
    comf = rf_com(-fit);
    r = min(sta_temp(:))/min(fit(:));
    ratio(cell_type_n, k) = para_temp(8)*(para_temp(6)/para_temp(3))^2;
    

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

subplot(2, 2, 2*(k-1));
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
title([cell_type{cell_type_n} '   ' condition{k-1}])

end



%%

figure
for cell_type_n = 1:n-1
for k = 3:-1:2
    
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
if k ~= 3
    if cell_type{cell_type_n}(2) == 'N'  
        r = max(sta_temp1(:))/max(sta_1(:));
    else
        r = min(sta_temp1(:))/min(sta_1(:));
    end
end
    
sta_temp = reshape(sta_1, 1, d^2);

if k == 3
    sta_temp1 = sta_temp;
end

if k == 3
    subplot(2, 3, cell_type_n)
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
xlim([0 d/2])
legend(condition{2}, condition{1})
title(cell_type{cell_type_n});
end

%%

figure
x = 1:n-1; y = ratio(1:n-1, 2:3); 
bar(x,y)
set(gca,'xticklabel',cell_type)
ylim([0 1.2])
legend('with melatonin', 'without melatonin')
ylabel('surround center volume ratio')





%% compare melatonin and non-melatonin

% ratio_m = zeros(n, 2);

cell_type_n = 6;

figure
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



subplot(2, 2, 2*k-1);
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
condition = {'melatonin', 'non-melatonin'};
title([cell_type{cell_type_n} '   ' condition{k}])




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
    ratio_m(cell_type_n, k) = para_temp(8)*(para_temp(6)/para_temp(3))^2;
    

else
    fit = fit2 - fit1;
    comf = rf_com(-fit);
    r = min(sta_temp(:))/min(fit(:));
    ratio_m(cell_type_n, k) = para_temp(8)*(para_temp(6)/para_temp(3))^2;
    

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

subplot(2, 2, 2*k);
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
title([cell_type{cell_type_n} '   ' condition{k}])

end



%%

figure
for cell_type_n = 1:n
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
    subplot(2, 3, cell_type_n)
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
xlim([0 d/2])
legend(condition{1}, condition{2})
title(cell_type{cell_type_n});
end

%%

figure
x = 1:n; y = ratio_m; 
bar(x,y)
set(gca,'xticklabel',cell_type)
ylim([0 1.2])
legend('with melatonin', 'without melatonin')
ylabel('surround center volume ratio')