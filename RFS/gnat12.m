%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun{1} = load_data('/Analysis/xyao/2013-05-28-0/data016/data016', opt);
datarun{2} = load_data('/Analysis/xyao/2013-05-29-0/data015-map/data015-map', opt);
datarun{3} = load_data('/Analysis/xyao/2013-05-30-0/data008/data008', opt);


cell_type = {'ON transient', 'ON brisk transient', 'ON sustained', 'OFF brisk transient', 'OFF transient'};
n = length(cell_type);
cell_id = cell(2, 1);
cell_idx = cell(2, 1);

for i = 1:n
    for j = 1:2
        [cell_idx{j}{i}, ~, ~] = get_cell_indices(datarun{j}, cell_type{i});
        cell_id{j}{i} = datarun{j}.cell_ids(cell_idx{j}{i});
    end
end

for i = [1 4 5]
        [cell_idx{3}{i}, ~, ~] = get_cell_indices(datarun{3}, cell_type{i});
        cell_id{3}{i} = datarun{3}.cell_ids(cell_idx{3}{i});
end

%% weighted addition of RF by SNR

d(1, :) = size(datarun{1}.stas.stas{1});
d(2, :) = size(datarun{2}.stas.stas{1});
d(3, :) = size(datarun{3}.stas.stas{1});

f(1) = 30;
f(2) = 30;
f(3) = 40;

sta_f_weighted = cell(n, 2);
sta_f_N = zeros(n, 2);

for k = 1:n;
    for g = 3:3
    index = cell_idx{g}{k};
    cell_numb = length(index);
    sta_f = zeros(f(g)*d(g, 1)/2, f(g)*d(g, 1)/2, 1, d(g, 4));
    
    sta_f_n = 0;
    
    if isempty(index) == 0
       for i = 1:cell_numb
           
           sta = datarun{g}.stas.stas{index(i)};
           if cell_type{k}(2) == 'N' 
              [~, b] = max(sta(:));
           else
              [~, b] = min(sta(:));
           end
           z = ceil(b/(d(g, 1)*d(g, 2)));
           sta_1 = sta(:, :, 1, z);
           noise = robust_std(sta_1(:), 3);
           sta_temp = zeros(f(g)*d(g, 1), f(g)*d(g, 2), d(g, 3), d(g, 4));
           for j = 1:d(g, 4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f(g), f(g)));
           end
           xy = floor(f(g)*rf_com(sta));
           if sum(xy<=f(g)*d(g, 1)/4) == 0 && sum(xy>=3*f(g)*d(g, 1)/4) == 0 && isempty(xy) == 0
           sta_f = sta_f + sta_temp(xy(2)-f(g)*d(g, 1)/4+1:xy(2)+f(g)*d(g, 1)/4, xy(1)-f(g)*d(g, 2)/4+1:xy(1)+f(g)*d(g, 2)/4, :, :)/noise;
           sta_f_n = sta_f_n + 1;
           end
           
       end       
    end
    
    sta_f_weighted{k, g} = sta_f;
   
    sta_f_N(k, g) = sta_f_n;
    end
end






%% mean RF for different sampling


sta_c = cell(n, 2);
for i = 1:n
    for g = 3:3
        sta = sta_f_weighted{i, g};
        sta_c_temp = zeros(d(g, 1), d(g, 2), d(g, 3), d(g, 4));
        for j = 1:d(g, 1)
            for k = 1:d(g, 2)
                for m = 1:d(g, 4)
                    sta_temp = sta(f(g)/2*(j-1)+1:f(g)/2*j, f(g)/2*(k-1)+1:f(g)/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c{i, g} = sta_c_temp;
    
    end
    
end



%% fit sta


fit_sta_c = cell(n, 2);
params_c = cell(n, 2);



for i = [1 4 5]
    for g = 3:3
    fit_temp = fit_sta_sequence(sta_c{i, g});
    fit_sta_c{i, g} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x*40/d(g, 1);
    params_temp(2) = fit_temp.center_point_y*40/d(g, 1);
    params_temp(3) = fit_temp.center_sd_x*40/d(g, 1);
    params_temp(4) = fit_temp.center_sd_y*40/d(g, 1);
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale*40/d(g, 1);
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale*40/d(g, 1);
    params_temp(8) = fit_temp.surround_amp_scale;
    params_c{i, g} = params_temp;
    
    end
    i
end


save('gnat12.mat', 'fit_sta_c', 'params_c', 'sta_f_weighted', 'sta_f_N', ...
    'sta_c', 'cell_type')











%% plot mean RF versus fitting
% ratio_all = zeros(n, 3);

dim = d(:, 1);
cell_type_n = 3;
figure;

for k = 1:2;
    

sta = sta_c{cell_type_n, k};
dd = dim(k);

if cell_type{cell_type_n}(2) == 'N'  
    [~, b] = max(sta(:));
else
    [~, b] = min(sta(:));
end

z = ceil(b/dd^2);
sta_1 = sta(:, :, 1, z);



subplot(2, 2, k);
im = norm_image(sta_1);
image(im)
hold on
clear params_temp
    params_temp = params_c{cell_type_n, k}(1:7)/40*dd;
    params_temp(5) = params_c{cell_type_n, k}(5);


[X1, Y1] = drawEllipse(params_temp(1:5));
[X2, Y2] = drawEllipse(params_temp([1, 2, 6, 7, 5]));
plot(X1, Y1, 'r')
plot(X2, Y2, 'b')
legend('center', 'surround');
line = {'Gnat2 KO', 'Gnat1 KO', 'Gnat1 KO'};
title([line{k} '   ' cell_type{cell_type_n} ])




if cell_type{cell_type_n}(2) == 'N'  
    com = rf_com(sta_1);
else
    com = rf_com(-sta_1);
end


pix = [1:dd] + 0.5;
[X Y] = meshgrid(pix);
xdis = com(1) - X;
ydis = com(2) - Y;

dis = zeros(dd, dd);
for i = 1:dd
    for j = 1:dd
        dis(i, j) = norm([xdis(i, j) ydis(i, j)]);
    end
end

dis1 = reshape(dis, 1, dd^2);
sta_temp = reshape(sta_1, 1, dd^2);

para_temp = params_c{cell_type_n, k};
fit_para1 = struct('center_point_x', para_temp(1)*dd/40, 'center_point_y', para_temp(2)*dd/40, ...
'sd_x', para_temp(3)*dd/40, 'sd_y', para_temp(4)*dd/40, 'rotation_angle', para_temp(5), ...
'x_dim', dd, 'y_dim', dd);

fit_para2 = struct('center_point_x', para_temp(1)*dd/40, 'center_point_y', para_temp(2)*dd/40, ...
'sd_x', para_temp(6)*dd/40, 'sd_y', para_temp(7)*dd/40, 'amp_scale', para_temp(8), ...
'rotation_angle', para_temp(5), 'x_dim', dd, 'y_dim', dd);


fit1 = make_Gaussian_two_d(fit_para1);
fit2 = make_Gaussian_two_d(fit_para2);
if cell_type{cell_type_n}(2) == 'N' 
    fit = fit1 - fit2;
    comf = rf_com(fit);
    r = max(sta_temp(:))/max(fit(:));
    ratio_all(cell_type_n, k) = para_temp(8)*(para_temp(6)/para_temp(3))^2;
else
    fit = fit2 - fit1;
    comf = rf_com(-fit);
    r = min(sta_temp(:))/min(fit(:));
    ratio_all(cell_type_n, k) = para_temp(8)*(para_temp(6)/para_temp(3))^2;


end



pix = [1:dd] + 0.5;
[X Y] = meshgrid(pix);
xdisf = comf(1) - X;
ydisf = comf(2) - Y;

disf = zeros(dd, dd);
for i = 1:dd
    for j = 1:dd
        disf(i, j) = norm([xdisf(i, j) ydisf(i, j)]);
    end
end

dis1f = reshape(disf, 1, dd^2);
fit1 = reshape(fit, 1, dd^2);

subplot(2, 2, k+2);
if cell_type{cell_type_n}(2) == 'N'
    plot(dis1, sta_temp, '.', 'MarkerSize', 7)
    hold on 
    plot(dis1f, fit1*r, '.', 'MarkerSize', 7, 'Color', 'r')
    plot(dis1, zeros(1, dd^2), 'm')
else
    plot(dis1, -sta_temp, '.', 'MarkerSize', 7)
    hold on
    plot(dis1f, -fit1*r, '.', 'MarkerSize', 7, 'Color', 'r')
    plot(dis1, zeros(1, dd^2), 'm')
end
xlim([0 20])
legend('data', 'fit');
title([line{k} '    ' cell_type{cell_type_n}]);
end









%% real data of mean RF

figure;

for cell_type_n = 1:5;

    if cell_type_n == 2 || cell_type_n == 3
        a = 1:2;
    else
        a = 1:3;
    end
for k = a
    
sta = sta_c{cell_type_n, k};
dd = dim(k);

if cell_type{cell_type_n}(2) == 'N'
    [~, b] = max(sta(:));
else
    [~, b] = min(sta(:));
end

z = ceil(b/dd^2);
sta_1 = sta(:, :, 1, z);


if cell_type{cell_type_n}(2) == 'N'
    com = rf_com(sta_1);
else
    com = rf_com(-sta_1);
end


pix = [1:dd] + 0.5;
[X Y] = meshgrid(pix);
xdis = com(1) - X;
ydis = com(2) - Y;

dis = zeros(dd, dd);
for i = 1:dd
    for j = 1:dd
        dis(i, j) = norm([xdis(i, j) ydis(i, j)]);
    end
end
dis = dis/dd*40;

dis1 = reshape(dis, 1, dd^2);

if k ~= 1
    if cell_type{cell_type_n}(2) == 'N'
        r = max(sta_temp1(:))/max(sta_1(:));
    else
        r = min(sta_temp1(:))/min(sta_1(:));
    end
end

sta_temp = reshape(sta_1, 1, dd^2);

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
plot(dis1, zeros(1, dd^2), 'm')

xlim([0 20])
if cell_type_n == 2 || cell_type_n == 3
    legend(line{1}, line{2});
else
    legend(line{1}, line{2}, line{3});
end

title(cell_type{cell_type_n});
end




%%
figure
x = [1:5]; y = ratio_all; 
bar(x,y)
set(gca,'xticklabel',cell_type(2:end))
legend('Gnat2 KO', 'Gnat1 KO', 'Gnat1 KO', 'location', 'northwest')