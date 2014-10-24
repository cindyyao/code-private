%% weighted addition of RF by SNR

d0 = size(datarun0.stas.stas{1});
d4 = size(datarun4.stas.stas{1});
f0 = 40;
f4 = 30;


sta_f_weighted = cell(n, 2);
for k = 1:n;
    index = cell_idx04{k};
    cell_numb = size(cell_id04{k}, 1);
    sta_f0 = zeros(f0*d0(1)/2, f0*d0(1)/2, 1, d0(4));
    sta_f4 = zeros(f4*d4(1)/2, f4*d4(1)/2, 1, d4(4));
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
           
    
           sta = datarun4.stas.stas{index(i, 2)};
           if k == 1 || k == 2
              [~, b] = max(sta(:));
           else
              [~, b] = min(sta(:));
           end
           z = ceil(b/(d4(1)*d4(2)));
           sta_1 = sta(:, :, 1, z);
           noise = robust_std(sta_1(:), 3);
           sta_temp = zeros(f4*d4(1), f4*d4(2), d4(3), d4(4));
           for j = 1:d4(4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f4, f4));
           end
           xy = floor(f4*rf_com(sta));
           if sum(xy<=f0*d0(1)/4) == 0 && sum(xy>=3*f0*d0(1)/4) == 0 && isempty(xy) == 0
           sta_f4 = sta_f4 + sta_temp(xy(2)-f4*d4(1)/4+1:xy(2)+f4*d4(1)/4, xy(1)-f4*d4(2)/4+1:xy(1)+f4*d4(2)/4, :, :)/noise;
           end
           
       end       
    end
    sta_f_weighted{k, 1} = sta_f0;
    sta_f_weighted{k, 2} = sta_f4;
end













%% mean RF for different sampling

% no significant stixel:  [4 1]

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
    sta_c_temp = zeros(d4(1)/2, d4(2)/2, d4(3), d4(4));
        for j = 1:d4(1)/2
            for k = 1:d4(2)/2
                for m = 1:d4(4)
                    sta_temp = sta(f4*(j-1)+1:f4*j, f4*(k-1)+1:f4*k, 1, m);
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
    sta_c_temp = zeros(d4(1), d4(2), d4(3), d4(4));
        for j = 1:d4(1)
            for k = 1:d4(2)
                for m = 1:d4(4)
                    sta_temp = sta(f4/2*(j-1)+1:f4/2*j, f4/2*(k-1)+1:f4/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c_30{i, 2} = sta_c_temp;
end


sta15 = sta_c_15{4, 1};
[~, b] = max(sta15(:));

z = ceil(b/225);
sta15 = sta15(:, :, 1, z);

sta30 = sta_c_30{4, 1};
[~, b] = max(sta30(:));

z = ceil(b/900);
sta30 = sta30(:, :, 1, z);

figure
hist(sta15(:), 20)
xlabel('stixel value')
ylabel('number of stixels')
title('15*15')
figure
hist(sta30(:), 20)
xlabel('stixel value')
ylabel('number of stixels')
title('30*30')












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

    a = clock;
    a(4:5)
    fit_temp = fit_sta_sequence(sta_c_15{i, 2}, 'verbose', true);
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
    a = clock;
    a(4:5)
    i
    

end


fit_sta_30 = cell(n, 2);
params_30 = cell(n, 2);

fit_sta_30_2 = cell(n, 2);
params_30_2 = cell(n, 2);


for i = 1:4
        
    fit_temp = fit_sta_sequence(sta_c_30{i, 1});
    fit_sta_30_2{i, 1} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x*4/3;
    params_temp(2) = fit_temp.center_point_y*4/3;
    params_temp(3) = fit_temp.center_sd_x*4/3;
    params_temp(4) = fit_temp.center_sd_y*4/3;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale*4/3;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale*4/3;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_30_2{i, 1} = params_temp;

    
    fit_temp = fit_sta_sequence(sta_c_30{i, 2});
    fit_sta_30_2{i, 2} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x;
    params_temp(2) = fit_temp.center_point_y;
    params_temp(3) = fit_temp.center_sd_x;
    params_temp(4) = fit_temp.center_sd_y;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_30_2{i, 2} = params_temp;
    
    i
end













%% plot mean RF versus fitting

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
if light_level == 1
    clear params_temp
    params_temp = params_30_2{cell_type_n, light_level}(1:7)/4*3;
    params_temp(5) = params_30_2{cell_type_n, light_level}(5);
else
    params_temp = params_30_2{cell_type_n, light_level};
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

para_temp = params_30_2{cell_type_n, light_level};
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









%% real data of mean RF

cell_type_n = 5;
dis = cell(2, 1);
sta = cell(2, 1)

for i = 1:2
    
sta = sta_c_15{cell_type_n, i};

if i == 1
    d = d0(1)/2;
else
    d = d4(1)/2;
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



figure
if cell_type_n == 1 || cell_type_n == 2 
    plot(dis1, sta_temp, '.', 'MarkerSize', 7)
    hold on 
    plot(dis1, zeros(1, d^2), 'm')
else
    plot(dis1, -sta_temp, '.', 'MarkerSize', 7)
    hold on
    plot(dis1, zeros(1, d^2), 'm')
end
xlim([0 10])
legend('data', 'fit')
title([cell_type{cell_type_n} '   NDF ' num2str(NDF(light_level))]);






%% check bad cell

% get index of bad cell
bad_cell = [];
for i = 1:n
    params_temp = params04{i};
    cell_numb = size(params_temp, 1);
    for j = 1:cell_numb
        for k = 1:2
            if params_temp(j, 6, k)/params_temp(j, 3, k) <= 1
                bad_cell = [bad_cell; i j k];
            end
        end
    end
end


% compare fit and data

cell_type_n = 5;
cell_numb = 7;
light_level = 1;
index = cell_idx04{cell_type_n}(cell_numb, light_level);
if light_level == 1
    sta = datarun0.stas.stas{index};
    d = d0(1);
else
    sta = datarun4.stas.stas{index};
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
% colormap gray
% imagesc(sta_1)
imshow(sta_1 + .5, 'InitialMagnification', 1000)
hold on
if light_level == 1
    clear params_temp
    params_temp = params04{cell_type_n}(cell_numb, 1:7, light_level)/4*3;
    params_temp(5) = params04{cell_type_n}(cell_numb, 5, light_level);
else
    params_temp = params04{cell_type_n}(cell_numb, :, light_level);
end

[X1, Y1] = drawEllipse(params_temp(1:5));
[X2, Y2] = drawEllipse(params_temp([1, 2, 6, 7, 5]));
plot(X1, Y1, 'r')
plot(X2, Y2, 'b')
legend('center', 'surround')
NDF = [4 0];
title([cell_type{cell_type_n} '  cell n = ' num2str(cell_numb) '  NDF = ' num2str(NDF(light_level))])

% amplitude vs distance

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

para_temp = params04{cell_type_n}(cell_numb, :, light_level);
fit_para1 = struct('center_point_x', para_temp(1), 'center_point_y', para_temp(2), ...
    'sd_x', para_temp(3), 'sd_y', para_temp(4), 'rotation_angle', para_temp(5), ...
    'x_dim', d, 'y_dim', d);

fit_para2 = struct('center_point_x', para_temp(1), 'center_point_y', para_temp(2), ...
    'sd_x', para_temp(6), 'sd_y', para_temp(7), 'amp_scale', para_temp(8), ...
    'rotation_angle', para_temp(5), 'x_dim', d, 'y_dim', d);


fit1 = make_Gaussian_two_d(fit_para1);
fit2 = make_Gaussian_two_d(fit_para2);

if cell_type_n == 1 || cell_type_n == 2
    fit = fit1 - fit2;
    comf = rf_com(fit);
    r = max(sta_temp(:))/max(fit(:));
    
else
    fit = fit2 - fit1;
    comf = rf_com(-fit);
    r = min(sta_temp(:))/min(fit(:));
end


% figure
% colormap gray
% imagesc(fit)

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

figure
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
title([cell_type{cell_type_n} '  cell n = ' num2str(cell_numb) '  NDF = ' num2str(NDF(light_level))])
legend('data', 'fit')


