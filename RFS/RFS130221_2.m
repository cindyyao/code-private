% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun0 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data000/data000', opt);
datarun2 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data002/data002', opt);
datarun4 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data004/data004', opt);

% map ei
[cell_list_map40, ~] = map_ei(datarun4, datarun0);
[cell_list_map42, ~] = map_ei(datarun4, datarun2);

% map data000 and data004
cell_type = {'ON transient', 'ON brisk transient', 'OFF brisk transient', 'OFF sustained'};
n = length(cell_type);
cell_id024 = cell(n, 1);
cell_idx024 = cell(n, 1);
cell_id_matched024 = [];
cell_idx_matched024 = [];

for i = 1:n
    [cell_idx4, ~, ~] = get_cell_indices(datarun4, cell_type{i});
    cell_id2 = cell_list_map42(cell_idx4);
    ept = 1-cellfun(@isempty, cell_id2);
    cell_idx4 = cell_idx4.*ept;
    cell_idx4(cell_idx4 == 0) = [];
    cell_id2 = cell2mat(cell_id2);
    cell_idx2 = get_cell_indices(datarun2, cell_id2);
    cell_id4 = datarun4.cell_ids(cell_idx4);
    
    cell_id0 = cell_list_map40(cell_idx4);
    ept = 1-cellfun(@isempty, cell_id0);
    cell_idx4 = cell_idx4.*ept;
    cell_idx2 = cell_idx2.*ept;
    cell_idx4(cell_idx4 == 0) = [];
    cell_idx2(cell_idx2 == 0) = [];
    cell_id0 = cell2mat(cell_id0);
    cell_idx0 = get_cell_indices(datarun0, cell_id0);
    cell_id4 = datarun4.cell_ids(cell_idx4);
    cell_id2 = datarun2.cell_ids(cell_idx2);
    
    id = [cell_id0' cell_id2' cell_id4'];
    idx = [cell_idx0' cell_idx2' cell_idx4'];
    cell_id024{i} = id;
    cell_idx024{i} = idx;
    cell_id_matched024 = [cell_id_matched024; id];
    cell_idx_matched024 = [cell_idx_matched024; idx];

    
   
end



% get parameters
datarun2 = compute_sta_fits_sequence(datarun2, cell_id_matched024(:, 1), 'verbose', true);
datarun4 = compute_sta_fits_sequence(datarun4, cell_id_matched024(:, 2), 'verbose', true);

save_sta_fits(datarun2);
save_sta_fits(datarun4);

datarun2 = load_sta_fits(datarun2);
datarun4 = load_sta_fits(datarun4);

params24 = cell(n, 1);
for i = 1:n
    index = cell_idx024{i};
    cell_numb = size(index, 1);
    params_temp = zeros(cell_numb, 8, 2);
    for j = 1:cell_numb
        params2 = datarun2.matlab.sta_fits{index(j, 1)};
        params4 = datarun4.matlab.sta_fits{index(j, 2)};
        if isempty(params2) == 1 || isempty(params4) == 1;
           params_temp(j, :, 1) = zeros(1, 8);
           params_temp(j, :, 2) = zeros(1, 8);

        else  
        params_temp(j, 1, 1) = params2.center_point_x;
        params_temp(j, 2, 1) = params2.center_point_y;
        params_temp(j, 3, 1) = params2.center_sd_x;
        params_temp(j, 4, 1) = params2.center_sd_y;
        params_temp(j, 5, 1) = params2.center_rotation_angle;
        params_temp(j, 6, 1) = params2.surround_sd_scale*params2.center_sd_x;
        params_temp(j, 7, 1) = params2.surround_sd_scale*params2.center_sd_y;
        params_temp(j, 8, 1) = params2.surround_amp_scale;
        
        params_temp(j, 1, 2) = params4.center_point_x;
        params_temp(j, 2, 2) = params4.center_point_y;
        params_temp(j, 3, 2) = params4.center_sd_x;
        params_temp(j, 4, 2) = params4.center_sd_y;
        params_temp(j, 5, 2) = params4.center_rotation_angle;
        params_temp(j, 6, 2) = params4.surround_sd_scale*params4.center_sd_x;
        params_temp(j, 7, 2) = params4.surround_sd_scale*params4.center_sd_y;
        params_temp(j, 8, 2) = params4.surround_amp_scale;
        end
        
    end
    params24{i} = params_temp;
end


%% plot mosaic

cell_type_n = 5;
params_temp = params24{cell_type_n};
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
d2 = size(datarun2.stas.stas{1});
d4 = size(datarun4.stas.stas{1});
f0 = 40;
f2 = 30;
f4 = 30;


sta_f_weighted = cell(n, 3);
sta_f_n = zeros(n, 3);

for k = 1:n;
    index = cell_idx024{k};
    cell_numb = size(cell_id024{k}, 1);
    sta_f0 = zeros(f0*d0(1)/2, f0*d0(1)/2, 1, d0(4));
    sta_f2 = zeros(f2*d2(1)/2, f2*d2(1)/2, 1, d2(4));
    sta_f4 = zeros(f4*d4(1)/2, f4*d4(1)/2, 1, d4(4));
    
    sta_f_n0 = 0;
    sta_f_n2 = 0;
    sta_f_n4 = 0;
    
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
           sta_f_n0 = sta_f_n0 + 1;
           end
           
           sta = datarun2.stas.stas{index(i, 2)};
           if k == 1 || k == 2 
              [~, b] = max(sta(:));
           else
              [~, b] = min(sta(:));
           end
           z = ceil(b/(d2(1)*d2(2)));
           sta_1 = sta(:, :, 1, z);
           noise = robust_std(sta_1(:), 3);
           sta_temp = zeros(f2*d2(1), f2*d2(2), d2(3), d2(4));
           for j = 1:d2(4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f2, f2));
           end
           xy = floor(f2*rf_com(sta));
           if sum(xy<=f2*d2(1)/4) == 0 && sum(xy>=3*f2*d2(1)/4) == 0 && isempty(xy) == 0
           sta_f2 = sta_f2 + sta_temp(xy(2)-f2*d2(1)/4+1:xy(2)+f2*d2(1)/4, xy(1)-f2*d2(2)/4+1:xy(1)+f2*d2(2)/4, :, :)/noise;
           sta_f_n2 = sta_f_n2 + 1;
           end
           
    
           sta = datarun4.stas.stas{index(i, 3)};
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
           if sum(xy<=f2*d2(1)/4) == 0 && sum(xy>=3*f2*d2(1)/4) == 0 && isempty(xy) == 0
           sta_f4 = sta_f4 + sta_temp(xy(2)-f4*d4(1)/4+1:xy(2)+f4*d4(1)/4, xy(1)-f4*d4(2)/4+1:xy(1)+f4*d4(2)/4, :, :)/noise;
           sta_f_n4 = sta_f_n4 + 1;
           end
           
       end       
    end
    
    sta_f_weighted{k, 1} = sta_f0;
    sta_f_weighted{k, 2} = sta_f2;
    sta_f_weighted{k, 3} = sta_f4;
    sta_f_n(k, 1) = sta_f_n0;
    sta_f_n(k, 2) = sta_f_n2;
    sta_f_n(k, 3) = sta_f_n4;
end













%% mean RF for different sampling


sta_c = cell(n, 3);
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
    sta_c{i, 1} = sta_c_temp;
    
    sta = sta_f_weighted{i, 2};
    sta_c_temp = zeros(d2(1), d2(2), d2(3), d2(4));
        for j = 1:d2(1)
            for k = 1:d2(2)
                for m = 1:d2(4)
                    sta_temp = sta(f2/2*(j-1)+1:f2/2*j, f2/2*(k-1)+1:f2/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c{i, 2} = sta_c_temp;
    
    sta = sta_f_weighted{i, 3};
    sta_c_temp = zeros(d4(1), d4(2), d4(3), d4(4));
        for j = 1:d4(1)
            for k = 1:d4(2)
                for m = 1:d4(4)
                    sta_temp = sta(f4/2*(j-1)+1:f4/2*j, f4/2*(k-1)+1:f4/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c{i, 3} = sta_c_temp;
end








%% fit sta


fit_sta_c = cell(n, 3);
params_c = cell(n, 3);



for i = 1:4
    
    fit_temp = fit_sta_sequence(sta_c{i, 1});
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


save('RFS130221_map3.mat', 'fit_sta_c', 'params_c', 'sta_f_weighted',...
    'sta_c', 'n', 'd0', 'd2', 'd4', 'f0', 'f2', 'f4', 'cell_type')











%% plot mean RF versus fitting

dim = [d0; d2; d4];
factor = [f0 f2 f4];
cell_type_n = 4;
ratio_all = zeros(n, 3);

for k = 1:3;
    
figure;
sta = sta_c{cell_type_n, k};
d = dim(k, 1);

if cell_type_n == 1 || cell_type_n == 2 
    [~, b] = max(sta(:));
else
    [~, b] = min(sta(:));
end

z = ceil(b/d^2);
sta_1 = sta(:, :, 1, z);



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
legend('center', 'surround');
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
legend('data', 'fit');
title([cell_type{cell_type_n} '   NDF ' num2str(NDF(k))]);
end









%% real data of mean RF

cell_type_n = 4;
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



%%
figure
load('RFS130214_map3.mat')
a = [-4 -2 0];
plot(a, ratio_all([1, 2, 3, 5, 4], :)', 'marker', 'x')
hold on

load('RFS130221_map3.mat')
plot(a, ratio_all', 'marker', 'x')

legend('ON transient', 'ON brisk transient', 'OFF brisk transient', 'OFF sustained', 'OFF transient', 'Location', 'northeastoutside')
xlabel('log(light intensity)')
ylabel('surround/center volume ratio')
