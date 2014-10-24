% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun21 = load_data('/Analysis/xyao/2013-04-09-0/data021/data021', opt);
datarun25 = load_data('/Analysis/xyao/2013-04-09-0/data025/data025', opt);

% map ei
[cell_list_map51, ~] = map_ei(datarun25, datarun21);

% map data000 and data006
cell_type = {'ON transient', 'OFF transient', 'OFF brisk transient', 'OFF sustained', 'OFF sustained2'};
n = length(cell_type);
cell_id51 = cell(n, 1);
cell_idx51 = cell(n, 1);
cell_id_matched51 = [];
cell_idx_matched51 = [];

for i = 1:n
    [cell_idx5, ~, ~] = get_cell_indices(datarun25, cell_type{i});
    cell_id1 = cell_list_map51(cell_idx5);
    ept = 1-cellfun(@isempty, cell_id1);
    cell_idx5 = cell_idx5.*ept;
    cell_idx5(cell_idx5 == 0) = [];
    cell_id1 = cell2mat(cell_id1);
    cell_idx1 = get_cell_indices(datarun21, cell_id1);
    cell_id5 = datarun25.cell_ids(cell_idx5);
    id = [cell_id1' cell_id5'];
    idx = [cell_idx1' cell_idx5'];
    cell_id51{i} = id;
    cell_idx51{i} = idx;
    cell_id_matched51 = [cell_id_matched51; id];
    cell_idx_matched51 = [cell_idx_matched51; idx];

end





%% weighted addition of RF by SNR

d1 = size(datarun21.stas.stas{1});
d5 = size(datarun25.stas.stas{1});
f1 = 30;
f5 = 30;


sta_f_weighted = cell(n, 2);
for k = 1:n;
    index = cell_idx51{k};
    cell_numb = size(cell_id51{k}, 1);
    sta_f1 = zeros(f1*d1(1)/2, f1*d1(1)/2, 1, d1(4));
    sta_f5 = zeros(f5*d5(1)/2, f5*d5(1)/2, 1, d5(4));
    if isempty(index) == 0
       for i = 1:cell_numb
           sta = datarun21.stas.stas{index(i, 1)};
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
           end
           
    
           sta = datarun25.stas.stas{index(i, 2)};
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
           if sum(xy<=f1*d1(1)/4) == 0 && sum(xy>=3*f1*d1(1)/4) == 0 && isempty(xy) == 0
           sta_f5 = sta_f5 + sta_temp(xy(2)-f5*d5(1)/4+1:xy(2)+f5*d5(1)/4, xy(1)-f5*d5(2)/4+1:xy(1)+f5*d5(2)/4, :, :)/noise;
           end
           
       end       
    end
    sta_f_weighted{k, 1} = sta_f1;
    sta_f_weighted{k, 2} = sta_f5;
end


%% mean RF for different sampling


sta_c = cell(n, 2);
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
    sta_c_temp = zeros(d5(1), d5(2), d5(3), d5(4));
        for j = 1:d5(1)
            for k = 1:d5(2)
                for m = 1:d5(4)
                    sta_temp = sta(f5/2*(j-1)+1:f5/2*j, f5/2*(k-1)+1:f5/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c{i, 2} = sta_c_temp;
end



%% fit sta



fit_sta_c = cell(n, 2);
params_c = cell(n, 2);

for i = 5:5
     
%     temp_marks_sta = significant_stixels(sta_c{i, 1}, 'thresh', 5, 'time', 'max');
%     fit_ins = struct('sig_stixels', temp_marks_sta);

    fit_temp = fit_sta_sequence(sta_c{i, 1}, 'fit_instructions', fit_ins);
    fit_sta_c{i, 1} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x;
    params_temp(2) = fit_temp.center_point_y;
    params_temp(3) = fit_temp.center_sd_x;
    params_temp(4) = fit_temp.center_sd_y;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_c{i, 1} = params_temp;


%     fit_temp = fit_sta_sequence(sta_c{i, 2});
%     fit_sta_c{i, 2} = fit_temp;
%     params_temp = zeros(1, 8);
%     params_temp(1) = fit_temp.center_point_x;
%     params_temp(2) = fit_temp.center_point_y;
%     params_temp(3) = fit_temp.center_sd_x;
%     params_temp(4) = fit_temp.center_sd_y;
%     params_temp(5) = fit_temp.center_rotation_angle;
%     params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale;
%     params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale;
%     params_temp(8) = fit_temp.surround_amp_scale;
%     params_c{i, 2} = params_temp;
%     
    i

end






%% plot mean RF versus fitting

% ratio_all = zeros(n, 2);
dim = [d1; d5];
factor = [f1 f5];
cell_type_n = 2;

figure
for k = 1:2
   
sta = sta_c{cell_type_n, k};
d = dim(k, 1);

if cell_type_n == 1  
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
clear params_temp


    params_temp = params_c{cell_type_n, k};


[X1, Y1] = drawEllipse(params_temp(1:5));
[X2, Y2] = drawEllipse(params_temp([1, 2, 6, 7, 5]));
plot(X1, Y1, 'r')
plot(X2, Y2, 'b')
legend('center', 'surround')
hepes = {'without HEPES', 'with hepes'};
title([cell_type{cell_type_n} '   ' hepes{k}])




if cell_type_n == 1 
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

    fit_para1 = struct('center_point_x', para_temp(1), 'center_point_y', para_temp(2), ...
    'sd_x', para_temp(3), 'sd_y', para_temp(4), 'rotation_angle', para_temp(5), ...
    'x_dim', d, 'y_dim', d);

    fit_para2 = struct('center_point_x', para_temp(1), 'center_point_y', para_temp(2), ...
    'sd_x', para_temp(6), 'sd_y', para_temp(7), 'amp_scale', para_temp(8), ...
    'rotation_angle', para_temp(5), 'x_dim', d, 'y_dim', d);


fit1 = make_Gaussian_two_d(fit_para1);
fit2 = make_Gaussian_two_d(fit_para2);
if cell_type_n == 1
    fit = fit1 - fit2;
    comf = rf_com(fit);
    r = max(sta_temp(:))/max(fit(:));
    
%     fit_center = fit;
%     fit_center(fit_center < 0) = 0;
%     center = sum(fit_center(:));
%     fit_surround = fit;
%     fit_surround(fit_surround > 0) = 0;
%     surround = sum(fit_surround(:));
%     ratio = -surround/center
    ratio = para_temp(8)*(para_temp(6)/para_temp(3))^2;    
    ratio_all(cell_type_n, k) = ratio;

else
    fit = fit2 - fit1;
    comf = rf_com(-fit);
    r = min(sta_temp(:))/min(fit(:));
    
%     fit_center = fit;
%     fit_center(fit_center > 0) = 0;
%     center = sum(fit_center(:));
%     fit_surround = fit;
%     fit_surround(fit_surround < 0) = 0;
%     surround = sum(fit_surround(:));
%     ratio = -surround/center
    ratio = para_temp(8)*(para_temp(6)/para_temp(3))^2;    
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

subplot(2, 2, 2*k);
if cell_type_n == 1 
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
title([cell_type{cell_type_n} '    ' hepes{k}]);
end
%% real data of mean RF

figure

for nn = 1:5
    ct_idx = [1 2 3 4 5];
cell_type_n = ct_idx(nn);
distance = cell(2, 1);
STA = cell(2, 1);

for k = 1:2
    
sta = sta_c{cell_type_n, k};
d = dim(k, 1);

if cell_type_n == 1
    [~, b] = max(sta(:));
else
    [~, b] = min(sta(:));
end

z = ceil(b/d^2);
sta_1 = sta(:, :, 1, z);


if cell_type_n == 1 
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
    if cell_type_n == 1
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
    subplot(2, 3, nn)
    if cell_type_n == 1 
        plot(dis1, sta_temp, '.', 'MarkerSize', 7)
        hold on 
    else
        plot(dis1, -sta_temp, '.', 'MarkerSize', 7)
        hold on
    end
else
    if cell_type_n == 1 
        plot(dis1, sta_temp*r, '.', 'MarkerSize', 7, 'color', 'r')
    else
        plot(dis1, -sta_temp*r, '.', 'MarkerSize', 7, 'color', 'r')
    end
end


end
plot(dis1, zeros(1, d^2), 'm')
xlim([0 20])
legend('without HEPES', 'with HEPES')
title([cell_type{cell_type_n} ]);
end


%%

figure
x = [1:5]; y = ratio_all; 
bar(x,y)
set(gca,'xticklabel',cell_type)
legend('without HEPES', 'with HEPES'); %, 'location', 'northwest')
ylabel('surround/center volume ratio')