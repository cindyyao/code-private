% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun1 = load_data('/Analysis/xyao/2013-02-13-0/data001/data001', opt);
datarun3 = load_data('/Analysis/xyao/2013-02-13-0/data003/data003', opt);
datarun5 = load_data('/Analysis/xyao/2013-02-13-0/data005/data005', opt);
datarun10 = load_data('/Analysis/xyao/2013-02-13-0/data010/data010', opt);


% map ei
[cell_list_map510, ~] = map_ei(datarun5, datarun10);
[cell_list_map51, ~] = map_ei(datarun5, datarun1);
[cell_list_map53, ~] = map_ei(datarun5, datarun3);


% map data000 and data004
cell_type = {'OFF brisk transient', 'OFF brisk transient large', 'OFF transient'};
n = length(cell_type);
cell_id510 = cell(n, 1);
cell_idx510 = cell(n, 1);
cell_id_matched510 = [];
cell_idx_matched510 = [];

for i = 1:n
    [cell_idx5, ~, ~] = get_cell_indices(datarun5, cell_type{i});
    cell_id10 = cell_list_map510(cell_idx5);
    ept = 1-cellfun(@isempty, cell_id10);
    cell_idx5 = cell_idx5.*ept;
    cell_idx5(cell_idx5 == 0) = [];
    cell_id10 = cell2mat(cell_id10);
    cell_idx10 = get_cell_indices(datarun10, cell_id10);
    cell_id5 = datarun5.cell_ids(cell_idx5);
    id = [cell_id10' cell_id5'];
    idx = [cell_idx10' cell_idx5'];
    cell_id510{i} = id;
    cell_idx510{i} = idx;
    cell_id_matched510 = [cell_id_matched510; id];
    cell_idx_matched510 = [cell_idx_matched510; idx];

end




%% weighted addition of RF by SNR

d5 = size(datarun5.stas.stas{1});
d10 = size(datarun10.stas.stas{1});
f5 = 30;
f10 = 30;


sta_f_weighted = cell(n, 2);
sta_f_n = zeros(n, 2);

for k = 1:n;
    index = cell_idx510{k};
    cell_numb = size(cell_id510{k}, 1);
    sta_f5 = zeros(f5*d5(1)/2, f5*d5(1)/2, 1, d5(4));
    sta_f10 = zeros(f10*d10(1)/2, f10*d10(1)/2, 1, d10(4));
    
    sta_f_n5 = 0;
    sta_f_n10 = 0;
    
    if isempty(index) == 0
       for i = 1:cell_numb
           sta5 = datarun5.stas.stas{index(i, 1)};
           
           [~, b5] = min(sta5(:));
           
           z5 = ceil(b5/(d5(1)*d5(2)));
           sta5_1 = sta5(:, :, 1, z);
           noise5 = robust_std(sta5_1(:), 3);
           sta5_temp = zeros(f5*d5(1), f5*d5(2), d5(3), d5(4));
           for j = 1:d5(4)
               sta5_temp(:, :, 1, j) = kron(sta5(:, :, 1, j), ones(f5, f5));
           end
           xy5 = floor(f5*rf_com(sta5));
           
           sta10 = datarun10.stas.stas{index(i, 2)};
          
           [~, b10] = min(sta10(:));
           z10 = ceil(b10/(d10(1)*d10(2)));
           sta10_1 = sta10(:, :, 1, z);
           noise10 = robust_std(sta10_1(:), 3);
           sta10_temp = zeros(f10*d10(1), f10*d10(2), d10(3), d10(4));
           for j = 1:d10(4)
               sta10_temp(:, :, 1, j) = kron(sta10(:, :, 1, j), ones(f10, f10));
           end
           xy10 = floor(f10*rf_com(sta10));
           
           
           
           
           
           if sum(xy5<=f5*d5(1)/4) == 0 && sum(xy5>=3*f5*d5(1)/4) == 0 && ...
                   isempty(xy5) == 0 && sum(xy10<=f10*d10(1)/4) == 0 && ...
                   sum(xy10>=3*f10*d10(1)/4) == 0 && isempty(xy10) == 0 
           sta_f5 = sta_f5 + sta5_temp(xy5(2)-f5*d5(1)/4+1:xy5(2)+f5*d5(1)/4, ...
               xy5(1)-f5*d5(2)/4+1:xy5(1)+f5*d5(2)/4, :, :)/noise5;
           sta_f_n5 = sta_f_n5 + 1;
           sta_f10 = sta_f10 + sta10_temp(xy10(2)-f10*d10(1)/4+1:xy10(2)+f10*d10(1)/4, ...
               xy10(1)-f10*d10(2)/4+1:xy10(1)+f10*d10(2)/4, :, :)/noise10;
           sta_f_n10 = sta_f_n10 + 1;
           end
           
       end       
    end
    
    sta_f_weighted{k, 1} = sta_f5;
    sta_f_weighted{k, 2} = sta_f10;
    sta_f_n(k, 1) = sta_f_n5;
    sta_f_n(k, 2) = sta_f_n10;
end













%% mean RF for different sampling


sta_c = cell(n, 2);
for i = 1:n
        sta = sta_f_weighted{i, 1};
        sta_c_temp = zeros(d5(1), d5(2), d5(3), d5(4));
        for j = 1:d5(1)
            for k = 1:d5(2)
                for m = 1:d5(4)
                    sta_temp = sta(f5/2*(j-1)+1:f5/2*j, f5/2*(k-1)+1:f5/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c{i, 1} = sta_c_temp;
    
    sta = sta_f_weighted{i, 2};
    sta_c_temp = zeros(d10(1), d10(2), d10(3), d10(4));
        for j = 1:d10(1)
            for k = 1:d10(2)
                for m = 1:d10(4)
                    sta_temp = sta(f10/2*(j-1)+1:f10/2*j, f10/2*(k-1)+1:f10/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c{i, 2} = sta_c_temp;
    
  
end








%% fit sta


fit_sta_c = cell(n, 2);
params_c = cell(n, 1);



for i = 1:n
    
    fit_temp = fit_sta_sequence(sta_c{i, 1});
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

    
    
    
    i
end


save('RFS130213_510.mat', 'fit_sta_c', 'params_c', 'sta_f_weighted',...
    'sta_c', 'n', 'd5', 'd10', 'f5', 'f10', 'cell_type')











%% plot mean RF versus fitting

dim = [d5; d10];
factor = [f5 f10];
cell_type_n = 2;
ratio_all = zeros(n, 2);

for k = 1:2;
    
figure;
sta = sta_c{cell_type_n, k};
d = dim(k, 1);


[~, b] = min(sta(:));


z = ceil(b/d^2);
sta_1 = sta(:, :, 1, z);



subplot(1, 2, 1);
im = norm_image(sta_1);
image(im)
hold on
clear params_temp

    params_temp = params_c{cell_type_n, k};


[X1, Y1] = drawEllipse(params_temp(1:5));
[X2, Y2] = drawEllipse(params_temp([1, 2, 6, 7, 5]));
plot(X1, Y1, 'r')
plot(X2, Y2, 'b')
legend('center', 'surround');
TP = [6 2];
title([cell_type{cell_type_n} ' TP = ' num2str(TP(k))])





    com = rf_com(-sta_1);



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

    fit = fit2 - fit1;
    [~, I] = min(fit(:));
    comf(1) = ceil(I/d);
    comf(2) = I - floor(I/d)*d;
    
    r = min(sta_temp(:))/min(fit(:));
    fit_center = fit;
    fit_center(fit_center > 0) = 0;
    center = sum(fit_center(:));
    fit_surround = fit;
    fit_surround(fit_surround < 0) = 0;
    surround = sum(fit_surround(:));
    ratio = -surround/center
    ratio_all(cell_type_n, k) = ratio;





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

    plot(dis1, -sta_temp, '.', 'MarkerSize', 7)
    hold on
    plot(dis1f, -fit1*r, '.', 'MarkerSize', 7, 'Color', 'r')
    plot(dis1, zeros(1, d^2), 'm')

xlim([0 20])
legend('data', 'fit');
title([cell_type{cell_type_n} '   TP ' num2str(TP(k))]);
end









%% real data of mean RF

cell_type_n = 3;
distance = cell(2, 1);
STA = cell(1, 1);

for k = 1:2
    
sta = sta_c{cell_type_n, k};
d = dim(k, 1);


    [~, b] = min(sta(:));


z = ceil(b/d^2);
sta_1 = sta(:, :, 1, z);



    com = rf_com(-sta_1);



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
   
        r = min(sta_temp1(:))/min(sta_1(:));
   
end

sta_temp = reshape(sta_1, 1, d^2);

if k == 1
    sta_temp1 = sta_temp;
end


if k == 1
    figure

        plot(dis1, -sta_temp, '.', 'MarkerSize', 7)
        hold on
    
else 
   
        plot(dis1, -sta_temp*r, '.', 'MarkerSize', 7, 'color', 'r')
   

        
end


end

plot(dis1, zeros(1, d^2), 'm')

xlim([0 20])
legend('TP 6', 'TP 2')
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
