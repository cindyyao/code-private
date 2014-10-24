% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun{1} = load_data('/Analysis/xyao/2013-05-30-0/data005/data005', opt);
datarun{2} = load_data('/Analysis/xyao/2013-05-30-0/data008/data008', opt);

% get cell ids and indices
cell_type = {'ON transient', 'OFF sustained2', 'OFF sustained', 'OFF slow', 'OFF transient'};
n = length(cell_type);
cell_id = cell(2, 1);
cell_idx = cell(2, 1);
for j = 1:2
    cell_id_temp = cell(n, 1); 
    cell_idx_temp = cell(n, 1);
    for i = 1:n
        cell_idx_temp{i} = get_cell_indices(datarun{j}, cell_type{i});
        cell_id_temp{i} = datarun{j}.cell_ids(cell_idx_temp{i});
    end
    cell_id{j} = cell_id_temp;
    cell_idx{j} = cell_idx_temp;
end



%% weighted addition of RF by SNR

sta_f_weighted = cell(2, 1);
sta_f_n = cell(2, 1);

for ll = 1:2
    d = size(datarun{ll}.stas.stas{1});
    f = 1200/d(1);
    sta_f_weighted_ll = cell(n, 1);
    sta_f_n_ll = zeros(n, 1);

for k = 1:n;
    index = cell_idx{ll}{k};
    cell_numb = size(cell_id{ll}{k}, 2);
    sta_f = zeros(f*d(1)/2, f*d(1)/2, 1, d(4));
    sta_f_n_temp = 0;
    
    if isempty(index) == 0
       for i = 1:cell_numb
           
           sta = datarun{ll}.stas.stas{index(i)};
           if cell_type{k}(2) == 'N' 
              [~, b] = max(sta(:));
           else
              [~, b] = min(sta(:));
           end
           z = ceil(b/(d(1)*d(2)));
           sta_1 = sta(:, :, 1, z);
           noise = robust_std(sta_1(:), 3);
           sta_temp = zeros(f*d(1), f*d(2), d(3), d(4));
           for j = 1:d(4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f, f));
           end
           xy = floor(f*rf_com(sta));
           if sum(xy<=f*d(1)/4) == 0 && sum(xy>=3*f*d(1)/4) == 0 && isempty(xy) == 0
           sta_f = sta_f + sta_temp(xy(2)-f*d(1)/4+1:xy(2)+f*d(1)/4, xy(1)-f*d(2)/4+1:xy(1)+f*d(2)/4, :, :)/noise;
           sta_f_n_temp = sta_f_n_temp + 1;
           end
           
       end       
    end
    
    sta_f_weighted_ll{k} = sta_f;
    sta_f_n_ll(k) = sta_f_n_temp;
end
sta_f_weighted{ll} = sta_f_weighted_ll;
sta_f_n{ll} = sta_f_n_ll;
end







%% mean RF for different sampling


sta_c = cell(2, 1);
for ll = 1:2
    d = size(datarun{ll}.stas.stas{1});
    f = 1200/d(1);

    sta_c_ll = cell(n, 1);
for i = 1:n
        sta = sta_f_weighted{ll}{i};
        sta_c_temp = zeros(d(1), d(2), d(3), d(4));
        for j = 1:d(1)
            for k = 1:d(2)
                for m = 1:d(4)
                    sta_temp = sta(f/2*(j-1)+1:f/2*j, f/2*(k-1)+1:f/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c_ll{i} = sta_c_temp;
end

sta_c{ll} = sta_c_ll;

end




%% fit sta


fit_sta_c = cell(2, 1);
params_c = cell(2, 1);

for ll = 1:2
    d = size(datarun{ll}.stas.stas{1});
    F = 40/d(1);

    fit_sta_c_ll = cell(n, 1);
    params_c_ll = cell(n, 1);
for i = 1:n
    
    temp_marks_sta = significant_stixels(sta_c{ll}{i}, 'thresh', 5, 'time', 'max');
    fit_ins = struct('sig_stixels', temp_marks_sta);

    fit_temp = fit_sta_sequence(sta_c{ll}{i}, 'fit_instructions', fit_ins);
    fit_sta_c_ll{i} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x*F;
    params_temp(2) = fit_temp.center_point_y*F;
    params_temp(3) = fit_temp.center_sd_x*F;
    params_temp(4) = fit_temp.center_sd_y*F;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale*F;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale*F;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_c_ll{i} = params_temp;
    
   
end
fit_sta_c{ll} = fit_sta_c_ll;
params_c{ll} = params_c_ll;

ll
end



save('gnat1_nomap.mat', 'fit_sta_c', 'params_c', 'sta_f_weighted', 'sta_c')











%% plot mean RF versus fitting

% ratio_all = zeros(n, 2);


cell_type_n = 3;

figure;
for k = 1:2;

    d = size(datarun{k}.stas.stas{1}, 1);
    F = 40/d;

sta = sta_c{k}{cell_type_n};

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
clear params_temp

params_temp = params_c{k}{cell_type_n}(1:7)/F;
params_temp(5) = params_c{k}{cell_type_n}(5);


[X1, Y1] = drawEllipse(params_temp(1:5));
[X2, Y2] = drawEllipse(params_temp([1, 2, 6, 7, 5]));
plot(X1, Y1, 'r')
plot(X2, Y2, 'b')
legend('center', 'surround');
NDF = [1 0];
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

dis1 = reshape(dis, 1, d^2);
sta_temp = reshape(sta_1, 1, d^2);

para_temp = params_c{k}{cell_type_n};
fit_para1 = struct('center_point_x', para_temp(1)/F, 'center_point_y', para_temp(2)/F, ...
'sd_x', para_temp(3)/F, 'sd_y', para_temp(4)/F, 'rotation_angle', para_temp(5), ...
'x_dim', d, 'y_dim', d);

fit_para2 = struct('center_point_x', para_temp(1)/F, 'center_point_y', para_temp(2)/F, ...
'sd_x', para_temp(6)/F, 'sd_y', para_temp(7)/F, 'amp_scale', para_temp(8), ...
'rotation_angle', para_temp(5), 'x_dim', d, 'y_dim', d);


fit1 = make_Gaussian_two_d(fit_para1);
fit2 = make_Gaussian_two_d(fit_para2);
if cell_type{cell_type_n}(2) == 'N' 
    fit = fit1 - fit2;
    comf = rf_com(fit);
    r = max(sta_temp(:))/max(fit(:));
    ratio = para_temp(8)*(para_temp(6)/para_temp(3))^2;
    ratio_all(cell_type_n, k) = ratio;
else
    fit = fit2 - fit1;
    comf = rf_com(-fit);
    r = min(sta_temp(:))/min(fit(:));
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
legend('data', 'fit');
title([cell_type{cell_type_n} '   NDF ' num2str(NDF(k))]);
end




%% real data of mean RF

figure

for cell_type_n = 1:n
distance = cell(2, 1);
STA = cell(2, 1);

for k = 1:2
    
sta = sta_c{k}{cell_type_n};
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
    subplot(2, 3, cell_type_n)
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
legend('NDF 1', 'NDF 0')
title([cell_type{cell_type_n} ]);
end

%%
figure
x = [1:5]; y = [ratio_all(:, 2) ratio_all(:, 1)];
bar(x,y)
set(gca,'xticklabel',cell_type)
legend('NDF 0', 'NDF 1', 'location', 'northwest')
