% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun0 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data000/data000', opt);
datarun4 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data004/data004', opt);

% map ei
[cell_list_map40, ~] = map_ei(datarun4, datarun0);

% map data000 and data004
cell_type = {'ON transient', 'ON brisk transient', 'OFF brisk transient', 'OFF sustained', 'DS'};
n = length(cell_type);
cell_id04 = cell(n, 1);
cell_idx04 = cell(n, 1);
cell_id_matched04 = [];
cell_idx_matched04 = [];

for i = 1:n
    [cell_idx4, ~, ~] = get_cell_indices(datarun4, cell_type{i});
    cell_id0 = cell_list_map40(cell_idx4);
    ept = 1-cellfun(@isempty, cell_id0);
    cell_idx4 = cell_idx4.*ept;
    cell_idx4(cell_idx4 == 0) = [];
    cell_id0 = cell2mat(cell_id0);
    cell_idx0 = get_cell_indices(datarun0, cell_id0);
    cell_id4 = datarun4.cell_ids(cell_idx4);
    id = [cell_id0' cell_id4'];
    idx = [cell_idx0' cell_idx4'];
    cell_id04{i} = id;
    cell_idx04{i} = idx;
    cell_id_matched04 = [cell_id_matched04; id];
    cell_idx_matched04 = [cell_idx_matched04; idx];

end


cell_ids04 = struct('ON_transient', cell_id04(1), 'ON_brisk_transient', cell_id04(2), ...
    'OFF_brisk_transient', cell_id04(3),'OFF_sustained', cell_id04(4),'DS', cell_id04(5));
cell_idxs04 = struct('ON_transient', cell_idx04(1), 'ON_brisk_transient', cell_idx04(2), ...
    'OFF_brisk_transient', cell_idx04(3),'OFF_sustained', cell_idx04(4),'DS', cell_idx04(5));

% get parameters
datarun0 = compute_sta_fits_sequence(datarun0, cell_id_matched04(:, 1), 'verbose', true);
datarun4 = compute_sta_fits_sequence(datarun4, cell_id_matched04(:, 2), 'verbose', true);

save_sta_fits(datarun0);
save_sta_fits(datarun4);

datarun0 = load_sta_fits(datarun0);
datarun4 = load_sta_fits(datarun4);

params04 = cell(n, 1);
for i = 1:n
    index = cell_idx04{i};
    cell_numb = size(index, 1);
    params_temp = zeros(cell_numb, 8, 2);
    for j = 1:cell_numb
        params0 = datarun0.matlab.sta_fits{index(j, 1)};
        params4 = datarun4.matlab.sta_fits{index(j, 2)};
        if isempty(params0) == 1 || isempty(params4) == 1;
           params_temp(j, :, 1) = zeros(1, 8);
           params_temp(j, :, 2) = zeros(1, 8);

        else  
        params_temp(j, 1, 1) = params0.center_point_x*4/3;
        params_temp(j, 2, 1) = params0.center_point_y*4/3;
        params_temp(j, 3, 1) = params0.center_sd_x*4/3;
        params_temp(j, 4, 1) = params0.center_sd_y*4/3;
        params_temp(j, 5, 1) = params0.center_rotation_angle;
        params_temp(j, 6, 1) = params0.surround_sd_scale*params0.center_sd_x*4/3;
        params_temp(j, 7, 1) = params0.surround_sd_scale*params0.center_sd_y*4/3;
        params_temp(j, 8, 1) = params0.surround_amp_scale;
        
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
    params04{i} = params_temp;
end

paramss04 = struct('ON_transient', params04(1), 'ON_brisk_transient', params04(2), ...
    'OFF_brisk_transient', params04(3),'OFF_sustained', params04(4),'DS', params04(5));

%% plot mosaic

cell_type_n = 5;
params_temp = params04{cell_type_n};
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





%% plot mean fitted RF

% % calculate mean RF
d0 = size(datarun0.stas.stas{1});
d4 = size(datarun4.stas.stas{1});
f0 = 40;
f4 = 30;
sta_f = cell(n, 2);
for k = 1:n;
    index = cell_idx04{k};
    cell_numb = size(cell_id04{k}, 1);
    sta_f0 = zeros(f0*d0(1)/2, f0*d0(1)/2, 1, d0(4));
    sta_f4 = zeros(f4*d4(1)/2, f4*d4(1)/2, 1, d4(4));
    cell_n0 = 0;
    cell_n4 = 0;
    if isempty(index) == 0
       for i = 1:cell_numb
           sta = datarun0.stas.stas{index(i, 1)};
           sta_temp = zeros(f0*d0(1), f0*d0(2), d0(3), d0(4));
           for j = 1:d0(4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f0, f0));
           end
           xy = floor(f0*rf_com(sta));
           if sum(xy<=f0*d0(1)/4) == 0 && sum(xy>=3*f0*d0(1)/4) == 0 && isempty(xy) == 0
           sta_f0 = sta_f0 + sta_temp(xy(2)-f0*d0(1)/4+1:xy(2)+f0*d0(1)/4, xy(1)-f0*d0(2)/4+1:xy(1)+f0*d0(2)/4, :, :);
           cell_n0 = cell_n0 + 1;
           end
           
    
           sta = datarun4.stas.stas{index(i, 2)};
           sta_temp = zeros(f4*d4(1), f4*d4(2), d4(3), d4(4));
           for j = 1:d4(4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f4, f4));
           end
           xy = floor(f4*rf_com(sta));
           if sum(xy<=f0*d0(1)/4) == 0 && sum(xy>=3*f0*d0(1)/4) == 0 && isempty(xy) == 0
           sta_f4 = sta_f4 + sta_temp(xy(2)-f4*d4(1)/4+1:xy(2)+f4*d4(1)/4, xy(1)-f4*d4(2)/4+1:xy(1)+f4*d4(2)/4, :, :);
           cell_n4 = cell_n4 + 1;
           end
           
       end       
    end
    sta_f0 = sta_f0/cell_n0;
    sta_f4 = sta_f4/cell_n4;
    sta_f{k, 1} = sta_f0;
    sta_f{k, 2} = sta_f4;
end




% % plot mean RF

sta = sta_f{1, 1};

sta = sta - min(sta(:));
sta = sta/max(sta(:));

figure
for i = 1:d0(4)
    colormap gray
    imagesc(sta(:, :, 1, i));
    pause
end


%% fit mean RF in coarse resolution
sta_c = cell(n, 2);
fit_sta_c = cell(n, 2);
params_c = cell(n, 2);
d = [40; 30];
for i = 1:n
        sta = sta_f{i, 1};
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
    fit_temp = fit_sta_sequence(sta_c_temp);
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


    sta = sta_f{i, 2};
    sta_c_temp = zeros(d4(1), d4(2), d4(3), d4(4));
        for j = 1:d4(1)
            for k = 1:d4(2)
                for m = 1:d4(4)
                    sta_temp = sta(f4/2*(j-1)+1:f4/2*j, f4/2*(k-1)+1:f4/2*k, 1, m);
                    sta_c_temp(j, k, 1, m) = mean(sta_temp(:));
                end
            end
        end
    sta_c{i, 2} = sta_c_temp;
    fit_temp = fit_sta_sequence(sta_c_temp);
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
    

end



%% plot RF
cell_type_n = 5;
params_temp0 = params_c{cell_type_n, 1};
params_temp4 = params_c{cell_type_n, 2};

figure;
[X1, Y1] = drawEllipse(params_temp0(1:5));
[X2, Y2] = drawEllipse(params_temp4(1:5));
[X3, Y3] = drawEllipse(params_temp0([1, 2, 6, 7, 5]));
[X4, Y4] = drawEllipse(params_temp4([1, 2, 6, 7, 5]));

h1 = plot(X1, Y1, 'r');
hold on
h2 = plot(X2, Y2, 'b');
plot(X3, Y3, 'r');
plot(X4, Y4, 'b');

axis([0 40 0 40]);
legend([h1 h2], 'NDF 4', 'NDF 0');
title([cell_type{cell_type_n}]);
