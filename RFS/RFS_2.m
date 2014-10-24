% load data
opt = struct('verbose',1,'load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun0 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data000/data000', opt);
datarun2 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data002/data002', opt);
datarun4 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data004/data004', opt);
datarun5 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data005/data005', opt);

% map ei
[cell_list_map40, ~] = map_ei(datarun4, datarun0);
[cell_list_map50, ~] = map_ei(datarun5, datarun0);

% map data000 and data004
cell_type = {'ON transient', 'ON brisk transient', 'OFF brisk transient', 'OFF sustained', 'DS'};
n = length(cell_type);
cell_id04 = cell(n, 1);
cell_idx04 = cell(n, 1);

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
end

cell_ids04 = struct('ON_transient', cell_id04(1), 'ON_brisk_transient', cell_id04(2), ...
    'OFF_brisk_transient', cell_id04(3),'OFF_sustained', cell_id04(4),'DS', cell_id04(5));
cell_idxs04 = struct('ON_transient', cell_idx04(1), 'ON_brisk_transient', cell_idx04(2), ...
    'OFF_brisk_transient', cell_idx04(3),'OFF_sustained', cell_idx04(4),'DS', cell_idx04(5));

% get parameters
load('RFS_TEMP.mat', 'datarun0', 'datarun4', 'datarun5');
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
        params_temp(j, 6, 2) = params4.surround_sd_scale*params0.center_sd_x;
        params_temp(j, 7, 2) = params4.surround_sd_scale*params0.center_sd_y;
        params_temp(j, 8, 2) = params4.surround_amp_scale;
        end
        
    end
    params04{i} = params_temp;
end

paramss04 = struct('ON_transient', params04(1), 'ON_brisk_transient', params04(2), ...
    'OFF_brisk_transient', params04(3),'OFF_sustained', params04(4),'DS', params04(5));

% plot mosaic

cell_type_n = 1;
params_temp = params04{cell_type_n};
cell_numb = size(params_temp, 1);
figure;
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
legend([h1 h2], 'NDF 0', 'NDF 4');
title([cell_type{cell_type_n} ' center']);

figure;
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
legend([h1 h2], 'NDF 0', 'NDF 4');
title([cell_type{cell_type_n} ' surround']);


% plot mean fitted RF

% % calculate mean RF
sta_f = cell(n, 2);
for k = 1:n;
    index = cell_idx04{k};
    cell_numb = size(cell_id04{k}, 1);
    d0 = size(datarun0.stas.stas{1});
    d4 = size(datarun4.stas.stas{1});
    f0 = 40;
    f4 = 30;
    sta_f0 = zeros(f0*d0(1)/2+1, f0*d0(1)/2+1, 1, d0(4));
    sta_f4 = zeros(f4*d4(1)/2+1, f4*d4(1)/2+1, 1, d4(4));
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
           sta_f0 = sta_f0 + sta_temp(xy(2)-f0*d0(1)/4:xy(2)+f0*d0(1)/4, xy(1)-f0*d0(2)/4:xy(1)+f0*d0(2)/4, :, :);
           cell_n0 = cell_n0 + 1;
           end
           
    
           sta = datarun4.stas.stas{index(i, 2)};
           sta_temp = zeros(f4*d4(1), f4*d4(2), d4(3), d4(4));
           for j = 1:d4(4)
               sta_temp(:, :, 1, j) = kron(sta(:, :, 1, j), ones(f4, f4));
           end
           xy = floor(f4*rf_com(sta));
           if sum(xy<=f0*d0(1)/4) == 0 && sum(xy>=3*f0*d0(1)/4) == 0 && isempty(xy) == 0
           sta_f4 = sta_f4 + sta_temp(xy(2)-f4*d4(1)/4:xy(2)+f4*d4(1)/4, xy(1)-f4*d4(2)/4:xy(1)+f4*d4(2)/4, :, :);
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
sta = sta_f{2, 2};
figure
for i = 1:d0(4)
    colormap gray
    imagesc(sta(:, :, 1, i));
    pause
end

% % fit mean RF
params04_f = cell(n, 2);
sta_fit = cell(n, 2);
fit_instructions = struct('fit_surround', true, 'fit_surround_sd_scale', true, 'fit_surround_amp_scale', ...
    true, 'fit_scale_one', false, 'fit_scale_two', false, 'fit_tau_one', false, 'fit_tau_two', false);
for i = 2:n
    for j = 1:2
        sta = sta_f{i, j};
        sta_temp = fit_sta(sta, fit_instructions);
        sta_fit{i, j} = sta_temp;
        params04f_temp = zeros(1, 8);
        params04f_temp(1) = sta_temp.center_point_x;
        params04f_temp(2) = sta_temp.center_point_y;
        params04f_temp(3) = sta_temp.center_sd_x;
        params04f_temp(4) = sta_temp.center_sd_y;
        params04f_temp(5) = sta_temp.center_rotation_angle;
        params04f_temp(6) = sta_temp.surround_sd_scale*sta_temp.center_sd_x;
        params04f_temp(7) = sta_temp.surround_sd_scale*sta_temp.center_sd_y;
        params04f_temp(8) = sta_temp.surround_amp_scale;
        params04_f{i, j} = params04f_temp;
        fprintf(['j = ' num2str(j)]);
    end
    fprintf(['i = ' num2str(i)])
end

% % plot fitted mean RF
cell_type_n = 2;
[X1, Y1] = drawEllipse(params04_f{cell_type_n, 1}(1:5));
[X2, Y2] = drawEllipse(params04_f{cell_type_n, 1}([1, 2, 6, 7, 5]));
[X3, Y3] = drawEllipse(params04_f{cell_type_n, 2}(1:5));
[X4, Y4] = drawEllipse(params04_f{cell_type_n, 2}([1, 2, 6, 7, 5]));



figure;
h1 = plot(X1, Y1, 'r');
hold on
plot(X2, Y2, 'r');
h2 = plot(X3, Y3, 'b');
plot(X4, Y4, 'b');

legend([h1 h2], 'NDF 0', 'NDF 4');


title(['cell id4 = ' num2str(cell_id4_ont(cell_type_n)) '   cell id2 = ' ...
    num2str(cell_id5_ont(cell_type_n)) '   cell n4 = ' num2str(cell_n4_ont(cell_type_n)) ...
    '   cell n2 = ' num2str(cell_n5_ont(cell_type_n))]);


% % fit mean RF in coarse resolution
sta_c = cell(n, 2);
for i = 1:n
    for j = 1:2
        sta_temp = sta_f{i, j};
        c_n = round(size(sta_f{1, 1}, 1)/20)*2-1;
        sta_c_temp = zeros(c_n, c_n, 1, d0(4));
        initial = ceil((size(sta_f{1, 1}, 1) - c_n*10)/2);
        for k = 1:c_n
            for m = 1:c_n
                for h = 1:d0(4)
                    sta_c_temp(k, m) = sta_temp((initial+(k-1)*10):(initial-1+k*10), (initial+(m-1)*10):(initial-1+m*10), 1, h)

