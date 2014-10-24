% load data
opt = struct('verbose',1,'load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun0 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data000/data000', opt);
datarun2 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data002/data002', opt);
datarun4 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data004/data004', opt);
datarun5 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data005/data005', opt);

% fit sta
fit_instructions = struct('fit_surround', true, 'fit_surround_sd_scale', true, 'fit_surround_amp_scale', ...
    true, 'fit_scale_one', false, 'fit_scale_two', false, 'fit_tau_one', false, 'fit_tau_two', false);
datarun0 = compute_sta_fits(datarun0, 'all', 'fit_instructions', fit_instructions, 'verbose', true);
datarun2 = compute_sta_fits(datarun2, 'all', 'fit_instructions', fit_instructions, 'verbose', true);
datarun4 = compute_sta_fits(datarun4, 'all', 'fit_instructions', fit_instructions, 'verbose', true);
datarun5 = compute_sta_fits(datarun5, 'all', 'fit_instructions', fit_instructions, 'verbose', true);

% map ei
[cell_list_map02, ~] = map_ei(datarun0, datarun2);
[cell_list_map04, ~] = map_ei(datarun0, datarun4);
[cell_list_map05, ~] = map_ei(datarun0, datarun5);
[cell_list_map24, ~] = map_ei(datarun2, datarun4);
[cell_list_map25, ~] = map_ei(datarun2, datarun5);
[cell_list_map45, ~] = map_ei(datarun4, datarun5);

% map data000 and data002
[cell_n0_offs, ~, ~] = get_cell_indices(datarun0, 'OFF sustained');
cell_id2_offs = cell_list_map02(cell_n0_offs);
ept = 1-cellfun(@isempty, cell_id2_offs);
cell_n0_offs = cell_n0_offs.*ept;
cell_n0_offs(cell_n0_offs == 0) = [];
cell_id2_offs = cell2mat(cell_id2_offs);
cell_n2_offs = get_cell_indices(datarun2, cell_id2_offs);
cell_id0_offs = datarun0.cell_ids(cell_n0_offs);

% map data000 and data004
cell_id4_offs = cell_list_map04(cell_n0_offs);
ept = 1-cellfun(@isempty, cell_id4_offs);
cell_n0_offs = cell_n0_offs.*ept;
cell_n0_offs(cell_n0_offs == 0) = [];
cell_n2_offs = cell_n2_offs.*ept;
cell_n2_offs(cell_n2_offs == 0) = [];
cell_id2_offs = cell_id2_offs.*ept;
cell_id2_offs(cell_id2_offs == 0) = [];
cell_id4_offs = cell2mat(cell_id4_offs);
cell_n4_offs = get_cell_indices(datarun4, cell_id4_offs);
cell_id0_offs = datarun0.cell_ids(cell_n0_offs);

% map data000 and data005
cell_id5_offs = cell_list_map05(cell_n0_offs);
ept = 1-cellfun(@isempty, cell_id5_offs);
cell_n0_offs = cell_n0_offs.*ept;
cell_n0_offs(cell_n0_offs == 0) = [];
cell_n2_offs = cell_n2_offs.*ept;
cell_n2_offs(cell_n2_offs == 0) = [];
cell_n4_offs = cell_n4_offs.*ept;
cell_n4_offs(cell_n4_offs == 0) = [];
cell_id2_offs = cell_id2_offs.*ept;
cell_id2_offs(cell_id2_offs == 0) = [];
cell_id4_offs = cell_id4_offs.*ept;
cell_id4_offs(cell_id4_offs == 0) = [];
cell_id5_offs = cell2mat(cell_id5_offs);
cell_n5_offs = get_cell_indices(datarun5, cell_id5_offs);
cell_id0_offs = datarun0.cell_ids(cell_n0_offs);

% compile cell id
cell_id_offs = zeros(length(cell_id0_offs), 4);
cell_id_offs(:, 1) = cell_id0_offs;
cell_id_offs(:, 2) = cell_id2_offs;
cell_id_offs(:, 3) = cell_id4_offs;
cell_id_offs(:, 4) = cell_id5_offs;

% compile cell index
cell_n_offs = zeros(length(cell_n0_offs), 4);
cell_n_offs(:, 1) = cell_n0_offs;
cell_n_offs(:, 2) = cell_n2_offs;
cell_n_offs(:, 3) = cell_n4_offs;
cell_n_offs(:, 4) = cell_n5_offs;

% get parameters 
params0_DS = zeros(length(cell_n0_offs), 8);
params2_DS = zeros(length(cell_n2_offs), 8);
params4_DS = zeros(length(cell_n4_offs), 8);
params5_DS = zeros(length(cell_n5_offs), 8);

for i = 1:length(cell_n0_offs)
    params0 = datarun0.matlab.sta_fits{cell_n0_offs(i)};
    if isempty(params0)
        params0_DS(i, :) = zeros(1, 8);
    else    
    params0_DS(i, 1) = params0.center_point_x*4/3;
    params0_DS(i, 2) = params0.center_point_y*4/3;
    params0_DS(i, 3) = params0.center_sd_x*4/3;
    params0_DS(i, 4) = params0.center_sd_y*4/3;
    params0_DS(i, 5) = params0.center_rotation_angle;
    params0_DS(i, 6) = params0.surround_sd_scale*params0.center_sd_x*4/3;
    params0_DS(i, 7) = params0.surround_sd_scale*params0.center_sd_y*4/3;
    params0_DS(i, 8) = params0.surround_amp_scale;
    end
    
    params2 = datarun2.matlab.sta_fits{cell_n2_offs(i)};
    if isempty(params2)
        params2_DS(i, :) = zeros(1, 8);
    else
    params2_DS(i, 1) = params2.center_point_x;
    params2_DS(i, 2) = params2.center_point_y;
    params2_DS(i, 3) = params2.center_sd_x;
    params2_DS(i, 4) = params2.center_sd_y;
    params2_DS(i, 5) = params2.center_rotation_angle;
    params2_DS(i, 6) = params2.surround_sd_scale*params2.center_sd_x;
    params2_DS(i, 7) = params2.surround_sd_scale*params2.center_sd_y;
    params2_DS(i, 8) = params2.surround_amp_scale;
    end
    
     params4 = datarun4.matlab.sta_fits{cell_n4_offs(i)};
    if isempty(params4)
        params4_DS(i, :) = zeros(1, 8);
    else
    params4_DS(i, 1) = params4.center_point_x;
    params4_DS(i, 2) = params4.center_point_y;
    params4_DS(i, 3) = params4.center_sd_x;
    params4_DS(i, 4) = params4.center_sd_y;
    params4_DS(i, 5) = params4.center_rotation_angle;
    params4_DS(i, 6) = params4.surround_sd_scale*params4.center_sd_x;
    params4_DS(i, 7) = params4.surround_sd_scale*params4.center_sd_y;
    params4_DS(i, 8) = params4.surround_amp_scale;
    end
    
     params5 = datarun5.matlab.sta_fits{cell_n5_offs(i)};
    if isempty(params5)
        params5_DS(i, :) = zeros(1, 8);
    else
    params5_DS(i, 1) = params5.center_point_x/3;
    params5_DS(i, 2) = params5.center_point_y/3;
    params5_DS(i, 3) = params5.center_sd_x/3;
    params5_DS(i, 4) = params5.center_sd_y/3;
    params5_DS(i, 5) = params5.center_rotation_angle;
    params5_DS(i, 6) = params5.surround_sd_scale*params5.center_sd_x/3;
    params5_DS(i, 7) = params5.surround_sd_scale*params5.center_sd_y/3;
    params5_DS(i, 8) = params5.surround_amp_scale;
    end
end

% compare one cell across 2 light levels

% pair_n = 10;
% [X1, Y1] = drawEllipse(params0_ont(pair_n, 1:5));
% [X2, Y2] = drawEllipse(params0_ont(pair_n, [1, 2, 6, 7, 5]));
% [X3, Y3] = drawEllipse(params2_ont(pair_n, 1:5));
% [X4, Y4] = drawEllipse(params2_ont(pair_n, [1, 2, 6, 7, 5]));
% [X5, Y5] = drawEllipse(params4_ont(pair_n, 1:5));
% [X6, Y6] = drawEllipse(params4_ont(pair_n, [1, 2, 6, 7, 5]));
% [X7, Y7] = drawEllipse(params5_ont(pair_n, 1:5));
% [X8, Y8] = drawEllipse(params5_ont(pair_n, [1, 2, 6, 7, 5]));
% 
% 
% 
% figure;
% h1 = plot(X5, Y5, 'r');
% hold on
% plot(X6, Y6, 'r');
% h2 = plot(X3, Y3, 'b');
% plot(X4, Y4, 'b');
% 
% legend([h1 h2], 'NDF 0', 'NDF 2');
% title(['cell id4 = ' num2str(cell_id4_ont(pair_n)) '   cell id2 = ' ...
%     num2str(cell_id5_ont(pair_n)) '   cell n4 = ' num2str(cell_n4_ont(pair_n)) ...
%     '   cell n2 = ' num2str(cell_n5_ont(pair_n))]);

% calculate the average RF size
n = size(params0_offbt, 1);
params_offbt = zeros(n, 10, 4);
params_offbt(:, 1:8, 1) = params0_offbt;
params_offbt(:, 1:8, 2) = params2_offbt;
params_offbt(:, 1:8, 3) = params4_offbt;
params_offbt(:, 1:8, 4) = params5_offbt;

temp = params_offbt;
n = size(temp, 1);
for i = 1:n
    for j = 1:4
        temp(i, 9, j) = temp(i, 3, j)*temp(i, 4, j)*pi;
        temp(i, 10, j) = temp(i, 6, j)*temp(i, 7, j)*pi;
    end
end
params_offbt = temp;


temp = params_offbt;
center_a = zeros(1, 4);
error_c = zeros(1, 4);
surround_a = zeros(1, 4);
error_s = zeros(1, 4);
for i = 1:4
    center_a(i) = mean(temp(:, 9, i));
    error_c(i) = std(temp(:, 9, i));
    surround_a(i) = mean(temp(:, 10, i));
    error_s(i) = std(temp(:, 10, i));
end
figure
errorbar(center_a, error_c, 'b')
hold on
errorbar(surround_a, error_s, 'r')
ylabel('area')
xlabel('background luminance')
title('DS')
legend('center', 'surround')

% plot mosaic
% % map data000 and data004
[cell_n0_offs, ~, ~] = get_cell_indices(datarun0, 'OFF sustained');
cell_id4_offs = cell_list_map02(cell_n0_offs);
ept = 1-cellfun(@isempty, cell_id4_offs);
cell_n0_offs = cell_n0_offs.*ept;
cell_n0_offs(cell_n0_offs == 0) = [];
cell_id4_offs = cell2mat(cell_id4_offs);
cell_n2_offs = get_cell_indices(datarun2, cell_id4_offs);
cell_id0_offs = datarun0.cell_ids(cell_n0_offs);
