opt = struct('verbose',1,'load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun0 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data000/data000', opt);
datarun2 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data002/data002', opt);
datarun4 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data004/data004', opt);
datarun5 = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data005/data005', opt);
[cell_list_map02, ~] = map_ei(datarun0, datarun2);
[cell_list_map04, ~] = map_ei(datarun0, datarun4);
[cell_list_map05, ~] = map_ei(datarun0, datarun5);
[cell_list_map24, ~] = map_ei(datarun2, datarun4);
[cell_list_map25, ~] = map_ei(datarun2, datarun5);
[cell_list_map45, ~] = map_ei(datarun4, datarun5);

fit_instructions = struct('fit_surround', true, 'fit_surround_sd_scale', true, 'fit_surround_amp_scale', true);
datarun0 = compute_sta_fits(datarun0, 'all', 'fit_instructions', fit_instructions, 'verbose', true);
datarun2 = compute_sta_fits(datarun2, 'all', 'fit_instructions', fit_instructions, 'verbose', true);
datarun4 = compute_sta_fits(datarun4, 'all', 'fit_instructions', fit_instructions, 'verbose', true);
datarun5 = compute_sta_fits(datarun5, 'all', 'fit_instructions', fit_instructions, 'verbose', true);

[cell_n0_ont, ~, ~] = get_cell_indices(datarun0, 'OFF sustained');
cell_id2_ont = cell_list_map02(cell_n0_ont);
ept = 1-cellfun(@isempty, cell_id2_ont);
cell_n0_ont = cell_n0_ont.*ept;
cell_n0_ont(cell_n0_ont == 0) = [];
cell_id2_ont = cell2mat(cell_id2_ont);
cell_n2_ont = get_cell_indices(datarun2, cell_id2_ont);
cell_id0_ont = datarun0.cell_ids(cell_n0_ont);


params0_ont = zeros(length(cell_n0_ont), 8);
params2_ont = zeros(length(cell_n2_ont), 8);
for i = 1:length(cell_n0_ont)
    params0 = datarun0.matlab.sta_fits{cell_n0_ont(i)};
    if isempty(params0)
        params0_ont(i, :) = [];
    else    
    params0_ont(i, 1) = params0.center_point_x*4/3;
    params0_ont(i, 2) = params0.center_point_y*4/3;
    params0_ont(i, 3) = params0.center_sd_x*4/3;
    params0_ont(i, 4) = params0.center_sd_y*4/3;
    params0_ont(i, 5) = params0.center_rotation_angle;
    params0_ont(i, 6) = params0.surround_sd_scale*params0.center_sd_x*4/3;
    params0_ont(i, 7) = params0.surround_sd_scale*params0.center_sd_y*4/3;
    params0_ont(i, 8) = params0.surround_amp_scale;
    params0_ont(i, 10) = params0.surround_sd_scale;
    end
    
    params2 = datarun2.matlab.sta_fits{cell_n2_ont(i)};
    if isempty(params2)
        params2_ont(i, :) = [];
    else
    params2_ont(i, 1) = params2.center_point_x;
    params2_ont(i, 2) = params2.center_point_y;
    params2_ont(i, 3) = params2.center_sd_x;
    params2_ont(i, 4) = params2.center_sd_y;
    params2_ont(i, 5) = params2.center_rotation_angle;
    params2_ont(i, 6) = params2.surround_sd_scale*params2.center_sd_x;
    params2_ont(i, 7) = params2.surround_sd_scale*params2.center_sd_y;
    params2_ont(i, 8) = params2.surround_amp_scale;
    end
end

pair_n = 5;
[X1, Y1] = drawEllipse(params0_ont(pair_n, 1:5));
[X2, Y2] = drawEllipse(params0_ont(pair_n, [1, 2, 6, 7, 5]));
[X3, Y3] = drawEllipse(params2_ont(pair_n, 1:5));
[X4, Y4] = drawEllipse(params2_ont(pair_n, [1, 2, 6, 7, 5]));

figure;
h1 = plot(X1, Y1, 'r');
hold on
plot(X2, Y2, 'r');
h2 = plot(X3, Y3, 'b');
plot(X4, Y4, 'b');

legend([h1 h2], 'NDF 0', 'NDF 2');
title(['cell id0 = ' num2str(cell_id0_ont(pair_n)) '   cell id2 = ' ...
    num2str(cell_id2_ont(pair_n)) '   cell n0 = ' num2str(cell_n0_ont(pair_n)) ...
    '   cell n2 = ' num2str(cell_n2_ont(pair_n))]);



opt = struct('fit_surround', true, 'fit_surround_sd_scale', true, 'fit_surround_amp_scale', true, 'fit_scale_one', false, 'fit_scale_two', false, 'fit_tau_one', false, 'fit_tau_two', false, 'verbose', true);

