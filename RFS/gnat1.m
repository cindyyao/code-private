
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun5 = load_data('/Analysis/xyao/2013-05-30-0/data005/data005', opt);
datarun8 = load_data('/Analysis/xyao/2013-05-30-0/data008/data008', opt);

cell_type = {'ON brisk transient', 'ON transient', 'OFF brisk sustained', 'OFF sustained', 'OFF slow', 'OFF transient'};

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun8, ...
    datarun5, cell_type);

[sta_f_weighted, sta_f_n, sta_c] = mean_rf2(datarun8, datarun5, cell_type, ...
    cell_id, cell_idx);


d0 = size(datarun8.stas.stas{1});
d1 = size(datarun5.stas.stas{1});
n = size(cell_type, 2);

fit_sta_c = cell(n, 2);
params_c = cell(n, 2);

for i = 4:n
    
    f = 40/d0(1);
%     temp_marks_sta = significant_stixels(sta_c{i, 1}, 'thresh', 3.5, 'time', 'max');
%     fit_ins = struct('sig_stixels', temp_marks_sta);
    fit_temp = fit_sta_sequence(sta_c{i, 1}); %, 'fit_instructions', fit_ins);
    fit_sta_c{i, 1} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x*f;
    params_temp(2) = fit_temp.center_point_y*f;
    params_temp(3) = fit_temp.center_sd_x*f;
    params_temp(4) = fit_temp.center_sd_y*f;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale*f;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale*f;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_c{i, 1} = params_temp;
    
    f = 40/d1(1);
%     temp_marks_sta = significant_stixels(sta_c{i, 2}, 'thresh', 5, 'time', 'max');
%     fit_ins = struct('sig_stixels', temp_marks_sta);
    fit_temp = fit_sta_sequence(sta_c{i, 2}); % 'fit_instructions', fit_ins);
    fit_sta_c{i, 2} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x*f;
    params_temp(2) = fit_temp.center_point_y*f;
    params_temp(3) = fit_temp.center_sd_x*f;
    params_temp(4) = fit_temp.center_sd_y*f;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale*f;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale*f;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_c{i, 2} = params_temp;


    i

end

save('gnat1.mat', 'sta_f_weighted', 'sta_f_n', 'sta_c', 'fit_sta_c', 'params_c')

ratio(2, :) = plot_mean_rf2(sta_c, cell_type, 2, params_c)

%%



figure
x = [1:5]; y = ratio(2:end, :); 
bar(x,y)
set(gca,'xticklabel',cell_type(2:end))
legend('NDF 0', 'NDF 1')