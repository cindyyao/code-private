%% 2013-05-28-0

opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun3 = load_data('/Analysis/xyao/2013-05-28-0/data003/data003', opt);
datarun11 = load_data('/Analysis/xyao/2013-05-28-0/data011/data011', opt);
datarun16 = load_data('/Analysis/xyao/2013-05-28-0/data016/data016', opt);

cell_type = {'ON brisk transient', 'ON transient', 'OFF brisk transient', 'OFF transient'};

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map3(datarun16, ...
    datarun11, datarun3, cell_type);

[sta_f_weighted, sta_f_n, sta_c] = mean_rf3(datarun16, datarun11, datarun3, ...
    cell_type, cell_id, cell_idx);






d0 = size(datarun16.stas.stas{1});
d1 = size(datarun11.stas.stas{1});
d2 = size(datarun3.stas.stas{1});
f0 = 1200/d0(1);
f1 = 1200/d1(1);
f2 = 1200/d2(1);
n = size(cell_type, 2);

fit_sta_c = cell(n, 3);
params_c = cell(n, 3);


for i = 1:n
   
        
    f = 40/d0(1);
    temp_marks_sta = significant_stixels(sta_c{i, 1}, 'thresh', 5, 'time', 'max');
    fit_ins = struct('sig_stixels', temp_marks_sta);
    fit_temp = fit_sta_sequence(sta_c{i, 1}, 'fit_instructions', fit_ins);
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
    temp_marks_sta = significant_stixels(sta_c{i, 2}, 'thresh', 5, 'time', 'max');
    fit_ins = struct('sig_stixels', temp_marks_sta);
    fit_temp = fit_sta_sequence(sta_c{i, 2}, 'fit_instructions', fit_ins);
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

    f = 40/d2(1);
    temp_marks_sta = significant_stixels(sta_c{i, 3}, 'thresh', 5, 'time', 'max');
    fit_ins = struct('sig_stixels', temp_marks_sta);
    fit_temp = fit_sta_sequence(sta_c{i, 3}, 'fit_instructions', fit_ins);
    fit_sta_c{i, 3} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x*f;
    params_temp(2) = fit_temp.center_point_y*f;
    params_temp(3) = fit_temp.center_sd_x*f;
    params_temp(4) = fit_temp.center_sd_y*f;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale*f;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale*f;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_c{i, 3} = params_temp;
    
    i

end


save('gnat2.mat', 'sta_f_weighted', 'sta_f_n', 'sta_c', 'fit_sta_c', 'params_c')

[ratio] = plot_mean_rf3(sta_c, cell_type, 3, params_c)

%% 2013-02-21-0

opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun4 = load_data('/lab/Analysis/2013-02-21-0/data004/data004', opt);
datarun2 = load_data('/lab/Analysis/2013-02-21-0/data002/data002', opt);
datarun0 = load_data('/Analysis/xyao/2013-02-21-0/data000/data000', opt);

cell_type = {'ON brisk transient', 'ON transient', 'OFF brisk transient', 'OFF transient'};

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map3(datarun4, ...
    datarun2, datarun0, cell_type);

[sta_f_weighted, sta_f_n, sta_c] = mean_rf3(datarun4, datarun2, datarun0, ...
    cell_type, cell_id, cell_idx);

d0 = size(datarun4.stas.stas{1});
d1 = size(datarun2.stas.stas{1});
d2 = size(datarun0.stas.stas{1});
f0 = 1200/d0(1);
f1 = 1200/d1(1);
f2 = 1200/d2(1);
n = size(cell_type, 2);

fit_sta_c = cell(n, 3);
params_c = cell(n, 3);


for i = 1:n
        
    f = 40/d0(1);
%     temp_marks_sta = significant_stixels(sta_c{i, 1}, 'thresh', 5, 'time', 'max');
%     fit_ins = struct('sig_stixels', temp_marks_sta);
    fit_temp = fit_sta_sequence(sta_c{i, 1}, 'verbose', true); %'fit_instructions', fit_ins);
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
    fit_temp = fit_sta_sequence(sta_c{i, 2},'verbose', true); % 'fit_instructions', fit_ins);
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

    f = 40/d2(1);
%     temp_marks_sta = significant_stixels(sta_c{i, 3}, 'thresh', 5, 'time', 'max');
%     fit_ins = struct('sig_stixels', temp_marks_sta);
    fit_temp = fit_sta_sequence(sta_c{i, 3}, 'verbose', true); %, 'fit_instructions', fit_ins);
    fit_sta_c{i, 3} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x*f;
    params_temp(2) = fit_temp.center_point_y*f;
    params_temp(3) = fit_temp.center_sd_x*f;
    params_temp(4) = fit_temp.center_sd_y*f;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale*f;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale*f;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_c{i, 3} = params_temp;
    
    i

end


save('map_ndf4_0221.mat', 'sta_f_weighted', 'sta_f_n', 'sta_c', 'fit_sta_c', 'params_c')

[ratio] = plot_mean_rf3(sta_c, cell_type, 1, params_c)

%% 2013-02-14-0

opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun6 = load_data('/Analysis/xyao/2013-02-14-0/data006/data006', opt);
datarun2 = load_data('/Analysis/xyao/2013-02-14-0/data002/data002', opt);
datarun0 = load_data('/Analysis/xyao/2013-02-14-0/data000/data000', opt);

cell_type = {'ON brisk transient', 'ON transient', 'OFF brisk transient', 'OFF transient'};

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map3(datarun6, ...
    datarun2, datarun0, cell_type);

[sta_f_weighted, sta_f_n, sta_c] = mean_rf3(datarun6, datarun2, datarun0, ...
    cell_type, cell_id, cell_idx);






d0 = size(datarun6.stas.stas{1});
d1 = size(datarun2.stas.stas{1});
d2 = size(datarun0.stas.stas{1});
f0 = 1200/d0(1);
f1 = 1200/d1(1);
f2 = 1200/d2(1);
n = size(cell_type, 2);

fit_sta_c = cell(n, 3);
params_c = cell(n, 3);


for i = 1:n
   
        
    f = 40/d0(1);
    temp_marks_sta = significant_stixels(sta_c{i, 1}, 'thresh', 5, 'time', 'max');
    fit_ins = struct('sig_stixels', temp_marks_sta);
    fit_temp = fit_sta_sequence(sta_c{i, 1}, 'fit_instructions', fit_ins);
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
    temp_marks_sta = significant_stixels(sta_c{i, 2}, 'thresh', 5, 'time', 'max');
    fit_ins = struct('sig_stixels', temp_marks_sta);
    fit_temp = fit_sta_sequence(sta_c{i, 2}, 'fit_instructions', fit_ins);
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

    f = 40/d2(1);
    temp_marks_sta = significant_stixels(sta_c{i, 3}, 'thresh', 5, 'time', 'max');
    fit_ins = struct('sig_stixels', temp_marks_sta);
    fit_temp = fit_sta_sequence(sta_c{i, 3}, 'fit_instructions', fit_ins);
    fit_sta_c{i, 3} = fit_temp;
    params_temp = zeros(1, 8);
    params_temp(1) = fit_temp.center_point_x*f;
    params_temp(2) = fit_temp.center_point_y*f;
    params_temp(3) = fit_temp.center_sd_x*f;
    params_temp(4) = fit_temp.center_sd_y*f;
    params_temp(5) = fit_temp.center_rotation_angle;
    params_temp(6) = fit_temp.center_sd_x*fit_temp.surround_sd_scale*f;
    params_temp(7) = fit_temp.center_sd_y*fit_temp.surround_sd_scale*f;
    params_temp(8) = fit_temp.surround_amp_scale;
    params_c{i, 3} = params_temp;
    
    i

end


save('map_ndf4_0214.mat', 'sta_f_weighted', 'sta_f_n', 'sta_c', 'fit_sta_c', 'params_c')

[ratio] = plot_mean_rf3(sta_c, cell_type, 4, params_c)

%%

cell_type = {'ON brisk transient', 'ON transient', 'OFF brisk transient', 'OFF sustained'};
NDF = [0 -2 -4];

figure
for i = 1:4
    subplot(2, 2, i);
    plot(NDF, ratio_gnat2(i, :))
    hold on
    plot(NDF, ratio_0221(i, :), 'color', 'r')
    plot(NDF, ratio_0214(i, :), 'color', 'r')
%     ylim([0 1])
    title(cell_type(i))
    xlabel('log(light intensity)')
    ylabel('surround/center volume ratio')
    if i == 1
        legend('Gnat2 -/-', 'WT')
    end
    

end

%%

opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun4 = load_data('/lab/Analysis/2013-02-21-0/data004/data004', opt);
datarun2 = load_data('/lab/Analysis/2013-02-21-0/data002/data002', opt);
datarun0 = load_data('/Analysis/xyao/2013-02-21-0/data000/data000', opt);

cell_type = {'ON brisk transient', 'ON transient', 'OFF brisk transient', 'OFF transient'};

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map3(datarun4, ...
    datarun2, datarun0, cell_type);

cell_spec = cell_id{4}(:, 1);
