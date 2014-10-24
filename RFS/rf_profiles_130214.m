datarun{1} = load_data('/Analysis/xyao/2013-02-21-0/data000/data000');
datarun{2} = load_data('/Analysis/xyao/2013-02-21-0/data004/data004');
 
for cn = 1:2
    datarun{cn} = load_neurons(datarun{cn});
    datarun{cn} = load_sta(datarun{cn}, 'load_sta', 'all');
    datarun{cn} = load_params(datarun{cn});
    datarun{cn} = load_ei(datarun{cn}, 'all');
    marks_params.thresh = 4.0;
    datarun{cn} = get_sta_summaries(datarun{cn}, 'all', 'marks_params', marks_params);
    datarun{cn} = get_sta_fits_from_vision(datarun{cn});
end
    
%% cell types comparison

cell_types = {'on brisk transient', 'on transient', 'off brisk transient', ...
     'off sustained'};
 
for i = 1:length(cell_types)
    temp_cell_type_nums{1}{i} = get_cell_type_nums(datarun{1}, cell_types{i});
    temp_cell_type_nums{2}{i} = get_cell_type_nums(datarun{2}, cell_types{i});
end


for i = 1:length(cell_types)
    [cell_id, cell_idx, ~, ~] = cell_map2(datarun{1}, datarun{2}, cell_types(i));
    datarun{1}.cell_types{temp_cell_type_nums{1}{i}}.cell_ids = cell_id{1}(:, 1);
    datarun{2}.cell_types{temp_cell_type_nums{2}{i}}.cell_ids = cell_id{1}(:, 2);
end


figure;
subplot(1, 2, 1)

% generate profile plots for a chosen cell type across light levels
profiles = interdig_compute_and_plot_neighbor_rf(datarun{1}, temp_cell_type_nums{1}, 'foa', -1,...
                        'extension_factor', 3.0, 'profile_points', 125, 'distance_cutoff_factor', 1.5, 'colors', 'krgb');
title('NDF 4')
 
subplot(1, 2, 2)
 
% generate profile plots for a chosen cell type across light levels
profiles = interdig_compute_and_plot_neighbor_rf(datarun{2},temp_cell_type_nums{2}, 'foa', -1,...
                        'extension_factor', 3.0, 'profile_points', 125, 'distance_cutoff_factor', 1.5, 'colors', 'krgb');
title('NDF 0')
                    
                    
%% comparison between light level

figure;

for tt = 1:4
    subplot(3, 4, tt)
    plot_rf_summaries(datarun{1}, cell_types{tt}, 'plot_fits', true)
    title('NDF 4')
    subplot(3, 4, tt+4)
    plot_rf_summaries(datarun{2}, cell_types{tt}, 'plot_fits', true)
    title('NDF 0')
    subplot(3, 4, tt+8)
    [~, h(1)] = interdig_compute_and_plot_neighbor_rf(datarun{1}, temp_cell_type_nums{1}(tt), 'foa', -1,...
                        'extension_factor', 3.0, 'profile_points', 125, 'distance_cutoff_factor', 1.5, 'colors', 'b');
    [~, h(2)] = interdig_compute_and_plot_neighbor_rf(datarun{2}, temp_cell_type_nums{2}(tt), 'foa', -1,...
                        'extension_factor', 3.0, 'profile_points', 125, 'distance_cutoff_factor', 1.5, 'colors', 'r');
                    legend(h, 'NDF 4', 'NDF 0')
                    ylim([-1 4])
                    title(cell_types{tt})
end

    

 
