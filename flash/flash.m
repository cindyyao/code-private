% load white noise data

opt = struct('verbose',1,'load_neurons',1, 'load_ei',1, 'load_sta', 1, 'load_params', 1);
datarun{1} = load_data('/Analysis/xyao/2013-05-30-0/data005/data005', opt);

% load flash data

opt = struct('verbose',1,'load_params',1,'load_neurons',1, 'load_ei',1);
datarun{2} = load_data('/Analysis/xyao/2013-05-30-0/data004/data004', opt);


% stimulus params

repeat = 30;
light_level = 10;
trial_duration = 2;

% load stimulus

datarun{2}.names.stimulus_path = '/Volumes/lab/Analysis/2013-05-30-0/stimulus/s04';
datarun{2} = load_stim(datarun{2}, 'user_defined_trigger_set', [1:2:repeat*light_level*2], 'find_trigger', 0);

% map ei

cell_type = {'on brisk transient', 'on transient', 'off sustained2', 'off sustained', 'off slow', 'off transient'};
% cell_type = {'on', 'off'};
n = length(cell_type);
cell_id = cell(n, 1);
cell_idx = cell(n, 1);
for i = 1:n
    [cell_list_map, ~] = map_ei(datarun{1}, datarun{2}, 'master_cell_type', cell_type{i});
    ept = 1-cellfun(@isempty, cell_list_map);
    master_id = datarun{1}.cell_ids.*ept;
    master_id(master_id == 0) = [];
    slave_id = cell2mat(cell_list_map);
    if isempty(slave_id) == 0
    cell_id{i}(1, :) = slave_id;
    cell_id{i}(2, :) = master_id;
    cell_idx{i}(1, :) = get_cell_indices(datarun{2}, slave_id);
    cell_idx{i}(2, :) = get_cell_indices(datarun{1}, master_id);
    end
end

%% classify trial index according to light level

triggers = datarun{2}.stimulus.triggers-trial_duration/2;
list = datarun{2}.stimulus.trial_list;

% if 1st trigger < 0, exclude the 1st repeat
if triggers(1)<0
    triggers = triggers(light_level+1:end);
    repeat = repeat - 1;
    list = list(light_level+1:end);
end


index = zeros(repeat, light_level);
params = zeros(1, light_level);
for i = 1:light_level
    index(:, i) = find(list == i);
    params(i) = datarun{2}.stimulus.trials(i).RGB(1);
end

[params, sequence] = sort(params);  
index_temp = index;
index = zeros(repeat, light_level);
for i = 1:light_level
    index(:, i) = index_temp(:, sequence(i));
end

spike_added = cell(n, 1);
raster_all_cells = cell(n, 1);

for i = 1:n
    cell_numb = size(cell_id{i}, 2);
    spike_added_ct = cell(cell_numb, 1);
    raster_ct = cell(cell_numb, 1);
    for rgc = 1:cell_numb
        spike = datarun{2}.spikes{cell_idx{i}(1, rgc)};
        raster = get_raster(spike, triggers, 'plot', false);
        spike_added_temp = cell(light_level, 1);
        raster_temp = cell(light_level, 1);
        for j = 1:light_level
            index_temp = index(:, j);
            raster_temp{j} = raster(index_temp);
            spike_added_temp{j} = cell2mat(raster(index_temp));
        end
    spike_added_ct{rgc} = spike_added_temp;
    raster_ct{rgc} = raster_temp;
    end
    spike_added{i} = spike_added_ct;
    raster_all_cells{i} = raster_ct;
end

%% plot raster

cell_type_n = 6;
for cell_n = 1:1;
raster = raster_all_cells{cell_type_n}{cell_n};
bin_size = 0.01;

N_total = [];

for i = 1:light_level
    t = spike_added{cell_type_n}{cell_n}{i};
    X = bin_size/2:bin_size:trial_duration-(bin_size/2); 
    N = hist(t, X);
    N_total = [N_total N];
end

figure;
for i = 1:light_level
    subplot(light_level, 2, 2*(i-1)+1);
    for j = 1:repeat
        SpikeTime = raster{i}{j};
        SpikeTime = SpikeTime';
        X = [SpikeTime; SpikeTime];
        Y = [ones(1, length(SpikeTime))*(j-0.9); ...
            ones(1, length(SpikeTime))*j];
        line(X, Y, 'color', 'b');
        xlim([0 trial_duration])
        ylim([0 repeat])
        hold on
    end
    ylabel([num2str(params(i))])
    if i == 1
        title(['cell id   ' num2str(cell_id{cell_type_n}(2, cell_n))]);
    end
    
    subplot(light_level, 2, 2*i);        
    t = spike_added{cell_type_n}{cell_n}{i};
    X = bin_size/2:bin_size:trial_duration-(bin_size/2); 
    hist(t, X);     
    xlim([0 trial_duration])  
    ylim([0,max(N_total)]);
   

end
end
%% plot sensitivity curve


sens_curve = cell(n, 1);
figure
for cell_type_n = 1:n;
bin_size = 0.01;
nn = size(cell_id{cell_type_n}, 2);
sens_curve_temp = zeros(nn, length(params));
subplot(2, 3, cell_type_n);
for cell_n = 1:nn;
sens= zeros(light_level, 1);

for j = 1:light_level
    t = spike_added{cell_type_n}{cell_n}{j};
    X = bin_size/2:bin_size:trial_duration-(bin_size/2); 
    N = hist(t, X);
    
    sens(j) = sqrt(mean((N(length(X)/2+1:length(X)/2+1/bin_size) - ...
        mean(N(length(X)/2-1/bin_size+1:length(X)/2))).^2)) - std(N(1:end/2));    
end


m_sens = max(sens(:));
sens = sens/m_sens;
sens_curve_temp(cell_n, :) = sens;
semilogx(params, sens);
xlim([params(1) params(end)])
hold on
end
% ylim([-0.2 1])
title(cell_type{cell_type_n});
xlabel('light level')
sens_curve{cell_type_n} = sens_curve_temp;
end
