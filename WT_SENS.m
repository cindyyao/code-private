% load data

opt = struct('verbose',1,'load_params',1,'load_neurons',1, 'load_ei',1);
datarun = load_data('/Analysis/xyao/2013-04-09-0/data001-020/data001-020', opt);
slave_datarun = datarun;
master_datarun = load_data('/Analysis/xyao/2013-04-09-0/data021/data021', opt);
[cell_list_map, failed_cells] = map_ei(master_datarun,slave_datarun);

% get triggers of every trial

n = length(datarun.triggers)/4;
trigger = zeros(n, 1);
cell_numb = length(datarun.cell_ids);
for i = 1:n
    trigger(i) = datarun.triggers(4*(i-1)+1);
end

% get psth and raster
n_level = 12;
spike_added = cell(cell_numb, n_level);
raster_all_cells = cell(cell_numb, 1);
hist_all_cells = cell(cell_numb, n_level);
for m = 1:cell_numb
    raster_temp = get_raster(datarun.spikes{m}, trigger, 'plot', false);
    raster_all_cells{m, 1} = raster_temp; 

    % reshape the raster cell
    raster = cell(49, n_level);
    % every row in cell raster is data from one light intensity of stimulus

    for i = 1:n_level
        if i <= 2
            for j = 1:49
               raster{j, i} = raster_temp{50*(i-1)+j};
            end
        else
            for j = 1:24
                raster{j, i} = raster_temp{25*(i-1)+j};
            end
        end
    end
    
  
    % get psth and sensitivity curve
    for i = 1:n_level
        t = [];
        if i <= 2
            for j = 1:49
                t = [t; raster{j, i}];
            end
            spike_added{m, i} = t;
            X = 0.05:0.1:11.95; % 120 bins, 100ms/bin
            N = hist(t, X);
            N = N/49;
            hist_all_cells{m, i} = N;
        else
            for j = 1:24
                t = [t; raster{j, i}];
            end
            spike_added{m, i} = t;
            X = 0.05:0.1:11.95; % 120 bins, 100ms/bin
            N = hist(t, X);
            N = N/24;
            hist_all_cells{m, i} = N;
        end
    end
    fprintf(num2str(m))
end

% % sensitivity curve
NDF = [5.6 5 4 3 2.6 2.3 2 1.6 1.3 1 0.3 0];

% off GC
cell_id = 316;
cell_n = find(datarun.cell_ids == cell_id);
sens_s_mean = zeros(n_level, 4);
for j = 1:n_level
    sens_s_mean(j, 2) = max(hist_all_cells{cell_n, j}(31:40))-mean(hist_all_cells{cell_n, j}(21:30));
    sens_s_mean(j, 3) = max(hist_all_cells{cell_n, j}(61:70))-mean(hist_all_cells{cell_n, j}(51:60));
    sens_s_mean(j, 1) = mean(hist_all_cells{cell_n, j}(111:120))-min(hist_all_cells{cell_n, j}(1:10));
    sens_s_mean(j, 4) = mean(hist_all_cells{cell_n, j}(81:90))-min(hist_all_cells{cell_n, j}(91:100));
end
m_sens = max(sens_s_mean(:));
sens_s_mean = sens_s_mean/m_sens;

figure;
plot(-NDF, sens_s_mean(:, 1), 'r');
hold on
plot(-NDF, sens_s_mean(:, 2), 'g');
plot(-NDF, sens_s_mean(:, 3), 'b');
plot(-NDF, sens_s_mean(:, 4), 'k');
title(['cell id = ' num2str(cell_id)]);
xlabel('light level')
legend('0s', '3s', '6s', '9s', 'location', 'NorthWest')



% on GC
cell_id = 3422;
cell_n = find(datarun.cell_ids == cell_id);
sens_s_mean = zeros(n_level, 4);
for j = 1:n_level
    sens_s_mean(j, 1) = max(hist_all_cells{cell_n, j}(1:10))-mean(hist_all_cells{cell_n, j}(111:120));
    sens_s_mean(j, 4) = max(hist_all_cells{cell_n, j}(91:100))-mean(hist_all_cells{cell_n, j}(81:90));
    sens_s_mean(j, 2) = mean(hist_all_cells{cell_n, j}(21:30))-min(hist_all_cells{cell_n, j}(31:40));
    sens_s_mean(j, 3) = mean(hist_all_cells{cell_n, j}(51:60))-min(hist_all_cells{cell_n, j}(61:70));
end
m_sens = max(sens_s_mean(:));
sens_s_mean = sens_s_mean/m_sens;

figure;
plot(-NDF, sens_s_mean(:, 1), 'r');
hold on
plot(-NDF, sens_s_mean(:, 2), 'g');
plot(-NDF, sens_s_mean(:, 3), 'b');
plot(-NDF, sens_s_mean(:, 4), 'k');
title(['cell id = ' num2str(cell_id)]);
xlabel('light level')
legend('0s', '3s', '6s', '9s', 'location', 'NorthWest')


% on off GC
cell_id = 4982;
cell_n = find(datarun.cell_ids == cell_id);
sens_s_mean = zeros(n_level, 4);
for j = 1:n_level
    sens_s_mean(j, 2) = max(hist_all_cells{cell_n, j}(31:40))-mean(hist_all_cells{cell_n, j}(21:30));
    sens_s_mean(j, 3) = max(hist_all_cells{cell_n, j}(61:70))-mean(hist_all_cells{cell_n, j}(51:60));
    sens_s_mean(j, 1) = max(hist_all_cells{cell_n, j}(1:10))-mean(hist_all_cells{cell_n, j}(111:120));
    sens_s_mean(j, 4) = max(hist_all_cells{cell_n, j}(91:100))-mean(hist_all_cells{cell_n, j}(81:90));
end
m_sens = max(sens_s_mean(:));
sens_s_mean = sens_s_mean/m_sens;

figure;
plot(-NDF, sens_s_mean(:, 1), 'r');
hold on
plot(-NDF, sens_s_mean(:, 2), 'g');
plot(-NDF, sens_s_mean(:, 3), 'b');
plot(-NDF, sens_s_mean(:, 4), 'k');
title(['cell id = ' num2str(cell_id)]);
xlabel('light level')
legend('0s', '3s', '6s', '9s', 'location', 'NorthWest')




% draw raster and psth
cell_id = 76;
cell_n = find(datarun.cell_ids == cell_id);
raster = raster_all_cells{cell_n, 1};
figure;
for i = 1:n_level
    subplot(n_level, 2, 2*(i-1)+1);
    if i <= 2
        for j = 1:49
        SpikeTime = raster{50*(i-1)+j, 1};
        SpikeTime = SpikeTime';
        X = [SpikeTime; SpikeTime];
        Y = [ones(1, length(SpikeTime))*(j-0.9); ...
            ones(1, length(SpikeTime))*j];
        line(X, Y, 'color', 'b');
        axis([0,12,0,50])
        hold on
        end
    else
        for j = 1:24
        SpikeTime = raster{25*(i-1)+j, 1};
        SpikeTime = SpikeTime';
        X = [SpikeTime; SpikeTime];
        Y = [ones(1, length(SpikeTime))*(j-0.9); ...
            ones(1, length(SpikeTime))*j];
        line(X, Y, 'color', 'b');
        axis([0,12,0,25])
        hold on
        end
    end
    ylabel(['NDF ' num2str(NDF(i))])
    if i == 1
        title(['cell id = ' num2str(cell_id)]);
    end
    
    subplot(n_level, 2, 2*i);        
    h = hist_all_cells{cell_n, i};
    bar(h)
    xlim([0 120])  
    
    M = max(h);
    if isempty(M) == 0;
        ylim([0,M+0.1]);
    end

end

