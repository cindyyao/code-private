% load data

opt = struct('verbose',1,'load_params',1,'load_neurons',1, 'load_ei',1);
datarun = load_data('/Analysis/xyao/2013-04-09-0/data001-data020/data001-data020', opt);
slave_datarun = datarun;
master_datarun = load_data('/Analysis/xyao/2013-03-11-0/data000/data000', opt);
[cell_list_map, failed_cells] = map_ei(master_datarun,slave_datarun);

% get triggers of every trial

n = length(datarun.triggers)/4;
trigger = zeros(n, 1);
cell_numb = length(datarun.cell_ids);
for i = 1:n
    trigger(i) = datarun.triggers(4*(i-1)+1);
end


% get psth and raster
spike_added = cell(cell_numb, 9);
raster_all_cells = cell(cell_numb, 1);
hist_all_cells = cell(cell_numb, 9);
for m = 1:cell_numb
    raster_temp = get_raster(datarun.spikes{m}, trigger, 'plot', true);
    raster_temp1 = raster_temp; % change sequence
    raster_temp1(51:100) = raster_temp(76:125);
    raster_temp1(101:125) = raster_temp(51:75);
    raster_temp = raster_temp1;
    raster_all_cells{m, 1} = raster_temp; 

    % reshape the raster cell
    raster = cell(24, 9);
    % every row in cell raster is data from one light intensity of stimulus

    for i = 1:24
        for j = 1:9
            raster{i, j} = raster_temp{25*(j-1)+i};
        end
    end


    % get psth and sensitivity curve
    for i = 1:9
        t = [];
        for j = 1:24
            t = [t; raster{j, i}];
        end
        spike_added{m, i} = t;
        X = 0.05:0.1:11.95; % 120 bins, 100ms/bin
        N = hist(t, X);
        hist_all_cells{m, i} = N;
    end
    fprintf(num2str(m))
end

% % sensitivity curve

% off GC
NDF = [4 3 2.6 2.3 2 1.6 1.3 1 0];
cell_id = 811;
cell_n = find(datarun.cell_ids == cell_id);
sens_s_mean = zeros(9, 4);
for j = 1:9
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
NDF = [4 3 2.6 2.3 2 1.6 1.3 1 0];
cell_id = 708;
cell_n = find(datarun.cell_ids == cell_id);
sens_s_mean = zeros(9, 4);
for j = 1:9
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
NDF = [4 3 2.6 2.3 2 1.6 1.3 1 0];
cell_id = 4982;
cell_n = find(datarun.cell_ids == cell_id);
sens_s_mean = zeros(9, 4);
for j = 1:9
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
cell_id = 1503;
cell_n = find(datarun.cell_ids == cell_id);
raster = raster_all_cells{cell_n, 1};
figure;
for i = 1:9
    subplot(9, 2, 2*(i-1)+1);
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
    ylabel(['NDF ' num2str(NDF(i))])
    if i == 1
        title(['cell id = ' num2str(cell_id)]);
    end
    
    subplot(9, 2, 2*i);        
    t = spike_added{cell_n, i};
    X = 0.05:0.1:11.95; % 120 bins, 100ms/bin
    N = hist(t, X);
    hist(t, X);     
    xlim([0 12])  
    
    M = max(N);
    if isempty(M) == 0;
        ylim([0,M+10]);
    end

end

