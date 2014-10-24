opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
datarun{1} = load_data('/Analysis/xyao/2014-01-14-0/data000-001map/data000-001map', opt);
datarun{2} = load_data('/Analysis/xyao/2014-01-16-0/data000-001map/data000-001map', opt);

opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun_wn{1} = load_data('/Analysis/xyao/2014-01-14-0/data008-map/data008-map', opt);
datarun_wn{2} = load_data('/Analysis/xyao/2014-01-16-0/data008-map/data008-map', opt);

% find the start time of stimulus 
for i = 1:2
    trigger_diff{i} = diff(datarun{i}.triggers); 
    start_I{i} = [1; find(trigger_diff{i} >3.1)+1]; % after turning on the output, the duration of each trial is 3sec.
    start_time{i} = datarun{i}.triggers(start_I{i}); % convert the index of triggers to time
end
start_time{1}(12) = [];
start_I{1}(12) = [];

% calculate stimuli strengths
load('LED_calibration.mat')
load('CNG_1401.mat')

% stimuli{2}(1:2, :) = [];
% start_time{2}(1:2) = [];
% start_I{2}(1:2) = [];

for i = 1:2
    for j = 1:size(stimuli{i}, 1)
        I = find(LED_calibration(:, 1) == stimuli{i}(j, 1));
        stimuli{i}(j, 1) = LED_calibration(I, 2);
    end
    stimuli{i} = stimuli{i}(:, 1).*stimuli{i}(:, 2);
    [stimuli{i}, seq{i}] = sort(stimuli{i});
    delta = diff(stimuli{i});
    delta(delta ~= 0) = 1;
    delta = [1; delta];
    k = 0;
    for j = 1:length(stimuli{i})
        if delta(j) == 1
            k = k+1;
            stimuli_nrpt{i}{k} = [];
            stimuli_nrpt{i}{k} = [stimuli_nrpt{i}{k} j];
        else
            stimuli_nrpt{i}{k} = [stimuli_nrpt{i}{k} j];
        end
    end
    start_time{i} = start_time{i}(seq{i});
    start_I{i} = start_I{i}(seq{i});
        
end


% get cell ids
flash_id{1} = datarun{1}.cell_ids;
flash_id{2} = datarun{2}.cell_ids;

% cell_type = {'ON 1', 'ON 2', 'OFF 1', 'OFF 2', 'OFF 3'};

% for i = 1:2
%     for j = 1:length(cell_type)
%         clear id_temp1 id_temp2 id_temp3 id_temp_all
%         temp = map_ei(datarun_wn{i}, datarun{i}, 'master_cell_type', cell_type{j}, 'slave_cell_type', flash_id{i});
%         if sum(cellfun(@isempty, temp)) ~= length(datarun_wn{i}.cell_ids)
%         id_temp1(:, 1) = datarun_wn{i}.cell_ids(cellfun(@isempty, temp) == 0);
%         id_temp1(:, 2) = cell2mat(temp);
%         else
%             id_temp1 = [];
%         end
%         id_temp2 = intersect(flash_id{i}, get_cell_ids(datarun_wn{i}, cell_type{j}));
%         id_temp2 = repmat(id_temp2', 1, 2);
%         id_temp_all = [id_temp1; id_temp2];
%         id_temp_all = unique(id_temp_all, 'rows');
%         cell_id{i}{j} = id_temp_all;
%     end
% end

% get raster

% by cell type

% for i = 1:2
%     for ct = 1:length(cell_type)
%         for cn = 1:size(cell_id{i}{ct}, 1)
%             idx = get_cell_indices(datarun{i}, cell_id{i}{ct}(cn, 2));
%             spike = datarun{i}.spikes{idx};
%             for st = 1:length(stimuli{i})
%                 trigger = datarun{i}.triggers(start_I{i}(st):start_I{i}(st)+39);
%                 raster{i}{ct}{cn}(:, st) = get_raster(spike, trigger, 'plot', false, 'stop', 3);
%                 raster_all{i}{ct}{cn}{st} = sort(cell2mat(raster{i}{ct}{cn}(:, st)));
%             end
%             
%             spike_number_mean_temp = mean(cellfun(@length, raster{i}{ct}{cn}));
%             for j = 1:length(stimuli_nrpt{i})
%                 spikes_number_mean{i}{ct}{cn}(j) = mean(spike_number_mean_temp(stimuli_nrpt{i}{j}));
%                 raster_all_nrpt{i}{ct}{cn}{j} = cell2mat(raster_all{i}{ct}{cn}(stimuli_nrpt{i}{j})');
%             end
%         end
%     end
% end

% all cells
for i = 1:2
    for cn = 1:length(flash_id{i})
        idx = get_cell_indices(datarun{i}, flash_id{i}(cn));
        spike = datarun{i}.spikes{idx};
        for st = 1:length(stimuli{i})
            trigger = datarun{i}.triggers(start_I{i}(st):start_I{i}(st)+39);
            raster{i}{cn}(:, st) = get_raster(spike, trigger, 'plot', false, 'stop', 3);
            raster_all{i}{cn}{st} = sort(cell2mat(raster{i}{cn}(:, st)));
        end
            
        spike_number_mean_temp = mean(cellfun(@length, raster{i}{cn}));
        for j = 1:length(stimuli_nrpt{i})
            spikes_number_mean{i}{cn}(j) = mean(spike_number_mean_temp(stimuli_nrpt{i}{j}));
            raster_all_nrpt{i}{cn}{j} = cell2mat(raster_all{i}{cn}(stimuli_nrpt{i}{j})');
        end
    end
end

%% plot single sensitivity curves

% n = 120;
% XX = 3/n/2:3/n:3-3/n/2;
% 
% condition = {'control', 'rescue'};
% for i = 1:5;
%     for j = 1:2;
%     FigHandle = figure;
%     set(FigHandle, 'Position', [0, 0, 1080, 1080]);
% 
%     cn = size(cell_id{j}{i}, 1);
%     for k = 1:cn
%         subplot(ceil(sqrt(cn)), ceil(sqrt(cn)), k)
% %         semilogx(unique(stimuli{j}), spikes_number_mean{j}{i}{k})
%         clear sens
%         for st = 1:length(stimuli_nrpt{j})
%             N = hist(raster_all_nrpt{j}{i}{k}{st}, XX)/length(stimuli_nrpt{j}{st});
%             sens(st) = sqrt(mean((N(1:2/3*length(XX)) - mean(N(2/3*length(XX)+1:end))).^2)) - std(N(2/3*length(XX)+1:end)); 
%         end
%         sensitivity{j}{i}{k} = sens;
%         semilogx(unique(stimuli{j}), sens)
%         ylim([0 max(sens)])
%         if k == 1
%             title([cell_type{i} ' ' condition{j}])
%         end
%     end
%     end
% end

n = 120;
XX = 3/n/2:3/n:3-3/n/2;
%     FigHandle = figure;
%     set(FigHandle, 'Position', [0, 0, 1080, 1080]);
% 
condition = {'control', 'rescue'};
    for j = 1:2;

    cn = size(flash_id{j}, 1);
    for k = 1:cn
%         subplot(ceil(sqrt(cn)), ceil(sqrt(cn)), k)
%         semilogx(unique(stimuli{j}), spikes_number_mean{j}{i}{k})
        clear sens
        for st = 1:length(stimuli_nrpt{j})
            N = hist(raster_all_nrpt{j}{k}{st}, XX)/length(stimuli_nrpt{j}{st});
            sens(st) = sqrt(mean((N(1:2/3*length(XX)) - mean(N(2/3*length(XX)+1:end))).^2)) - std(N(2/3*length(XX)+1:end)); 
        end
        sensitivity_all{j}(k, :) = sens;
%         semilogx(unique(stimuli{j}), sens)
%         ylim([0 max(sens)])
%         if k == 1
%             title([cell_type{i} ' ' condition{j}])
%         end
    end
    end



%% plot raster and histgram


% for cd = 1:2;
% a = ceil(length(stimuli{cd})/2);
% 
% for ct = 1:5;
% 
% for cn =1:size(cell_id{cd}{ct}, 1);
% st = 1;
% 
% b = 1;
% 
%     FigHandle = figure;
%     set(FigHandle, 'Position', [0, 0, 1920, 1080]);
% 
% 
% while st <= length(stimuli{cd})
%     if  st == a+1
%         b = 3;
%     end
%     subplot(a, 4, b)
%     for j = 1:40
%     SpikeTime = raster{cd}{ct}{cn}{j, st};
%     SpikeTime = SpikeTime';
%     X = [SpikeTime; SpikeTime];
%     Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
%     line(X, Y, 'color', 'b');
%     axis([0, 3, 0, 40]);
%     hold on
%     end
%     if b == 1
%         title(num2str(cell_id{cd}{ct}(cn, 2)))
%     end
%     b = b+1;
%     subplot(a, 4, b)
%     hist(raster_all{cd}{ct}{cn}{st}, XX)
%     xlim([0 3])
%        
%     b = b+3;
%     st = st+1;
% end
% 
% screen_size = [24 12];
% set(figure(1), 'paperpositionmode', 'auto');
% % set(figure(1), 'PaperPosition', [-0.5 -0.25 22 10]);
% set(gcf, 'PaperUnits', 'inch');
% set(figure(1), 'PaperSize', screen_size);
% print(figure(1), '-dpdf', [cell_type{ct} '_' condition{cd} '_' num2str(cn)])
% close
% 
% 
% 
% end
% end
% end

%% plot sensitivity curves of cell type

% figure
% for i = 1:5
%     color = 'br';
%     subplot(2, 3, i)
%     for cd = 1:2
%     for j = 1:length(avg_i{cd}{i})
%     h(cd) = semilogx(unique(stimuli{cd}), sensitivity{cd}{i}{avg_i{cd}{i}(j)}/max(sensitivity{cd}{i}{avg_i{cd}{i}(j)}), color(cd));
%     %
%     hold on
%     end
%     
%     end
%     legend(h, condition{1}, condition{2}, 'location', 'northwest')
%     title(cell_type{i})
% end
    

figure
    color = 'br';
    for cd = 1:2
    for j = 1:length(flash_id{cd})
    h(cd) = semilogx(unique(stimuli{cd}), sensitivity_all{cd}(j, :), color(cd));
    %
    hold on
    end
    
    end
    legend(h, condition{1}, condition{2}, 'location', 'northwest')
    
figure
semilogx(unique(stimuli{1}), mean(sensitivity_all{1}), 'b')
hold on 
semilogx(unique(stimuli{2}), mean(sensitivity_all{2}), 'r')
legend('control', 'rescue', 'location', 'northwest')
xlabel('flash strength')
ylabel('response')
