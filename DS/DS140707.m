%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
datarun{1} = load_data('/Analysis/xyao/2014-07-07-0/data001-map/data001-map', opt);
datarun{1}.names.stimulus_path = '/Analysis/xyao/2014-07-07-0/stimuli/s01';
datarun{1} = load_stim(datarun{1}, 'user_defined_trigger_interval', 10);

datarun{2} = load_data('/Analysis/xyao/2014-07-07-0/data002-map/data002-map', opt);
datarun{2}.names.stimulus_path = '/Analysis/xyao/2014-07-07-0/stimuli/s02';
datarun{2} = load_stim(datarun{2}, 'user_defined_trigger_interval', 10);
% datarun{2} = get_autocorrelations(datarun{2}, id);

datarun{3} = load_data('/Analysis/xyao/2014-07-07-0/data004-map/data004-map', opt);
datarun{3}.names.stimulus_path = '/Analysis/xyao/2014-07-07-0/stimuli/s04';
datarun{3} = load_stim(datarun{3}, 'user_defined_trigger_interval', 10);

datarun{4} = load_data('/Analysis/xyao/2014-07-07-0/data005-map/data005-map', opt);
datarun{4}.names.stimulus_path = '/Analysis/xyao/2014-07-07-0/stimuli/s05';
datarun{4} = load_stim(datarun{4}, 'user_defined_trigger_interval', 10);

opt = struct('load_params', 1,'load_neurons', 1);

datarun{5} = load_data('/Analysis/xyao/2014-07-07-0/data003-map/data003-map', opt);
datarun{6} = load_data('/Analysis/xyao/2014-07-07-0/data006-map/data006-map', opt);
datarun{7} = load_data('/Analysis/xyao/2014-07-07-0/data000-map/data000-map', opt);
datarun{8} = load_data('/Analysis/xyao/2014-07-07-0/data000/data000', opt);

datarun{9} = load_data('/Analysis/xyao/2014-07-07-0/data008-map/data008-map', opt);
datarun{9}.names.stimulus_path = '/Analysis/xyao/2014-07-07-0/stimuli/s08';
datarun{9} = load_stim(datarun{9}, 'user_defined_trigger_interval', 10);

% datarun{10} = load_data('/Analysis/xyao/2014-07-07-0/data007-map/data007-map', opt);
% datarun{10} = load_sta(datarun{10});
% datarun{10} = get_autocorrelations(datarun{10}, intersect(id, datarun{10}.cell_ids));
%% 

[NumSpikesCell, StimComb] = get_spikescellstim(datarun{2},datarun{2}.cell_ids,0);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);

% pull out DS cells

figure
plot(ds_struct.mag{2, 1}, ds_struct.mag{3, 1}, 'o')
title('data002 after mapping')
xlabel('TP 24')
ylabel('TP 60')
hold on
[x, y] = ginput;
plot(x, y);

IN = inpolygon(ds_struct.mag{2, 1}, ds_struct.mag{3, 1}, x, y);
[~, I] = find(IN == 1);
id = datarun{2}.cell_ids(I);
idx = get_cell_indices(datarun{2}, id);
I = zeros(length(id), 4);
for i = 1:4
    I(:, i) = ismember(id, datarun{i}.cell_ids);
end
I = sum(I');
id4 = id(I == 4); % DS cells that can be found in all light levels
idx4 = arrayfun(@(x) find(id == x), id4);

[raster, raster_p_sum, p_idx] = deal(cell(4, 1));
for i = 1:4    
    [NumSpikesCell, StimComb] = get_spikescellstim(datarun{i},id,0);
    DS{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster{i} = get_ds_raster(datarun{i}, id, 'MG', 'stop', 8);
end

param_p = 5; % choose which params to use to calculate prefer direction indices 
MAG_all_norm = cell(4, 1);
max_r = cell(4, 1);
norm_max_r = cell(4, 1);

for i = 1:4
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster{i}, DS{i}.angle{param_p});
    MAG_all_norm{i} = normalize_MAG(DS{i});
    rep = datarun{i}.stimulus.repetitions;
    max_r{i} = max_firing_rate(raster_p_sum{i}, 0.2, 8)/rep;
    norm_max_r{i} = max_r{i}./repmat(max(max_r{i}, [], 2), 1, size(max_r{i}, 2));
end

ll = {'NDF3 SP60', 'NDF3 SP240', 'NDF0 SP60', 'NDF0 SP240'};


%% plot cell summary

for cc = 31:31 %length(id)
    plot_ds_raster(DS, raster, cc, id(cc), ll, 2, 2, 0)
end

%% plot single cell tuning curve

% use unnormalized vector sum as response
figure
for i = 1:4
    v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    subplot(2, 2, i)
    semilogx(v, exciseColumn(MAG_all_norm{i}), 'b')
    xlabel('speed')
    ylabel('Response')
    title(ll{i})
    xlim([v(end) v(1)])
end

% use maximum firing rate of preferred direction as response

figure
for i = 1:4
    v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    subplot(2, 2, i)
    semilogx(v, exciseColumn(norm_max_r{i}'), 'b')
    xlabel('speed')
    ylabel('Response')
    title(ll{i})
    xlim([v(end) v(1)])
end

%% classification based on speed tunning
%% pca
mag_pca = MAG_all_norm{2};
% mag_pca = MAG_all_norm{1}(:, idx4);
mag_pca = mag_pca';
[id_sub, idx_sub] = deal(cell(3, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(mag_pca);
pc1 = 1; pc2 = 2;
plot(scores(:, pc1), scores(:, pc2), 'o')
hold on
for i = 1:3
    [x, y] = ginput;
    plot(x, y)
    IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
    [~, idx_sub{i}] = find(IN' == 1);
    id_sub{i} = id(idx_sub{i});
end
xlabel('1st Principal Component')
ylabel('2nd Principal Component')

for j = 1:3
    figure
    for i = 1:4
        v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
        subplot(2, 2, i)
        semilogx(v, exciseColumn(MAG_all_norm{i}(:, idx_sub{j})), 'b')
        xlabel('speed')
        ylabel('Response')
        title(ll{i})
        xlim([v(end) v(1)])
    end
end

%% peak_area

% % peak location
% v = datarun{1}.stimulus.params.SPATIAL_PERIOD./datarun{1}.stimulus.params.TEMPORAL_PERIOD;
% mag = MAG_all_norm{1};
% mag(:, any(isnan(mag))) = [];
% tpn = size(mag, 1);
% peak = v(mod(find(mag == 1), tpn));
% 
% % area under curve
% area = sum(mag) - mag(1, :)/2 - mag(end, :)/2;
% 
% figure
% scatter(peak, area)
% set(gca, 'XScale', 'log')
% xlabel('peak location')
% ylabel('area under curve')
% title('NDF3 SP60')
% 
%% concatenated tuning curves
% mag_all = cell2mat(MAG_all_norm);
% mag_all = mag_all(:, idx4);
% mag_all = mag_all';
% 
% [id_sub, idx_sub] = deal(cell(2, 1));
% 
% FigHandle = figure;
% set(FigHandle, 'Position', [1 1 380 400])
% 
% [~,scores,~,~] = princomp(mag_all);
% pc1 = 1; pc2 = 2;
% plot(scores(:, pc1), scores(:, pc2), 'o')
% hold on
% for i = 1:2
%     [x, y] = ginput;
%     plot(x, y)
%     IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
%     [~, idx] = find(IN' == 1);
%     idx_sub{i} = idx4(idx);
%     id_sub{i} = id4(idx);
% end
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
% 
% for j = 1:2
%     figure
%     for i = 1:4
%         v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
%         subplot(2, 2, i)
%         semilogx(v, exciseColumn(MAG_all_norm{i}(:, idx_sub{j})), 'b')
%         xlabel('speed')
%         ylabel('Response')
%         title(ll{i})
%         xlim([v(end) v(1)])
%     end
% end
% 
% 

%% plot average tunning curve
color = 'brk';
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:4
    v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    subplot(2, 2, i)
    for ct = 1:2
        mag_temp = exciseColumn(MAG_all_norm{i}(:, idx_sub{ct}));
        tuning_avg{i}(:, ct) = mean(mag_temp, 2);
        tuning_ste{i}(:, ct) = std(mag_temp, [], 2)/sqrt(size(mag_temp, 2));
        errorbar(v, tuning_avg{i}(:, ct), tuning_ste{i}(:, ct), color(ct))
        hold on
    end
    set(gca, 'XScale', 'log')
    title(ll{i})
    xlim([min(v) max(v)])
    xlabel('speed')
    ylabel('response')
end
legend('on-off DSGC', 'on DSGC', 'location', 'southeast')

%% vector sum plot of individual type
for ct = 1:2
    [NumSpikesCell, StimComb] = get_spikescellstim(datarun{2},id_sub{ct},0);
    dscellanalysis(NumSpikesCell, StimComb);
end

%% compare across light level
SP = {'SP 60', 'SP 240'};
CT = {'on-off DSGC', 'on-DSGC'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:2
    v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    for j = 1:2
        subplot(2, 2, i+2*(j-1))
        errorbar(v, tuning_avg{i}(:, j), tuning_ste{i}(:, j), 'b')
        hold on
        errorbar(v, tuning_avg{i+2}(:, j), tuning_ste{i+2}(:, j), 'r')
        title([SP{i} '  ' CT{j}])
        set(gca, 'XScale', 'log')
        xlim([min(v) max(v)])
        xlabel('speed')
        ylabel('response')

    end
end
legend('NDF 3', 'NDF 0', 'location', 'northeast')


%% full field pulses

[raster_ff, raster_ff_all] = deal(cell(2, 1));
for d = 5:6
    [raster_ff{d-4}, raster_ff_all{d-4}] = get_ffp_raster(datarun{d}, id, 3);
end

for i = 7:7 %length(id) 
    if ~isempty(raster_ff{1}{i}) || ~isempty(raster_ff{2}{i})
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 800 800])
        for d = 1:2
        subplot(1, 2, d)
        plot_ffp(raster_ff{d}, raster_ff_all{d}, i, 3)
        title([num2str(id(i)) ' ' ll{2*d}(1:4)])
        end
        
%         print_close(1, [12, 12], num2str(id(i)))
    end
end

%% ffp psth
raster_ff_all{2}{33} = [];
raster_ff_all{2}{6} = [];
raster_ff_all{1}{6} = [];
raster_ff{2}{33} = [];
raster_ff{2}{6} = [];
raster_ff{1}{6} = [];

LL = {'NDF 3', 'NDF 0'};

XX = 0.05:0.05:11.95;
figure
for i = 1:4
    idx_temp = idx_dir_oo{i};
    for cd = 1:2
        psth{cd}{i} = [];
        subplot(4, 2, 2*(i-1)+cd)
        a = 0;
        for cc = 1:length(idx_temp)
            if ~isempty(raster_ff_all{cd}{idx_temp(cc)})
                h = hist(raster_ff_all{cd}{idx_temp(cc)}, XX);
                plot(XX, h/30, color(cd+1))
                xlabel('time/sec')
                ylabel(ctype{i})
                if i == 1
                    title(LL{cd})
                end
                hold on
                a = a+1;
                psth{cd}{i}(a, :) = h/30;
            end
        end
        xlim([0 12])    
    end
end

% for i = 1:4
%     for cd = 1:2
%         psth_temp = psth{cd}{i};
%         if ~isempty(psth_temp)
%             subplot(4, 2, 2*i)
%             errorbar(XX, mean(psth_temp, 1), std(psth_temp, [], 1)/sqrt(size(psth_temp, 1)), color(cd+1));
%             hold on
%         end
%     end
% end
%                 
%% 

trial_n = 25;
for cd = 1:2
[FFP.burst_spike_n{cd}, FFP.burst_spike_n_avg{cd}, FFP.burst_spike_n_ste{cd}, FFP.burst_start_time{cd}, FFP.burst_start_time_avg{cd}, FFP.burst_start_time_ste{cd}] = deal(cell(4, 1));
    for ct = 1:4
        cell_idx = idx_dir_oo{ct};
        for cc = 1:length(cell_idx)
            if ~isempty(raster_ff{cd}{cell_idx(cc)})
                for trial = 1:trial_n
                    raster_temp = raster_ff{cd}{cell_idx(cc)}{trial};
                    raster_temp = raster_temp(raster_temp > 9);
                    if ~isempty(raster_temp)
                        spike_interval = find(diff(raster_temp) > 0.2);
                        start_time_i = [1; spike_interval+1];
                        start_time = raster_temp(start_time_i);
                        spike_n = diff([0; spike_interval; length(raster_temp)]);
                        FFP.time{cd}{ct}{cc}{trial} = start_time;
                        FFP.spike_n{cd}{ct}{cc}{trial} = spike_n;
                    end
                end

                burst_n = cellfun(@(x)(length(x)), FFP.time{cd}{ct}{cc});
                max_burst_n = max(burst_n);
                for burst = 1:max_burst_n
                    FFP.burst_start_time{cd}{ct}{cc}{burst} = cellfun(@(x)(x(burst)), FFP.time{cd}{ct}{cc}(burst_n >= burst));
                    FFP.burst_spike_n{cd}{ct}{cc}{burst} = cellfun(@(x)(x(burst)), FFP.spike_n{cd}{ct}{cc}(burst_n >= burst));
                    FFP.burst_start_time_avg{cd}{ct}{cc}(burst) = mean(FFP.burst_start_time{cd}{ct}{cc}{burst});
                    FFP.burst_spike_n_avg{cd}{ct}{cc}(burst) = mean(FFP.burst_spike_n{cd}{ct}{cc}{burst});
                    FFP.burst_start_time_ste{cd}{ct}{cc}(burst) = std(FFP.burst_start_time{cd}{ct}{cc}{burst})/sqrt(length(FFP.burst_start_time{cd}{ct}{cc}{burst}));
                    FFP.burst_spike_n_ste{cd}{ct}{cc}(burst) = std(FFP.burst_spike_n{cd}{ct}{cc}{burst})/sqrt(length(FFP.burst_spike_n{cd}{ct}{cc}{burst}));
                end


            end
        end
    end
end


dir = {'Inferior', 'A(P)', 'Superior', 'P(A)'};
condition = {'NDF 3', 'NDF 0'};
color = 'brgmckybrgmcky';
figure
for cd = 1:2
    h = subplot(1, 2, cd);
    for ct = 1:4
        if ~isempty(FFP.burst_start_time_avg{cd}{ct})
            
            for cc = 1:length(FFP.burst_start_time_avg{cd}{ct});
            %     x = FFP.burst_start_time_avg{cd}{cc};
            %     x_bar = FFP.burst_start_time_ste{cd}{cc};
            %     y = FFP.burst_spike_n_avg{cd}{cc};
            %     y_bar = FFP.burst_spike_n_ste{cd}{cc};


            %     errorbar(x, y, y_bar, y_bar, 's-')
            %     hold on
                if ~isempty(FFP.burst_start_time_avg{cd}{ct}{cc})
                    x = ct;
                    y = FFP.burst_start_time_avg{cd}{ct}{cc};
                    y_bar = FFP.burst_start_time_ste{cd}{ct}{cc};

                   for i = 1:length(y)
                       errorbar(x, y(i), y_bar(i), y_bar(i), 'Marker', 's', 'color', color(i))
                       hold on
                   end
                end
            end
        end
    end
    xlim([0.5 4.5])
    set(h,'XTicklabel',dir)                
    ylabel('time (second)')  
    title(condition{cd})
end
%% autocorrelation

id_wn = intersect(id, datarun{10}.cell_ids);
figure
plot_autocorrelograms(datarun{10}, intersect(id, datarun{10}.cell_ids), 'foa', -1, 'normalize', true, 'line_color', 'b')
autoc = [];
for i = 1:length(datarun{10}.cell_ids)
    if ~isempty(datarun{10}.autocorrelation{i})
        autoc = [autoc; datarun{10}.autocorrelation{i}.probabilities];
    end
end

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(autoc);
pc1 = 1; pc2 = 2;
plot(scores(:, pc1), scores(:, pc2), 'o')
hold on
for i = 1:2
    [x, y] = ginput;
    plot(x, y)
    IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
    [~, idx_auto{i}] = find(IN' == 1);
    id_auto{i} = id_wn(idx_auto{i});
end

figure
for i = 1:3
    subplot(2, 2, i)
    plot_autocorrelograms(datarun{10}, id_auto{i}, 'foa', -1, 'normalize', true)
end
%% frequency analysis
    
duration = 8;
bin_rate = 10000;
hist_spikes = cell(4, 1);
signal_length = duration*bin_rate;                
NFFT = 2^nextpow2(signal_length);
f = bin_rate/2*linspace(0,1,NFFT/2+1);
fft_spikes = cell(4, 1);
[DC, F1, F2] = deal(cell(4, 1));

for i = 1:4
    tp = datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    [DC{i}, F1{i}, F2{i}] = deal(zeros(length(id), length(tp)));
    for rgc = 1:length(id)
        if ~isempty(raster{i}{rgc})
        for time = 1:length(tp)
            spikes = floor(raster_p_sum{i}{rgc}{time}*bin_rate);
            tmp_binned_spikes = zeros(1, signal_length);
            tmp_binned_spikes(spikes) = 1;
            hist_spikes{i}{rgc}(time, :) = tmp_binned_spikes;
            
            f1 = 60/tp(time); %Hz
            f2 = f1*2;
            f_diff1 = f - f1;
            f_diff2 = f - f2;
            [~,f1_index] = min(abs(f_diff1));
            [~,f2_index] = min(abs(f_diff2));
            tmp_fft = fft(tmp_binned_spikes, NFFT)./ signal_length;
            fft_spikes{i}{rgc}(time, :) = 2*abs(tmp_fft(1:NFFT/2+1));
            if f1_index > 1
                fund_power(time) = sum(fft_spikes{i}{rgc}(time, f1_index-1:f1_index+1)); % f1_index+2???
                sec_power(time) = sum(fft_spikes{i}{rgc}(time, f2_index-1:f2_index+1));
            else
                fund_power(time) = sum(fft_spikes{i}{rgc}(time, f1_index:f1_index+2)); % f1_index+2???
                sec_power(time) = sum(fft_spikes{i}{rgc}(time, f2_index:f2_index+2));
            end
            DC_power(time) = fft_spikes{i}{rgc}(time, 1);
        end
    % stores info for this cell into the matrix tuning curves
        F1{i}(rgc,:) = fund_power ./ max(DC_power);
        F2{i}(rgc,:) = sec_power ./ max(DC_power);
        DC{i}(rgc,:) = DC_power ./ max(DC_power);
        
        clear fund_power sec_power DC_power

        end
        
    end
    ratio{i} = F1{i}./F2{i};
end

% plot 

figure
for i = 1:4
    speed = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    subplot(2, 2, i) 
    for rgc = 1:length(id)
    if ~isempty(raster{i}{rgc})
    semilogx(speed, DC{i}(rgc,:), 'b')
    hold on
    end
    end
    xlabel('speed')
    xlim([speed(end) speed(1)])
    title(ll{i})
end

figure
for i = 1:4
    speed = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    subplot(2, 2, i) 
    for rgc = 1:length(id)
    if ~isempty(raster{i}{rgc})
    semilogx(speed, F1{i}(rgc,:), 'r')
    hold on
    end
    end
    xlabel('speed')
    xlim([speed(end) speed(1)])
    title(ll{i})

end

figure
for i = 1:4
    speed = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    subplot(2, 2, i) 
    for rgc = 1:length(id)
    if ~isempty(raster{i}{rgc})
    semilogx(speed, F2{i}(rgc,:), 'k')
    hold on
    end
    end
    xlabel('speed')
    xlim([speed(end) speed(1)])
    title(ll{i})

end

      
for i = 1:6
    figure
    for rgc = 1:length(id)
    if ~isempty(raster{i}{rgc})
        subplot(7, 8, rgc)
        semilogx(speed, DC{i}(rgc,:), 'b')
        hold on
        semilogx(speed, F1{i}(rgc,:), 'r')
        semilogx(speed, F2{i}(rgc,:), 'k')
        title(num2str(id(rgc)))
    end
    end
end
    
for i = 1:6
    figure
    for rgc = 1:length(id_sub)
    if ~isempty(raster{i}{idx_sub(rgc)})
        subplot(6, 8, rgc)
        semilogx(speed, DC{i}(idx_sub(rgc),:), 'b')
        hold on
        semilogx(speed, F1{i}(idx_sub(rgc),:), 'r')
        semilogx(speed, F2{i}(idx_sub(rgc),:), 'k')
        title(num2str(id_sub(rgc)))
    end
    end
end

for i = 1:6
    figure
    for rgc = 1:length(id2)
    if ~isempty(raster{i}{idx2(rgc)})
        subplot(3, 3, rgc)
        semilogx(speed, DC{i}(idx2(rgc),:), 'b')
        hold on
        semilogx(speed, F1{i}(idx2(rgc),:), 'r')
        semilogx(speed, F2{i}(idx2(rgc),:), 'k')
        title(num2str(id2(rgc)))
    end
    end
end

% type 1
[DC_avg_1, DC_ste_1, f1_avg_1, f1_ste_1,f2_avg_1, f2_ste_1, ...
    DC_avg_2, DC_ste_2, f1_avg_2, f1_ste_2,f2_avg_2, f2_ste_2] = deal(zeros(4, length(tp)));
for i = 1:4
    DC_tuning_1 = DC{i}(idx_sub, :);
    temp = sum(DC_tuning_1, 2);
    [~, idx] = find(temp' ~= 0);
    DC_tuning_1 = DC_tuning_1(idx, :);
    DC_avg_1(i, 1:length(v{i})) = mean(DC_tuning_1);
    DC_ste_1(i, 1:length(v{i})) = std(DC_tuning_1)/sqrt(size(DC_tuning_1, 1));
    
    f1_tuning_1 = F1{i}(idx_sub, :);
    temp = sum(f1_tuning_1, 2);
    [~, idx] = find(temp' ~= 0);
    f1_tuning_1 = f1_tuning_1(idx, :);
    f1_avg_1(i, 1:length(v{i})) = mean(f1_tuning_1);
    f1_ste_1(i, 1:length(v{i})) = std(f1_tuning_1)/sqrt(size(f1_tuning_1, 1));
    
    f2_tuning_1 = F2{i}(idx_sub, :);
    temp = sum(f2_tuning_1, 2);
    [~, idx] = find(temp' ~= 0);
    f2_tuning_1 = f2_tuning_1(idx, :);
    f2_avg_1(i, 1:length(v{i})) = mean(f2_tuning_1);
    f2_ste_1(i, 1:length(v{i})) = std(f2_tuning_1)/sqrt(size(f2_tuning_1, 1));
end
    

% type 2
for i = 1:4
    DC_tuning_2 = DC{i}(idx2, :);
    temp = sum(DC_tuning_2, 2);
    [~, idx] = find(temp' ~= 0);
    DC_tuning_2 = DC_tuning_2(idx, :);
    DC_avg_2(i, 1:length(v{i})) = mean(DC_tuning_2);
    DC_ste_2(i, 1:length(v{i})) = std(DC_tuning_2)/sqrt(size(DC_tuning_2, 1));
    
    f1_tuning_2 = F1{i}(idx2, :);
    temp = sum(f1_tuning_2, 2);
    [~, idx] = find(temp' ~= 0);
    f1_tuning_2 = f1_tuning_2(idx, :);
    f1_avg_2(i, 1:length(v{i})) = mean(f1_tuning_2);
    f1_ste_2(i, 1:length(v{i})) = std(f1_tuning_2)/sqrt(size(f1_tuning_2, 1));
    
    f2_tuning_2 = F2{i}(idx2, :);
    temp = sum(f2_tuning_2, 2);
    [~, idx] = find(temp' ~= 0);
    f2_tuning_2 = f2_tuning_2(idx, :);
    f2_avg_2(i, 1:length(v{i})) = mean(f2_tuning_2);
    f2_ste_2(i, 1:length(v{i})) = std(f2_tuning_2)/sqrt(size(f2_tuning_2, 1));
end

% plot
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:4
    subplot(2, 2, i)
    h1 = errorbar(speed{i}, DC_avg_1(i, 1:length(v{i})), DC_ste_1(i,1:length(v{i})), 'b');
    set(get(h1,'Parent'), 'XScale', 'log')
    hold on
    h2 = errorbar(speed{i}, f1_avg_1(i, 1:length(v{i})), f1_ste_1(i, 1:length(v{i})), 'r');
    set(get(h2,'Parent'), 'XScale', 'log')
    h3 = errorbar(speed{i}, f2_avg_1(i, 1:length(v{i})), f2_ste_1(i, 1:length(v{i})), 'k');
    set(get(h3,'Parent'), 'XScale', 'log')
    xlabel('speed')
    xlim([0.1 20])
    title(ll{i})
    legend('DC', 'F1', 'F2', 'location', 'northwest')
end

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:4
    subplot(2, 2, i)
    h1 = errorbar(speed{i}, DC_avg_2(i, 1:length(v{i})), DC_ste_2(i, 1:length(v{i})), 'b');
    set(get(h1,'Parent'), 'XScale', 'log')
    hold on
    h2 = errorbar(speed{i}, f1_avg_2(i, 1:length(v{i})), f1_ste_2(i, 1:length(v{i})), 'r');
    set(get(h2,'Parent'), 'XScale', 'log')
    h3 = errorbar(speed{i}, f2_avg_2(i, 1:length(v{i})), f2_ste_2(i, 1:length(v{i})), 'k');
    set(get(h3,'Parent'), 'XScale', 'log')
    xlabel('speed')
    xlim([0.1 20])
    title(ll{i})
    legend('DC', 'F1', 'F2')
end

%% DSI

for i = 1:4
    DS{i}.dsindex = cell2mat(DS{i}.dsindex');
end

% compare across cell types

FigHandle = figure;
set(FigHandle, 'Position', [1 1 900 900])
for i = 1:4
    subplot(2, 2, i)
    h1 = semilogx(v{i}, DS{i}.dsindex(idx4(idx_sub), :), 'b');
    h(1) = h1(1);
    hold on
    h2 = semilogx(v{i}, DS{i}.dsindex(idx4(idx2), :), 'r');
    h(2) = h2(1);
    xlim([v{2}(end) v{2}(1)])
    if mod(i, 2) == 0
        ylim([.5 1])
    end
    title(ll{i})
    xlabel('speed')
    ylabel('DSI')
    legend(h, 'on-off DSGC', 'on-DSGC', 'location', 'southeast')
end

FigHandle = figure;
set(FigHandle, 'Position', [1 1 900 900])
set(FigHandle, 'DefaultLineLineWidth', 1.5)

for i = 1:4
    subplot(2, 2, i)
    h1 = errorbar(v{i}, mean(DS{i}.dsindex(idx4(idx_sub), :)), std(DS{i}.dsindex(idx4(idx_sub), :))/sqrt(length(id_sub)), 'b');
    set(get(h1,'Parent'), 'XScale', 'log')
    hold on
    h2 = errorbar(v{i}, mean(DS{i}.dsindex(idx4(idx2), :)), std(DS{i}.dsindex(idx4(idx2), :))/sqrt(length(id2)), 'r');
    set(get(h2,'Parent'), 'XScale', 'log')
    xlim([v{2}(end) v{2}(1)])
    if mod(i, 2) == 1
        ylim([-.5 1])
    else
        ylim([.5 1])
    end
    title(ll{i})
    xlabel('speed')
    ylabel('DSI')
    legend([h1 h2], 'on-off DSGC', 'on-DSGC', 'location', 'southeast')
    set([h1 h2], 'DefaultLineLineWidth', 1.5)

end

% compare across light levels

FigHandle = figure;
set(FigHandle, 'Position', [1 1 900 900])
for i = 1:2
    subplot(2, 2, i)
    clear h
    h1 = semilogx(v{i}, DS{i}.dsindex(idx4(idx_sub), :), 'b');
    h(1) = h1(1);
    hold on
    h2 = semilogx(v{i+2}, DS{i+2}.dsindex(idx4(idx_sub), :), 'r');
    h(2) = h2(1);
    xlim([v{2}(end) v{2}(1)])
    if i == 1
        ylim([-.5 1])
    else
        ylim([.5 1])
    end
    title([SP{i} '  ' CT{1}])    
    xlabel('speed')
    ylabel('DSI')
    legend(h, 'NDF 3', 'NDF 0', 'location', 'southeast')
    
    subplot(2, 2, i+2)
    clear h
    h1 = semilogx(v{i}, DS{i}.dsindex(idx4(idx2), :), 'b');
    h(1) = h1(1);
    hold on
    h2 = semilogx(v{i+2}, DS{i+2}.dsindex(idx4(idx2), :), 'r');
    h(2) = h2(1);
    xlim([v{2}(end) v{2}(1)])
    if i == 1
        ylim([-.5 1])
    else
        ylim([.5 1])
    end
    title([SP{i} '  ' CT{2}])
    xlabel('speed')
    ylabel('DSI')
    legend(h, 'NDF 3', 'NDF 0', 'location', 'southeast')
end

FigHandle = figure;
set(FigHandle, 'Position', [1 1 900 900])
set(FigHandle, 'DefaultLineLineWidth', 1.5)

for i = 1:2
    subplot(2, 2, i)
    clear h
    h(1) = errorbar(v{i}, mean(DS{i}.dsindex(idx4(idx_sub), :)), std(DS{i}.dsindex(idx4(idx_sub), :))/sqrt(length(idx_sub)), 'b');
    set(get(h(1),'Parent'), 'XScale', 'log')
    hold on
    h(2) = errorbar(v{i+2}, mean(DS{i+2}.dsindex(idx4(idx_sub), :)), std(DS{i+2}.dsindex(idx4(idx_sub), :))/sqrt(length(idx_sub)), 'r');
    set(get(h(2),'Parent'), 'XScale', 'log')
    xlim([v{2}(end) v{2}(1)])
    if i == 1
        ylim([-.5 1])
    else
        ylim([.5 1])
    end
    title([SP{i} '  ' CT{1}])    
    xlabel('speed')
    ylabel('DSI')
    legend(h, 'NDF 3', 'NDF 0', 'location', 'southeast')
    
    subplot(2, 2, i+2)
    clear h
    h(1) = errorbar(v{i}, mean(DS{i}.dsindex(idx4(idx2), :)), std(DS{i}.dsindex(idx4(idx2), :))/sqrt(length(idx2)), 'b');
    set(get(h(1),'Parent'), 'XScale', 'log')
    hold on
    h(2) = errorbar(v{i+2}, mean(DS{i+2}.dsindex(idx4(idx2), :)), std(DS{i+2}.dsindex(idx4(idx2), :))/sqrt(length(idx2)), 'r');
    set(get(h(2),'Parent'), 'XScale', 'log')
    xlim([v{2}(end) v{2}(1)])
    if i == 1
        ylim([-.5 1])
    else
        ylim([.5 1])
    end
    title([SP{i} '  ' CT{2}])
    xlabel('speed')
    ylabel('DSI')
    legend(h, 'NDF 3', 'NDF 0', 'location', 'southeast')
end

%% flash responses

id_flash = intersect(id, datarun{7}.cell_ids);

idx_flash = arrayfun(@(x) find(datarun{7}.cell_ids == x), id_flash);

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:length(id_flash)
    semilogx(unique(stimuli), Pc{idx_flash(i)}, 'b')
    hold on
    pause
end
xlabel('flash strength')
ylabel('performance')
% title('flash only')
% legend('KO', 'rescue', 'WT', 'location', 'northwest')

id1_dg = id_ds_dg([1 5]);
id2_dg = id_ds_dg([2 3 4 6 7]);
idx1_flash = idx_flash([1 5]);
idx2_flash = idx_flash([2 3 4 6 7]);


raster_ds = get_ds_raster(datarun{2}, id_ds_dg);

%% compare ds with non ds

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:length(idx1_flash)
    h(1) = semilogx(unique(stimuli{3}), Pc{3}{idx1_flash(i)}, 'b');
    hold on
%     pause
end
for i = 1:length(idx2_flash)
    h(2) = semilogx(unique(stimuli{3}), Pc{3}{idx2_flash(i)}, 'c');
    hold on
%     pause
end

color = 'rgkym';
for i = 1:5
    for j = 1:length(nds_id{i}(:, 1))
    h(i+2) = semilogx(unique(stimuli{3}), Pc{3}{nds_idx{i}(j, 1)}, 'color', color(i));
    hold on
%     pause
    end
end


xlabel('flash strength')
ylabel('performance')
legend(h, 'DS 1', 'DS 2', 'ON 1', 'ON 2', 'OFF 1', 'OFF 2', 'OFF 3', 'location', 'northwest')

%% average pc
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
semilogx(unique(stimuli{3}), mean(cell2mat(Pc{3}(idx1_flash)')), 'b');
hold on
semilogx(unique(stimuli{3}), mean(cell2mat(Pc{3}(idx2_flash)')), 'c');

color = 'rgkym';
for i = 1:5
    semilogx(unique(stimuli{3}), mean(cell2mat(Pc{3}(nds_idx{i}(:, 1))'), 1), 'color', color(i));
end


xlabel('flash strength')
ylabel('performance')
legend('DS 1', 'DS 2', 'ON 1', 'ON 2', 'OFF 1', 'OFF 2', 'OFF 3', 'location', 'northwest')


%% plot ds raster
tt = datarun{2}.stimulus.params.DIRECTION*pi/180;
v = datarun{2}.stimulus.params.SPATIAL_PERIOD./datarun{2}.stimulus.params.TEMPORAL_PERIOD;
x = 3;
y = 6;
a = [8 9 3 2 1 7 13 14 15; 11 12 6 5 4 10 16 17 18];


for cc = 1:length(id_ds_dg)
    FigHandle = figure;
    set(FigHandle, 'Position', get(0, 'ScreenSize'))
    for t = 1:length(v)
        subplot(x, y, a(t, 1)) ; polar(tt, DS.rho{t}(idx_dg(cc), :))
        for i = 2:9
            subplot(x, y, a(t, i)); plot_raster(squeeze(raster_ds{cc}(1, t, i-1, :)), 0, 8)
            if i == 4
                title(id_ds_dg(cc))
            end
        end
%         name = [num2str(id_nDS_all(idx_nds_4(cc))) '_' num2str(t)];
%         screen_size = [24 12];
%         set(figure(1), 'paperpositionmode', 'auto');
%         % set(figure(1), 'PaperPosition', [-0.5 -0.25 22 10]);
%         set(gcf, 'PaperUnits', 'inch');
%         set(figure(1), 'PaperSize', screen_size);
%         print(figure(1), '-dpdf', name)
%         close

    end
    
end
%% flash summary
bin_n = 10;
n = 60;
XX = 3/n/2:3/n:3-3/n/2;

a = ceil(length((stimuli)+1)/2);
for cc =6:6 %length(id_flash(:, 1));
    cn = idx_flash(cc);
    st = 1;
    trial_n = length(raster_dark{cn});
    b = 1;
    FigHandle = figure;
    set(FigHandle, 'Position', [0, 0, 1920, 1080]);
    while st <= length(stimuli)+1
        if  st == a+1
            b = 4;
        end
        subplot(a, 6, b)
        if st > 1
            trial_n = length(raster{cn}{st-1});
        end
        for j = 1:trial_n
            if st == 1
                SpikeTime = raster_dark{cn}{j};
            else
                SpikeTime = raster{cn}{st-1}{j};
            end
            SpikeTime = SpikeTime';
            X = [SpikeTime; SpikeTime];
            Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
            line(X, Y, 'color', 'b');
            axis([0, 3, 0, trial_n]);
            hold on
        end
        if b == 1
            title(num2str(id_flash(cc)))
        end
        b = b+1;
        subplot(a, 6, b)
        if st > 1
            hist(raster_all{cn}{st-1}, XX)
        else
            hist(cell2mat(raster_dark{cn}), XX)
        end
        xlim([0 3])

        b = b+1;
        if st > 1
            corrd = corr_dark{cn}{st-1};
            corrf = corr_flash{cn}{st-1};
            xx = linspace(min([corrd corrf]), max([corrd corrf]), bin_n);
            clear theta
            theta(:, 1) = hist(corrd, xx);
            theta(:, 2) = hist(corrf, xx);
            subplot(a, 6, b)
            bar(xx, theta)
        else
            subplot(a, 6, b)
            semilogx(unique(stimuli), Pc{cn})
            ylim([0.5 1])
        end


        b = b+4;
        st = st+1;
    end

%     screen_size = [24 12];
%     set(figure(1), 'paperpositionmode', 'auto');
%     % set(figure(1), 'PaperPosition', [-0.5 -0.25 22 10]);
%     set(gcf, 'PaperUnits', 'inch');
%     set(figure(1), 'PaperSize', screen_size);
%     print(figure(1), '-dpdf', ['Flash_' num2str(nds_id_matched(cc, 1)) '_WN_' num2str(nds_id_matched(cc, 2))]);
%     close



end


%% plot sinewave grating response

x = 3;
y = 3;
j = 9;
a = [5 6 3 2 1 4 7 8 9];

for cc = 1:length(id)
    if ~isempty(raster{j}{cc})
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 1080 1080])
        h = subplot(x, y, a(1)); polar(tt{j}, DS{j}.rho{1}(cc, :));
        polar_theta_off(h)
        for i = 2:9
            subplot(x, y, a(i)); plot_raster(squeeze(raster{j}{cc}(1, 1, i-1, :)), 0, 8)
            if i == 4
                title(['cell id: ' num2str(id(cc)) '  sp:240 tp:480 sinewave grating'])
            end 
            if mod(a(i), y) == 1
                ylabel('trial number')
            end
            if a(i) > y*(x-1)
                xlabel('time (s)')
            end
        end

        
        name = num2str(id(cc));
        screen_size = [13 13];
        set(figure(1), 'paperpositionmode', 'auto');
        % set(figure(1), 'PaperPosition', [-0.5 -0.25 22 10]);
        set(gcf, 'PaperUnits', 'inch');
        set(figure(1), 'PaperSize', screen_size);
        print(figure(1), '-dpdf', name)
        close

    end
end

%% classify DSGC into subtypes (directions)
d = 2;
t = 6;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
subplot(1, 2, 1)
idx_temp = idx_sub{1};
compass(DS{d}.U{t}(idx_temp), DS{d}.V{t}(idx_temp))
color = 'brgk';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DS{d}.U{t}(idx_temp), DS{d}.V{t}(idx_temp), x, y);
    [~, I] = find(IN == 1);
    idx_dir_on{i} = idx_temp(I);
end
title('on-off DSGC')

p_direction = DS{d}.angle{t}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;


subplot(1, 2, 2)
for i = 1:dirn
    for cc = 1:length(idx_dir_oo{i})
        [xsort, seq] = sort(xx(idx_dir_oo{i}(cc), :));
        y_temp = DS{d}.rho{t}(idx_dir_oo{i}(cc), :);
        plot(xsort, y_temp(seq), color(i))
        hold on
    end
end
xlabel('direction (rad)')
ylabel('normalized response')
xlim([-pi pi])
    
FigHandle = figure;
set(FigHandle, 'Position', [1 1 900 900])
for i = 1:4
    subplot(2, 2, i)
    for j = 1:dirn
        h = semilogx(v{i}, DS{i}.dsindex(idx_dir_oo{j}, :), color(j));
        H(j) = h(1);
        hold on
    end
    xlim([v{2}(end) v{2}(1)])
    if mod(i, 2) == 0
        ylim([.5 1])
    end
    if i == 1
        legend(H, 'dorsal', 'N/T', 'ventral', 'N/T', 'location', 'southeast')
    end
    title(ll{i})
    xlabel('speed')
    ylabel('DSI')
end
FigHandle = figure;
set(FigHandle, 'Position', [1 1 900 900])
for i = 1:4
    subplot(2, 2, i)
    for j = 1:dirn
        errorbar(v{i}, mean(DS{i}.dsindex(idx_dir_oo{j}, :)), std(DS{i}.dsindex(idx_dir_oo{j}, :))/sqrt(length(idx_dir_oo{j})), color(j));
        set(gca, 'XScale', 'log')
        hold on
    end
    xlim([v{2}(end) v{2}(1)])
    if mod(i, 2) == 0
        ylim([.5 1])
    end
    if i == 1
        legend('dorsal', 'N/T', 'ventral', 'N/T', 'location', 'southeast')
    end
    title(ll{i})
    xlabel('speed')
    ylabel('DSI')
end

%% ei mosaic
% idx_dir_oo{3} = sort([idx_dir_oo{3} 38]);
ctype = {'Inferior', 'P(A)', 'Superior', 'A(P)'};
% ctype = {'Inferior', 'P(A)', 'Superior'};

FigHandle = figure;
set(FigHandle, 'Position', [1 1 900 900])
color = 'ybrgkmc';
for i = 1:4
    idx_temp = idx(idx_dir_oo{i});
    subplot(2, 2, i)
    for cc = 1:length(idx_temp)
        plot_ei_(datarun{2}.ei.eis{idx_temp(cc)}, datarun{2}.ei.position, 0, 'cutoff', 0.25, 'neg_color', color(mod(cc, 7)+1), 'scale', 2);
        hold on
    end
    title(ctype{i})
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end
%% speed tuning of different direction preferance 

figure
for i = 1:4
    subplot(2, 2, i)
    for dir = 1:4
        h1 = semilogx(v{i}, MAG_all_norm{i}(:, idx_dir_oo{dir}), color(dir));
        h(dir) = h1(1);
        hold on
    end
    xlim([v{i}(end) v{i}(1)])
    if i == 4
        legend(h, 'ventral', 'N/T', 'dorsal', 'N/T', 'location', 'southeast')
    end
    title(ll{i})
    xlabel('speed')
    ylabel('normalized response')
end

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:4
    subplot(2, 2, i)
    for dir = 1:4
        tuning = MAG_all_norm{i}(:, idx_dir_oo{dir});
        h1 = errorbar(v{i}, mean(tuning, 2), std(tuning,[], 2)/sqrt(size(tuning, 2)), color(dir));
        h(dir) = h1(1);
        set(gca, 'XScale', 'log')
        hold on
    end
    xlim([v{i}(end) v{i}(1)])
    if i == 4
        legend(h, 'ventral', 'N/T', 'dorsal', 'N/T', 'location', 'southeast')
    end
    title(ll{i})
    xlabel('speed')
    ylabel('normalized response')

end

%% ei mapping

datarun_a{1} = load_data('/Analysis/xyao/2014-07-07-0/data003/data003', opt);
datarun_a{2} = load_data('/Analysis/xyao/2014-07-07-0/data006/data006', opt);
datarun{2}.cell_types{9}.name = 'DS';
datarun{2}.cell_types{9}.cell_ids = id;
[cell_list_map, failed_cells] = map_ei(datarun{2},datarun_a{2},'master_cell_type', 'DS');
idx = find(~cellfun(@isempty,cell_list_map)); 

%% cross-correlation of nearby DS cells
% full field pulses
% superior: 4742_5146, 5896_6181, 992_1607
% A(P): 3437_3631
% P(A): 2867_3302, 7084_7382, 7668_7741, 6302_6484, 6484_7247, 7247_7668
% superior & inferior: 5896_5416, 1607_724, 6527_724, 992_724
% inferior & P(A): 724_7382
% P(A) & superior: 4069_4742, 4981_4742, 4832_5146
a = 3; b = 4;

bin_size = 0.001; % unit: s

for cc1 = 1:length(idx_dir_oo{a})
    cell_id(1) = id(idx_dir_oo{a}(cc1));
    for cc2 = 1:length(idx_dir_oo{b})
        cell_id(2) = id(idx_dir_oo{b}(cc2));


        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 1500 400])
        for i = 1:2
            if ~isempty(raster_ff{i}{id == cell_id(1)}) && ~isempty(raster_ff{i}{id == cell_id(2)})
            shift_predictor_corr = shift_predictor_ffp(datarun{i+4}, cell_id(1), cell_id(2), bin_size);
            x = (length(shift_predictor_corr)-1)/2*bin_size;
            subplot(1, 3, i)
            bar(-x:bin_size:x, shift_predictor_corr)
            xlabel('time (s)')
            ylabel('spikes/s')
            title([num2str(cell_id(1)) ' & ' num2str(cell_id(2)) '  ' LL{i}])
            xlim([-x x])
            end
        end

        subplot(1, 3, 3)

        for cc = 1:2
            idx_temp = get_cell_indices(datarun{2}, cell_id(cc));
            plot_ei_(datarun{2}.ei.eis{idx_temp}, datarun{2}.ei.position, 0, 'cutoff', 0.25, 'neg_color', color(mod(cc, 7)+1), 'scale', 2);
            hold on
        end
        
        screen_size = [18 6];
        set(figure(1), 'paperpositionmode', 'auto');
        set(gcf, 'PaperUnits', 'inch');
        set(figure(1), 'PaperSize', screen_size);
        print(figure(1), '-dpdf', [num2str(cell_id(1)) '_' num2str(cell_id(2))])
        close


    end
end

%%
figure
ct = [1 2];
for i = ct;
    idx_temp = idx(idx_dir_oo{i});
    for cc = 1:length(idx_temp)
        plot_ei_(datarun{2}.ei.eis{idx_temp(cc)}, datarun{2}.ei.position, 0, 'cutoff', 0.25, 'neg_color', color(i+1), 'scale', 2);
        hold on
    end
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end
title([ctype{ct(1)} ' & ' ctype{ct(2)}])
%% cross-correlation of nearby DS cells
% drifting grating
LL = {'NDF3', 'NDF0'};
ct = 'AIS';
bin_size = 0.001; % unit: s
maxlag = 0.1;

idx_dir = idx_dir_oo;
cell_id = [5146 5896];
for a = 1:3
    for b = a:3
        
        for cc1 = 1:length(idx_dir{a})
            cell_id(1) = id(idx_dir{a}(cc1));
            if a == b
                cc2_start = cc1+1;
            else
                cc2_start = 1;
            end
            for cc2 = cc2_start:length(idx_dir{b})
                cell_id(2) = id(idx_dir{b}(cc2));


                FigHandle = figure;
                set(FigHandle, 'Position', [1 1 1500 400])
                for i = 1:2
                    if ~isempty(raster{i*2}{id == cell_id(1)}) && ~isempty(raster{i*2}{id == cell_id(2)}) 
                        dir(1) = DS{2}.angle{5}(id == cell_id(1));
                        dir(2) = DS{2}.angle{5}(id == cell_id(2));
                        dir(dir<0) = dir(dir<0) + 2*pi;
                        dir = round(dir*4/pi)+1;
                        dir(dir == 9) = 1;

                        shift_predictor_corr(1, :) = shift_predictor_dg(datarun{2*i}, cell_id(1), cell_id(2), bin_size, maxlag, 0);
                        shift_predictor_corr(2, :) = shift_predictor_dg(datarun{2*i}, cell_id(1), cell_id(2), bin_size, maxlag, 0);
                        shift_predictor_corr = mean(shift_predictor_corr);
                        x = (length(shift_predictor_corr)-1)/2*bin_size;
                        subplot(1, 3, i)
                        bar(-x:bin_size:x, shift_predictor_corr)
                        xlabel('time (s)')
                        ylabel('# of pairs')
                        title([num2str(cell_id(1)) ' & ' num2str(cell_id(2)) '  ' LL{i}])
                        xlim([-x x])
                    end
                end

                subplot(1, 3, 3)

                for cc = 1:2
%                     subplot(1, 2, cc)
                    idx_temp = get_cell_indices(datarun{2}, cell_id(cc));
                    plot_ei_(datarun{2}.ei.eis{idx_temp}, datarun{2}.ei.position, 0, 'cutoff', 0, 'neg_color', color(mod(cc, 7)+1), 'scale', 2);
%                     hold on
                    set(gca,'xtick',[])
                    set(gca,'ytick',[])
                end




                screen_size = [18 6];
                set(figure(1), 'paperpositionmode', 'auto');
                set(gcf, 'PaperUnits', 'inch');
                set(figure(1), 'PaperSize', screen_size);
                print(figure(1), '-dpdf', [ct([a b]) '_' num2str(cell_id(1)) '_' num2str(cell_id(2))])
                close


            end
        end
    end
end

%% correlation coefficient

for i = 1:2
    for cc = 1:length(raster{2*i})
        if ~isempty(raster{2*i}{cc})
            spike_count_temp = cellfun('length', raster{2*i}{cc});
            spike_count_pspeed_temp = spike_count_temp(1, 7, :, :);
            spike_count_pspeed{i}{cc} = squeeze(mean(spike_count_pspeed_temp, 4));
            raster_pspeed_temp = squeeze(raster{2*i}{cc}(1, 7, :, :));
            
        end
    end
end