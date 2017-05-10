%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
datarun{1} = load_data('/Analysis/xyao/2013-12-09-0/data001-map/data001-map', opt);
datarun{1}.names.stimulus_path = '/Analysis/xyao/2013-12-09-0/stimuli/s01';
datarun{1} = load_stim(datarun{1}, 'user_defined_trigger_interval', 10);

datarun{2} = load_data('/Analysis/xyao/2013-12-09-0/data002-map/data002-map', opt);
datarun{2}.names.stimulus_path = '/Analysis/xyao/2013-12-09-0/stimuli/s02';
datarun{2} = load_stim(datarun{2}, 'user_defined_trigger_interval', 10);

datarun{3} = load_data('/Volumes/Suk/Analysis/xyao/2013-12-09-0/data003-map/data003-map', opt);
datarun{3}.names.stimulus_path = '/Volumes/Suk/Analysis/xyao/2013-12-09-0/stimuli/s03';
datarun{3} = load_stim(datarun{3}, 'user_defined_trigger_interval', 10);

datarun{4} = load_data('/Analysis/xyao/2013-12-09-0/data005-map/data005-map', opt);
datarun{4}.names.stimulus_path = '/Analysis/xyao/2013-12-09-0/stimuli/s05';
datarun{4} = load_stim(datarun{4}, 'user_defined_trigger_interval', 10);

datarun{5} = load_data('/Analysis/xyao/2013-12-09-0/data006-map/data006-map', opt);

opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);

datarun{6} = load_data('/Analysis/xyao/2013-12-09-0/data007-map/data007-map', opt);
datarun{6} = load_java_movie(datarun{6}, '/Volumes/lab/acquisition/movie-xml/BW-15-1-0.48-11111-40x40-60.35.xml');
datarun{6} = get_sta_summaries(datarun{6}, 'all');
datarun{6} = get_autocorrelations(datarun{6}, 'all');

%% get DSGC ids

[NumSpikesCell, ~,StimComb] = get_spikescellstim(datarun{3},datarun{3}.cell_ids,0, 1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);

% pull out DS cells

figure
plot(ds_struct.mag{3, 1}, ds_struct.mag{4, 1}, 'o')
hold on
[x, y] = ginput;
plot(x, y);
xlabel('TP 30')
ylabel('TP 60')
title('Vector sum plot')

IN = inpolygon(ds_struct.mag{3, 1}, ds_struct.mag{4, 1}, x, y);
[~, I] = find(IN == 1);
id = datarun{3}.cell_ids(I);

% DSGC only 
% [NumSpikesCell, StimComb] = get_spikescellstim(datarun{3},id,0);
% [mag dsindex magmax magave angle rho RHO theta num U V] = dscellanalysis(NumSpikesCell, StimComb);

I = zeros(length(id), length(datarun));
for i = 1:4
    I(:, i) = ismember(id, datarun{i}.cell_ids);
end

I = sum(I');
id = id(I >= 2); % DS cells that can be found in all light level


[raster, raster_p_sum, p_idx] = deal(cell(4, 1));
for i = 1:4    
    [NumSpikesCell, StimComb] = get_spikescellstim(datarun{i},id,0);
    DS{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster{i} = get_ds_raster(datarun{i}, id);
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

ll = {'NDF4', 'NDF3', 'NDF2', 'NDF0'};






%% Fourier transform
duration = 8;
bin_rate = 10000;
hist_spikes = cell(4, 1);
signal_length = duration*bin_rate;                
NFFT = 2^nextpow2(signal_length);
f = bin_rate/2*linspace(0,1,NFFT/2+1);
fft_spikes = cell(4, 1);
[DC, F1, F2] = deal(cell(4, 1));

for i = 1:4
    [DC{i}, F1{i}, F2{i}] = deal(zeros(length(id), length(tp)));
    for rgc = 1:length(id)
        if ~isempty(raster{i}{rgc})
        for t = 1:length(tp)
            spikes = floor(raster_sum{i}{rgc}{t}*bin_rate);
            tmp_binned_spikes = zeros(1, signal_length);
            tmp_binned_spikes(spikes) = 1;
            hist_spikes{i}{rgc}(t, :) = tmp_binned_spikes;
            
            f1 = 60/tp(t); %Hz
            f2 = f1*2;
            f_diff1 = f - f1;
            f_diff2 = f - f2;
            [~,f1_index] = min(abs(f_diff1));
            [~,f2_index] = min(abs(f_diff2));
            tmp_fft = fft(tmp_binned_spikes, NFFT)./ signal_length;
            fft_spikes{i}{rgc}(t, :) = 2*abs(tmp_fft(1:NFFT/2+1));
            fund_power(t) = sum(fft_spikes{i}{rgc}(t, f1_index:f1_index+2)); % f1_index+2???
            sec_power(t) = sum(fft_spikes{i}{rgc}(t, f2_index:f2_index+2));
            DC_power(t) = fft_spikes{i}{rgc}(t, 1);
        end
    % stores info for this cell into the matrix tuning curves
        F1{i}(rgc,:) = fund_power ./ max(fund_power);
        F2{i}(rgc,:) = sec_power ./ max(sec_power);
        DC{i}(rgc,:) = DC_power ./ max(DC_power);
        end
        
    end
end

% plot 

speed = sp./tp;
figure
for i = 1:4
    subplot(2, 2, i) 
    for rgc = 1:length(id)
    if ~isempty(raster{i}{rgc})
    semilogx(speed, DC{i}(rgc,:), 'b')
    hold on
    end
    end
    xlabel('speed')
    title(ll{i})
end

figure

for i = 1:4
    subplot(2, 2, i) 
    for rgc = 1:length(id)
    if ~isempty(raster{i}{rgc})
    semilogx(speed, F1{i}(rgc,:), 'r')
    hold on
    end
    end
    xlabel('speed')
    title(ll{i})

end

figure

for i = 1:4
    subplot(2, 2, i) 
    for rgc = 1:length(id)
    if ~isempty(raster{i}{rgc})
    semilogx(speed, F2{i}(rgc,:), 'k')
    hold on
    end
    end
    xlabel('speed')
    title(ll{i})

end

      
    for i = 1:4
        figure
        for rgc = 1:length(id)
        if ~isempty(raster{i}{rgc})
            subplot(5, 8, rgc)
            semilogx(speed, DC{i}(rgc,:), 'b')
            hold on
            semilogx(speed, F1{i}(rgc,:), 'r')
            semilogx(speed, F2{i}(rgc,:), 'k')
            title(num2str(id(rgc)))
        end
        end
end
    
    for i = 1:4
        figure
        for rgc = 1:length(id1)
        if ~isempty(raster{i}{idx1(rgc)})
            subplot(5, 6, rgc)
            semilogx(speed, DC{i}(idx1(rgc),:), 'b')
            hold on
            semilogx(speed, F1{i}(idx1(rgc),:), 'r')
            semilogx(speed, F2{i}(idx1(rgc),:), 'k')
            title(num2str(id1(rgc)))
        end
        end
    end
    
    for i = 1:4
        figure
        for rgc = 1:length(id2)
        if ~isempty(raster{i}{idx2(rgc)})
            subplot(3, 5, rgc)
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
    DC_tuning_1 = DC{i}(idx1, :);
    temp = sum(DC_tuning_1, 2);
    [~, idx] = find(temp' ~= 0);
    DC_tuning_1 = DC_tuning_1(idx, :);
    DC_avg_1(i, :) = mean(DC_tuning_1);
    DC_ste_1(i, :) = std(DC_tuning_1)/sqrt(size(DC_tuning_1, 1));
    
    f1_tuning_1 = F1{i}(idx1, :);
    temp = sum(f1_tuning_1, 2);
    [~, idx] = find(temp' ~= 0);
    f1_tuning_1 = f1_tuning_1(idx, :);
    f1_avg_1(i, :) = mean(f1_tuning_1);
    f1_ste_1(i, :) = std(f1_tuning_1)/sqrt(size(f1_tuning_1, 1));
    
    f2_tuning_1 = F2{i}(idx1, :);
    temp = sum(f2_tuning_1, 2);
    [~, idx] = find(temp' ~= 0);
    f2_tuning_1 = f2_tuning_1(idx, :);
    f2_avg_1(i, :) = mean(f2_tuning_1);
    f2_ste_1(i, :) = std(f2_tuning_1)/sqrt(size(f2_tuning_1, 1));
end
    

% type 2
for i = 1:4
    DC_tuning_2 = DC{i}(idx2, :);
    temp = sum(DC_tuning_2, 2);
    [~, idx] = find(temp' ~= 0);
    DC_tuning_2 = DC_tuning_2(idx, :);
    DC_avg_2(i, :) = mean(DC_tuning_2);
    DC_ste_2(i, :) = std(DC_tuning_2)/sqrt(size(DC_tuning_2, 1));
    
    f1_tuning_2 = F1{i}(idx2, :);
    temp = sum(f1_tuning_2, 2);
    [~, idx] = find(temp' ~= 0);
    f1_tuning_2 = f1_tuning_2(idx, :);
    f1_avg_2(i, :) = mean(f1_tuning_2);
    f1_ste_2(i, :) = std(f1_tuning_2)/sqrt(size(f1_tuning_2, 1));
    
    f2_tuning_2 = F2{i}(idx2, :);
    temp = sum(f2_tuning_2, 2);
    [~, idx] = find(temp' ~= 0);
    f2_tuning_2 = f2_tuning_2(idx, :);
    f2_avg_2(i, :) = mean(f2_tuning_2);
    f2_ste_2(i, :) = std(f2_tuning_2)/sqrt(size(f2_tuning_2, 1));
end

% plot
figure
for i = 1:4
    subplot(2, 2, i)
    h1 = errorbar(speed, DC_avg_1(i, :), DC_ste_1(i, :), 'b');
    set(get(h1,'Parent'), 'XScale', 'log')
    hold on
    h2 = errorbar(speed, f1_avg_1(i, :), f1_ste_1(i, :), 'r');
    set(get(h2,'Parent'), 'XScale', 'log')
    h3 = errorbar(speed, f2_avg_1(i, :), f2_ste_1(i, :), 'k');
    set(get(h3,'Parent'), 'XScale', 'log')
    xlabel('speed')
    xlim([0.1 10])
    title(ll{i})
    legend('DC', 'F1', 'F2', 'location', 'northwest')
end

figure
for i = 1:4
    subplot(2, 2, i)
    h1 = errorbar(speed, DC_avg_2(i, :), DC_ste_2(i, :), 'b');
    set(get(h1,'Parent'), 'XScale', 'log')
    hold on
    h2 = errorbar(speed, f1_avg_2(i, :), f1_ste_2(i, :), 'r');
    set(get(h2,'Parent'), 'XScale', 'log')
    h3 = errorbar(speed, f2_avg_2(i, :), f2_ste_2(i, :), 'k');
    set(get(h3,'Parent'), 'XScale', 'log')
    xlabel('speed')
    xlim([0.1 10])
    title(ll{i})
    legend('DC', 'F1', 'F2')
end

%% plot raster and psth

ll = {'NDF 4', 'NDF 3', 'NDF 2', 'NDF 0'};
bin = 320;
spike_rate = cell(4, 1);

for rgc = 1:length(id)
    for i = 1:4
        if isempty(raster{i}{rgc}) == 0
%         FigHandle = figure;
%         set(FigHandle, 'Position', get(0, 'ScreenSize'))

        for t = 1:length(tp)
%             subplot(length(tp), 2, 2*t-1)
%             plot_raster(squeeze(raster{i}{rgc}(1, t, pdirection(rgc), :)), 0, 8)
%             subplot(length(tp), 2, 2*t)
%             hist(raster_sum{i}{rgc}{t}, bin);
            a = hist(raster_sum{i}{rgc}{t}, bin);
            spike_rate{i}{rgc}(t, :) = a(round(bin/8):end);
%             xlim([0 8])
%             if t == 1
%                 title(ll{i})
%             end

        end
        end
    end
end

%% calculate speed tuning curves

speed = sp./tp;

tuning_curve = cell(4, 1);
for i = 1:4
    for rgc = 1:length(id)
        if isempty(raster{i}{rgc}) == 0
        for t = 1:length(tp)
            n = tp(t)/60*bin/8;
            r = n;
            m = [];
            while r <= size(spike_rate{1}{1}, 2)
                m = [m max(spike_rate{i}{rgc}(t, r-n+1:r))];
                r = r+n;
            end
            tuning_curve{i}{rgc}(t) = mean(m);
        end
        tuning_curve{i}{rgc} = tuning_curve{i}{rgc}/max(tuning_curve{i}{rgc});
        end
    end
end

figure
for i = 1:4
    subplot(2, 2, i)
    for rgc = 1:length(tuning_curve{i})
        if isempty(tuning_curve{i}{rgc}) == 0
            semilogx(tp, tuning_curve{i}{rgc})
            hold on
        end
    end
    xlabel('temporal period')
    ylabel('normalized max firing rate')
    xlim([min(tp) max(tp)])
    title(ll{i})
end


FigHandle = figure;
set(FigHandle, 'Position', [1 1 800 800])

[~,scores,~,~] = princomp(cell2mat(tuning_curve{3}'));
subplot(2, 2, 1); 
plot(scores(:, 1), scores(:, 2), 'o')
hold on
[x, y] = ginput;
plot(x, y)

xlabel('1st Principal Component')
ylabel('2nd Principal Component')
subplot(2, 2, 2); 
plot(scores(:, 1), scores(:, 3), 'o')
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
subplot(2, 2, 3); 
plot(scores(:, 2), scores(:, 3), 'o')
xlabel('2nd Principal Component')
ylabel('3rd Principal Component')
IN = inpolygon(scores(:, 1), scores(:, 2), x, y);
[~, idx1] = find(IN' == 1);
id1 = id(idx1);
[~, idx2] = find(IN' == 0);
id2 = id(idx2);


% plot speed tunning curves of cells in same type

figure
for i = 1:4
    subplot(2, 2, i)
    for cc = 1:length(id1)
        if ~isempty(raster{i}{idx1(cc)})
            m = tuning_curve{i}{idx1(cc)};
            semilogx(tp, m)
            hold on
        end
    end
    xlabel('Temporal Period')
    ylabel('normalized max firing rate')
    title(ll{i})
    xlim([6 360])
end

figure
for i = 1:4
    subplot(2, 2, i)
    for cc = 1:length(id2)
        if ~isempty(raster{i}{idx2(cc)})
            m = tuning_curve{i}{idx2(cc)};
            semilogx(tp, m)
            hold on
        end
    end
    xlabel('Temporal Period')
    ylabel('normalized max firing rate')
    title(ll{i})
    xlim([6 360])
end

% plot average tunning curve
for i = 1:4
    mag_temp = [];
    for cc = 1:length(id1)
        if ~isempty(raster{i}{idx1(cc)})
            mag_temp = [mag_temp; tuning_curve{i}{idx1(cc)}];
        end
    end
    tuning_avg(:, 1, i) = mean(mag_temp);
    
    mag_temp = [];
    for cc = 1:length(id2)
        if ~isempty(raster{i}{idx2(cc)})
            mag_temp = [mag_temp; tuning_curve{i}{idx2(cc)}];
        end
    end
    tuning_avg(:, 2, i) = mean(mag_temp);
end

FigHandle = figure;
set(FigHandle, 'Position', [1 1 1800 500])
subplot(1, 3, 1)
semilogx(tp, squeeze(tuning_avg(:, 1, :)))
xlabel('temporal period')
ylabel('normalized response')
xlim([6 360])
legend(ll, 'location', 'northwest')
title('type 1')

subplot(1, 3, 2)
semilogx(tp, squeeze(tuning_avg(:, 2, :)))
xlabel('temporal period')
ylabel('normalized response')
xlim([6 360])
legend(ll, 'location', 'northwest')
title('type 2')

subplot(1, 3, 3)
semilogx(tp, squeeze(tuning_avg(:, 1, :)))
hold on
semilogx(tp, squeeze(tuning_avg(:, 2, :)), 'x-')

xlabel('temporal period')
ylabel('normalized response')
xlim([6 360])
legend(ll, 'location', 'northwest')

% vector sum plot of individual type
[NumSpikesCell, StimComb] = get_spikescellstim(datarun{3},id2,0);
[mag dsindex magmax magave angle rho RHO theta num U V] = dscellanalysis(NumSpikesCell, StimComb);


%% plot cell summary

x = 6;
y = 6;
name = cell(length(id)*length(tp), 1);
for cc = 15:15
    for t = 1:length(tp)
        FigHandle = figure;
        set(FigHandle, 'Position', get(0, 'ScreenSize'))
        if ~isempty(raster{1}{cc})
        subplot(x, y, 8); polar(tt, DS{1}.rho{t}(cc, :))
        subplot(x, y, 9); plot_raster(squeeze(raster{1}{cc}(1, t, 1, :)), 0, 8)
        subplot(x, y, 3); plot_raster(squeeze(raster{1}{cc}(1, t, 2, :)), 0, 8)
        subplot(x, y, 2); plot_raster(squeeze(raster{1}{cc}(1, t, 3, :)), 0, 8); title('NDF 4')
        subplot(x, y, 1); plot_raster(squeeze(raster{1}{cc}(1, t, 4, :)), 0, 8)
        subplot(x, y, 7); plot_raster(squeeze(raster{1}{cc}(1, t, 5, :)), 0, 8)
        subplot(x, y, 13); plot_raster(squeeze(raster{1}{cc}(1, t, 6, :)), 0, 8)
        subplot(x, y, 14); plot_raster(squeeze(raster{1}{cc}(1, t, 7, :)), 0, 8)
        subplot(x, y, 15); plot_raster(squeeze(raster{1}{cc}(1, t, 8, :)), 0, 8)
        end
        
        if ~isempty(raster{2}{cc})
        subplot(x, y, 11); polar(tt, DS{2}.rho{t}(cc, :))
        subplot(x, y, 12); plot_raster(squeeze(raster{2}{cc}(1, t, 1, :)), 0, 8)
        subplot(x, y, 6); plot_raster(squeeze(raster{2}{cc}(1, t, 2, :)), 0, 8)
        subplot(x, y, 5); plot_raster(squeeze(raster{2}{cc}(1, t, 3, :)), 0, 8); title('NDF 3')
        subplot(x, y, 4); plot_raster(squeeze(raster{2}{cc}(1, t, 4, :)), 0, 8)
        subplot(x, y, 10); plot_raster(squeeze(raster{2}{cc}(1, t, 5, :)), 0, 8)
        subplot(x, y, 16); plot_raster(squeeze(raster{2}{cc}(1, t, 6, :)), 0, 8)
        subplot(x, y, 17); plot_raster(squeeze(raster{2}{cc}(1, t, 7, :)), 0, 8)
        subplot(x, y, 18); plot_raster(squeeze(raster{2}{cc}(1, t, 8, :)), 0, 8)
        end
        
        if ~isempty(raster{3}{cc})
        subplot(x, y, 26); polar(tt, DS{3}.rho{t}(cc, :))
        subplot(x, y, 27); plot_raster(squeeze(raster{3}{cc}(1, t, 1, :)), 0, 8)
        subplot(x, y, 21); plot_raster(squeeze(raster{3}{cc}(1, t, 2, :)), 0, 8)
        subplot(x, y, 20); plot_raster(squeeze(raster{3}{cc}(1, t, 3, :)), 0, 8); title('NDF 2')
        subplot(x, y, 19); plot_raster(squeeze(raster{3}{cc}(1, t, 4, :)), 0, 8)
        subplot(x, y, 25); plot_raster(squeeze(raster{3}{cc}(1, t, 5, :)), 0, 8)
        subplot(x, y, 31); plot_raster(squeeze(raster{3}{cc}(1, t, 6, :)), 0, 8)
        subplot(x, y, 32); plot_raster(squeeze(raster{3}{cc}(1, t, 7, :)), 0, 8)
        subplot(x, y, 33); plot_raster(squeeze(raster{3}{cc}(1, t, 8, :)), 0, 8)
        end
        
        if ~isempty(raster{4}{cc})
        subplot(x, y, 29); polar(tt, DS{4}.rho{t}(cc, :))
        subplot(x, y, 30); plot_raster(squeeze(raster{4}{cc}(1, t, 1, :)), 0, 8)
        subplot(x, y, 24); plot_raster(squeeze(raster{4}{cc}(1, t, 2, :)), 0, 8)
        subplot(x, y, 23); plot_raster(squeeze(raster{4}{cc}(1, t, 3, :)), 0, 8); title('NDF 0')
        subplot(x, y, 22); plot_raster(squeeze(raster{4}{cc}(1, t, 4, :)), 0, 8)
        subplot(x, y, 28); plot_raster(squeeze(raster{4}{cc}(1, t, 5, :)), 0, 8)
        subplot(x, y, 34); plot_raster(squeeze(raster{4}{cc}(1, t, 6, :)), 0, 8)
        subplot(x, y, 35); plot_raster(squeeze(raster{4}{cc}(1, t, 7, :)), 0, 8)
        subplot(x, y, 36); plot_raster(squeeze(raster{4}{cc}(1, t, 8, :)), 0, 8)
        end
        
        
        name{length(tp)*(cc-1)+t} = [num2str(id(cc)) '_' num2str(tp(t))];
    end
end

%% plot single cell tuning curve
% choose the preferred direction under NDF 2, TP 30, as the preferred
% direction for plotting


p_direction = DS{3}.angle{3}';
I = repmat(1:8, length(id), 1);
I = I + repmat(round(p_direction*4/pi)+4, 1, length(tt));
I = mod(I, length(tt));
I(I == 0) = length(tt);

RHO_all = cell(4, 1);
MAG_all = cell(4, 1);
% pn = repmat(tp', 1, length(id));
% pn = 8*60.35./pn;

for i = 1:4
    RHO_temp = reshape(cell2mat(DS{i}.RHO'), length(id), length(tt), length(tp));
    for cc = 1:length(id)
        r = RHO_temp(cc, :, :);
        RHO_all{i}(cc, :, :) = r(:, I(cc, :), :);
    end
    MAG_all{i} = cell2mat(DS{i}.MAG); %./pn;
end


x = tt-3/4*pi;
y = log10(tp);
[XX YY] = meshgrid(x, y);

scr_size = get(0, 'ScreenSize');
name = cell(length(id), 1);
for cc = 1:24 %length(id)
%         FigHandle = figure;
%         set(FigHandle, 'Position', scr_size([1 2 4 4]))
%     if ~isempty(raster{1}{cc})
%         subplot(2, 2, 1)
%         mesh(XX, YY, squeeze(RHO_all{1}(cc, :, :))')
%         hidden off
%         xlabel('direction/rad')
%         ylabel('log(temporal period)')
%         zlabel('spike number')
%         title('NDF 4')
%     end
%     
%     if ~isempty(raster{2}{cc})
%         subplot(2, 2, 2)
%         mesh(XX, YY, squeeze(RHO_all{2}(cc, :, :))')
%         hidden off
%         title('NDF 3')
%     end
%     
%     if ~isempty(raster{3}{cc})
%         subplot(2, 2, 3)
%         mesh(XX, YY, squeeze(RHO_all{3}(cc, :, :))')
%         hidden off
%         title('NDF 2')
%     end
% 
%     if ~isempty(raster{4}{cc})
%         subplot(2, 2, 4)
%         mesh(XX, YY, squeeze(RHO_all{4}(cc, :, :))')
%         hidden off
%         title('NDF 0')
%     end
% 

        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 1520 1080])
    if ~isempty(raster{1}{cc})
        subplot(2, 3, 1)
        for i = 1:length(tp)
            plot(x, RHO_all{1}(cc, :, i), 'color', [1/length(tp)*i 0 1-1/length(tp)*i])
            hold on            
        end
        title('NDF 4')
        xlabel('direction/rad')
        ylabel('spike number')
        legend('TP 360', 'location', 'northwest')
    end
    
    if ~isempty(raster{2}{cc})
        subplot(2, 3, 2)
        for i = 1:length(tp)
            plot(x, RHO_all{2}(cc, :, i), 'color', [1/length(tp)*i 0 1-1/length(tp)*i])
            hold on
        end
        title('NDF 3')
        xlabel('direction/rad')
        ylabel('spike number')
        legend('TP 360', 'location', 'northwest')
    end
    
    if ~isempty(raster{3}{cc})
        subplot(2, 3, 4)
        for i = 1:length(tp)
            plot(x, RHO_all{3}(cc, :, i), 'color', [1/length(tp)*i 0 1-1/length(tp)*i])
            hold on
        end
        title('NDF 2')
        xlabel('direction/rad')
        ylabel('spike number')
        legend('TP 360', 'location', 'northwest')
    end
    
    if ~isempty(raster{4}{cc})
        subplot(2, 3, 5)
        for i = 1:length(tp)
            plot(x, RHO_all{4}(cc, :, i), 'color', [1/length(tp)*i 0 1-1/length(tp)*i])
            hold on
        end
        title('NDF 0')
        xlabel('direction/rad')
        ylabel('spike number')
        legend('TP 360', 'location', 'northwest')
    end
    
    subplot(2, 3, 3)
    color = {'b', 'r', 'k', 'g'};
    lgd = cell(4, 1);
    h = [];
    for i = 1:4
        if ~isempty(raster{i}{cc})
            h = [h semilogx(tp, MAG_all{i}(:, cc), color{i})];
            lgd{i} = ll{i};
            hold on
        end
    end
        lgd = lgd(~cellfun('isempty',lgd));
        xlabel('Temporal Period')
        ylabel('response')
        legend(h, lgd, 'location', 'northwest')
        
%     name{2*cc-1} = [num2str(id(cc)) '_3d_tunning'];
    name{cc} = [num2str(id(cc)) '_2d_tunning'];

end

% plot speed tunning curves of all cells in one figure

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

% classify cells into two groups according to speed tunning
% pca
mag_pca = MAG_all_norm{3};
% mag_pca = MAG_all_norm{1}(:, idx4);
mag_pca = mag_pca';

FigHandle = figure;
set(FigHandle, 'Position', [1 1 800 800])

[~,scores,~,~] = princomp(mag_pca);
subplot(2, 2, 1); 
plot(scores(:, 1), scores(:, 2), 'o')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
subplot(2, 2, 2); 
plot(scores(:, 1), scores(:, 3), 'o')
hold on
[x, y] = ginput;
plot(x, y)
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
subplot(2, 2, 3); 
plot(scores(:, 2), scores(:, 3), 'o')
xlabel('2nd Principal Component')
ylabel('3rd Principal Component')


IN = inpolygon(scores(:, 1), scores(:, 3), x, y);
[~, idx1] = find(IN' == 1);
id1 = id(idx1);
[~, idx2] = find(IN' == 0);
id2 = id(idx2);

% plot speed tunning curves of cells in same type
speed = datarun{1}.stimulus.params.SPATIAL_PERIOD./datarun{1}.stimulus.params.TEMPORAL_PERIOD;
figure
for i = 1:4
    subplot(2, 2, i)
    for cc = 1:length(id1)
        if ~isempty(raster{i}{idx1(cc)})
            m = MAG_all_norm{i}(:, idx1(cc));
            semilogx(speed, m)
            hold on
        end
    end
    xlabel('speed')
    ylabel('Response')
    title(ll{i})
    xlim([min(speed) max(speed)])
end

figure
for i = 1:4
    subplot(2, 2, i)
    for cc = 1:length(id2)
        if ~isempty(raster{i}{idx2(cc)})
            m = MAG_all_norm{i}(:, idx2(cc));
            semilogx(speed, m)
            hold on
        end
    end
    xlabel('speed')
    ylabel('Response')
    title(ll{i})
    xlim([min(speed) max(speed)])
end

% plot average tunning curve
for i = 1:4
    mag_temp = [];
    for cc = 1:length(id1)
        if ~isempty(raster{i}{idx1(cc)})
            mag_temp = [mag_temp MAG_all_norm{i}(:, idx1(cc))];
        end
    end
    tunning_avg(:, 1, i) = mean(mag_temp, 2);
    tunning_ste(:, 1, i) = std(mag_temp')/sqrt(length(id1));
    mag_temp = [];
    for cc = 1:length(id2)
        if ~isempty(raster{i}{idx2(cc)})
            mag_temp = [mag_temp MAG_all_norm{i}(:, idx2(cc))];
        end
    end
    tunning_avg(:, 2, i) = mean(mag_temp, 2);
    tunning_ste(:, 2, i) = std(mag_temp')/sqrt(length(id2));

end


color = 'brk';
figure
set(gcf, 'Position', [1 1 800 800])
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:4
    subplot(2, 2, i)
    for j = 1:2
%     semilogx(v, tunning_avg(:, j, i)/max(tunning_avg(:, j, i)), color(j))
    errorbar(speed, tunning_avg(:, j, i), tunning_ste(:, j, i), color(j))
    hold on
    end
    set(gca, 'XScale', 'log')
    title(ll{i})
    xlim([min(speed) max(speed)])
    xlabel('Speed')
    ylabel('Normalized Response')
end
legend('type 1', 'type 2')


%
figure;
set(gcf, 'Position', [1 1 600 500])
set(gcf, 'DefaultLineLineWidth', 1.5)
% subplot(1, 3, 1)
errorbar(repmat(speed', 1, 4), squeeze(tunning_avg(:, 1, :)), squeeze(tunning_ste(:, 1, :)))
set(gca, 'XScale', 'log')
xlabel('speed')
ylabel('normalized response')
xlim([min(speed) max(speed)])
legend(ll)
title('on-off DSGC')

% subplot(1, 3, 2)
figure;
set(gcf, 'Position', [1 1 600 500])
set(gcf, 'DefaultLineLineWidth', 1.5)

errorbar(repmat(speed', 1, 4), squeeze(tunning_avg(:, 2, :)), squeeze(tunning_ste(:, 2, :)))
set(gca, 'XScale', 'log')
xlabel('speed')
ylabel('normalized response')
xlim([min(speed) max(speed)])
legend(ll)
title('on DSGC')

subplot(1, 3, 3)
semilogx(tp, squeeze(tunning_avg(:, 1, :)))
hold on
semilogx(tp, squeeze(tunning_avg(:, 2, :)), 'x-')

xlabel('temporal period')
ylabel('normalized response')
xlim([6 360])
legend(ll, 'location', 'northwest')

% vector sum plot of individual type
[NumSpikesCell, StimComb] = get_spikescellstim(datarun{3},id2,0);
[mag dsindex magmax magave angle rho RHO theta num U V] = dscellanalysis(NumSpikesCell, StimComb);

%% full field pulses

n = length(datarun{5}.triggers)/4;
trigger = zeros(n, 1);
for i = 1:n
    trigger(i) = datarun{5}.triggers(4*(i-1)+1);
end

trigger = trigger - 1.5;

raster_ff = cell(length(id), 1);
for i = 1:length(id)
    if ~isempty(intersect(datarun{5}.cell_ids, id(i)))
        idx = get_cell_indices(datarun{5}, id(i));
        raster_ff{i, 1} = get_raster(datarun{5}.spikes{idx}, trigger, 'plot', false);
    end
end

%% plot wn & full field flashes

a = [0 1.5 1.5 4.5 4.5 7.5 7.5 10.5 10.5 12]';
b = [1 1 2 2 1 1 0 0 1 1]';
for i = 36:length(id) 
    FigHandle = figure;
    set(FigHandle, 'Position', [1 1 800 800])

    subplot(4, 2, [1 3])
    plot_rf(datarun{6}, id(i), 'title', false, 'scale', 8)
    subplot(4, 2, 5)
    plot_time_course(datarun{6}, id(i), 'clear', false, 'figure', -1)
    subplot(4, 2, 7)
    plot_autocorrelograms(datarun{6}, id(i), 'foa', -1)
    if ~isempty(intersect(datarun{5}.cell_ids, id(i)))
    h = subplot(4, 2, [2 4 6 8]);
    for j = 1:length(raster_ff{i})
        SpikeTime = raster_ff{i}{j};
        SpikeTime = SpikeTime';
        X = [SpikeTime; SpikeTime];
        Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
        line(X, Y, 'color', 'b');
        xlim([0 12])
        hold on
    end   
    line(a, b+27)
    xlabel('time/sec')
    set(h,'yticklabel',[]);
    end
end


id1wn = get_cell_ids(datarun{6}, 'DS type1');
id2wn = get_cell_ids(datarun{6}, 'DS type2');

figure;

for i = 1:length(id1wn)
    plot_autocorrelograms(datarun{6}, id1wn(i), 'foa', -1, 'clear_fig', false, 'normalize', true, 'line_color', 'b')
    hold on
end

for i = 1:length(id2wn)
    plot_autocorrelograms(datarun{6}, id2wn(i), 'foa', -1, 'clear_fig', false, 'normalize', true, 'line_color', 'r')
    hold on
end

%% DSI

for i = 1:4
    DS{i}.dsindex = cell2mat(DS{i}.dsindex');
    dsi_temp = DS{i}.dsindex(idx1, :);
    z = zeros(size(dsi_temp));
    i_out = find(sum(dsi_temp == z, 2) == length(tp));
    dsi_temp(i_out, :) = [];
    sample_size(i, 1) = size(dsi_temp, 1);
    dsi{i}{1} = dsi_temp;
    
    dsi_temp = DS{i}.dsindex(idx2, :);
    z = zeros(size(dsi_temp));
    i_out = find(sum(dsi_temp == z, 2) == length(tp));
    dsi_temp(i_out, :) = [];
    sample_size(i, 2) = size(dsi_temp, 1);
    dsi{i}{2} = dsi_temp;
end

% compare across cell types

FigHandle = figure;
set(FigHandle, 'Position', [1 1 900 1080])
for i = 1:4
    subplot(2, 2, i)
    h1 = semilogx(speed, dsi{i}{1}', 'b');
    h(1) = h1(1);
    hold on
    h2 = semilogx(speed, dsi{i}{2}', 'r');
    h(2) = h2(1);
    xlim([speed(end) speed(1)])
%     if mod(i, 2) == 0
%         ylim([0 1])
%     end
    title(ll{i})
    xlabel('speed')
    ylabel('DSI')
    legend(h, 'on-off DSGC', 'on-DSGC', 'location', 'southeast')
end

FigHandle = figure;
set(FigHandle, 'Position', [1 1 900 1080])
set(FigHandle, 'DefaultLineLineWidth', 1.5)

for i = 1:4
    subplot(2, 2, i)
    h1 = errorbar(speed, mean(dsi{i}{1}), std(dsi{i}{1})/sqrt(sample_size(i, 1)), 'b');
    set(get(h1,'Parent'), 'XScale', 'log')
    hold on
    h2 = errorbar(speed, mean(dsi{i}{2}), std(dsi{i}{2})/sqrt(sample_size(i, 2)), 'r');
    set(get(h2,'Parent'), 'XScale', 'log')
    xlim([speed(end) speed(1)])
    if i >= 3
        ylim([0 1])
    end
    title(ll{i})
    xlabel('speed')
    ylabel('DSI')
    legend([h1 h2], 'on-off DSGC', 'on-DSGC', 'location', 'southeast')
    set([h1 h2], 'DefaultLineLineWidth', 1.5)

end

% compare across light levels

FigHandle = figure;
set(FigHandle, 'Position', [1 1 1480 540])
set(FigHandle, 'DefaultLineLineWidth', 1.5)

subplot(1, 2, 1)
clear h
h(1) = errorbar(speed, mean(dsi{1}{1}), std(dsi{1}{1})/sqrt(sample_size(1, 1)), 'b');
set(get(h(1),'Parent'), 'XScale', 'log')
hold on
h(2) = errorbar(speed, mean(dsi{2}{1}), std(dsi{2}{1})/sqrt(sample_size(2, 1)), 'r');
set(get(h(2),'Parent'), 'XScale', 'log')
h(3) = errorbar(speed, mean(dsi{3}{1}), std(dsi{3}{1})/sqrt(sample_size(3, 1)), 'k');
set(get(h(3),'Parent'), 'XScale', 'log')
h(4) = errorbar(speed, mean(dsi{4}{1}), std(dsi{4}{1})/sqrt(sample_size(4, 1)), 'c');
set(get(h(4),'Parent'), 'XScale', 'log')
xlim([speed(end) speed(1)])
%     if i == 1
%         ylim([-.5 1])
%     else
%         ylim([.5 1])
%     end
title('on-off DSGC')    
xlabel('speed')
ylabel('DSI')
legend(h, 'NDF 4', 'NDF 3', 'NDF 2', 'NDF 0', 'location', 'southeast')

subplot(1, 2, 2)
clear h
h(1) = errorbar(speed, mean(dsi{1}{2}), std(dsi{1}{2})/sqrt(sample_size(1, 2)), 'b');
set(get(h(1),'Parent'), 'XScale', 'log')
hold on
h(2) = errorbar(speed, mean(dsi{2}{2}), std(dsi{2}{2})/sqrt(sample_size(2, 2)), 'r');
set(get(h(2),'Parent'), 'XScale', 'log')
h(3) = errorbar(speed, mean(dsi{3}{2}), std(dsi{3}{2})/sqrt(sample_size(3, 2)), 'k');
set(get(h(3),'Parent'), 'XScale', 'log')
h(4) = errorbar(speed, mean(dsi{4}{2}), std(dsi{4}{2})/sqrt(sample_size(4, 2)), 'c');
set(get(h(4),'Parent'), 'XScale', 'log')
xlim([speed(end) speed(1)])
%     if i == 1
%         ylim([-.5 1])
%     else
%         ylim([.5 1])
%     end
title('on-DSGC')
xlabel('speed')
ylabel('DSI')
legend(h, 'NDF 4', 'NDF 3', 'NDF 2', 'NDF 0', 'location', 'southeast')


