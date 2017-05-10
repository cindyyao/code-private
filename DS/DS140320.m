%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
datarun{1} = load_data('/Analysis/xyao/2014-03-20-0/data000/data000', opt);
datarun{1}.names.stimulus_path = '/Analysis/xyao/2014-03-20-0/stimuli/s00';
datarun{1} = load_stim(datarun{1}, 'user_defined_trigger_interval', 10);

datarun{2} = load_data('/Analysis/xyao/2014-03-20-0/data001-map/data001-map', opt);
datarun{2}.names.stimulus_path = '/Analysis/xyao/2014-03-20-0/stimuli/s01';
datarun{2} = load_stim(datarun{2}, 'user_defined_trigger_interval', 10);

datarun{3} = load_data('/Analysis/xyao/2014-03-20-0/data002-map/data002-map', opt);
datarun{3}.names.stimulus_path = '/Analysis/xyao/2014-03-20-0/stimuli/s02';
datarun{3} = load_stim(datarun{3}, 'user_defined_trigger_interval', 10);

datarun{4} = load_data('/Analysis/xyao/2014-03-20-0/data003-map/data003-map', opt);
datarun{4}.names.stimulus_path = '/Analysis/xyao/2014-03-20-0/stimuli/s03';
datarun{4} = load_stim(datarun{4}, 'user_defined_trigger_interval', 10);

datarun{5} = load_data('/Analysis/xyao/2014-03-20-0/data004-map/data004-map', opt);
datarun{5}.names.stimulus_path = '/Analysis/xyao/2014-03-20-0/stimuli/s04';
datarun{5} = load_stim(datarun{5}, 'user_defined_trigger_interval', 10);

datarun{6} = load_data('/Analysis/xyao/2014-03-20-0/data005-map/data005-map', opt);
datarun{6}.names.stimulus_path = '/Analysis/xyao/2014-03-20-0/stimuli/s05';
datarun{6} = load_stim(datarun{6}, 'user_defined_trigger_interval', 10);

% datarun{7} = load_data('/Analysis/xyao/2014-03-20-0/data007-map/data007-map', opt);
% datarun{7} = load_sta(datarun{7});
% datarun{7} = get_autocorrelations(datarun{7}, intersect(id, datarun{7}.cell_ids));
% 
% 
datarun{8} = load_data('/Analysis/xyao/2014-03-20-0/data006-map/data006-map', opt);





[NumSpikesCell, StimComb] = get_spikescellstim(datarun{6},datarun{6}.cell_ids,0);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);

% pull out DS cells

figure
plot(ds_struct.mag{2, 1}, ds_struct.mag{3, 1}, 'o')
title('data005 after mapping')
xlabel('TP 60')
ylabel('TP 120')

hold on
[x, y] = ginput;
plot(x, y);

IN = inpolygon(ds_struct.mag{2, 1}, ds_struct.mag{3, 1}, x, y);
[~, I] = find(IN == 1);
id = datarun{6}.cell_ids(I);


I = zeros(length(id), 6);
for i = 1:6
    I(:, i) = ismember(id, datarun{i}.cell_ids);
end

idx3th = find(I(:, 3) == 1);

I = sum(I');
id3 = id(I >= 3); % DS cells that can be found in more than 3 light levels
id6 = id(I == 6); % DS cells that can be found in all light levels

idx6 = arrayfun(@(x) find(id == x), id6);

[raster, raster_p_sum, p_idx] = deal(cell(6, 1));
for i = 1:6    
    [NumSpikesCell, StimComb] = get_spikescellstim(datarun{i},id,0);
    DS{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster{i} = get_ds_raster(datarun{i}, id);
end

%% 
param_p = 5; % choose which params to use to calculate prefer direction indices 
MAG_all_norm = cell(6, 1);
max_r = cell(6, 1);
norm_max_r = cell(6, 1);

for i = 1:6
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster{i}, DS{i}.angle{param_p});
    MAG_all_norm{i} = normalize_MAG(DS{i});
    rep = datarun{i}.stimulus.repetitions;
    max_r{i} = max_firing_rate(raster_p_sum{i}, 0.2, 8)/rep;
    norm_max_r{i} = max_r{i}./repmat(max(max_r{i}, [], 2), 1, size(max_r{i}, 2));
end

ll = {'NDF3 SP60', 'NDF3 SP120', 'NDF3 SP240', 'NDF0 SP60', 'NDF3 SP120', 'NDF0 SP240'};

%% non-DS types

id_nDS = cell(8, 1);
cell_type = {'ON transient', 'ON brist transient', 'OFF transient1', 'OFF transient2', ...
    'OFF slow', 'OFF sustained', 'OFF slow sustained', 'OFF slow sustained2'};
for i = 1:8
    id_nDS{i} = get_cell_ids(datarun{7}, cell_type{i});
end

id_nDS_all = cell2mat(id_nDS');
I = zeros(length(id_nDS_all), 6);
for i = 1:6
    I(:, i) = ismember(id_nDS_all, datarun{i}.cell_ids);
end
% I = sum(I');
% idx_nds_6 = find(I == 6);
I = sum(I(:, 3:6)');
idx_nds_4 = find(I == 4);


for i = 1:6    
    [NumSpikesCell, StimComb] = get_spikescellstim(datarun{i}, id_nDS_all,0);
    [mag MAG dsindex magmax magave angle rho RHO theta num U V] = dscellanalysis(NumSpikesCell, StimComb);
    nDS{i} = v2struct(mag, MAG, dsindex, magmax, magave, angle, rho, RHO, theta, num, U, V);
end

% sort the sequence of directions
for i = 1:6
    for j = 1:9
        [theta_seq, I_seq] = sort(nDS{i}.theta{j}(1, :));
        r = nDS{i}.rho{j};
        R = nDS{i}.RHO{j};
        rho_seq = r(:, I_seq);
        RHO_seq = R(:, I_seq);
        nDS{i}.rho{j} = rho_seq;
        nDS{i}.RHO{j} = RHO_seq;
    end
end


% get raster
for i = 1:6
    raster_nDS{i} = get_ds_raster(datarun{i}, cell2mat(id_nDS'));
end
%% plot

for cc = 1:1 %length(id)
    plot_ds_raster(DS, raster, cc, id(cc), ll, 2, 3, 1)
end


%% plot single cell tuning curve
% choose the preferred direction under NDF 2, TP 30, as the preferred
% direction for plotting

tt = datarun{1}.stimulus.params.DIRECTION*pi/180;
v = datarun{1}.stimulus.params.SPATIAL_PERIOD./datarun{1}.stimulus.params.TEMPORAL_PERIOD;

p_direction = DS{6}.angle{3}';
I = repmat(1:length(tt), length(id), 1);
I = I + repmat(round(p_direction*4/pi)+4, 1, length(tt));
I = mod(I, length(tt));
I(I == 0) = length(tt);

RHO_all = cell(6, 1);
MAG_all = cell(6, 1);
ll = {'NDF3 SP60', 'NDF3 SP120', 'NDF3 SP240', 'NDF0 SP60', 'NDF0 SP120', 'NDF0 SP240'};
% pn = repmat(tp', 1, length(id));
% pn = 8*60./pn;

for i = 1:6
    RHO_temp = reshape(cell2mat(DS{i}.RHO'), length(id), length(tt), length(v));
    for cc = 1:length(id)
        r = RHO_temp(cc, :, :);
        RHO_all{i}(cc, :, :) = r(:, I(cc, :), :);
    end
    MAG_all{i} = cell2mat(DS{i}.MAG)%./pn;
end


x = tt-3/4*pi;
y = log10(v);
[XX YY] = meshgrid(x, y);

scr_size = get(0, 'ScreenSize');
name = cell(length(id), 1);
for cc = 1:5 %length(id)
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
        subplot(2, 4, 1)
        for i = 1:length(v)
            plot(x, RHO_all{1}(cc, :, i), 'color', [1/length(v)*i 0 1-1/length(v)*i])
            hold on            
        end
        title('NDF3 SP60')
        xlabel('direction/rad')
        ylabel('spike number')
        legend('fast', 'location', 'northwest')
    end
    
    if ~isempty(raster{2}{cc})
        subplot(2, 4, 2)
        for i = 1:length(v)
            plot(x, RHO_all{2}(cc, :, i), 'color', [1/length(v)*i 0 1-1/length(v)*i])
            hold on
        end
        title('NDF3 SP120')
        xlabel('direction/rad')
        ylabel('spike number')
        legend('fast', 'location', 'northwest')
    end
    
    if ~isempty(raster{3}{cc})
        subplot(2, 4, 3)
        for i = 1:length(v)
            plot(x, RHO_all{3}(cc, :, i), 'color', [1/length(v)*i 0 1-1/length(v)*i])
            hold on
        end
        title('NDF3 SP240')
        xlabel('direction/rad')
        ylabel('spike number')
        legend('fast', 'location', 'northwest')
    end
    
    if ~isempty(raster{4}{cc})
        subplot(2, 4, 5)
        for i = 1:length(v)
            plot(x, RHO_all{4}(cc, :, i), 'color', [1/length(v)*i 0 1-1/length(v)*i])
            hold on
        end
        title('NDF0 SP60')
        xlabel('direction/rad')
        ylabel('spike number')
        legend('fast', 'location', 'northwest')
    end
    
    if ~isempty(raster{5}{cc})
        subplot(2, 4, 6)
        for i = 1:length(v)
            plot(x, RHO_all{5}(cc, :, i), 'color', [1/length(v)*i 0 1-1/length(v)*i])
            hold on
        end
        title('NDF0 SP120')
        xlabel('direction/rad')
        ylabel('spike number')
        legend('fast', 'location', 'northwest')
    end
    
    if ~isempty(raster{6}{cc})
        subplot(2, 4, 7)
        for i = 1:length(v)
            plot(x, RHO_all{6}(cc, :, i), 'color', [1/length(v)*i 0 1-1/length(v)*i])
            hold on
        end
        title('NDF0 SP240')
        xlabel('direction/rad')
        ylabel('spike number')
        legend('fast', 'location', 'northwest')
    end
    
    subplot(2, 4, 4)
    color = {'b', 'r', 'k', 'g', 'c', 'y'};
    lgd = cell(6, 1);
    h = [];
    for i = 1:6
        if ~isempty(raster{i}{cc})
            h = [h semilogx(v, MAG_all{i}(:, cc), color{i})];
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

MAG_all_norm = cell(6, 1);
figure
for i = 1:6
    v = datarun{i}.stimulus.params.SPATIAL_PERIOD./datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    subplot(2, 3, i)
    MAG_all_norm{i} = zeros(length(v), length(id));
    for cc = 1:length(id)
        if ~isempty(raster{i}{cc})
            m = MAG_all{i}(:, cc)/max(MAG_all{i}(:, cc));
            semilogx(v, m)
            MAG_all_norm{i}(:, cc) = m;
            hold on
        end
    end
    xlabel('speed')
    ylabel('Response')
    title(ll{i})
    xlim([v(end) v(1)])
end

% classify cells into two groups according to speed tunning
% pca
mag_pca = zeros(length(id), length(v));
for cc = 1:length(id)
    mag_temp1 = MAG_all{6}(:, cc)/max(MAG_all{6}(:, cc));
    mag_pca(cc, :) = mag_temp1;
end

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 800])

[~,scores,~,~] = princomp(mag_pca);
subplot(2, 1, 1); 
plot(scores(:, 1), scores(:, 2), 'o')
hold on
[x, y] = ginput;
plot(x, y)
IN = inpolygon(scores(:, 1), scores(:, 2), x, y);
[~, idx1] = find(IN' == 1);
id1 = id(idx1);
[x, y] = ginput;
plot(x, y)
IN = inpolygon(scores(:, 1), scores(:, 2), x, y);
[~, idx2] = find(IN' == 1);
id2 = id(idx2);

% [x, y] = ginput;
% plot(x, y)
% IN = inpolygon(scores(:, 1), scores(:, 2), x, y);
% [~, idx3] = find(IN' == 1);
% id3 = id(idx3);


xlabel('1st Principal Component')
ylabel('2nd Principal Component')
subplot(2, 1, 2); 
plot(scores(:, 1), scores(:, 3), 'o')
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
ylim([-0.4 0.3])
IN = inpolygon(scores(:, 1), scores(:, 3), x, y);
[~, idx1] = find(IN' == 1);
id1 = id(idx1);
[~, idx2] = find(IN' == 0);
id2 = id(idx2);

% plot speed tunning curves of cells in same type

figure
for i = 1:6
    subplot(2, 3, i)
    for cc = 1:length(id1)
        if ~isempty(raster{i}{idx1(cc)})
            m = MAG_all_norm{i}(:, idx1(cc));
            semilogx(v, m)
            hold on
        end
    end
    xlabel('speed')
    ylabel('Response')
    title(ll{i})
    xlim([v(end) v(1)])
end

figure
for i = 1:6
    subplot(2, 3, i)
    for cc = 1:length(id2)
        if ~isempty(raster{i}{idx2(cc)})
            m = MAG_all_norm{i}(:, idx2(cc));
            semilogx(v, m)
            hold on
        end
    end
    xlabel('Temporal Period')
    ylabel('Response')
    title(ll{i})
    xlim([v(end) v(1)])
end

figure
for i = 1:6
    subplot(2, 3, i)
    for cc = 1:length(id3)
        if ~isempty(raster{i}{idx3(cc)})
            m = MAG_all_norm{i}(:, idx3(cc));
            semilogx(v, m)
            hold on
        end
    end
    xlabel('Temporal Period')
    ylabel('Response')
    title(ll{i})
    xlim([v(end) v(1)])
end


% plot average tunning curve
for i = 1:6
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

    
%     mag_temp = [];
%     for cc = 1:length(id3)
%         if ~isempty(raster{i}{idx3(cc)})
%             mag_temp = [mag_temp MAG_all_norm{i}(:, idx3(cc))];
%         end
%     end
%     tunning_avg(:, 3, i) = mean(mag_temp, 2);
end

for i = 1:3
FigHandle = figure;
% set(FigHandle, 'Position', [1 1 1800 500])
subplot(2, 2, i)
semilogx(v, squeeze(tunning_avg(:, i, :)))
xlabel('temporal period')
ylabel('normalized response')
xlim([v(end) v(1)])
legend(ll, 'location', 'northwest')
title(['type ' num2str(i)])
end

subplot(2, 2, 4)
semilogx(v, squeeze(tunning_avg(:, 1, :)))
hold on
semilogx(v, squeeze(tunning_avg(:, 2, :)), 'x-')
semilogx(v, squeeze(tunning_avg(:, 3, :)), 'o-')
xlabel('speed')
ylabel('normalized response')
xlim([v(end) v(1)])
legend(ll, 'location', 'northwest')

%% average tuning curve
color = 'brk';
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
ii = [2 3 5 6];
for i = 1:4
    subplot(2, 2, i)
    for j = 1:2
%     semilogx(v, tunning_avg(:, j, i)/max(tunning_avg(:, j, i)), color(j))
    errorbar(v, tunning_avg(:, j, ii(i)), tunning_ste(:, j, ii(i)), color(j))
    hold on
    end
    set(gca, 'XScale', 'log')
    title(ll{ii(i)})
    xlim([min(v) max(v)])
end
legend('on-off DSGC', 'on DSGC', 'location', 'southeast')

% vector sum plot of individual type
[NumSpikesCell, StimComb] = get_spikescellstim(datarun{6},id1,0);
[mag dsindex magmax magave angle rho RHO theta num U V] = dscellanalysis(NumSpikesCell, StimComb);

%% compare across light level
SP = {'SP 60', 'SP 120', 'SP 240'};
CT = {'on-off DSGC', 'on DSGC'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 2:3
    for j = 1:2
        subplot(2, 2, i+2*(j-1)-1)
        errorbar(v, tunning_avg(:, j, i), tunning_ste(:, j, i), 'b')
        hold on
        errorbar(v, tunning_avg(:, j, i+3), tunning_ste(:, j, i+3), 'r')
        title([SP{i} '  ' CT{j}])
        set(gca, 'XScale', 'log')
        xlim([min(v) max(v)])
    end
end
legend('NDF 3', 'NDF 0', 'location', 'southeast')


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

%%
n_ffp = 1;

[raster_ff, raster_ff_all] = deal(cell(n_ffp, 1));
[raster_ff{1}, raster_ff_all{1}] = get_ffp_raster(datarun{8}, id, 3);


for i = 36:36 %length(id) 
%     if ~isempty(raster_ff{1}{i}) || ~isempty(raster_ff{2}{i})
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 400 400])
        for d = 1:n_ffp
        subplot(1, n_ffp, d)
        plot_ffp(raster_ff{d}, raster_ff_all{d}, i, 3)
        title([num2str(id(i))]) % ' ' ll{d}])
        end
        
%         print_close(1, [24, 12], num2str(ds_id(i)))
%     end
end

%% maxmium firing rate
% psth

ll = {'NDF 4', 'NDF 3', 'NDF 2', 'NDF 0'};
bin = 320;
spike_rate = cell(6, 1);

for rgc = 1:length(id)
    for i = 1:6
        if isempty(raster{i}{rgc}) == 0
%         FigHandle = figure;
%         set(FigHandle, 'Position', get(0, 'ScreenSize'))

        for time = 1:length(v)
%             subplot(length(tp), 2, 2*t-1)
%             plot_raster(squeeze(raster{i}{rgc}(1, t, pdirection(rgc), :)), 0, 8)
%             subplot(length(tp), 2, 2*t)
%             hist(raster_sum{i}{rgc}{t}, bin);
            a = hist(raster_sum{i}{rgc}{time}, bin);
            spike_rate{i}{rgc}(time, :) = a;
%             spike_rate{i}{rgc}(t, :) = a(round(bin/8):end);
%             xlim([0 8])
%             if t == 1
%                 title(ll{i})
%             end

        end
        end
    end
end

% calculate speed tuning curves
tuning_curve = cell(6, 1);
for i = 1:6
    for rgc = 1:length(id)
        if isempty(raster{i}{rgc}) == 0
        for time = 1:length(tp)
            n = tp(time)/60*bin/8;
            r = n;
            m = [];
            while r <= size(spike_rate{1}{1}, 2)
                m = [m max(spike_rate{i}{rgc}(time, r-n+1:r))];
                r = r+n;
            end
            tuning_curve{i}{rgc}(time) = mean(m);
        end
        tuning_curve{i}{rgc} = tuning_curve{i}{rgc}/max(tuning_curve{i}{rgc});
        end
    end
end

figure
for i = 1:6
    subplot(2, 3, i)
    for rgc = 1:length(tuning_curve{i})
        if isempty(tuning_curve{i}{rgc}) == 0
            semilogx(v, tuning_curve{i}{rgc})
            hold on
        end
    end
    xlabel('temporal period')
    ylabel('normalized max firing rate')
    xlim([min(v) max(v)])
%     title(ll{i})
end


% FigHandle = figure;
% set(FigHandle, 'Position', [1 1 800 800])
% 
% [~,scores,~,~] = princomp(cell2mat(tuning_curve{3}'));
% subplot(2, 2, 1); 
% plot(scores(:, 1), scores(:, 2), 'o')
% hold on
% [x, y] = ginput;
% plot(x, y)
% 
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
% subplot(2, 2, 2); 
% plot(scores(:, 1), scores(:, 3), 'o')
% xlabel('1st Principal Component')
% ylabel('3rd Principal Component')
% subplot(2, 2, 3); 
% plot(scores(:, 2), scores(:, 3), 'o')
% xlabel('2nd Principal Component')
% ylabel('3rd Principal Component')
% IN = inpolygon(scores(:, 1), scores(:, 2), x, y);
% [~, idx1] = find(IN' == 1);
% id1 = id(idx1);
% [~, idx2] = find(IN' == 0);
% id2 = id(idx2);
% 
% 
% % plot speed tunning curves of cells in same type
% 
% figure
% for i = 1:4
%     subplot(2, 2, i)
%     for cc = 1:length(id1)
%         if ~isempty(raster{i}{idx1(cc)})
%             m = tuning_curve{i}{idx1(cc)};
%             semilogx(tp, m)
%             hold on
%         end
%     end
%     xlabel('Temporal Period')
%     ylabel('normalized max firing rate')
%     title(ll{i})
%     xlim([6 360])
% end
% 
% figure
% for i = 1:4
%     subplot(2, 2, i)
%     for cc = 1:length(id2)
%         if ~isempty(raster{i}{idx2(cc)})
%             m = tuning_curve{i}{idx2(cc)};
%             semilogx(tp, m)
%             hold on
%         end
%     end
%     xlabel('Temporal Period')
%     ylabel('normalized max firing rate')
%     title(ll{i})
%     xlim([6 360])
% end
% 
% % plot average tunning curve
% for i = 1:4
%     mag_temp = [];
%     for cc = 1:length(id1)
%         if ~isempty(raster{i}{idx1(cc)})
%             mag_temp = [mag_temp; tuning_curve{i}{idx1(cc)}];
%         end
%     end
%     tuning_avg(:, 1, i) = mean(mag_temp);
%     
%     mag_temp = [];
%     for cc = 1:length(id2)
%         if ~isempty(raster{i}{idx2(cc)})
%             mag_temp = [mag_temp; tuning_curve{i}{idx2(cc)}];
%         end
%     end
%     tuning_avg(:, 2, i) = mean(mag_temp);
% end
% 
% FigHandle = figure;
% set(FigHandle, 'Position', [1 1 1800 500])
% subplot(1, 3, 1)
% semilogx(tp, squeeze(tuning_avg(:, 1, :)))
% xlabel('temporal period')
% ylabel('normalized response')
% xlim([6 360])
% legend(ll, 'location', 'northwest')
% title('type 1')
% 
% subplot(1, 3, 2)
% semilogx(tp, squeeze(tuning_avg(:, 2, :)))
% xlabel('temporal period')
% ylabel('normalized response')
% xlim([6 360])
% legend(ll, 'location', 'northwest')
% title('type 2')
% 
% subplot(1, 3, 3)
% semilogx(tp, squeeze(tuning_avg(:, 1, :)))
% hold on
% semilogx(tp, squeeze(tuning_avg(:, 2, :)), 'x-')
% 
% xlabel('temporal period')
% ylabel('normalized response')
% xlim([6 360])
% legend(ll, 'location', 'northwest')
% 
% % vector sum plot of individual type
% [NumSpikesCell, StimComb] = get_spikescellstim(datarun{3},id2,0);
% [mag dsindex magmax magave angle rho RHO theta num U V] = dscellanalysis(NumSpikesCell, StimComb);


%% frequency analysis

duration = 8;
bin_rate = 10000;
hist_spikes = cell(6, 1);
signal_length = duration*bin_rate;
NFFT = 2^nextpow2(signal_length);
f = bin_rate/2*linspace(0,1,NFFT/2+1);
fft_spikes = cell(6, 1);
[DC, F1, F2] = deal(cell(6, 1));

for i = 1:6
    tp = datarun{i}.stimulus.params.TEMPORAL_PERIOD;
    [DC{i}, F1{i}, F2{i}] = deal(zeros(length(id), length(tp)-4));
    for rgc = 1:length(id)
        if ~isempty(raster{i}{rgc})
        for time = 1:length(tp)-4
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
%         F1{i}(rgc,:) = fund_power ./ max(fund_power);
%         F2{i}(rgc,:) = sec_power ./ max(sec_power);
%         DC{i}(rgc,:) = DC_power ./ max(DC_power);
        
        F1{i}(rgc,:) = fund_power ./ max(DC_power);
        F2{i}(rgc,:) = sec_power ./ max(DC_power);
        DC{i}(rgc,:) = DC_power ./ max(DC_power);
        clear fund_power sec_power DC_power
        end
        
    end
    ratio{i} = F2{i}./F1{i};
    for ct = 1:4
        ratio_dir{ct}{i} = ratio{i}(idx_dir_oo{ct}, :);
        ratio_dir{ct}{i} = exciseRows_empty(ratio_dir{ct}{i});
        ratio_dir_mean(i, ct, :) = mean(ratio_dir{ct}{i}, 1);
        ratio_dir_ste(i, ct, :) = std(ratio_dir{ct}{i}, [], 1)/sqrt(size(ratio_dir{ct}{i}, 1));
    end
    for ct = 1:2
        ratio_dir_on{ct}{i} = ratio{i}(idx_dir_on{ct}, :);
        ratio_dir_on{ct}{i} = exciseRows_empty(ratio_dir_on{ct}{i});
        ratio_dir_mean_on(i, ct, :) = mean(ratio_dir_on{ct}{i}, 1);
        ratio_dir_ste_on(i, ct, :) = std(ratio_dir_on{ct}{i}, [], 1)/sqrt(size(ratio_dir_on{ct}{i}, 1));
    end
%     ratio{i} = exciseRows_empty(ratio{i});
%     ratio_mean(i, :) = mean(ratio{i});
%     ratio_ste(i, :) = std(ratio{i})/sqrt(size(ratio{i}, 1));
end

% plot 

speed = v;
figure
for i = 1:6
    subplot(2, 3, i) 
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
for i = 1:6
    subplot(2, 3, i) 
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
for i = 1:6
    subplot(2, 3, i) 
    for rgc = 1:length(id)
    if ~isempty(raster{i}{rgc})
    semilogx(speed, F2{i}(rgc,:), 'k')
    hold on
    end
    end
    xlabel('speed')
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
    for rgc = 1:length(id1)
    if ~isempty(raster{i}{idx1(rgc)})
        subplot(6, 8, rgc)
        semilogx(speed, DC{i}(idx1(rgc),:), 'b')
        hold on
        semilogx(speed, F1{i}(idx1(rgc),:), 'r')
        semilogx(speed, F2{i}(idx1(rgc),:), 'k')
        title(num2str(id1(rgc)))
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
for i = 1:6
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
for i = 1:6
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
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:6
    subplot(2, 3, i)
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
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:6
    subplot(2, 3, i)
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

%% F1/F2 plot
figure
set(gcf, 'DefaultLineLineWidth', 1.5)

% type 1
% subplot(2, 2, 1)
% ratio_temp = ratio{2}(idx1, :);
% avg_temp = nanmean(ratio_temp);
% ste_temp = nanstd(ratio_temp)/sqrt(length(idx1));
% errorbar(speed, avg_temp, ste_temp, 'b')
% hold on
% ratio_temp = ratio{5}(idx1, :);
% avg_temp = nanmean(ratio_temp);
% ste_temp = nanstd(ratio_temp)/sqrt(length(idx1));
% errorbar(speed, avg_temp, ste_temp, 'r')
% set(gca, 'XScale', 'log')
% xlabel('speed')
% ylabel('F2/F1')
% title('SP120  on-off DSGC')
% xlim([min(speed) max(speed)])
% legend('NDF 3', 'NDF 0')
speed = speed(1:5);
subplot(1, 2, 1)
ratio_temp = exciseRows_empty(ratio{3}(idx1, :));
avg_temp = mean(ratio_temp);
ste_temp = std(ratio_temp)/sqrt(size(ratio_temp, 1));
errorbar(speed, avg_temp, ste_temp, 'b')
hold on
ratio_temp = exciseRows_empty(ratio{6}(idx1, :));
avg_temp = mean(ratio_temp);
ste_temp = std(ratio_temp)/sqrt(size(ratio_temp, 1));
errorbar(speed, avg_temp, ste_temp, 'r')
set(gca, 'XScale', 'log')
xlabel('speed')
ylabel('F2/F1')
title('SP240  on-off DSGC')
xlim([min(speed)-1 max(speed)+1])


% type 2
% subplot(2, 2, 3)
% ratio_temp = ratio{2}(idx2, :);
% avg_temp = nanmean(ratio_temp);
% ste_temp = nanstd(ratio_temp)/sqrt(length(idx2));
% errorbar(speed, avg_temp, ste_temp, 'b')
% hold on
% ratio_temp = ratio{5}(idx2, :);
% avg_temp = nanmean(ratio_temp);
% ste_temp = nanstd(ratio_temp)/sqrt(length(idx2));
% errorbar(speed, avg_temp, ste_temp, 'r')
% set(gca, 'XScale', 'log')
% xlabel('speed')
% ylabel('F2/F1')
% title('SP120  on DSGC')
% xlim([min(speed) max(speed)])

subplot(1, 2, 2)
ratio_temp = exciseRows_empty(ratio{3}(idx2, :));
avg_temp = mean(ratio_temp);
ste_temp = std(ratio_temp)/sqrt(size(ratio_temp, 1));
errorbar(speed, avg_temp, ste_temp, 'b')
hold on
ratio_temp = exciseRows_empty(ratio{6}(idx2, :));
avg_temp = mean(ratio_temp);
ste_temp = std(ratio_temp)/sqrt(size(ratio_temp, 1));
errorbar(speed, avg_temp, ste_temp, 'r')
set(gca, 'XScale', 'log')
xlabel('speed')
ylabel('F2/F1')
title('SP240  on DSGC')
xlim([min(speed)-1 max(speed)+1])

% individual direction
% on-off
figure
for ct = 1:4
    subplot(2, 2, ct)
    errorbar(v2, ratio_dir_mean(3, ct, :), ratio_dir_ste(3, ct, :), 'b')
    hold on
    errorbar(v2, ratio_dir_mean(6, ct, :), ratio_dir_ste(6, ct, :), 'r')
    set(gca, 'XScale', 'log')
end
    
figure
for ct = 1:2
    subplot(1, 2, ct)
    errorbar(v2, ratio_dir_mean_on(3, ct, :), ratio_dir_ste_on(3, ct, :), 'b')
    hold on
    errorbar(v2, ratio_dir_mean_on(6, ct, :), ratio_dir_ste_on(6, ct, :), 'r')
    set(gca, 'XScale', 'log')
end
%% full field pulses

n = length(datarun{8}.triggers)/4;
trigger = zeros(n, 1);
for i = 1:n
    trigger(i) = datarun{8}.triggers(4*(i-1)+1);
end

trigger = trigger - 1.5;

raster_ff = cell(length(id), 1);
for i = 1:length(id)
    if ~isempty(intersect(datarun{8}.cell_ids, id(i)))
        idx = get_cell_indices(datarun{8}, id(i));
        raster_ff{i, 1} = get_raster(datarun{8}.spikes{idx}, trigger, 'plot', false);
    end
end


%% plot wn & full field flashes

a = [0 1.5 1.5 4.5 4.5 7.5 7.5 10.5 10.5 12]';
b = [1 1 2 2 1 1 0 0 1 1]';
id1ff = intersect(id1, datarun{8}.cell_ids);
id2ff = intersect(id2, datarun{8}.cell_ids);
idx1ff = arrayfun(@(x) find(id == x), id1ff);
idx2ff = arrayfun(@(x) find(id == x), id2ff);

for i = 1:length(id1ff) 
    FigHandle = figure;
    set(FigHandle, 'Position', [1 1 400 800])
%     set(FigHandle,'yticklabel',[]);

    for j = 1:length(raster_ff{idx1ff(i)})
        SpikeTime = raster_ff{idx1ff(i)}{j};
        SpikeTime = SpikeTime';
        X = [SpikeTime; SpikeTime];
        Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
        line(X, Y, 'color', 'b');
        xlim([0 12])
        hold on
    end   
    line(a, b+27)
    xlabel('time/sec')
end

for i = 1:length(id2ff) 
    FigHandle = figure;
    set(FigHandle, 'Position', [1 1 400 800])
%     set(FigHandle,'yticklabel',[]);

    for j = 1:length(raster_ff{idx2ff(i)})
        SpikeTime = raster_ff{idx2ff(i)}{j};
        SpikeTime = SpikeTime';
        X = [SpikeTime; SpikeTime];
        Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
        line(X, Y, 'color', 'b');
        xlim([0 12])
        hold on
    end   
    line(a, b+27)
    xlabel('time/sec')
end


id1wn = intersect(id1, datarun{7}.cell_ids);
id2wn = intersect(id2, datarun{7}.cell_ids);

figure;

for i = 1:length(id1wn)
    plot_autocorrelograms(datarun{7}, id1wn(i), 'foa', -1, 'clear_fig', false, 'normalize', true, 'line_color', 'b')
    hold on
end

% figure

for i = 1:length(id2wn)
    plot_autocorrelograms(datarun{7}, id2wn(i), 'foa', -1, 'clear_fig', false, 'normalize', true, 'line_color', 'r')
    hold on
end


    

%% 3D plot
speed = v;
MAG_all_3{1}(:, 1, :) = MAG_all{1};
MAG_all_3{1}(:, 2, :) = MAG_all{2};
MAG_all_3{1}(:, 3, :) = MAG_all{3};
MAG_all_3{2}(:, 1, :) = MAG_all{4};
MAG_all_3{2}(:, 2, :) = MAG_all{5};
MAG_all_3{2}(:, 3, :) = MAG_all{6};

for i = 1:2
    for cc = 1:length(id)
        MAG_all_norm_3{i}(:, :, cc) = MAG_all_3{i}(:, :, cc)/max(max(MAG_all_3{i}(:, :, cc)));
    end
end

[XX YY] = meshgrid([60 120 240], speed);
mesh(XX, YY, MAG_all_norm_3{2}(:, :, 2))


I = zeros(length(id), 6);
for i = 1:6
    I(:, i) = ismember(id, datarun{i}.cell_ids);
end

% % NDF 0

I0 = I(:, 4:6);
I0 = sum(I0');
idx01 = intersect(find(I0 == 3), idx1);
MAG_avg_3{2}{1} = mean(MAG_all_norm_3{2}(:, :, idx01), 3);
% mesh
figure
mesh(XX, YY, MAG_avg_3{2}{1})
xlabel('SP')
ylabel('speed')
zlabel('normalized response')
title('NDF 0  on-off DSGC')

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
subplot(1, 2, 1)
semilogx(speed, MAG_avg_3{2}{1})
legend('SP 60', 'SP 120', 'SP 240', 'location', 'northwest')
xlabel('speed')
ylabel('response')
title('NDF 0  on-off DSGC')

subplot(1, 2, 2)
semilogx(speed, MAG_avg_3{2}{1}./repmat(max(MAG_avg_3{2}{1}), length(speed), 1))
xlabel('speed')
ylabel('normalized response')
title('NDF 0  on-off DSGC')

idx02 = intersect(find(I0 == 3), idx2);
MAG_avg_3{2}{2} = mean(MAG_all_norm_3{2}(:, :, idx02), 3);
% mesh
figure
mesh(XX, YY, MAG_avg_3{2}{2})
xlabel('SP')
ylabel('speed')
zlabel('normalized response')
title('NDF 0  ON DSGC')

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
subplot(1, 2, 1)
semilogx(speed, MAG_avg_3{2}{2})
legend('SP 60', 'SP 120', 'SP 240')
xlabel('speed')
ylabel('response')
title('NDF 0  ON DSGC')

subplot(1, 2, 2)
semilogx(speed, MAG_avg_3{2}{2}./repmat(max(MAG_avg_3{2}{2}), length(speed), 1))
xlabel('speed')
ylabel('normalized response')
title('NDF 0  ON DSGC')

% % NDF 3
I3 = I(:, 2:3);
I3 = sum(I3');
idx31 = intersect(find(I3 == 2), idx1);
MAG_avg_3{1}{1} = mean(MAG_all_norm_3{1}(:, 2:3, idx31), 3);
[XX YY] = meshgrid([120 240], speed);

% mesh
figure
mesh(XX, YY, MAG_avg_3{1}{1})
xlabel('SP')
ylabel('speed')
zlabel('normalized response')
title('NDF 3  on-off DSGC')

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
subplot(1, 2, 1)
semilogx(speed, MAG_avg_3{1}{1})
legend('SP 120', 'SP 240', 'location', 'northwest')
xlabel('speed')
ylabel('response')
title('NDF 3  on-off DSGC')

subplot(1, 2, 2)
semilogx(speed, MAG_avg_3{1}{1}./repmat(max(MAG_avg_3{1}{1}), length(speed), 1))
xlabel('speed')
ylabel('normalized response')
title('NDF 3  on-off DSGC')

% on
idx32 = intersect(find(I3 == 2), idx2);
MAG_avg_3{1}{2} = mean(MAG_all_norm_3{1}(:, 2:3, idx32), 3);
[XX YY] = meshgrid([120 240], speed);

% mesh
figure
mesh(XX, YY, MAG_avg_3{1}{2})
xlabel('SP')
ylabel('speed')
zlabel('normalized response')
title('NDF 3  on DSGC')

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
subplot(1, 2, 1)
semilogx(speed, MAG_avg_3{1}{2})
legend('SP 120', 'SP 240', 'location', 'northwest')
xlabel('speed')
ylabel('response')
title('NDF 3  on DSGC')

subplot(1, 2, 2)
semilogx(speed, MAG_avg_3{1}{2}./repmat(max(MAG_avg_3{1}{2}), length(speed), 1))
xlabel('speed')
ylabel('normalized response')
title('NDF 3  on DSGC')

%% phase

%psth for ndf 3, sp 240, tp 480, prefer direction

XX = 0.25:0.5:7.75;

for i = 1:length(id)
    temp = raster_sum{3}{i};
    if ~isempty(temp)
        psth{i} = hist(temp{5}, XX);
    end
end

figure
for i = 1:length(idx3th)
    subplot(5, 8, i)
    hist(raster_sum{3}{idx3th(i)}{5}, XX)
end

% calculate time that start firing        
for i = 1:length(id)
    temp = psth{i};
    if ~isempty(temp)
        time(i) = 0.5*(find(temp > max(temp)/2, 1) - 1);
%         time(i) = 0.5*(find(temp == max(temp), 1) - 1);
    end
end

% get rf location
for i = 1:length(id6)
    if ~isempty(intersect(datarun{7}.cell_ids, id6(i)))
        plot_rf(datarun{7}, id6(i), 'foa', 0)
        [x, y] = ginput;
        rf{i}(:, 1) = 15*x;
        rf{i}(:, 2) = 15*y;
    end
end

% for i = 1:length(id6)
%     if ~isempty(intersect(datarun{7}.cell_ids, id6(i)))
%         rf{i}(:, 1) = rf{i}(:, 1)+100;
%     end
% end

% get stimulus frame
spec.x_start = 0; spec.x_end = 800;
spec.y_start = 0; spec.y_end = 600;
spec.spatial_period = 240;
spec.temporal_period = 480;

figure
for i = 1:length(id6)
    if ~isempty(intersect(datarun{7}.cell_ids, id6(i)))
        spec.orientation = (pdirection(idx6(i)) - 1) * 45;
        [framesin, framecos, tdrift] = calc_drifting_grating_frame_intensities_jf(spec);
        ph = tdrift(mod(time(idx6(i))*60, spec.temporal_period)+1);
        im = framesin*cos(ph) + framecos*sin(ph);
        im = (sign(im(:, 101:700)) + 1)/2;
%         im = im(:, end:-1:1);
        subplot(3, 4, i)
        colormap gray
        imagesc(im)
        set(gca,'YDir','normal')
        hold on 
        plot_rf_fit(datarun{7}, id6(i), 'scale', 15, 'edge_color', [1 0 0], 'reverse_y', 0)
        title(num2str(spec.orientation));
    end
end
        
        
figure
for i = 1:length(id6)
    if ~isempty(intersect(datarun{7}.cell_ids, id6(i)))
        subplot(3, 4, i)
        plot_rf(datarun{7}, id6(i), 'title', false)
    end
end

%% DSI

for i = 1:6
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
set(FigHandle, 'Position', [1 1 1480 1080])
for i = 1:6
    subplot(2, 3, i)
    h1 = semilogx(v, dsi{i}{1}', 'b');
    h(1) = h1(1);
    hold on
    h2 = semilogx(v, dsi{i}{2}', 'r');
    h(2) = h2(1);
    xlim([v(end) v(1)])
    if mod(i, 3) == 0
        ylim([0 1])
    end
    title(ll{i})
    xlabel('speed')
    ylabel('DSI')
    legend(h, 'on-off DSGC', 'on-DSGC', 'location', 'southeast')
end

FigHandle = figure;
set(FigHandle, 'Position', [1 1 1480 1080])
set(FigHandle, 'DefaultLineLineWidth', 1.5)

for i = 1:6
    subplot(2, 3, i)
    h1 = errorbar(v, mean(dsi{i}{1}), std(dsi{i}{1})/sqrt(sample_size(i, 1)), 'b');
    set(get(h1,'Parent'), 'XScale', 'log')
    hold on
    h2 = errorbar(v, mean(dsi{i}{2}), std(dsi{i}{2})/sqrt(sample_size(i, 2)), 'r');
    set(get(h2,'Parent'), 'XScale', 'log')
    xlim([v(end) v(1)])
%     if mod(i, 2) == 1
%         ylim([-.5 1])
%     else
%         ylim([.5 1])
%     end
    title(ll{i})
    xlabel('speed')
    ylabel('DSI')
    legend([h1 h2], 'on-off DSGC', 'on-DSGC', 'location', 'southeast')
    set([h1 h2], 'DefaultLineLineWidth', 1.5)

end

% compare across light levels

FigHandle = figure;
set(FigHandle, 'Position', [1 1 1480 1080])

for i = 1:3
    subplot(2, 3, i)
    clear h
    h1 = semilogx(v, dsi{i}{1}', 'b');
    h(1) = h1(1);
    hold on
    h2 = semilogx(v, dsi{i+3}{1}', 'r');
    h(2) = h2(1);
    xlim([v(end) v(1)])
    if i > 1
        ylim([0 1])
    end
    title([SP{i} '  ' CT{1}])    
    xlabel('speed')
    ylabel('DSI')
    legend(h, 'NDF 3', 'NDF 0', 'location', 'southeast')
    
    subplot(2, 3, i+3)
    clear h
    h1 = semilogx(v, dsi{i}{2}', 'b');
    h(1) = h1(1);
    hold on
    h2 = semilogx(v, dsi{i+3}{2}', 'r');
    h(2) = h2(1);
    xlim([v(end) v(1)])
    if i >= 2
        ylim([0 1])
    end
    title([SP{i} '  ' CT{2}])
    xlabel('speed')
    ylabel('DSI')
    legend(h, 'NDF 3', 'NDF 0', 'location', 'southeast')
end

FigHandle = figure;
set(FigHandle, 'Position', [1 1 1480 1080])
set(FigHandle, 'DefaultLineLineWidth', 1.5)

for i = 1:3
    subplot(2, 3, i)
    clear h
    h(1) = errorbar(v, mean(dsi{i}{1}), std(dsi{i}{1})/sqrt(sample_size(i, 1)), 'b');
    set(get(h(1),'Parent'), 'XScale', 'log')
    hold on
    h(2) = errorbar(v, mean(dsi{i+3}{1}), std(dsi{i+3}{1})/sqrt(sample_size(i+3, 1)), 'r');
    set(get(h(2),'Parent'), 'XScale', 'log')
    xlim([v(end) v(1)])
%     if i == 1
%         ylim([-.5 1])
%     else
%         ylim([.5 1])
%     end
    title([SP{i} '  ' CT{1}])    
    xlabel('speed')
    ylabel('DSI')
    legend(h, 'NDF 3', 'NDF 0', 'location', 'southeast')
    
    subplot(2, 3, i+3)
    clear h
    h(1) = errorbar(v, mean(dsi{i}{2}), std(dsi{i}{2})/sqrt(sample_size(i, 2)), 'b');
    set(get(h(1),'Parent'), 'XScale', 'log')
    hold on
    h(2) = errorbar(v, mean(dsi{i+3}{2}), std(dsi{i+3}{2})/sqrt(sample_size(i+3, 2)), 'r');
    set(get(h(2),'Parent'), 'XScale', 'log')
    xlim([v(end) v(1)])
%     if i == 1
%         ylim([-.5 1])
%     else
%         ylim([.5 1])
%     end
    title([SP{i} '  ' CT{2}])
    xlabel('speed')
    ylabel('DSI')
    legend(h, 'NDF 3', 'NDF 0', 'location', 'southeast')
end

%% direction

d = 6;
t = 3;
h = figure;
dirn = 4;
ct = 1;
set(h, 'Position', [1 1 1080 500])
subplot(1, 2, 1)
idx_temp = idx2;
compass(DS{d}.U{t}(idx_temp), DS{d}.V{t}(idx_temp))
color = 'brgk';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DS{d}.U{t}(idx_temp), DS{d}.V{t}(idx_temp), x, y);
    [~, I] = find(IN == 1);
    idx_dir_oo{i} = idx_temp(I);
    id_dir_oo{i} = id(idx_dir_oo{i});
end

d = 6;
t = 3;
h = figure;
dirn = 2;
ct = 1;
set(h, 'Position', [1 1 1080 500])
subplot(1, 2, 1)
idx_temp = idx1;
compass(DS{d}.U{t}(idx_temp), DS{d}.V{t}(idx_temp))
color = 'brgk';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DS{d}.U{t}(idx_temp), DS{d}.V{t}(idx_temp), x, y);
    [~, I] = find(IN == 1);
    idx_dir_on{i} = idx_temp(I);
    id_dir_on{i} = id(idx_dir_on{i});
end

%
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
set(FigHandle, 'Position', [1 1 1300 900])
for i = 1:6
    subplot(2, 3, i)
    for j = 1:dirn
        dsi_temp = DS{i}.dsindex(idx_dir_oo{j}, :);
        i_out = find(sum(dsi_temp, 2) == 0);
        dsi_temp(i_out, :) = [];
        if ~isempty(dsi_temp)
            h = semilogx(v, dsi_temp, color(j));
            H(j) = h(1);
            hold on
        end
    end
    xlim([v(end) v(1)])
%     if mod(i, 2) == 0
%         ylim([.5 1])
%     end
    if i == 6
        legend(H, 'dorsal', 'N/T', 'ventral', 'N/T', 'location', 'southeast')
    end
    title(ll{i})
    xlabel('speed')
    ylabel('DSI')
end



FigHandle = figure;
set(FigHandle, 'Position', [1 1 1300 900])
for i = 1:6
    subplot(2, 3, i)
    for j = 1:dirn
        dsi_temp = DS{i}.dsindex(idx_dir_oo{j}, :);
        i_out = find(sum(dsi_temp, 2) == 0);
        dsi_temp(i_out, :) = [];
        if ~isempty(dsi_temp)
            errorbar(v, mean(dsi_temp), std(dsi_temp)/sqrt(size(dsi_temp, 1)), color(j));
            set(gca, 'XScale', 'log')
            hold on
        end
    end
    xlim([v(end) v(1)])
%     if mod(i, 2) == 0
%         ylim([.5 1])
%     end
    if i == 6
        legend('dorsal', 'N/T', 'ventral', 'N/T', 'location', 'southeast')
    end
    title(ll{i})
    xlabel('speed')
    ylabel('DSI')
end

%% DS tuning curves (drifting grating)
% all ds cells
% t = 2;
dirn = 4;
D = 6;
T = 4;

data_idx = [3 6];
for d = 1:2
    p_direction = DS{D}.angle{T}';
    xx = 0:pi/4:7*pi/4;
    xx = repmat(xx, length(id), 1) - repmat(p_direction, 1, 8);
    xx(xx>pi) = xx(xx>pi)-2*pi;
    xx(xx<-pi) = xx(xx<-pi)+2*pi;


    subplot(1, 2, d)
    for i = 1:dirn
        for cc = 1:length(id)
            if ~dg_idx(cc, d)
            [xsort, seq] = sort(xx(cc, :));
            y_temp = DS{data_idx(d)}.rho{T}(cc, :);
            plot(xsort, y_temp(seq), 'b')
            hold on
            end
        end
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    title(ll{data_idx(d)})
    xlim([-pi pi])
end

%subtypes
clear rho_dg_mean rho_dg_ste dsi_dg_mean dsi_dg_ste
for d = 1:2
    subplot(1, 2, d)
    for i = 1:dirn
        rho_dg{d}{i} = [];
        dsi_dg{d}{i} = [];
        for cc = 1:length(idx_dir_oo{i})
            if ~dg_idx(idx_dir_oo{i}(cc), d) && sum(DS{data_idx(d)}.rho{T}(idx_dir_oo{i}(cc), :))>0
            [xsort, seq] = sort(xx(idx_dir_oo{i}(cc), :));
            y_temp = DS{data_idx(d)}.rho{T}(idx_dir_oo{i}(cc), :);
            plot(xsort, y_temp(seq), color(i))
            ylim([0 1])
%             pause
            hold on
            rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
            dsi_dg{d}{i} = [dsi_dg{d}{i}; DS{data_idx(d)}.dsindex{T}(idx_dir_oo{i}(cc))];
            end
        end
        rho_dg_mean{d}(i, :) = mean(rho_dg{d}{i});
        rho_dg_ste{d}(i, :) = std(rho_dg{d}{i})/sqrt(size(rho_dg{d}{i}, 1));
        dsi_dg_mean{d}(i) = mean(dsi_dg{d}{i});
        dsi_dg_ste{d}(i) = std(dsi_dg{d}{i})/sqrt(length(dsi_dg{d}{i}));
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    title(ll{d})
    xlim([-pi pi])
end
dsi_dg_mean = cell2mat(dsi_dg_mean');
dsi_dg_ste = cell2mat(dsi_dg_ste');
% plot average (cell type)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for d = 1:2
    subplot(1, 2, d)
    for i = 1:dirn
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(i));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
    title(ll{data_idx(d)});
end
% legend(ct)

% plot average (light level)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:dirn
    subplot(2, 2, i)
    for d = 1:5
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(d));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
    title(ct{i})
end
legend(ll)
% DSI
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = ct;
model_series = [dsi_dg_mean(1,1) dsi_dg_mean(2,1) dsi_dg_mean(3,1) dsi_dg_mean(4,1) dsi_dg_mean(5,1); dsi_dg_mean(1,2) dsi_dg_mean(2,2) dsi_dg_mean(3,2) dsi_dg_mean(4,2) dsi_dg_mean(5,2);dsi_dg_mean(1,3) dsi_dg_mean(2,3) dsi_dg_mean(3,3) dsi_dg_mean(4,3) dsi_dg_mean(5,3); dsi_dg_mean(1,4) dsi_dg_mean(2,4) dsi_dg_mean(3,4) dsi_dg_mean(4,4) dsi_dg_mean(5,4)];   
model_error = [dsi_dg_ste(1,1) dsi_dg_ste(2,1) dsi_dg_ste(3,1) dsi_dg_ste(4,1) dsi_dg_ste(5,1); dsi_dg_ste(1,2) dsi_dg_ste(2,2) dsi_dg_ste(3,2) dsi_dg_ste(4,2) dsi_dg_ste(5,2);dsi_dg_ste(1,3) dsi_dg_ste(2,3) dsi_dg_ste(3,3) dsi_dg_ste(4,3) dsi_dg_ste(5,3); dsi_dg_ste(1,4) dsi_dg_ste(2,4) dsi_dg_ste(3,4) dsi_dg_ste(4,4) dsi_dg_ste(5,4)];
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('DSI')
legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

