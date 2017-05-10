%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load drifting grating data
datadg = load_data('/Volumes/lab/analysis/2017-03-10-0/data002-sorted/data002-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2017-03-10-0/stimuli/s02.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

% load moving bar data
datarun = load_data('/Volumes/lab/analysis/2017-03-10-0/data003-005-map/data003-005-map', opt);
time_points = [1909 3809];
datamb(1:3) = split_datarun(datarun, time_points);

datamb{1}.names.stimulus_path = '/Volumes/lab/analysis/2017-03-10-0/stimuli/s03.txt';
datamb{1} = load_stim(datamb{1}, 'user_defined_trigger_set', [1:2:1120]);
datamb{1}.stimulus.triggers = datamb{1}.stimulus.triggers';
datamb{2}.names.stimulus_path = '/Volumes/lab/analysis/2017-03-10-0/stimuli/s04.txt';
datamb{2} = load_stim(datamb{2}, 'user_defined_trigger_set', [1:2:1120]);
datamb{2}.stimulus.triggers = datamb{2}.stimulus.triggers';
datamb{3}.names.stimulus_path = '/Volumes/lab/analysis/2017-03-10-0/stimuli/s05.txt';
datamb{3} = load_stim(datamb{3}, 'user_defined_trigger_set', [1:2:1120]);
datamb{3}.stimulus.triggers = datamb{3}.stimulus.triggers';
load('DS170310.mat')
%%
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load drifting grating data
datadg = load_data('/Volumes/janacek/Analysis/2017-03-10-0/data002-sorted/data002-sorted', opt);
datadg.names.stimulus_path = '/Volumes/janacek/Analysis/2017-03-10-0/stimuli/s02.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

% load moving bar data
datarun = load_data('/Volumes/janacek/Analysis/2017-03-10-0/data003-005-map/data003-005-map', opt);
time_points = [1909 3809];
datamb(1:3) = split_datarun(datarun, time_points);

datamb{1}.names.stimulus_path = '/Volumes/janacek/Analysis/2017-03-10-0/stimuli/s03.txt';
datamb{1} = load_stim(datamb{1}, 'user_defined_trigger_set', [1:2:1120]);
datamb{1}.stimulus.triggers = datamb{1}.stimulus.triggers';
datamb{2}.names.stimulus_path = '/Volumes/janacek/Analysis/2017-03-10-0/stimuli/s04.txt';
datamb{2} = load_stim(datamb{2}, 'user_defined_trigger_set', [1:2:1120]);
datamb{2}.stimulus.triggers = datamb{2}.stimulus.triggers';
datamb{3}.names.stimulus_path = '/Volumes/janacek/Analysis/2017-03-10-0/stimuli/s05.txt';
datamb{3} = load_stim(datamb{3}, 'user_defined_trigger_set', [1:2:1120]);
datamb{3}.stimulus.triggers = datamb{3}.stimulus.triggers';


% identify DS cells
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [4 3]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

%% dg
% ds_id = datadg.cell_ids;
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);

delta_p = 4; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(DG{i});
rep = datadg.stimulus.repetitions;

%% mb
n = 3;
for i = 1:n
    duration(i) = datamb{i}.triggers(end);
end
bin_size = 0.025; %sec
[raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx, raster_mb_all] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},ds_id,duration(i),bin_size);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb{i} = get_mb_raster(datamb{i}, ds_id, duration(i));
    raster_mb_all{i} = combine_repeats(raster_mb{i});
    for j = 1:length(raster_mb{i})
        if(mb_idx(j))
            raster_mb{i}{j} = [];
            raster_mb_all{i}{j} = [];
        end
    end
    trial_dur{i} = get_mb_trial_dur(datamb{i}, 400, 400, 0.5);
end

ctr_p = 7; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_mb = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}, raster_p_sum_all{i}] = get_pdirection_raster(raster_mb{i}, MB{1}.angle{ctr_p});
    [raster_n_sum{i}, n_idx{i}, raster_n_sum_all{i}] = get_ndirection_raster(raster_mb{i}, MB{1}.angle{ctr_p});
    MAG_all_norm_mb{i} = normalize_MAG(MB{i});
    rep = datamb{i}.stimulus.repetitions;
end
close all

%% classification based on speed tunning
L = 1;
mag_pca = MAG_all_norm_dg{L};
mag_pca = mag_pca';
[id_sub, idx_sub] = deal(cell(2, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(mag_pca);
pc1 = 1; pc2 = 2;
plot(scores(:, pc1), scores(:, pc2), 'o')
hold on
for i = 1:2
    [x, y] = ginput;
    plot(x, y)
    IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
    [~, idx_sub{i}] = find(IN' == 1);
    id_sub{i} = ds_id(idx_sub{i});
end
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
title('NDF 0')

figure
plot(scores(idx_sub{1}, pc1), scores(idx_sub{1}, pc2), 'ro', scores(idx_sub{2}, pc1), scores(idx_sub{2}, pc2), 'bo')
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
title('NDF 0')

v = 4*datadg.stimulus.params.SPATIAL_PERIOD./datadg.stimulus.params.TEMPORAL_PERIOD;
figure
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
xlim([v(end) v(1)])

t = 3;
figure
compass(DG{1}.U{t}(idx_sub{1}), DG{1}.V{t}(idx_sub{1}), 'r')
hold on
compass(DG{1}.U{t}(idx_sub{2}), DG{1}.V{t}(idx_sub{2}), 'b')

%%
d = 1;
t = 5;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(DG{d}.U{t}(idx_sub{2}), DG{d}.V{t}(idx_sub{2}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}(idx_sub{2}), DG{d}.V{t}(idx_sub{2}), x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = idx_sub{2}(I);
    id_dir{i} = ds_id(idx_dir{i});
end

d = 1;
t = 3;
h = figure;
dirn = 3;
set(h, 'Position', [1 1 1080 500])
compass(DG{d}.U{t}(idx_sub{1}), DG{d}.V{t}(idx_sub{1}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}(idx_sub{1}), DG{d}.V{t}(idx_sub{1}), x, y);
    [~, I] = find(IN == 1);
    idx_dir_on{i} = idx_sub{1}(I);
    id_dir_on{i} = ds_id(idx_dir_on{i});
end

%% plot cell summary
for dir = 3:3
    for cc = 1:1%length(id_dir{dir})
        plot_mb_raster_ctr(MB, raster_mb, trial_dur, idx_dir{dir}(cc), id_dir{dir}(cc), '', 3, 7, 1)
    end
end

%% 
for dir = 1:4
    id_dir_mb{dir} = id_dir{dir}(~mb_idx(idx_dir{dir}));
    idx_dir_mb{dir} = idx_dir{dir}(~mb_idx(idx_dir{dir}));
end
for dir = 1:3
    id_dir_on_mb{dir} = id_dir_on{dir}(~mb_idx(idx_dir_on{dir}));
    idx_dir_on_mb{dir} = idx_dir_on{dir}(~mb_idx(idx_dir_on{dir}));
end
%% get spontaneous activity
for drug = 1:3
    for dir = 1:4
        for cc = 1:length(id_dir_mb{dir})
            idx = get_cell_indices(datamb{drug},id_dir_mb{dir}(cc));
            spikes_temp = datamb{drug}.spikes{idx};
            bgnd_firing{dir}(drug, cc) = length(spikes_temp(spikes_temp < 1830 & spikes_temp > 1780))/50;
        end
    end
    for cc = 1:length(id_sub{1})
        if ~mb_idx(idx_sub{1}(cc))
            idx = get_cell_indices(datamb{drug},id_sub{1}(cc));
            spikes_temp = datamb{drug}.spikes{idx};
            bgnd_firing_on(drug, cc) = length(spikes_temp(spikes_temp < 1830 & spikes_temp > 1780))/50;
        end
    end
%     for cc = 1:length(id_ot_mb)
%         idx = get_cell_indices(datamb{drug},id_ot_mb(cc));
%         spikes_temp = datamb{drug}.spikes{idx};
%         bgnd_firing_ot(drug, cc) = length(spikes_temp(spikes_temp < 1850 & spikes_temp > 1800))/50;
%     end

end

%% max window
% on-off DSGC
clear Max_i
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'wash'};
step_size = 30;
for drug = 1:3
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                for ctr = 7:-1:1
                    a = raster_p_sum{drug}{idx_dir{dir}(cc)}{ctr};
                    hist_temp = hist(a, xx);
                    if drug == 1 && ctr == 7
                        [max_p, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
                        Max_i{dir}(CC) = max_i;
                    else
                        max_p = conv(hist_temp, ones(1,step_size), 'valid');
                        max_p = max_p(Max_i{dir}(CC));
                    end
                    response_pmax{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
                end
                CC = CC + 1;
            end
        end
        response_pmax_norm{drug}{dir} = response_pmax{drug}{dir}./repmat(max(response_pmax{1}{dir}, [], 2), 1, size(response_pmax{drug}{dir},2));
%         response_pmax_norm{drug}{dir} = response_pmax{drug}{dir}./repmat(max(response_pmax{drug}{dir}, [], 2), 1, size(response_pmax{drug}{dir},2));
    end
    response_pmax_all{drug} = cell2mat(response_pmax{drug}(1)');
    response_pmax_norm_all{drug} = cell2mat(response_pmax_norm{drug}(1)');
%     response_pmax_all{drug} = cell2mat(response_pmax{drug}');
%     response_pmax_norm_all{drug} = cell2mat(response_pmax_norm{drug}');

end


ctr_x = [10 20 40 80 150 300 400];
color = 'brgkc';
figure
set(gcf, 'Position', [1 1 900 800])

for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:3
        errorbar(ctr_x, mean(response_pmax_norm{drug}{dir}), std(response_pmax_norm{drug}{dir})/sqrt(size(response_pmax_norm{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
%     xlim([3 400])
end

%% contrast response function (spike count)

% DS cell
ctr_x = [10 20 40 80 150 300 400];
color = 'brgkc';
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'wash'};
for drug = 1:3
    for cc = 1:length(raster_p_sum{1})
        if ~isempty(raster_p_sum{drug}{cc})
            pd_spikes{drug}(cc,:) = cellfun('length', raster_p_sum{drug}{cc})/datamb{drug}.stimulus.repetitions;
        end
    end
    
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                spike_temp = cellfun('length', raster_p_sum{drug}{idx_dir{dir}(cc)})/datamb{drug}.stimulus.repetitions;
                if sum(spike_temp) ~= 0
                    pd_dir_spikes{drug}{dir}(CC,:) = cellfun('length', raster_p_sum{drug}{idx_dir{dir}(cc)})/datamb{drug}.stimulus.repetitions;
%                     pd_dir_spikes_nor{drug}{dir}(CC,:) = pd_dir_spikes{drug}{dir}(CC,:)/max(pd_dir_spikes{1}{dir}(CC,:));
                    pd_dir_spikes_nor{drug}{dir}(CC,:) = pd_dir_spikes{drug}{dir}(CC,:)/max(pd_dir_spikes{drug}{dir}(CC,:));
                    CC = CC + 1;
                end
            end
        end
        pd_dir_spikes_nor{drug}{dir} = nan2empty(pd_dir_spikes_nor{drug}{dir});
    end
    for dir = 1:2
        CC = 1;
        for cc = 1:length(idx_dir_on{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir_on{dir}(cc)})
                spike_temp = cellfun('length', raster_p_sum{drug}{idx_dir_on{dir}(cc)})/datamb{drug}.stimulus.repetitions;
                if sum(spike_temp) ~= 0
                    pd_dir_on_spikes{drug}{dir}(CC,:) = cellfun('length', raster_p_sum{drug}{idx_dir_on{dir}(cc)})/datamb{drug}.stimulus.repetitions;
%                     pd_dir_on_spikes_nor{drug}{dir}(CC,:) = pd_dir_on_spikes{drug}{dir}(CC,:)/max(pd_dir_on_spikes{1}{dir}(CC,:));
                    pd_dir_on_spikes_nor{drug}{dir}(CC,:) = pd_dir_on_spikes{drug}{dir}(CC,:)/max(pd_dir_on_spikes{drug}{dir}(CC,:));
                    CC = CC + 1;
                end
            end
        end
        pd_dir_on_spikes_nor{drug}{dir} = nan2empty(pd_dir_on_spikes_nor{drug}{dir});
    end
    response_p_all{drug} = cell2mat(pd_dir_spikes{drug}');
    response_p_norm_all{drug} = cell2mat(pd_dir_spikes{drug}');

end


figure
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:3
        errorbar(ctr_x, mean(pd_dir_spikes{drug}{dir}), std(pd_dir_spikes{drug}{dir})/sqrt(size(pd_dir_spikes{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition)
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike count')
    title(dscell_type{dir})
    xlim([3 500])
end

%% fit ds
for dir = 1:4
    for drug = 1:3
        figure
        set(gcf, 'Position', [1 1 900 800])
        for cc = 1:size(response_pmax{drug}{dir}, 1)

            ydata = response_pmax{drug}{dir}(cc, :);
            xdata = log10(ctr_x);
            [f, G] = fit_nr(xdata, ydata, 'upper', [100, 100, log10(400), 0]);
            fit_all{drug}{dir}{cc} = f;
            G_all{drug}{dir}{cc} = G;

            x = linspace(min(xdata), max(xdata), 100);
            y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;

            subplot(5,6,cc)
            plot(xdata, ydata)
            hold on
            plot(x, y)

            sigma{drug}{dir}(cc) = f.sigma;
            ymax{drug}{dir}(cc) = f.ymax;
            aa{drug}{dir}(cc) = f.a;
            bb{drug}{dir}(cc) = f.b;
        end
    end
end

fitting = struct();
fitting.sigma = sigma{1};
fitting.ymax = ymax{1};
fitting.aa = aa{1};
fitting.bb = bb{1};
%%

for dir = 1:4
    for drug = 1:3
        for cc = 1:size(response_pmax_norm{drug}{dir}, 1)
            x = linspace(sigma{dir}(1,cc)-0.9, sigma{dir}(1, cc)+1.1, 10);
            y = ymax{dir}(drug, cc)*x.^aa{dir}(drug, cc)./(x.^aa{dir}(drug, cc) + sigma{dir}(drug, cc)^aa{dir}(drug, cc))+bb{dir}(drug, cc);
            pn_y{drug}{dir}(cc, :) = y;
%             pause
        end
    end
end
        
figure
set(gcf, 'Position', [1 1 900 800])
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4 
        errorbar(linspace(0.1,2.1,10), mean(pn_y{drug}{dir}), std(pn_y{drug}{dir})/sqrt(size(pn_y{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
%     set(gca, 'Xscale', 'log')
    xlabel('normalized contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
%     xlim([3 400])
end

figure
set(gcf, 'Position', [1 1 900 800])
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:3
        ydata = mean(response_pmax_norm{drug}{dir});
        yste = std(response_pmax_norm{drug}{dir})/sqrt(size(response_pmax_norm{drug}{dir}, 1));
        xdata = log10(ctr_x);
        [f, G] = fit_nr(xdata, ydata, 'upper', [2, 100, log10(400), 0]);
        sigma_avg(dir, drug) = f.sigma;
        x = linspace(0, max(xdata)+0.5, 100);
        y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a)+f.b;
        h(drug) = plot(10.^x,y,'color', color(drug));
        hold on
        errorbar(10.^xdata, ydata, yste, 'color', color(drug), 'marker', 'o', 'linestyle', 'none')
    end
    legend([h(1), h(2), h(3)], condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
end

%% separate ON and OFF responses for CRF
% ctr = 7;
% bin_size = 0.01;
% xx = bin_size/2:bin_size:trial_dur-bin_size/2;
% for drug = 1:1
%     for cc = 1:length(ds_id)
%         if ~mb_idx(cc)
%             h = figure(1);
%             set(h, 'position', [1 1 1000 1000])
%             for i = 7:-1:1
%                 subplot(7,1,i)
%                 temp = hist(raster_p_sum{drug}{cc}{ctr-7+i}, xx);
%                 plot(xx, temp)
%                 hold on
%             end
% %             [x,~] = ginput;
% %             BreakIndices{drug}(cc, :) = round(x'/bin_size);
% %             close(1)
% %         else
% %             BreakIndices{drug}(cc, :) = nan(1, ctr);
% 
% pause
% close(1)
%         
%         end
%     end
% end
% 
% ctr = 7;
% drug = 1;
% for dir = 1:1
%     for cc = 1:length(id_dir{dir})
%         if ~mb_idx(idx_dir{dir}(cc))
%             h = figure(1);
%             set(h, 'position', [1 1 1000 1000])
%             for i = 7:-1:1
%                 subplot(7,1,i)
%                 temp = hist(raster_p_sum{drug}{idx_dir{dir}(cc)}{ctr-7+i}, xx);
%                 plot(xx, temp)
%                 hold on
%             end
%             [x,~] = ginput;
%             BreakIndices(idx_dir{dir}(cc), :) = round(x'/bin_size);
%             close(1)
%         else
%             BreakIndices(idx_dir{dir}(cc), :) = nan(1, ctr);
% % pause
% % close(1)
%         end
%     end
% end
% 
% %% max window
% % on-off DSGC
% clear Max_i
% trial_dur = mean(diff(datamb{1}.stimulus.triggers));
% bin_size = 0.01;
% xx = bin_size/2:bin_size:trial_dur-bin_size/2;
% dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
% condition = {'control', 'AP5', 'wash'};
% step_size = 60;
% for drug = 1:3
%     for dir = 1:4
%         CC = 1;
%         for cc = 1:length(idx_dir{dir})
%             if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
%                 for ctr = 7:-1:1
%                     a = raster_p_sum{drug}{idx_dir{dir}(cc)}{ctr};
%                     hist_temp = hist(a, xx);
%                     hist_temp(BreakIndices(idx_dir{dir}(cc), ctr):end) = zeros(1, length(hist_temp)-BreakIndices(idx_dir{dir}(cc), ctr)+1);
% %                     if drug == 1 && ctr == 7
%                         [max_p, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
% %                         Max_i{dir}(CC) = max_i;
% %                     else
% %                         max_p = conv(hist_temp, ones(1,step_size), 'valid');
% %                         max_p = max_p(Max_i{dir}(CC));
% %                     end
%                     response_pmax{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
%                 end
%                 CC = CC + 1;
%             end
%         end
%         response_pmax_norm{drug}{dir} = response_pmax{drug}{dir}./repmat(max(response_pmax{1}{dir}, [], 2), 1, size(response_pmax{drug}{dir},2));
% %         response_pmax_norm{drug}{dir} = response_pmax{drug}{dir}./repmat(max(response_pmax{drug}{dir}, [], 2), 1, size(response_pmax{drug}{dir},2));
%     end
% end
% 
% 
% ctr_x = [10 20 40 80 150 300 400];
% color = 'brgkc';
% figure
% set(gcf, 'Position', [1 1 900 800])
% 
% for dir = 1:4
%     subplot(2,2,dir)
%     for drug = 1:3
%         errorbar(ctr_x, mean(response_pmax_norm{drug}{dir}), std(response_pmax_norm{drug}{dir})/sqrt(size(response_pmax_norm{drug}{dir}, 1)), 'color', color(drug));
%         hold on
%     end
%     legend(condition, 'location', 'northwest')
%     set(gca, 'Xscale', 'log')
%     xlabel('% contrast')
%     ylabel('spike rate')
%     title(dscell_type{dir})
% %     xlim([3 400])
% end
% 
%%
% B
Drug = {'control', 'AP5', 'wash'};
color = {[0,0,0]/255,[255,20,60]/255,[0,191,255]/255,  [0.5 0.5 0.5]};

h = figure;
set(h, 'Position', [1 1 1520,1080])
for cl = 1:7
    subplot(1, 7, cl)
    for drug = 1:3
        errorbar(xsort/pi*180, rho_mb_all_mean{cl}{drug}, rho_mb_all_ste{cl}{drug}, 'color', color{drug});
        hold on
    end
    xlabel('degrees')

    if cl == 1
        ylabel(ct{i})
    end

end
legend(Drug)

% C
ctr_i = 4;
x = 0:45:315;
temp = RHO_mb_all{1}{ctr_i} - RHO_mb_all{2}{ctr_i};
temp_mean = mean(temp, 1);
temp_ste = std(temp, [], 1)/sqrt(size(temp, 1));
figure
errorbar(x,temp_mean, temp_ste)
xlabel('degree')
ylabel('absolute suppression (spike #)')


temp = (RHO_mb_all{1}{ctr_i} - RHO_mb_all{2}{ctr_i})./RHO_mb_all{1}{ctr_i};
temp(isinf(temp)) = 0;
temp(isnan(temp)) = 0;
temp(temp < -1) = 0;
temp_mean = nanmean(temp, 1);
temp_ste = nanstd(temp, [], 1)/sqrt(size(temp, 1));
figure
errorbar(x, temp_mean, temp_ste)
ylim([0 1])
xlabel('degree')
ylabel('percentage suppression')

figure
plot(repmat(x', 1, size(temp, 1)), temp', 'color', [.5 .5 .5]) 
xlabel('degree')
ylabel('absolute suppression (spike #)')


% G
ctr_x = [10 20 40 80 150 300 400];
temp = response_pmax_all{1} - response_pmax_all{2};
temp_mean = mean(temp, 1);
temp_ste = std(temp, [], 1)/sqrt(size(temp, 1));
figure
errorbar(ctr_x, temp_mean, temp_ste)
set(gca, 'xscal', 'log')
xlabel('% contrast')
ylabel('absolute suppression (spike #)')

figure
plot(repmat(ctr_x', 1, size(temp, 1)), temp', 'color', [.5 .5 .5]) 
set(gca, 'xscal', 'log')
xlabel('% contrast')
ylabel('absolute suppression (spike #)')


temp = (response_pmax_all{1} - response_pmax_all{2})./response_pmax_all{1};
temp(isinf(temp)) = nan;
temp(temp > 1) = 1;

temp_mean = nanmean(temp, 1);
temp_ste = nanstd(temp, [], 1)/sqrt(size(temp, 1));
figure
errorbar(ctr_x, temp_mean, temp_ste)
ylim([0 1])
set(gca, 'xscal', 'log')
xlabel('% contrast')
ylabel('percentage suppression')


figure
plot(repmat(ctr_x', 1, size(temp, 1)), temp', 'color', [.5 .5 .5]) 
set(gca, 'xscal', 'log')
xlabel('% contrast')
ylabel('percentage suppression (spike #)')
ylim([0 1])

%% direction tuning curves
% all ds cells
ctr = {'10%', '20%', '40%', '80%', '150%', '300%', '400%'};
color = 'brgkcmy';
dirn = 4;
D = 1;
T = 1;
BW = 1;
CL = 7;

p_direction = MB{D}.angle{T,BW,CL}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;


%subtypes
clear rho_mb_all RHO_mb_all rho_mb_mean rho_mb_ste dsi_mb_mean dsi_mb_ste rho_mb dsi_mb dsi_mb_mean_all dsi_mb_ste_all dsi_mb_mean_all_on dsi_mb_ste_all_on

for drug = 1:3
    for cl = 1:7
%         subplot(3, 7, (drug-1)*7+cl)
        rho_mb_all{drug}{cl} = [];
        RHO_mb_all{drug}{cl} = [];
        for i = 1:4
            rho_mb{drug}{i}{cl} = [];
            RHO_mb{drug}{i}{cl} = [];
            dsi_mb{drug}{i}{cl} = [];
            for cc = 1:length(idx_dir{i})
                if ~mb_idx(idx_dir{i}(cc))% && sum(MB{drug}.RHO{T, BW,cl}(idx_dir{i}(cc), :))>0
                [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
                Y_temp = MB{drug}.RHO{T,BW,cl}(idx_dir{i}(cc), :);
                y_temp = MB{drug}.rho{T,BW,cl}(idx_dir{i}(cc), :);
%                 y_temp = MB{drug}.RHO{T,BW,cl}(idx_dir{i}(cc), :)/bin_size;
%                 plot(xsort, y_temp(seq), color(i))
%                     ylim([0 1])
    %             pause
%                 hold on
                rho_mb{drug}{i}{cl} = [rho_mb{drug}{i}{cl}; y_temp(seq)];
                RHO_mb{drug}{i}{cl} = [RHO_mb{drug}{i}{cl}; Y_temp(seq)];
                dsi_mb{drug}{i}{cl} = [dsi_mb{drug}{i}{cl}; MB{drug}.dsindex{T,BW,cl}(idx_dir{i}(cc))];
                rho_mb_all{drug}{cl} = [rho_mb_all{drug}{cl}; y_temp(seq)];
                RHO_mb_all{drug}{cl} = [RHO_mb_all{drug}{cl}; Y_temp(seq)];
                end
            end
            if ~isempty(rho_mb{drug}{i}{cl})
                rho_mb_mean{cl}{drug}(i, :) = mean(rho_mb{drug}{i}{cl}, 1);
                rho_mb_ste{cl}{drug}(i, :) = std(rho_mb{drug}{i}{cl}, [], 1)/sqrt(size(rho_mb{drug}{i}{cl}, 1));
                RHO_mb_mean{cl}{drug}(i, :) = mean(RHO_mb{drug}{i}{cl}, 1);
                RHO_mb_ste{cl}{drug}(i, :) = std(RHO_mb{drug}{i}{cl}, [], 1)/sqrt(size(RHO_mb{drug}{i}{cl}, 1));
                dsi_mb_mean{cl}{drug}(i) = mean(dsi_mb{drug}{i}{cl});
                dsi_mb_ste{cl}{drug}(i) = std(dsi_mb{drug}{i}{cl})/sqrt(length(dsi_mb{drug}{i}{cl}));
                rho_mb_all_mean{cl}{drug} = mean(rho_mb_all{drug}{cl});
                rho_mb_all_ste{cl}{drug} = std(rho_mb_all{drug}{cl}, [], 1)/sqrt(size(rho_mb_all{drug}{cl}, 1));
                RHO_mb_all_mean{cl}{drug} = mean(RHO_mb_all{drug}{cl});
                RHO_mb_all_ste{cl}{drug} = std(RHO_mb_all{drug}{cl}, [], 1)/sqrt(size(RHO_mb_all{drug}{cl}, 1));
            else
                rho_mb_mean{cl}{drug}(i, :) = nan;
                rho_mb_ste{cl}{drug}(i, :) = nan;
                RHO_mb_mean{cl}{drug}(i, :) = nan;
                RHO_mb_ste{cl}{drug}(i, :) = nan;
                dsi_mb_mean{cl}{drug}(i) = nan;
                dsi_mb_ste{cl}{drug}(i) = nan;
                rho_mb_all_mean{cl}{drug} = nan;
                rho_mb_all_ste{cl}{drug} = nan;
                RHO_mb_all_mean{cl}{drug} = nan;
                RHO_mb_all_ste{cl}{drug} = nan;
            end
        end
%         xlabel('direction (rad)')
%         ylabel('spike number')
%             title(ll{d})
%         xlim([-pi pi])
        dsi_mb_mean_all{cl} = cell2mat(dsi_mb_mean{cl}');
        dsi_mb_ste_all{cl} = cell2mat(dsi_mb_ste{cl}'); 

        for i = 1:2
            rho_mb_on{drug}{i}{cl} = [];
            RHO_mb_on{drug}{i}{cl} = [];
            dsi_mb_on{drug}{i}{cl} = [];
            for cc = 1:length(idx_dir_on{i})
                if ~mb_idx(idx_dir_on{i}(cc))% && sum(MB_NDF{drug, d}.RHO{T, BW,cl}(idx_dir{i}(cc), :))>0
                [xsort, seq] = sort(xx(idx_dir_on{i}(cc), :));
                Y_temp = MB{drug}.RHO{T,BW,cl}(idx_dir_on{i}(cc), :);
                y_temp = MB{drug}.rho{T,BW,cl}(idx_dir{i}(cc), :);
%                 y_temp = MB{drug}.RHO{T,BW,cl}(idx_dir{i}(cc), :)/bin_size;
%                 plot(xsort, y_temp(seq), color(i))
%                     ylim([0 1])
    %             pause
%                 hold on
                rho_mb_on{drug}{i}{cl} = [rho_mb_on{drug}{i}{cl}; y_temp(seq)];
                RHO_mb_on{drug}{i}{cl} = [RHO_mb_on{drug}{i}{cl}; Y_temp(seq)];
                dsi_mb_on{drug}{i}{cl} = [dsi_mb_on{drug}{i}{cl}; MB{drug}.dsindex{T,BW,cl}(idx_dir_on{i}(cc))];
                end
            end
            if ~isempty(rho_mb_on{drug}{i}{cl})
                rho_mb_mean_on{cl}{drug}(i, :) = mean(rho_mb_on{drug}{i}{cl},1);
                rho_mb_ste_on{cl}{drug}(i, :) = std(rho_mb_on{drug}{i}{cl},[],1)/sqrt(size(rho_mb_on{drug}{i}{cl}, 1));
                RHO_mb_mean_on{cl}{drug}(i, :) = mean(RHO_mb_on{drug}{i}{cl},1);
                RHO_mb_ste_on{cl}{drug}(i, :) = std(RHO_mb_on{drug}{i}{cl},[],1)/sqrt(size(RHO_mb_on{drug}{i}{cl}, 1));
                dsi_mb_mean_on{cl}{drug}(i) = mean(dsi_mb_on{drug}{i}{cl});
                dsi_mb_ste_on{cl}{drug}(i) = std(dsi_mb_on{drug}{i}{cl})/sqrt(length(dsi_mb_on{drug}{i}{cl}));
            end
        end
%         xlabel('direction (rad)')
%         ylabel('spike number')
%             title(ll{d})
%         xlim([-pi pi])
        dsi_mb_mean_all_on{cl} = cell2mat(dsi_mb_mean_on{cl}');
        dsi_mb_ste_all_on{cl} = cell2mat(dsi_mb_ste_on{cl}');            


    end
end
dsi_mb_mean_all = reshape(cell2mat(dsi_mb_mean_all), size(dsi_mb_mean_all{1}, 1), size(dsi_mb_mean_all{1}, 2), length(dsi_mb_mean_all));
dsi_mb_ste_all = reshape(cell2mat(dsi_mb_ste_all), size(dsi_mb_ste_all{1}, 1), size(dsi_mb_ste_all{1}, 2), length(dsi_mb_ste_all));
dsi_mb_mean_all_on = reshape(cell2mat(dsi_mb_mean_all_on), size(dsi_mb_mean_all_on{1}, 1), size(dsi_mb_mean_all_on{1}, 2), length(dsi_mb_mean_all_on));
dsi_mb_ste_all_on = reshape(cell2mat(dsi_mb_ste_all_on), size(dsi_mb_ste_all_on{1}, 1), size(dsi_mb_ste_all_on{1}, 2), length(dsi_mb_ste_all_on));

% plot average (cell type)
ct = {'superior', 'anterior', 'inferior', 'posterior'};
h = figure;
set(h, 'Position', [1 1 1520,1080])
for drug = 1:3
    for cl = 1:7
        subplot(3, 7, (drug-1)*7+cl)
        for i = 1:dirn
            if i<=size(rho_mb_mean{cl}{drug},1)
                errorbar(xsort/pi*180, rho_mb_mean{cl}{drug}(i, :), rho_mb_ste{cl}{drug}(i, :), color(i));
                hold on
            end
        end
        xlabel('degrees')
        ylabel('spike number')
%                 title(ll{d});

    end
end
legend(ct)

h = figure;
set(h, 'Position', [1 1 1520,1080])
for drug = 1:3
    for cl = 1:7
        subplot(3, 7, (drug-1)*7+cl)
        for i = 1:2
            if i<=size(rho_mb_mean_on{cl}{drug},1)
                errorbar(xsort/pi*180, rho_mb_mean_on{cl}{drug}(i, :), rho_mb_ste_on{cl}{drug}(i, :), color(i));
                hold on
            end
        end
        xlabel('degrees')
        ylabel('spike number')
%                 title(ll{d});

    end
end
legend(ct)

% plot average (cell type)
% Drug = {'control', 'AP5', 'AP5+HEX', 'wash'};
Drug = {'control', 'AP5', 'wash'};

h = figure;
set(h, 'Position', [1 1 1520,1080])
for i = 1:4
    for cl = 1:7
        subplot(4, 7, (i-1)*7+cl)
        for drug = 1:3
            if i<=size(rho_mb_mean{cl}{drug},1)
                errorbar(xsort/pi*180, rho_mb_mean{cl}{drug}(i, :), rho_mb_ste{cl}{drug}(i, :), color(drug));
                hold on
            end
        end
        xlabel('degrees')

        if cl == 1
            ylabel(ct{i})
        end
        if i == 1
            title(ctr{cl})
        end
    end
    if i == 1
        legend(Drug)
    end
end

h = figure;
set(h, 'Position', [1 1 1520,1080])
for i = 1:2
    for cl = 1:7
        subplot(4, 7, (i-1)*7+cl)
        for drug = 1:3
            if i<=size(rho_mb_mean_on{cl}{drug},1)
                errorbar(xsort/pi*180, rho_mb_mean_on{cl}{drug}(i, :), rho_mb_ste_on{cl}{drug}(i, :), color(drug));
                hold on
            end
        end
        xlabel('degrees')

        if cl == 1
            title(ct{i})
        end
    end
    if i == 1
        legend(Drug)
    end
end

% control
% ctr = {'5%', '10%', '20%', '40%', '80%', '150%', '300%'};
ct = {'superior', 'anterior', 'inferior', 'posterior'};
h = figure;
set(h, 'Position', [1 1 1520,1080])
drug = 1;
for i = 1:4
    subplot(2, 2, i)
    for cl = 1:7
        if i<=size(rho_mb_mean{cl}{drug},1)
            errorbar(xsort/pi*180, rho_mb_mean{cl}{drug}(i, :)/max(rho_mb_mean{cl}{drug}(i, :)), rho_mb_ste{cl}{drug}(i, :), color(cl));
            hold on
        end
    end
    xlabel('degrees')
    ylabel('spike number')
    title(ct{i})
%                 title(ll{d});

end
legend(ctr)
% subplot(3,2,5)

h = figure;
set(h, 'Position', [1 1 1520,1080])
drug = 1;
for i = 1:2
    subplot(2, 2, i)
    for cl = 1:7
        if i<=size(rho_mb_mean_on{cl}{drug},1)
            errorbar(xsort/pi*180, rho_mb_mean_on{cl}{drug}(i, :), rho_mb_ste_on{cl}{drug}(i, :), color(cl));
            hold on
        end
    end
    xlabel('degrees')
    ylabel('spike number')
    title(ct{i})
%                 title(ll{d});

end
legend(ctr)
subplot(2,2,3)

% individual cell
figure
drug = 1;
for i = 7:7%size(RHO_mb{drug}{2}{1}, 1)
%     subplot(5,6,i)
    for cl = 1:7
        plot(xsort/pi*180, RHO_mb{drug}{2}{cl}(i,:))
        hold on
    end
end
xlim([-200 200])
xlabel('degrees')
ylabel('spike number')

% DSI against contrast
drug = 1;
ctr_x = [300 150 80 40 20 10 5];
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
figure;
for i = 1:4
    errorbar(ctr_x, dsi_mb_mean_all(drug, i, :), dsi_mb_ste_all(drug, i, :))
    hold on
end
set(gca, 'Xscale', 'log')
legend(dscell_type, 'location', 'NorthEastOutside')
xlabel('contrast %')
ylabel('DSI')

%% direction tuning (sliding window)

clear Max_i
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'wash'};
step_size = 60;
for drug = 1:3
    for dir = 1:4
        for theta = 1:8
            CC = 1;
            for cc = 1:length(idx_dir{dir})
                if ~isempty(raster_mb_all{drug}{idx_dir{dir}(cc)})
                    for ctr = 7:-1:1
                        a = raster_mb_all{drug}{idx_dir{dir}(cc)}{1,1,theta,ctr};
                        hist_temp = hist(a, xx);
                        if drug == 1 && ctr == 7
                            [max_fr, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
                            Max_i{dir}(theta, CC) = max_i;
                        else
                            max_fr = conv(hist_temp, ones(1,step_size), 'valid');
                            max_fr = max_fr(Max_i{dir}(theta, CC));
                        end
                        response_max{drug}{dir}(CC, theta, ctr) = max_fr/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
                    end
                    CC = CC + 1;
                end
            end
            response_max_norm{drug}{dir}(:, theta, :) = squeeze(response_max{drug}{dir}(:, theta, :))./repmat(squeeze(max(response_max{1}{dir}(:, theta, :), [], 1)), 1, size(response_max{drug}{dir},1))';
        end
    end
end


for dir = 1:4
    for cc = 1:length(idx_dir_mb{dir})
        pindex{dir}(cc) = get_pindex(response_max{1}{dir}(cc, :, 7));
        for drug = 1:3
            for ctr = 1:7
                response_max{drug}{dir}(cc, :, ctr) = circshift(response_max{drug}{dir}(cc, :, ctr), [0, 5-pindex{dir}(cc)]);
                response_max_norm{drug}{dir}(cc, :, ctr) = response_max{drug}{dir}(cc, :, ctr)/max(response_max{drug}{dir}(cc, :, ctr));
            end
        end
    end
end

for drug = 1:3
    for dir = 1:4
        for ctr = 1:7
            response_max_mean{drug}{dir}{ctr} = mean(response_max{drug}{dir}(:, :, ctr));
            response_max_ste{drug}{dir}{ctr} = std(response_max{drug}{dir}(:, :, ctr))./sqrt(size(response_max{drug}{dir}, 1));
            response_max_norm_mean{drug}{dir}{ctr} = nanmean(response_max_norm{drug}{dir}(:, :, ctr));
            response_max_norm_ste{drug}{dir}{ctr} = nanstd(response_max_norm{drug}{dir}(:, :, ctr))./sqrt(size(response_max_norm{drug}{dir}, 1));
        end
    end
end

theta = linspace(-pi, pi, 9);
theta = theta(1:8);
color = 'brgkcmy';
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
contrast = {'10%', '20%', '40%', '80%', '150%', '300%', '400%'};

figure
for dir = 1:4
    for ctr = 1:7
        subplot(4, 7, 7 * (dir - 1) + ctr);
        for drug = 1:3
            errorbar(theta, response_max_norm_mean{drug}{dir}{ctr}, response_max_norm_ste{drug}{dir}{ctr}, 'color', color(drug));
            hold on
        end
        xlim([-3.5 3])
        if dir == 1
            title(contrast{ctr})
            if ctr == 7
                legend('control', 'AP5', 'wash')
            end
        end
        if (ctr == 1)
            ylabel(dscell_type{dir});
        end
        if (dir == 4)
            xlabel('direction')
        end
    end
end


clear RHO_thresh RHO_thresh_mean RHO_thresh_ste
for drug = 1:3
    for dir = 1:4
        CC = 1;
        RHO_thresh{drug}{dir} = [];
        for cc = 1:length(id_dir_mb{dir})
            if ~mb_idx(idx_dir_mb{dir}(cc))
                for ctr = 7:-1:1
                    spikeNum = length(raster_p_sum{drug}{idx_dir_mb{dir}(cc)}{ctr});
                    if spikeNum > 0 && spikeNum < 50
                        RHO_thresh{drug}{dir}(CC,:) = RHO_mb{drug}{dir}{ctr}(cc,:);
                        CC = CC + 1;
                        break;
                    end
                end
            end
        end
        if (~isempty(RHO_thresh{drug}{dir}))
            RHO_thresh_mean{drug}(dir, :) = mean(RHO_thresh{drug}{dir}, 1);
            RHO_thresh_ste{drug}(dir, :) = std(RHO_thresh{drug}{dir}, [], 1)/sqrt(size(RHO_thresh{drug}{dir}, 1));
        end
    end
end

figure
for dir = 1:4
    subplot(2,2,dir)
    errorbar(theta, RHO_thresh_mean{1}(dir, :), RHO_thresh_ste{1}(dir, :), 'b');
    hold on
    errorbar(theta, RHO_thresh_mean{2}(dir, :), RHO_thresh_ste{2}(dir, :), 'r');
    xlabel('direction')
    ylabel('spike #')
    title(dscell_type{dir})
    legend('control', 'AP5')
end
    
figure
RHO_thresh_temp = cell2mat(RHO_thresh{1}(2:4)');
for cc = 1:size(RHO_thresh_temp, 1)
    subplot(6,6,cc)
    plot(theta, RHO_thresh_temp(cc, :))
end

figure
RHO_thresh_temp = cell2mat(RHO_thresh{2}(2:4)');
for cc = 1:size(RHO_thresh_temp, 1)
    subplot(6,6,cc)
    plot(theta, RHO_thresh_temp(cc, :))
end

%% 
for drug = 1:2
    CC = 1;
    for dir = 2:4
    for cc = 1:length(idx_dir{dir})
        if ~mb_idx(idx_dir{dir}(cc))
            for ctr = 1:7
                dsiVector{drug}(CC, ctr) = MB{drug}.dsindex{ctr}(idx_dir{dir}(cc));
                spikeCount{drug}(CC, ctr) = length(raster_p_sum{drug}{idx_dir{dir}(cc)}{ctr})/rep;
            end
            CC = CC + 1;
        end
    end
    end
end


figure
[spikeCountAvg{1}, dsiVectorAvg{1}, disError{1}] = curve_from_binning(spikeCount{1}(:), dsiVector{1}(:), 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:2:40);
h1 = errorbar(spikeCountAvg{1}, dsiVectorAvg{1}, disError{1}, 'b');
hold on
scatter(spikeCount{1}(:), dsiVector{1}(:), [], [0.7 0.7 1]);
[spikeCountAvg{2}, dsiVectorAvg{2}, disError{2}] = curve_from_binning(spikeCount{2}(:), dsiVector{2}(:), 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:2:40);
h2 = errorbar(spikeCountAvg{2}, dsiVectorAvg{2}, disError{2}, 'r');
scatter(spikeCount{2}(:), dsiVector{2}(:), [], [1 0.7 0.7]);


xlabel('spike# at PD')
ylabel('DSI')
title('2017-03-10-0')
legend([h1, h2],'control', 'AP5', 'location', 'southeast')
ylim([0 1])

