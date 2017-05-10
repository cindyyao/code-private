cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);


datadg = load_data('/Volumes/lab/analysis/2016-06-17-0/data005-sorted/data005-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-06-17-0/stimuli/s05.mat';
datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [5 6]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

datarun = load_data('/Volumes/lab/analysis/2016-06-17-0/data002-004-map/data002-004-map', opt);
time_points = [1600 3100];
datamb(1:3) = split_datarun(datarun, time_points);
datamb{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-06-17-0/stimuli/s02.mat';
datamb{1} = load_stim_matlab(datamb{1});
datamb{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-06-17-0/stimuli/s03.mat';
datamb{2} = load_stim_matlab(datamb{2});
datamb{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-06-17-0/stimuli/s04.mat';
datamb{3} = load_stim_matlab(datamb{3});

datawn = load_data('/Volumes/lab/analysis/2016-06-17-0/data000-map/data000-map', opt);
datawn = load_sta(datawn);
datawn = get_rfs(datawn, 'all');
%% dg
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
load('DS160617.mat')
n = 3;
duration = [1406.5 1406.5 1406.7]; %sec
bin_size = 0.025; %sec
[raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},ds_id,duration(i),bin_size);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb{i} = get_mb_raster(datamb{i}, ds_id, duration(i));
    for j = 1:length(raster_mb{i})
        if(mb_idx(j))
            raster_mb{i}{j} = [];
        end
    end
    trial_dur{i} = get_mb_trial_dur(datamb{i}, 400, 400, 0.5);
end

ctr_p = 1; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_mb = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}, raster_p_sum_all{i}] = get_pdirection_raster(raster_mb{i}, MB{1}.angle{ctr_p});
    MAG_all_norm_mb{i} = normalize_MAG(MB{i});
    rep = datamb{i}.stimulus.repetitions;
end

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

t = 4;
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
for cc = 1:length(ds_id)
    plot_mb_raster_ctr(MB, raster_mb, trial_dur, cc, ds_id(cc), 'NDF0', 3, 7, 1)
end

%% Direction tuning (moving bar)
% all ds cells

color = 'brgkcmy';
dirn = 4;
D = 1;
T = 1;
BW = 1;
CL = 1;

p_direction = MB{D}.angle{T,BW,CL}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;


%subtypes
clear rho_mb_mean rho_mb_ste dsi_mb_mean dsi_mb_ste rho_mb dsi_mb dsi_mb_mean_all dsi_mb_ste_all dsi_mb_mean_all_on dsi_mb_ste_all_on

for drug = 1:3
    for cl = 1:7
%         subplot(3, 7, (drug-1)*7+cl)
        for i = 1:4
            rho_mb{drug}{i}{cl} = [];
            RHO_mb{drug}{i}{cl} = [];
            dsi_mb{drug}{i}{cl} = [];
            for cc = 1:length(idx_dir{i})
                if ~mb_idx(idx_dir{i}(cc)) && sum(MB{drug}.RHO{T, BW,cl}(idx_dir{i}(cc), :))>0
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
                end
            end
            if ~isempty(rho_mb{drug}{i}{cl})
                rho_mb_mean{cl}{drug}(i, :) = mean(rho_mb{drug}{i}{cl},1);
                rho_mb_ste{cl}{drug}(i, :) = std(rho_mb{drug}{i}{cl},[],1)/sqrt(size(rho_mb{drug}{i}{cl}, 1));
                RHO_mb_mean{cl}{drug}(i, :) = mean(RHO_mb{drug}{i}{cl},1);
                RHO_mb_ste{cl}{drug}(i, :) = std(RHO_mb{drug}{i}{cl},[],1)/sqrt(size(RHO_mb{drug}{i}{cl}, 1));
                dsi_mb_mean{cl}{drug}(i) = mean(dsi_mb{drug}{i}{cl});
                dsi_mb_ste{cl}{drug}(i) = std(dsi_mb{drug}{i}{cl})/sqrt(length(dsi_mb{drug}{i}{cl}));
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
Drug = {'control', 'drug', 'wash'};
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
            title(ct{i})
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
ctr = {'300%', '150%', '80%', '40%', '20%', '10%', '5%'};
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
%% contrast response function (spike count)

% DS cell
ctr_x = [300 150 80 40 20 10 5];
color = 'brgkc';
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'Hex', 'wash'};
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
                    pd_dir_spikes_nor{drug}{dir}(CC,:) = pd_dir_spikes{drug}{dir}(CC,:)/max(pd_dir_spikes{1}{dir}(CC,:));
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
                    pd_dir_on_spikes_nor{drug}{dir}(CC,:) = pd_dir_on_spikes{drug}{dir}(CC,:)/max(pd_dir_on_spikes{1}{dir}(CC,:));
                    CC = CC + 1;
                end
            end
        end
        pd_dir_on_spikes_nor{drug}{dir} = nan2empty(pd_dir_on_spikes_nor{drug}{dir});
    end

end

for drug = 1:3
    figure
    for dir = 1:4
        errorbar(ctr_x, mean(pd_dir_spikes_nor{drug}{dir}), std(pd_dir_spikes_nor{drug}{dir})/sqrt(size(pd_dir_spikes_nor{drug}{dir}, 1)), 'color', color(dir));
        hold on
    end
    legend(dscell_type)
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('normalized response')
    title(condition{drug})
end

for drug = 1:3
    figure
    for dir = 1:2
        errorbar(ctr_x, mean(pd_dir_on_spikes_nor{drug}{dir}), std(pd_dir_on_spikes_nor{drug}{dir})/sqrt(size(pd_dir_on_spikes_nor{drug}{dir}, 1)), 'color', color(dir));
        hold on
    end
    legend(dscell_type)
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('normalized response')
    title(condition{drug})
end
%% CRF (Peak firing rate)
bin_size = 0.25;
xx = bin_size/2:bin_size:2.6-bin_size/2;
ctr_x = [300 150 80 40 20 10 5];
color = 'brgkc';
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'Hex', 'wash'};
for drug = 1:3    
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                a = raster_p_sum{drug}{idx_dir{dir}(cc)};
                hist_temp = cellfun(@(a) hist(a, xx), a, 'UniformOutput', false);
                hist_temp = cell2mat(squeeze(cellfun(@(hist_temp) max(hist_temp), hist_temp, 'UniformOutput', false)))/datamb{drug}.stimulus.repetitions/bin_size;
                if sum(hist_temp) ~= 0
                    pd_dir_spikes{drug}{dir}(CC,:) = hist_temp;
                    pd_dir_spikes_nor{drug}{dir}(CC,:) = pd_dir_spikes{drug}{dir}(CC,:)/max(pd_dir_spikes{1}{dir}(CC,:));
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
                a = raster_p_sum{drug}{idx_dir_on{dir}(cc)};
                hist_temp = cellfun(@(a) hist(a, xx), a, 'UniformOutput', false);
                hist_temp = cell2mat(squeeze(cellfun(@(hist_temp) max(hist_temp), hist_temp, 'UniformOutput', false)))/datamb{drug}.stimulus.repetitions/bin_size;
                if sum(hist_temp) ~= 0
                    pd_dir_on_spikes{drug}{dir}(CC,:) = hist_temp;
                    pd_dir_on_spikes_nor{drug}{dir}(CC,:) = pd_dir_on_spikes{drug}{dir}(CC,:)/max(pd_dir_on_spikes{1}{dir}(CC,:));
                    CC = CC + 1;
                end
            end
        end
        pd_dir_on_spikes_nor{drug}{dir} = nan2empty(pd_dir_on_spikes_nor{drug}{dir});
    end

end

for drug = 1:3
    figure
    for dir = 1:4
        errorbar(ctr_x, mean(pd_dir_spikes_nor{drug}{dir}), std(pd_dir_spikes_nor{drug}{dir})/sqrt(size(pd_dir_spikes_nor{drug}{dir}, 1)), 'color', color(dir));
        hold on
    end
    legend(dscell_type)
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('normalized response')
    title(condition{drug})
end

for drug = 1:3
    figure
    for dir = 1:2
        errorbar(ctr_x, mean(pd_dir_on_spikes_nor{drug}{dir}), std(pd_dir_on_spikes_nor{drug}{dir})/sqrt(size(pd_dir_on_spikes_nor{drug}{dir}, 1)), 'color', color(dir));
        hold on
    end
    legend(dscell_type)
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('normalized response')
    title(condition{drug})
end

figure
set(gcf, 'Position', [1 1 900 800])
for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:3
        errorbar(ctr_x, mean(pd_dir_spikes{drug}{dir}), std(pd_dir_spikes{drug}{dir})/sqrt(size(pd_dir_spikes{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('firing rate')
    title(dscell_type{dir})
    xlim([3 400])
end

figure
set(gcf, 'Position', [1 1 900 800])
for dir = 1:2
    subplot(2,2,dir)
    for drug = 1:3
        errorbar(ctr_x, mean(pd_dir_on_spikes{drug}{dir}), std(pd_dir_on_spikes{drug}{dir})/sqrt(size(pd_dir_on_spikes{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('firing rate')
    title(dscell_type{dir})
    xlim([3 400])
end

%% non-DSRGC
clear nonds_id
nonds_type = {'ON transient', 'ON sustained', 'OFF transient slow', 'OFF transient large', 'OFF transient'};
for i = 1:length(nonds_type)
    nds_temp = get_cell_ids(datawn, nonds_type{i});
    nonds_id{i} = intersect(nds_temp, datarun.cell_ids);
end
% nonds_id{2}([2 4 5 6]) = [];

nonds_id_all = sort(cell2mat(nonds_id));
for i = 1:length(nonds_type)
    [~,nonds_idx{i}] = intersect(nonds_id_all, nonds_id{i});
end

for i = 1:n
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},nonds_id_all,duration(i),bin_size);
    MB_n{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb_nds{i} = get_mb_raster(datamb{i}, nonds_id_all, duration(i));
end

ctr_p = 1; % choose which params to use to calculate prefer direction indices 

for i = 1:n
    [raster_p_sum_nds{i}, p_idx{i}, raster_p_sum_nds_all{i}] = get_pdirection_raster(raster_mb_nds{i}, MB_n{i}.angle{ctr_p});
end

%% spike count
for drug = 1:3
    
    for cc = 1:length(raster_p_sum_nds{1})
        if ~isempty(raster_p_sum_nds{drug}{cc})
            pd_spikes_nds{drug}(cc,:) = cellfun('length', raster_p_sum_nds{drug}{cc})/datamb{drug}.stimulus.repetitions;
        end
    end

    
    for ct = 1:length(nonds_type)
        CC = 1; 
        for cc = 1:length(nonds_id{ct})
            if ~isempty(raster_p_sum_nds{drug}{nonds_idx{ct}(cc)})
                spike_temp = cellfun('length', raster_p_sum_nds{drug}{nonds_idx{ct}(cc)})/datamb{drug}.stimulus.repetitions;
                if sum(spike_temp) ~= 0
                    pd_ct_spikes{drug}{ct}(CC,:) = cellfun('length', raster_p_sum_nds{drug}{nonds_idx{ct}(cc)})/datamb{drug}.stimulus.repetitions;
                    pd_ct_spikes_nor{drug}{ct}(CC,:) = pd_ct_spikes{drug}{ct}(CC,:)/max(pd_ct_spikes{1}{ct}(CC,:));
                    CC = CC + 1;
                end
            end         
        end
        pd_ct_spikes_nor{drug}{ct} = nan2empty(pd_ct_spikes_nor{drug}{ct});
    end
end

drug = 1;

figure
for ct = 1:5
    errorbar(ctr_x, mean(pd_ct_spikes_nor{drug}{ct}), std(pd_ct_spikes_nor{drug}{ct})/sqrt(size(pd_ct_spikes_nor{drug}{ct}, 1)));
    hold on
%     pause
end
set(gca, 'Xscale', 'log')
legend(nonds_type)
xlabel('% contrast')
ylabel('normalized response')
title(condition{drug})

figure
for ct = 1:5
    h1 = errorbar(ctr_x, mean(pd_ct_spikes_nor{drug}{ct}), std(pd_ct_spikes_nor{drug}{ct})/sqrt(size(pd_ct_spikes_nor{drug}{ct}, 1)), 'b');
    hold on
end
for ct = 1:4
    h2 = errorbar(ctr_x, mean(pd_dir_spikes_nor{drug}{ct}), std(pd_dir_spikes_nor{drug}{ct})/sqrt(size(pd_dir_spikes_nor{drug}{ct}, 1)), 'r');
    hold on
end
legend([h1 h2], {'nDS cell', 'DS cell'})
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('normalized response')
title(condition{drug})

%% peak firing rate
for drug = 1:3
    for ct = 1:length(nonds_type)
        CC = 1; 
        for cc = 1:length(nonds_id{ct})
            if ~isempty(raster_p_sum_nds{drug}{nonds_idx{ct}(cc)})
                a = raster_p_sum_nds{drug}{nonds_idx{ct}(cc)};
                hist_temp = cellfun(@(a) hist(a, xx), a, 'UniformOutput', false);
                hist_temp = cell2mat(squeeze(cellfun(@(hist_temp) max(hist_temp), hist_temp, 'UniformOutput', false)))/datamb{drug}.stimulus.repetitions/bin_size;
                if sum(hist_temp) ~= 0
                    pd_ct_spikes{drug}{ct}(CC,:) = hist_temp;
                    pd_ct_spikes_nor{drug}{ct}(CC,:) = pd_ct_spikes{drug}{ct}(CC,:)/max(pd_ct_spikes{1}{ct}(CC,:));
                    CC = CC + 1;
                end
            end         
        end
        pd_ct_spikes_nor{drug}{ct} = nan2empty(pd_ct_spikes_nor{drug}{ct});
    end
end
%% CRF
marker = 'sodx';
for drug = 1:3
    pd_ds_nor{drug} = cell2mat(pd_dir_spikes_nor{drug}');
    pd_ds{drug} = cell2mat(pd_dir_spikes{drug}');
end

drug = 1;
figure
for ct = 1:3
    errorbar(ctr_x, mean(pd_ct_spikes_nor{drug}{ct}), std(pd_ct_spikes_nor{drug}{ct})/sqrt(size(pd_ct_spikes_nor{drug}{ct}, 1)), 'color', color(ct));
    hold on
end
errorbar(ctr_x, mean(pd_ds_nor{drug}), std(pd_ds_nor{drug})/sqrt(size(pd_ds_nor{drug}, 1)), 'color', color(4));
legend('ON type 1', 'ON type 2', 'OFF type 1', 'DS')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('normalized response')
xlim([3 500])

figure
for drug = 1:3
    errorbar(ctr_x, mean(pd_ds{drug}), std(pd_ds{drug})/sqrt(size(pd_ds{drug}, 1)), 'color', 'k', 'Marker', marker(drug), 'MarkerSize', 6);
    hold on
end
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike count')
xlim([3 500])
% 
% figure
for drug = 1:3
    errorbar(ctr_x, mean(pd_ct_spikes{drug}{2}), std(pd_ct_spikes{drug}{2})/sqrt(size(pd_ct_spikes_nor{drug}{2}, 1)), 'color', 'r', 'Marker', marker(drug), 'MarkerSize', 6);
    hold on
end
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike count')
xlim([3 500])
legend('DS control', 'DS hex','DS wash','ON 2 control',  'ON 2 hex',  'ON 2 wash');

%% spike timing precision
drug = 1;
for ctr = 1:7
    for dir = 1:4
        spike_1st_std{ctr}{dir} = [];
        for cc = 1:length(id_dir{dir});
            idx = idx_dir{dir}(cc);
            
            if ~isempty(raster_p_sum_all{drug}{idx})
            
                spike = squeeze(raster_p_sum_all{drug}{idx}(1,1,ctr,:));
                spike_1st = [];
                for trial = 1:length(spike)
                    if ~isempty(spike{trial})
                        spike_1st = [spike_1st spike{trial}(1)];
                    end
                end
                if length(spike_1st)>=3
                    spike_1st_all{dir}{cc}{ctr} = spike_1st;
                    spike_1st_std{ctr}{dir} = [spike_1st_std{ctr}{dir} robust_std(spike_1st)];
                end
            end
        end
    end

    spike_1st_std_mean(:, ctr) = cellfun(@mean, spike_1st_std{ctr})';
    spike_1st_std_ste(:, ctr) = cellfun(@std,spike_1st_std{ctr})'/sqrt(length(spike_1st_std{ctr}{1}));

%     spike_1st_std{dir}(isnan(spike_1st_std{dir})) = 0;
%     spike_1st_std{dir}(sum(spike_1st_std{dir}, 2) == 0, :) = [];
    spike_1st_std{ctr}{1} = [];
end

spike_1st_std_all = cellfun(@cell2mat, spike_1st_std, 'UniformOutput', false);
spike_1st_std_all_mean = cellfun(@mean, spike_1st_std_all);
spike_1st_std_all_ste = cellfun(@std,spike_1st_std_all)./sqrt(cellfun(@length, spike_1st_std_all));
figure
errorbar(ctr_x, spike_1st_std_all_mean, spike_1st_std_all_ste) 
set(gca, 'Xscale', 'log')
xlim([3 400])

marker = 'sodx';
figure
for dir = 1:4
    errorbar(ctr_x, spike_1st_std_mean(dir, :), spike_1st_std_ste(dir, :), 'Marker', marker(dir), 'color', 'k', 'MarkerSize', 8)
    hold on
end

set(gca, 'Xscale', 'log')
legend(dscell_type)
xlabel('contrast %')
ylabel('std(second)')
title('spike timing precision')
% ylim([0 0.5])
%% plot single cell
for ct = 4:4
    for cc = 1:length(id_dir{ct});
        plot_mb_raster_ctr_one(MB(1), raster_mb(1), trial_dur(1), idx_dir{ct}(cc), id_dir{ct}(cc), 'NDF0', 1, 1, 1)
    end
end
%% tuning width
marker = 'sodx';
drug = 1;
for ctr = 1:7
    C = 1;
    for ct = 1:4
        figure
        CC = 1;
        for cc = 1:size(rho_mb{1}{ct}{ctr}, 1)
            if ~(sum(rho_mb{1}{ct}{ctr}(cc, :)) == 0)
                [f{ct}{ctr}{CC}, G{ct}{ctr}{CC}] = fit_gaussian(xsort, rho_mb{1}{ct}{ctr}(cc, :));
                xfit = linspace(xsort(1), xsort(end), 100);
                yfit = f{ct}{ctr}{CC}.a*exp(-((xfit-f{ct}{ctr}{CC}.b)/f{ct}{ctr}{CC}.c).^2)+f{ct}{ctr}{CC}.d;
                subplot(5,7,cc)
                plot(xsort, rho_mb{1}{ct}{ctr}(cc, :), 'b')
                hold on
                plot(xfit, yfit, 'r')
                tuning_width{ct}{ctr}(CC) = f{ct}{ctr}{CC}.c*0.83*180*2/pi;
                tuning_width_all{ctr}(C) = f{ct}{ctr}{CC}.c*0.83*180*2/pi;
                dc{ct}{ctr}(CC) = f{ct}{ctr}{CC}.d;
                CC = CC + 1;
                C = C + 1;
            end
        end
        tuning_width_mean{ct}(ctr) = mean(tuning_width{ct}{ctr});
        tuning_width_ste{ct}(ctr) = std(tuning_width{ct}{ctr})/sqrt(length(tuning_width{ct}{ctr}));
        dc_mean{ct}(ctr) = mean(dc{ct}{ctr});
        dc_ste{ct}(ctr) = std(dc{ct}{ctr})/sqrt(length(dc{ct}{ctr}));
    end
    tuning_width_mean_all(ctr) = mean(tuning_width_all{ctr});
    tuning_width_ste_all(ctr) = std(tuning_width_all{ctr})/sqrt(length(tuning_width_all{ctr}));

end

aa = 7;
figure
for ct = 1:4
    errorbar(ctr_x(1:aa), tuning_width_mean{ct}(1:aa), tuning_width_ste{ct}(1:aa), 'Marker', marker(ct), 'color', 'k', 'MarkerSize', 8);
    hold on
end
set(gca, 'Xscale', 'log')
ylim([0 200])
xlim([3 400])
legend(dscell_type)
xlabel('contrast %')
ylabel('degree')


figure
errorbar(ctr_x, tuning_width_mean_all, tuning_width_ste_all)
set(gca, 'Xscale', 'log')
xlabel('contrast %')
ylabel('degree')
title('tuning width')
ylim([0 200])

%% CRF
for drug = 1:3
    pd_ds_nor{drug} = cell2mat(pd_dir_spikes_nor{drug}');
    pd_ds{drug} = cell2mat(pd_dir_spikes{drug}');
end

drug = 1;
figure
for ct = 1:3
    errorbar(ctr_x, mean(pd_ct_spikes_nor{drug}{ct}), std(pd_ct_spikes_nor{drug}{ct})/sqrt(size(pd_ct_spikes_nor{drug}{ct}, 1)), 'color', color(ct));
    hold on
end
errorbar(ctr_x, mean(pd_ds_nor{drug}), std(pd_ds_nor{drug})/sqrt(size(pd_ds_nor{drug}, 1)), 'color', color(4));
legend('ON type 1', 'ON type 2', 'OFF type 1', 'DS')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('normalized response')
xlim([3 500])

figure
for drug = 1:3
    errorbar(ctr_x, mean(pd_ds{drug}), std(pd_ds{drug})/sqrt(size(pd_ds{drug}, 1)), 'color', 'k', 'Marker', marker(drug), 'MarkerSize', 6);
    hold on
end
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike count')
xlim([3 500])
% 
% figure
for drug = 1:3
    errorbar(ctr_x, mean(pd_ct_spikes{drug}{2}), std(pd_ct_spikes{drug}{2})/sqrt(size(pd_ct_spikes_nor{drug}{2}, 1)), 'color', 'r', 'Marker', marker(drug), 'MarkerSize', 6);
    hold on
end
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike count')
xlim([3 500])
legend('DS control', 'DS hex','DS wash','ON 2 control',  'ON 2 hex',  'ON 2 wash');

%%
for cc = 2:2%length(nonds_id{1})
    id = nonds_id{2}(cc);
    CC = find(nonds_id_all == id);
    plot_mb_raster_ctr_one(MB_n(1), raster_mb_nds(1), trial_dur(1), CC, id, 'NDF0', 1, 1, 0)
end
%% map camera and display coordinates
cd /Volumes/lab/Experiments/Array/Images/2016-06-17-0/
im_s = imread('WN.jpg'); % load stimulus picture taken by camera
im_array = imread('array.jpg'); % load array image taken by camera
stixel_size = 30; % frame shown in WN.jpg
movie_path = '/Volumes/lab/acquisition/movie-xml/BW-30-6-0.48-11111-20x20-60.35.xml';

cd /Users/xyao/matlab/code-private/DS_new/
mov = get_movie(movie_path, 0, 1);
mov_frame = matrix_scaled_up(squeeze(mov(:,:,1)), stixel_size);
clear movingPoints fixedPoints
cpselect(im_s, mov_frame) % select 4 control points


%% register two images
tform = fitgeotrans(movingPoints, fixedPoints, 'projective');
registered = imwarp(im_s, tform,'OutputView',imref2d(size(mov_frame)));
figure 
imshow(registered);
figure
imshowpair(mov_frame,registered,'blend');

% transform array image into display coordinates
registered_array = imwarp(im_array, tform, 'OutputView', imref2d(size(mov_frame)));
figure
imshow(registered_array);

% get array location in display coordinates

%                 EI                               DISPLAY
%
%               195 (1)                         386(5)  264(6)
%                 / \                               ______
%               /     \                            /      \
%   264 (6)    |       |    126 (2)               /        \
%   386 (5)    |       |    4   (3)       455(4)  \        / 195(1)
%               \     /                            \      /
%                 \ /                               ------
%                455 (4)                          4(3)  126(2)
array_location_display = ginput;
load('electrode_position.mat')
% get array location in ei coordinates
elec_corner = [195 126 4 455];
array_location_ei = position(elec_corner,:);
Tform = maketform('projective', array_location_ei, array_location_display);
test = tformfwd(Tform, array_location_ei)-array_location_display % should be equal or close to zeros

%%
id = [2628 4970 3336 4097 4096 5088 3556 4488];
figure
for i = 1:length(id)
    plot_rf_summaries(datawn, id(i), 'clear', false, 'fit_color', color(mod(i, 7)+1))
end
plot(array_location_display(:,1)/10, array_location_display(:,2)/10, 'k')

%%
repeats = 8;
figure
for i = 1:4
    idx = find(ds_id == id(i));
    subplot(3,3,i)
    for ctr = 1:4
        plot_raster(squeeze(raster_p_sum_all{1}{idx}(1,1,ctr+3,:)), 0, trial_dur{1}, 'color', color(i), 'first_trial', repeats*(ctr-1)+1)
        hold on
    end
    plot([0 2], [repeats repeats], 'k')
    plot([0 2], 2*[repeats repeats], 'k')
    plot([0 2], 3*[repeats repeats], 'k')
end

subplot(3,3,5)
for i = 1:length(id)
    plot_rf_summaries(datawn, id(i), 'clear', false, 'fit_color', color(i))
end
plot(array_location_display(:,1)/10, array_location_display(:,2)/10, 'k')

for i = 6:9
    idx = find(nonds_id_all == id(i-1));
    subplot(3,3,i)
    for ctr = 1:4
        plot_raster(squeeze(raster_p_sum_nds_all{1}{idx}(1,1,ctr+3,:)), 0, trial_dur{1}, 'color', color(i-4), 'first_trial', repeats*(ctr-1)+1)
        hold on
    end
    plot([0 2], [repeats repeats], 'k')
    plot([0 2], 2*[repeats repeats], 'k')
    plot([0 2], 3*[repeats repeats], 'k')
end

        
%% spike timing precision
marker = 'osdx';
drug = 1;
for dir = 1:4
    for cc = 1:length(id_dir{dir});
        idx = idx_dir{dir}(cc);
        if ~isempty(raster_p_sum_all{drug}{idx})
            for ctr = 1:6
                spike = squeeze(raster_p_sum_all{drug}{idx}(1,1,ctr,:));
                spike_1st = [];
                for trial = 1:length(spike)
                    if ~isempty(spike{trial})
                        spike_1st = [spike_1st spike{trial}(1)];
                    end
                end
                
                spike_1st_all{dir}{cc}{ctr} = spike_1st;
                spike_1st_std{dir}(cc, ctr) = std(spike_1st);
            end
        end
    end
    spike_1st_std{dir}(isnan(spike_1st_std{dir})) = 0;
    spike_1st_std{dir}(sum(spike_1st_std{dir}, 2) == 0, :) = [];
    spike_1st_std_mean{dir} = mean(spike_1st_std{dir});
    spike_1st_std_ste{dir} = std(spike_1st_std{dir})/sqrt(size(spike_1st_std{dir}, 1));
end


figure
for dir = 1:4
    errorbar(ctr_x, spike_1st_std_mean{dir}, spike_1st_std_ste{dir}, 'Marker', marker(dir), 'color', 'k', 'MarkerSize', 8)
    hold on
end
set(gca, 'Xscale', 'log')
legend(ct)
xlabel('contrast %')
ylabel('second')
title('spike timing precision')


%% spike timing precision (nds cells)
drug = 1;
for ctr = 1:7
    for dir = 1:1
        spike_1st_std{ctr}{dir} = [];
        for cc = 1:length(nonds_id{3});
            idx = nonds_idx{3}(cc);
            
            if ~isempty(raster_p_sum_nds_all{drug}{idx})
            
                spike = squeeze(raster_p_sum_nds_all{drug}{idx}(1,1,ctr,:));
                spike_1st = [];
                for trial = 1:length(spike)
                    if ~isempty(spike{trial})
                        spike_1st = [spike_1st spike{trial}(1)];
                    end
                end
                if length(spike_1st)>=3
                    spike_1st_all{dir}{cc}{ctr} = spike_1st;
                    spike_1st_std{ctr}{dir} = [spike_1st_std{ctr}{dir} robust_std(spike_1st)];
                end
            end
        end
    end

    spike_1st_std_mean(:, ctr) = cellfun(@mean, spike_1st_std{ctr})';
    spike_1st_std_ste(:, ctr) = cellfun(@std,spike_1st_std{ctr})'/sqrt(length(spike_1st_std{ctr}{1}));

%     spike_1st_std{dir}(isnan(spike_1st_std{dir})) = 0;
%     spike_1st_std{dir}(sum(spike_1st_std{dir}, 2) == 0, :) = [];
%     spike_1st_std{ctr}{1} = [];
end

spike_1st_std_all = cellfun(@cell2mat, spike_1st_std, 'UniformOutput', false);
spike_1st_std_all_mean = cellfun(@mean, spike_1st_std_all);
spike_1st_std_all_ste = cellfun(@std,spike_1st_std_all)./sqrt(cellfun(@length, spike_1st_std_all));
figure
errorbar(ctr_x, spike_1st_std_all_mean, spike_1st_std_all_ste) 
set(gca, 'Xscale', 'log')
xlim([3 400])

marker = 'sodx';
figure
for dir = 1:4
    errorbar(ctr_x, spike_1st_std_mean(dir, :), spike_1st_std_ste(dir, :), 'Marker', marker(dir), 'color', 'k', 'MarkerSize', 8)
    hold on
end

set(gca, 'Xscale', 'log')
legend(dscell_type)
xlabel('contrast %')
ylabel('std(second)')
title('spike timing precision')
% ylim([0 0.5])

%%
alpha = [0 45 90 135 180 225 270 315]/180*pi;
for ctr = 1:7
    t_width_all{ctr} = [];
    for dir = 1:4
        t_width{dir, ctr} = [];
        for cc = 1:size(rho_mb{1}{dir}{ctr}, 1)
            if sum(rho_mb{1}{dir}{ctr}(cc,:)) ~= 0
                [~, s0] = circ_std(alpha', rho_mb{1}{dir}{ctr}(cc,:)');
                t_width{dir, ctr} = [t_width{dir, ctr} s0];
                t_width_all{ctr} = [t_width_all{ctr} s0];
            end
        end
    end
end

t_width_mean = cellfun(@mean, t_width);
t_width_ste = cellfun(@std, t_width)./sqrt(cellfun(@length, t_width));


figure
for dir = 1:4
    errorbar(ctr_x, t_width_mean(dir, :), t_width_ste(dir, :))
    hold on
end
set(gca, 'Xscale', 'log')

t_width_mean_all = cellfun(@mean, t_width_all);
t_width_ste_all = cellfun(@std, t_width_all)./sqrt(cellfun(@length, t_width_all));


figure
errorbar(ctr_x, t_width_mean_all, t_width_ste_all)
set(gca, 'Xscale', 'log')
xlabel('contrast %')
ylabel('circular std')
title('tuning width')
ylim([0 2])

%% 

idx = find(ds_id == 7295);
figure

raster_high = squeeze(raster_mb{1}{idx}(1,1,:,1,:));
raster_low = squeeze(raster_mb{1}{idx}(1,1,:,6,:));

for dir = 1:8
    subplot(3, 8, dir+8) 
    if isempty(min(cell2mat(raster_high(dir, :)')))
        start = 0;
    else
        start = max(min(cell2mat(raster_high(dir, :)'))-0.2, 0);
    end

    plot_raster(raster_high(dir, :), start, min(trial_dur{1}, start+1.5))
    axis off
%     set(gca, 'XTickLabel', [], 'YTickLabel', [])
    subplot(3, 8, dir) 
    plot_raster(raster_low(dir, :), start, min(trial_dur{1}, start+1.5))
    axis off
%     set(gca, 'XTickLabel', [], 'YTickLabel', [])
end
subplot(3,8,17)
%%
color = 'brgck';
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
for drug = 1:1
    for dir = 1:4
        figure
        for i = 1:size(pd_dir_spikes_nor{drug}{dir}, 1)
            plot(ctr_x, pd_dir_spikes_nor{drug}{dir}(i, :), 'color', color(dir));
            hold on
        end
        errorbar(ctr_x, mean(pd_dir_spikes_nor{drug}{dir}), std(pd_dir_spikes_nor{drug}{dir})/sqrt(size(pd_dir_spikes_nor{drug}{dir}, 1)), 'k');
        set(gca, 'Xscale', 'log')
        xlabel('% contrast')
        ylabel('normalized response')
        title(dscell_type{dir})
    end
end
