%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);
load('DS170524.mat')
datadgshort = load_data('/Volumes/lab/analysis/2017-05-24-0/data012-sorted/data012-sorted', opt);
datadgshort.names.stimulus_path = '/Volumes/lab/analysis/2017-05-24-0/stimuli/s12.txt';
datadgshort = load_stim(datadgshort, 'user_defined_trigger_interval', 10);

datadg = load_data('/Volumes/lab/analysis/2017-05-24-0/data011-map/data011-map', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2017-05-24-0/stimuli/s11.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

% identify DS cells
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadgshort,datadgshort.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadgshort);
params_idx = [1 2]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadgshort, ds_struct, params_idx);
ds_id = ds_id(dg_classify_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

datarun = load_data('/Volumes/lab/analysis/2017-05-24-0/datadg-map/datadg-map', opt);
time_points = [650 1300 1950];
dataDG(1:4) = split_datarun(datarun, time_points);
dataDG{1} = datadgshort;

dataDG{2}.names.stimulus_path = '/Volumes/lab/analysis/2017-05-24-0/stimuli/s02.txt';
dataDG{2} = load_stim(dataDG{2}, 'user_defined_trigger_interval', 10);
dataDG{3}.names.stimulus_path = '/Volumes/lab/analysis/2017-05-24-0/stimuli/s05.txt';
dataDG{3} = load_stim(dataDG{3}, 'user_defined_trigger_interval', 10);
dataDG{4}.names.stimulus_path = '/Volumes/lab/analysis/2017-05-24-0/stimuli/s08.txt';
dataDG{4} = load_stim(dataDG{4}, 'user_defined_trigger_interval', 10);

datarun = load_data('/Volumes/lab/analysis/2017-05-24-0/dataflash-map/dataflash-map', opt);
time_points = [3379 6424];
dataflash(1:3) = split_datarun(datarun, time_points);

dataflash{1}.DfParams.NDF =   [4,4,4,3,3,3,3,3,2,2,2,1,1,0] ; % on filter turret 
dataflash{1}.DfParams.Ftime = [2,4,8,2,3,4,6,8,2,4,8,2,8,2] ; % ms
dataflash{1}.DfParams.interFlashInt = [3] ; % sec

dataflash{2}.DfParams.NDF =   [4,4,4,3,3,3,3,3,2,2,2,1,1] ; % on filter turret 
dataflash{2}.DfParams.Ftime = [2,4,8,2,3,4,6,8,2,4,8,2,8] ; % ms
dataflash{2}.DfParams.interFlashInt = [3] ; % sec

dataflash{3}.DfParams.NDF =   [4,4,4,3,3,3,3,3,2,2,2,1,1,0] ; % on filter turret 
dataflash{3}.DfParams.Ftime = [2,4,8,2,3,4,6,8,2,4,8,2,8,2] ; % ms
dataflash{3}.DfParams.interFlashInt = [3] ; % sec

load('DS170524.mat')
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

t = 2;
figure
compass(DG{1}.U{t}(idx_sub{1}), DG{1}.V{t}(idx_sub{1}), 'r')
hold on
compass(DG{1}.U{t}(idx_sub{2}), DG{1}.V{t}(idx_sub{2}), 'b')

%%
d = 1;
t = 2;
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

%% DG
n = 4;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, ~, StimComb] = get_spikescellstim(dataDG{i},ds_id,0,1);
    DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,dataDG{i}));
    raster_dg{i} = get_ds_raster(dataDG{i}, ds_id);
    for j = 1:length(raster_dg{i})
        if(dg_idx(j, i))
            raster_dg{i}{j} = [];
        end
    end
end

delta_p = 2; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{1}.angle{delta_p});
    [raster_n_sum{i}, n_idx{i}] = get_ndirection_raster(raster_dg{i}, DG{1}.angle{delta_p});
    MAG_all_norm_dg{i} = normalize_MAG(DG{i});
    rep = dataDG{i}.stimulus.repetitions;
end

ll = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};

%% plot cell summary
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [4 5]);
title = {'NDF 0', 'NDF 4 ctr', 'NDF 4 SR', 'NDF 4 wash'};
for dir = 3:4
    for cc = 1:length(id_dir{dir})
        plot_ds_raster(DG, raster_dg, idx_dir{dir}(cc), id_dir{dir}(cc), title, 2, 2, 1)
    end
end

%% DS tuning curves (drifting grating)
% all ds cells
color = 'brgkc';
% LL = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};
ct = {'superior', 'anterior', 'inferior', 'posterior'};
% t = 2;
dirn = 4;
D = 1;
T = 1;

p_direction = DG{D}.angle{T}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

%subtypes
clear rho_dg_mean rho_dg_ste dsi_dg_mean dsi_dg_ste
for d = 1:4
    subplot(2, 3, d)
    for i = 1:dirn
        rho_dg{d}{i} = [];
        dsi_dg{d}{i} = [];
        for cc = 1:length(idx_dir{i})
            if ~dg_idx(idx_dir{i}(cc), d) && sum(DG{d}.rho{T}(idx_dir{i}(cc), :))>0
            [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
            y_temp = DG{d}.rho{T}(idx_dir{i}(cc), :);
            plot(xsort, y_temp(seq), color(i))
            ylim([0 1])
%             pause
            hold on
            rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
            dsi_dg{d}{i} = [dsi_dg{d}{i}; DG{d}.dsindex{T}(idx_dir{i}(cc))];
            end
        end
        rho_dg_mean{d}(i, :) = mean(rho_dg{d}{i});
        rho_dg_ste{d}(i, :) = std(rho_dg{d}{i})/sqrt(size(rho_dg{d}{i}, 1));
        dsi_dg_mean{d}(i) = mean(dsi_dg{d}{i});
        dsi_dg_ste{d}(i) = std(dsi_dg{d}{i})/sqrt(length(dsi_dg{d}{i}));
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
%     title(ll{d})
    xlim([-pi pi])
end
dsi_dg_mean = cell2mat(dsi_dg_mean');
dsi_dg_ste = cell2mat(dsi_dg_ste');
% plot average (cell type)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for d = 1:4
    subplot(2, 2, d)
    for i = 1:dirn
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(i));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
%     title(ll{d});
end
legend(ct)

% plot average (light level)
color = 'rbgkc';
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:dirn
    subplot(2, 2, i)
    for d = 1:4
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(d));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
    title(ct{i})
end
% legend(LL)
% DSI
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = ct;
model_series = dsi_dg_mean';
model_error = dsi_dg_ste';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('DSI')
% legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

% dsi curve
marker = 'xo*d';
figure
for i = 1:4
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
errorbar(0:3, model_series(i,:), model_error(i,:), 'color', color(i));
hold on
end
xlim([-0.5 4.5])
ylim([0 1.1])
legend(ct)
xlabel('log(light intensity)')
ylabel('DSI')


for dir = 1:4
    figure
    for ll = 1:4
        subplot(2, 3, ll)
        hist(dsi_dg{ll}{dir}, 0:0.05:1)
        xlim([0 1])
%         title(LL{ll})
        xlabel('DSI')
        ylabel('cell#')
    end
end

%%
T = 1;
for drug = 1:2
    temp = cellfun(@length, vertcat(raster_n_sum{drug + 1}{idx_dir{1}}));
    nullspikes{1}(:, drug) = temp(:, T)/32;
    temp = cellfun(@length, vertcat(raster_n_sum{drug + 1}{cell2mat(idx_dir(2:4))}));
    nullspikes{2}(:, drug) = temp(:, T)/32;
end

figure
errorbar(mean(nullspikes{2}(:, 1)), mean(nullspikes{1}(:, 1)), std(nullspikes{1}(:, 1))/sqrt(size(nullspikes{1}(:, 1), 1)), 'bo');
hold on
errorbar(mean(nullspikes{2}(:, 2)), mean(nullspikes{1}(:, 2)), std(nullspikes{1}(:, 2))/sqrt(size(nullspikes{1}(:, 2), 1)), 'ro');
herrorbar(mean(nullspikes{2}(:, 1)), mean(nullspikes{1}(:, 1)), std(nullspikes{2}(:, 1))/sqrt(size(nullspikes{2}(:, 1), 1)), 'bo');
herrorbar(mean(nullspikes{2}(:, 2)), mean(nullspikes{1}(:, 2)), std(nullspikes{2}(:, 2))/sqrt(size(nullspikes{2}(:, 2), 1)), 'ro');
plot([0 8], [0 8], 'k--')
xlabel('others (Hz)')
ylabel('superior (Hz)')
legend('control', 'SR')
title('ND firing rate')

%%
T = 1;
for drug = 1:2
    temp = cellfun(@length, vertcat(raster_p_sum{drug + 1}{idx_dir{1}}));
    preferspikes{1}(:, drug) = temp(:, T)/32;
    temp = cellfun(@length, vertcat(raster_p_sum{drug + 1}{cell2mat(idx_dir(2:4))}));
    preferspikes{2}(:, drug) = temp(:, T)/32;
end

figure
errorbar(mean(preferspikes{2}(:, 1)), mean(preferspikes{1}(:, 1)), std(preferspikes{1}(:, 1))/sqrt(size(preferspikes{1}(:, 1), 1)), 'bo');
hold on
errorbar(mean(preferspikes{2}(:, 2)), mean(preferspikes{1}(:, 2)), std(preferspikes{1}(:, 2))/sqrt(size(preferspikes{1}(:, 2), 1)), 'ro');
herrorbar(mean(preferspikes{2}(:, 1)), mean(preferspikes{1}(:, 1)), std(preferspikes{2}(:, 1))/sqrt(size(preferspikes{2}(:, 1), 1)), 'bo');
herrorbar(mean(preferspikes{2}(:, 2)), mean(preferspikes{1}(:, 2)), std(preferspikes{2}(:, 2))/sqrt(size(preferspikes{2}(:, 2), 1)), 'ro');
plot([0 8], [0 8], 'k--')
xlabel('others (Hz)')
ylabel('superior (Hz)')
legend('control', 'SR')
title('PD firing rate')


%%
ds_id_flash = ds_id(~flash_idx);
for ct = 1:4
    [id_dir_flash{ct}, idx_dir_flash{ct}] = intersect(ds_id_flash, id_dir{ct});
end

%% parameters
interFlashIntVar = 0.005; % (sec) expected precision of trigger intervals

%% find sets of flash stimuli
% trigger = {[1:3:300], [601:3:900], [601:3:900], [1:3:300]};
for ds = 1:3
    clear trigger_set_i
    ts = 1;
    trigger_set_i{1} = [] ;
    for a=1:length(dataflash{ds}.triggers) ; % for each trigger
        trigger_set_i{ts} = [trigger_set_i{ts},a] ; % put it in a set
        if a<length(dataflash{ds}.triggers) ; % if its not the last trigger
            if sum(abs(dataflash{ds}.triggers(a+1)-dataflash{ds}.triggers(a)-dataflash{ds}.DfParams.interFlashInt)<interFlashIntVar)>0 ; % next trigger is has the right interval
                ts = ts ; % keep it in the set
            else
                ts = ts+1 ; % put it in a new set
                trigger_set_i{ts} = [] ;
            end
        end
    end
    
    i = 1;
    while(i<=length(trigger_set_i))
        if length(trigger_set_i{i}) < 60
            trigger_set_i(i) = [];
        else
            i = i+1;
        end
    end
    
    if ds == 1
        trigger_set_i(4:6) = [];
    end
    %%
    % interval = dataflash{ds}.DfParams.interFlashInt;
    interval = 1;
    ds_idx_flash = get_cell_indices(dataflash{ds}, ds_id_flash);
    bin_size = 0.01;
    start = 0;
    XX = start+bin_size/2:bin_size:dataflash{ds}.DfParams.interFlashInt-bin_size/2;
    load('DS161208.mat', 'NdfCalibration')
    % NdfCalibration(2, :) = NdfCalibration(2, :)/10;
    Irel = (dataflash{ds}.DfParams.Ftime/1000).*NdfCalibration(2,dataflash{ds}.DfParams.NDF+1) ;
    if ds == 1
        Irel(4:6) = [];
    end
    [Irel, i] = sort(Irel);

    
    ds_flash_raster = cell(length(ds_id_flash), 1);
    ds_flash_hist_trial = cell(length(ds_id_flash), 1);
    ds_flash_hist_mean = cell(length(ds_id_flash), 1);
    for cc = 1:length(ds_id_flash)
        raster_1cell = cell(length(trigger_set_i), 1);
        raster_1cell_all = cell(length(trigger_set_i), 1);
        hist_1cell = cell(length(trigger_set_i), 1);
        hist_1cell_all = cell(length(trigger_set_i), 1);
        for ts = 1:length(trigger_set_i)
            raster_1cell{ts} = get_raster(dataflash{ds}.spikes{ds_idx_flash(cc)}, ...
                dataflash{ds}.triggers(trigger_set_i{ts}), 'start', start, 'stop', ...
                interval, 'plot', 0);
            for t = 1:length(raster_1cell{ts})
                hist_1cell{ts}{t} = hist(raster_1cell{ts}{t}, XX);
            end
            raster_1cell_all{ts} = sort(cell2mat(raster_1cell{ts}));
            hist_1cell_all{ts} = hist(raster_1cell_all{ts}, XX)/length(trigger_set_i{ts});
        end
        ds_flash_raster{cc} = raster_1cell(i, :);
        ds_flash_hist_trial{cc} = hist_1cell(i, :);
        ds_flash_hist_mean{cc} = hist_1cell_all(i, :);
    end
    
    %% get raster and psth for dark trials
    % use 400-580 second (roughly) as dark trials
    
    trigger = 1:3:300;
    ds_dark_raster = cell(length(ds_id_flash), 1);
    ds_dark_raster_all = cell(length(ds_id_flash), 1);
    ds_dark_hist_trial = cell(length(ds_id_flash), 1);
    ds_dark_hist_mean = cell(length(ds_id_flash), 1);
    for cc = 1:length(ds_id_flash)
        ds_dark_raster{cc} = get_raster(dataflash{ds}.spikes{ds_idx_flash(cc)}, ...
            trigger, 'stop', interval, ...
            'plot', 0);
        for t = 1:length(ds_dark_raster{cc})
            ds_dark_hist_trial{cc}{t} = hist(ds_dark_raster{cc}{t}, XX);
        end
        ds_dark_raster_all{cc} = sort(cell2mat(ds_dark_raster{cc}));
        ds_dark_hist_mean{cc} = hist(ds_dark_raster_all{cc}, XX)/length(trigger);
    end
    DS_dark_raster_all(:, ds) = ds_dark_raster_all;
    
    %% 2 alternative forced choice test
    % use 400-580 second (roughly) as dark trials
    clear Pc Pc_dir Pc_dir_mean Pc_dir_ste
    % Irel = (dataflash{ds}.DfParams.Ftime/1000).*NdfCalibration(2,dataflash{ds}.DfParams.NDF+1) ;
    % [Irel, i] = sort(Irel);
    
    window = 3;
    bin_n = window/bin_size;
    
    Pc = zeros(length(ds_id_flash), length(trigger_set_i));
    for cc = 1:length(ds_id_flash)
        for ts = 1:length(trigger_set_i)
            %         trial_n = length(ds_flash_raster{1}{ts});
            trial_n = 60;
            corr_flash = zeros(trial_n, 1);
            corr_dark = zeros(trial_n, 1);
            temp = ds_dark_hist_trial{cc}(1:trial_n);
            ds_dark_hist_sum = sum(cell2mat(temp'));
            ds_flash_hist_sum = sum(cell2mat(ds_flash_hist_trial{cc}{ts}'));
            for t = 1:trial_n
                template_flash = ds_flash_hist_sum - ds_flash_hist_trial{cc}{ts}{t};
                template_flash = template_flash/norm(template_flash);
                template_flash(isnan(template_flash)) = 0;
                template_dark = ds_dark_hist_sum - ds_dark_hist_trial{cc}{t};
                template_dark = template_dark/norm(template_dark);
                template_dark(isnan(template_dark)) = 0;
                DV = template_flash - template_dark;
                corr_flash(t) = ds_flash_hist_trial{cc}{ts}{t}(1:bin_n) * DV(1:bin_n)';
                corr_dark(t) = ds_dark_hist_trial{cc}{t}(1:bin_n) * DV(1:bin_n)';
            end
            Pc(cc, ts) = (sum(corr_flash > 0) + sum(corr_dark <= 0))/(trial_n*2);
        end
    end
    
    
    
    
    % color = 'rbgk';
    % figure(1)
    % subplot(2,2,ds)
    % for ct = 1:4
    %     plot(log10(Irel), Pc(idx_dir_flash{ct}, :)', 'color', color(ct))
    %     hold on
    % end
    % xlabel('log(R*/rod)')
    % ylabel('probability')
    % title('wash')
    % compare across cell type
    for ct = 1:4
        Pc_dir{ct} = Pc(idx_dir_flash{ct}, :);
        Pc_dir_mean(ct, :) = mean(Pc_dir{ct}, 1);
        Pc_dir_ste(ct, :) = std(Pc_dir{ct}, [], 1)/sqrt(length(idx_dir_flash{ct}));
        %     if ds == 1
        %         Pc_dir_mean(ct, 4:6) = [];
        %         Pc_dir_ste(ct, 4:6) = [];
        %     end
    end
    
    
    figure(2)
    subplot(2,2,ds)
    
    for ct = 1:4
        errorbar(log10(Irel), Pc_dir_mean(ct, :), Pc_dir_ste(ct, :), 'color', color(ct));
        hold on
    end
    legend('superior', 'Anterior', 'inferior', 'posterior')
    xlabel('log(R*/rod)')
    ylabel('probability')
    title('wash')
    
end
%% raster plot

a = ceil((length(trigger_set_i)+1)/2);
for cc =1:length(ds_id_flash);
    ts = 1;
    b = 1;
    trial_n = length(ds_dark_raster{1});
    FigHandle = figure(10);
    set(FigHandle, 'Position', [0, 0, 1920, 1080]);
    while ts <= length(trigger_set_i)+1
        if  ts == a+1
            b = 3;
        end
        subplot(a, 4, b)
        if ts > 1
            trial_n = length(ds_flash_raster{cc}{ts-1});
        end
        for j = 1:trial_n
            if ts == 1
                SpikeTime = ds_dark_raster{cc}{j};
            else
                SpikeTime = ds_flash_raster{cc}{ts-1}{j};
            end
            SpikeTime = SpikeTime';
            X = [SpikeTime; SpikeTime];
            Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
            line(X, Y, 'color', 'b');
            axis([0, 3, 0, trial_n]);
            hold on
        end
        if b == 1
            title(num2str(ds_id_flash(cc)))
        end
        b = b+1;
        subplot(a, 4, b)
        if ts > 1
            bar(XX, ds_flash_hist_mean{cc}{ts-1}, 1)
        else
            bar(XX, ds_dark_hist_mean{cc}, 1)
        end
        xlim([0 3])

        b = b+3;
        ts = ts+1;
    end

    print_close(10, [24 12], [num2str(ds_id_flash(cc)) '_' num2str(ds)])

end

end