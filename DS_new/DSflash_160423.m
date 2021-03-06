cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);


datadg = load_data('/Volumes/lab/analysis/2016-04-23-0/data002-DSsorted/data002-DSsorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-04-23-0/stimuli/s02.mat';
datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);

dataflash = load_data('/Volumes/lab/analysis/2016-04-23-0/data000-map/data000-map', opt);
dataflash.DfParams.NDF =   [5,5,5,4,4,4,4,4,3,3,3,2,2,1,0] ; % on filter turret 
dataflash.DfParams.Ftime = [2,4,8,2,3,4,6,8,2,4,8,2,4,2,2] ; % ms
dataflash.DfParams.interFlashInt = [3] ; % sec


[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [4 5]; % which parameters to use for classification

[ds_id, nonds_id, id_init] = classify_ds(datadg, ds_struct, params_idx);

load('DS160423.mat')
ds_id_flash = ds_id(flash_idx);
ds_idx_flash = get_cell_indices(dataflash, ds_id_flash);

%% 
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg,ds_id,0,1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);

delta_p = 4; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(DG{i});
rep = datadg.stimulus.repetitions;

%%
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [4 5]);
id = 7098;
cc = get_cell_indices(datadg, id);
plot_ds_raster(DG_cut, raster_dg_cut, cc, id, '', 1, 1, 0)

%%
figure
v = datadg.stimulus.params.SPATIAL_PERIOD./datadg.stimulus.params.TEMPORAL_PERIOD*4;
semilogx(v, exciseColumn(MAG_all_norm_dg{1}), 'b')
xlabel('micron/second')
ylabel('Response')
% title(ll{2})
xlim([v(end) v(1)])



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
ylabel('3rd Principal Component')
title('NDF 0')


figure
plot(scores(idx_sub{1}, pc1), scores(idx_sub{1}, pc2), 'ro', scores(idx_sub{2}, pc1), scores(idx_sub{2}, pc2), 'bo')
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
title('NDF 0')

v = 4*datadg.stimulus.params.SPATIAL_PERIOD./datadg.stimulus.params.TEMPORAL_PERIOD;
% subplot(1, 2, 2)
figure
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
% title(ll{5})
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

%%
for ct = 1:4
    [id_dir_flash{ct}, idx_dir_flash{ct}] = intersect(ds_id_flash, id_dir{ct});
end
for ct = 1:3
    [id_dir_on_flash{ct}, idx_dir_on_flash{ct}] = intersect(ds_id_flash, id_dir_on{ct});
end

%% parameters
interFlashIntVar = 0.005; % (sec) expected precision of trigger intervals

%% find sets of flash stimuli
ts=1 ; 
trigger_set_i{1} = [] ;
for a=1:length(dataflash.triggers) ; % for each trigger
    trigger_set_i{ts} = [trigger_set_i{ts},a] ; % put it in a set   
    if a<length(dataflash.triggers) ; % if its not the last trigger
        if sum(abs(dataflash.triggers(a+1)-dataflash.triggers(a)-dataflash.DfParams.interFlashInt)<interFlashIntVar)>0 ; % next trigger is has the right interval
            ts = ts ; % keep it in the set
        else
            ts = ts+1 ; % put it in a new set
            trigger_set_i{ts} = [] ;
        end
    end
end
trigger_set_i([1 3 18]) = [];

%% get raster and psth for stimulus trials
bin_size = 0.05; 
start = 0;
% tau = 2*bin_size;
% tt = -3*tau:bin_size:3*tau;
% filter = exp(-tt.^2/(2*tau^2));
% circ = (length(filter)-1)/2;

XX = start+bin_size/2:bin_size:dataflash.DfParams.interFlashInt-bin_size/2;
% load /Volumes/lab/Experiments/Calibration/NdfCalibration
Irel = (dataflash.DfParams.Ftime/1000).*NdfCalibration(2,dataflash.DfParams.NDF+1) ;
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
        raster_1cell{ts} = get_raster(dataflash.spikes{ds_idx_flash(cc)}, ...
            dataflash.triggers(trigger_set_i{ts}), 'start', start, 'stop', ...
            dataflash.DfParams.interFlashInt, 'plot', 0);
        for t = 1:length(raster_1cell{ts})
%             psth_temp = hist(raster_1cell{ts}{t}, XX);
%             psth_temp = [psth_temp(end-circ+1:end) psth_temp psth_temp(1:circ)];
%             hist_1cell{ts}{t} = conv(psth_temp, filter, 'valid');
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

trigger = 0:3:180;
ds_dark_raster = cell(length(ds_id_flash), 1);
ds_dark_raster_all = cell(length(ds_id_flash), 1);
ds_dark_hist_trial = cell(length(ds_id_flash), 1);
ds_dark_hist_mean = cell(length(ds_id_flash), 1);
for cc = 1:length(ds_id_flash)
    ds_dark_raster{cc} = get_raster(dataflash.spikes{ds_idx_flash(cc)}, ...
        trigger, 'stop', dataflash.DfParams.interFlashInt, ...
        'plot', 0);
    for t = 1:length(ds_dark_raster{cc})
%         psth_temp = hist(ds_dark_raster{cc}{t}, XX);
%         psth_temp = [psth_temp(end-circ+1:end) psth_temp psth_temp(1:circ)];
%         ds_dark_hist_trial{cc}{t} = conv(psth_temp, filter, 'valid');
        ds_dark_hist_trial{cc}{t} = hist(ds_dark_raster{cc}{t}, XX);
    end
    ds_dark_raster_all{cc} = sort(cell2mat(ds_dark_raster{cc}));
    ds_dark_hist_mean{cc} = hist(ds_dark_raster_all{cc}, XX)/length(trigger);
end

%% 2 alternative forced choice test
% use 400-580 second (roughly) as dark trials
window = 1;
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
%             template_dark = sqrt(1/300)*ones(1, 300); % use flat psth as dark template
            DV = template_flash - template_dark;
            corr_flash(t) = ds_flash_hist_trial{cc}{ts}{t}(1:bin_n) * DV(1:bin_n)';
            corr_dark(t) = ds_dark_hist_trial{cc}{t}(1:bin_n) * DV(1:bin_n)';
        end
        Pc(cc, ts) = (sum(corr_flash > 0) + sum(corr_dark <= 0))/(trial_n*2);
%         Pc(cc, ts) = sum(corr_flash > 0)/trial_n;
    end
end


color = 'rbgk';
figure
for ct = 1:4
    for cc = 1:length(idx_dir_flash{ct})
        plot(log10(Irel), Pc(idx_dir_flash{ct}(cc), :)', 'color', color(ct))
        hold on
%         pause
    end
end
xlabel('log(R*/rod)')
ylabel('probability')

% compare across cell type
for ct = 1:4
    Pc_dir{ct} = Pc(idx_dir_flash{ct}, :);
    Pc_dir_mean(ct, :) = mean(Pc_dir{ct}, 1);
    Pc_dir_ste(ct, :) = std(Pc_dir{ct}, [], 1)/sqrt(length(idx_dir_flash{ct}));
end

figure
for ct = 1:4
    errorbar(log10(Irel), Pc_dir_mean(ct, :), Pc_dir_ste(ct, :), 'color', color(ct));
    hold on
end
legend('superior', 'Anterior', 'inferior', 'posterior')
xlabel('log(R*/rod)')
ylabel('probability')

%% raster plot

a = ceil((length(trigger_set_i)+1)/2);
for cc =21:21%length(ds_id_flash);
    ts = 1;
    b = 1;
    trial_n = length(ds_dark_raster{1});
    FigHandle = figure;
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

    print_close(1, [24 12], num2str(ds_id_flash(cc)))

end

%% fit
Pc_temp = Pc;
% Pc_temp(:, 4) = [];
Irel_temp = Irel;
% Irel_temp(4) = [];

% superior
ct = 1;
idx = idx_dir_flash{ct};
figure
for cc = 1:length(idx)

    ydata = Pc_temp(idx(cc), 1:10)-0.5;
    xdata = log10(Irel_temp(1:10))+4;
    [f, G] = fit_mm(xdata, ydata);
    fit_all{ct}{cc} = f;
    G_all{ct}{cc} = G;

    x = linspace(min(xdata), max(xdata), 100);
    y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a);

    subplot(4,4,cc)
    plot(xdata, ydata)
    hold on
    plot(x, y)
    if cc == 1
        title('superior')
    end
end

% anterior
ct = 2;
idx = idx_dir_flash{ct};
figure
for cc = 1:length(idx)

    ydata = Pc_temp(idx(cc), 1:15)-0.5;
    xdata = log10(Irel_temp(1:15))+4;
    [f, G] = fit_mm(xdata, ydata);
    fit_all{ct}{cc} = f;
    G_all{ct}{cc} = G;

    x = linspace(min(xdata), max(xdata), 100);
    y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a);

    subplot(3,3,cc)
    plot(xdata, ydata)
    hold on
    plot(x, y)
    if cc == 1
        title('anteriror')
    end

end

% inferior
ct = 3;
idx = idx_dir_flash{ct};
figure
for cc = 1:length(idx)

    ydata = Pc_temp(idx(cc), 1:14)-0.5;
    xdata = log10(Irel_temp(1:14))+4;
    [f, G] = fit_mm(xdata, ydata);
    fit_all{ct}{cc} = f;
    G_all{ct}{cc} = G;

    x = linspace(min(xdata), max(xdata), 100);
    y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a);

    subplot(3,3,cc)
    plot(xdata, ydata)
    hold on
    plot(x, y)
    if cc == 1
        title('inferior')
    end

end

% posterior
ct = 4;
idx = idx_dir_flash{ct};
figure
for cc = 1:length(idx)

    ydata = Pc_temp(idx(cc), 1:15)-0.5;
    xdata = log10(Irel_temp(1:15))+4;
    [f, G] = fit_mm(xdata, ydata);
    fit_all{ct}{cc} = f;
    G_all{ct}{cc} = G;

    x = linspace(min(xdata), max(xdata), 100);
    y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a);

    subplot(2,2,cc)
    plot(xdata, ydata)
    hold on
    plot(x, y)
    if cc == 1
        title('posterior')
    end

end

% close all

%%
xdata = log10(Irel_temp)+4;
x = linspace(min(xdata), max(xdata), 1000);

threshold = 0.34;
for ct = 1:4
    idx = idx_dir_flash{ct};
    for cc = 1:length(idx)
        f = fit_all{ct}{cc};
        y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a);
        [~, i] = min(abs(y-threshold));
%         xthreshold{ct}(cc) = 10^(x(i)-4);
        xthreshold{ct}(cc) = x(i)-4;
        rmse{ct}(cc) = G_all{ct}{cc}.rmse;
    end
end

% xthreshold{4}(1) = [];
color = 'brgk';
figure
for ct = 1:4
    plot(ct*ones(1, length(xthreshold{ct})), xthreshold{ct}, [color(ct) 'o'])
    hold on
    errorbar(ct+0.2, mean(xthreshold{ct}), std(xthreshold{ct})/sqrt(length(xthreshold{ct})), [color(ct) 'd']);
end
xlim([0.5 4.5])
ylabel('log(R*/rod)')
title('Pc = 0.83')
xtick = {'superior'; 'anterior'; 'inferior'; 'posterior'};
set(gca,'XTicklabel',xtick)

