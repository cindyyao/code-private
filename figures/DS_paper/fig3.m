cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);


datadg = load_data('/Volumes/lab/analysis/2016-03-24-0/data003-sorted/data003-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-03-24-0/stimuli/s03.mat';
datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);

dataflash = load_data('/Volumes/lab/analysis/2016-03-24-0/data000-map/data000-map', opt);
dataflash.DfParams.NDF =   [5,5,5,4,4,4,4,4,3,3,3,2,2,1,0] ; % on filter turret 
dataflash.DfParams.Ftime = [2,4,8,2,3,4,6,8,2,4,8,8,2,2,2] ; % ms
dataflash.DfParams.interFlashInt = [3] ; % sec


[NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [4 5]; % which parameters to use for classification

[ds_id, nonds_id, id_init] = classify_ds(datadg, ds_struct, params_idx);

load('DS160324.mat')
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

load('DS160324.mat')

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

%% get raster and psth for stimulus trials
bin_size = 0.02; 
start = 0;

tau = 4*bin_size;
tt = -3*tau:bin_size:3*tau;
filter = exp(-tt.^2/(2*tau^2));
circ = (length(filter)-1)/2;

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
            psth_temp = hist(raster_1cell{ts}{t}, XX);
            psth_temp = [psth_temp(end-circ+1:end) psth_temp psth_temp(1:circ)];
            hist_1cell{ts}{t} = conv(psth_temp, filter, 'valid');

%             hist_1cell{ts}{t} = hist(raster_1cell{ts}{t}, XX);
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
        psth_temp = hist(ds_dark_raster{cc}{t}, XX);
        psth_temp = [psth_temp(end-circ+1:end) psth_temp psth_temp(1:circ)];
        ds_dark_hist_trial{cc}{t} = conv(psth_temp, filter, 'valid');

%         ds_dark_hist_trial{cc}{t} = hist(ds_dark_raster{cc}{t}, XX);
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
            DV = template_flash - template_dark;
            corr_flash(t) = ds_flash_hist_trial{cc}{ts}{t}(1:bin_n) * DV(1:bin_n)';
            corr_dark(t) = ds_dark_hist_trial{cc}{t}(1:bin_n) * DV(1:bin_n)';
        end
        Pc(cc, ts) = (sum(corr_flash > 0) + sum(corr_dark <= 0))/(trial_n*2);
    end
end

for ct = 1:4
    Pc_dir{ct} = Pc(idx_dir_flash{ct}, :);
end
load('DS160304.mat', 'Pc_dir_0304')

for i = 1:3
    Pc_dir{i} = [Pc_dir{i}; Pc_dir_0304{i}];
end

for ct = 1:4
    Pc_dir_mean(ct, :) = mean(Pc_dir{ct}, 1);
    Pc_dir_ste(ct, :) = std(Pc_dir{ct}, [], 1)/sqrt(size(Pc_dir{ct}, 1));
end

%% fit
color = 'brgkc';
Pc_temp = Pc_dir_mean;
Pc_ste_temp = Pc_dir_ste;
Irel_temp = Irel;
a = 13;
figure
for ct = 1:4
    ydata = Pc_temp(ct, :)-0.5;
    xdata = log10(Irel_temp)+4;
    [f, G] = fit_mm(xdata, ydata, 'Startpoints', [0.5 2 0]);
    fit_avg{ct} = f;
    G_avg{ct} = G;

    x = linspace(min(xdata(1:a)), max(xdata(1:a)), 100);
    y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a);

    errorbar(log10(Irel(1:a)), Pc_dir_mean(ct, 1:a), Pc_dir_ste(ct, 1:a), [color(ct) 'o']);
    hold on
    h(ct) = plot(x-4, y+0.5, 'color', color(ct));
end
xlim([-inf 0.5])

load('DS160304.mat')

Pc_temp = mean(nds_Pc);
Pc_ste_temp = std(nds_Pc)/sqrt(size(nds_Pc, 1));
Irel_temp = Irel;
a = 13;
% for ct = 1:4
    ydata = Pc_temp-0.5;
    xdata = log10(Irel_temp)+4;
    [f, G] = fit_mm(xdata, ydata, 'Startpoints', [0.5 2 0]);
    fit_avg_nds = f;
    G_avg_nds = G;

    x = linspace(min(xdata(1:a)), max(xdata(1:a)), 100);
    y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a);

    errorbar(log10(Irel(1:a)), Pc_temp(1:a), Pc_ste_temp(1:a), 'LineStyle', 'none', 'Marker', 'o', 'color', [.5 .5 .5]);
    hold on
    h(5) = plot(x-4, y+0.5, 'color', [.5 .5 .5]);
% end
xlabel('R*/rod')
ylabel('Pc')
legend([h(1), h(2), h(3), h(4), h(5)], 'superior', 'Anterior', 'inferior', 'posterior', 'ON alpha');

%% raster plot
% id_all = [5883 7476 4295 4911];
% a = ceil((length(trigger_set_i)+1)/2);
% FigHandle = figure;
% set(FigHandle, 'Position', [0, 0, 1920, 1080]);
% for i = 1:4
%     id = id_all(i);
%     cc = find(ds_id_flash == id);
%     ts = 1;
%     b = (i-1)*2 + 1;
%     trial_n = length(ds_dark_raster{1});
%     while ts <= length(trigger_set_i)+1
% %         if  ts == a+1
% %             b = 3;
% %         end
%         subplot(a, 8, b)
%         if ts > 1
%             trial_n = length(ds_flash_raster{cc}{ts-1});
%         end
%         for j = 1:trial_n
%             if ts == 1
%                 SpikeTime = ds_dark_raster{cc}{j};
%                 SpikeTime = SpikeTime(SpikeTime < window);
%             else
%                 SpikeTime = ds_flash_raster{cc}{ts-1}{j};
%                 SpikeTime = SpikeTime(SpikeTime < window);
%             end
%             SpikeTime = SpikeTime';
%             X = [SpikeTime; SpikeTime];
%             Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
%             line(X, Y, 'color', 'b');
%             axis([0, window, 0, trial_n]);
%             hold on
%         end
%         if b == 1
%             title(num2str(ds_id_flash(cc)))
%         end
%         b = b+1;
%         subplot(a, 8, b)
%         if ts > 1
%             plot(XX(1:bin_n), ds_flash_hist_mean{cc}{ts-1}(1:bin_n)/bin_size)
%         else
%             plot(XX(1:bin_n), ds_dark_hist_mean{cc}(1:bin_n)/bin_size)
%         end
%         xlim([0 window])
% %         ylim([0 4])
% 
%         b = b+7;
%         ts = ts+2;
%     end
% 
% %     print_close(1, [24 12], num2str(ds_id_flash(cc)))
% 
% end

%%
id_all = [5883 7476 4295 4911];
a = ceil((length(trigger_set_i)+1)/2);
FigHandle = figure;
set(FigHandle, 'Position', [0, 0, 1920, 1080]);
for i = 1:4
    id = id_all(i);
    cc = find(ds_id_flash == id);
    ts = 1;
    b = i;
    trial_n = length(ds_dark_raster{1});
    while ts <= length(trigger_set_i)+1
        h = subplot(a, 5, b);
        if ts > 1
            trial_n = length(ds_flash_raster{cc}{ts-1});
        end
        for j = 1:trial_n
            if ts == 1
                SpikeTime = ds_dark_raster{cc}{j};
                SpikeTime = SpikeTime(SpikeTime < window);
            else
                SpikeTime = ds_flash_raster{cc}{ts-1}{j};
                SpikeTime = SpikeTime(SpikeTime < window);
            end
            SpikeTime = SpikeTime';
            X = [SpikeTime; SpikeTime];
            Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
            line(X, Y, 'color', 'b');
            axis([0, window, 0, trial_n]);
            set(h, 'xticklabel',[]);
            set(h, 'yticklabel',[]);
            hold on
        end
        if b == 1
            title(num2str(ds_id_flash(cc)))
        end
        b = b+5;
        ts = ts+2;
    end
end

cc = 10;
ts = 1;
b = 5;
trial_n = length(nds_dark_raster{2});
    while ts <= length(trigger_set_i)+1
        h = subplot(a, 5, b);
        if ts > 1
            trial_n = length(nds_flash_raster{cc}{ts-1});
        end
        for j = 1:trial_n
            if ts == 1
                SpikeTime = nds_dark_raster{cc}{j};
                SpikeTime = SpikeTime(SpikeTime < window);
            else
                SpikeTime = nds_flash_raster{cc}{ts-1}{j};
                SpikeTime = SpikeTime(SpikeTime < window);
            end
            SpikeTime = SpikeTime';
            X = [SpikeTime; SpikeTime];
            Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
            line(X, Y, 'color', 'b');
            axis([0, window, 0, trial_n]);
            set(h, 'xticklabel',[]);
            set(h, 'yticklabel',[]);
            hold on
        end
        if b == 1
            title(num2str(ds_id_flash(cc)))
        end
        b = b+5;
        ts = ts+2;
    end

%% fit
Pc_temp = Pc;
Irel_temp = Irel;

% superior
ct = 1;
idx = idx_dir_flash{ct};
for cc = 1:length(idx)

    ydata = Pc_temp(idx(cc), 1:10)-0.5;
    xdata = log10(Irel_temp(1:10))+4;
    [f, G] = fit_mm(xdata, ydata);
    fit_all{ct}{cc} = f;
    G_all{ct}{cc} = G;
end

% anterior
ct = 2;
idx = idx_dir_flash{ct};
for cc = 1:length(idx)

    ydata = Pc_temp(idx(cc), 1:11)-0.5;
    xdata = log10(Irel_temp(1:11))+4;
    [f, G] = fit_mm(xdata, ydata);
    fit_all{ct}{cc} = f;
    G_all{ct}{cc} = G;
end

% inferior
ct = 3;
idx = idx_dir_flash{ct};
for cc = 1:length(idx)

    ydata = Pc_temp(idx(cc), 1:14)-0.5;
    xdata = log10(Irel_temp(1:14))+4;
    [f, G] = fit_mm(xdata, ydata);
    fit_all{ct}{cc} = f;
    G_all{ct}{cc} = G;
end

% posterior
ct = 4;
idx = idx_dir_flash{ct};
for cc = 1:length(idx)

    ydata = Pc_temp(idx(cc), 1:15)-0.5;
    xdata = log10(Irel_temp(1:15))+4;
    [f, G] = fit_mm(xdata, ydata);
    fit_all{ct}{cc} = f;
    G_all{ct}{cc} = G;
end

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

load('DS160304.mat', 'xthreshold_0304')
for ct = 1:3
    xthreshold{ct} = [xthreshold{ct} xthreshold_0304{ct}];
end

% exclude outliers
for ct = 1:4
    notdone = 1;
    xthreshold_temp = xthreshold{ct};
    while notdone
        a = length(xthreshold_temp);
        xthreshold_temp(xthreshold_temp > std(xthreshold_temp)*2+mean(xthreshold_temp)) = [];
        b = length(xthreshold_temp);
        if a == b
            notdone = 0;
            xthreshold{ct} = xthreshold_temp;
        end
    end
end

% xthreshold{4}(1) = [];
color = 'brgk';
figure
for ct = 1:4
    plot(ct*ones(1, length(xthreshold{ct})), xthreshold{ct}, [color(ct) 'o'])
    hold on
    errorbar(ct+0.2, mean(xthreshold{ct}), std(xthreshold{ct})/sqrt(length(xthreshold{ct})), [color(ct) 'd'], 'MarkerSize', 3);
end
xlim([0.5 4.5])
% ylim([-3 0])
ylabel('log(R*/rod)')
title('Pc = 0.83')
xtick = {'superior'; 'anterior'; 'inferior'; 'posterior'};
set(gca,'XTicklabel',xtick)

%
ct = {'superior', 'anterior', 'inferior', 'posterior'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = ct;
model_series = cellfun(@mean, xthreshold)';
model_error = cellfun(@std, xthreshold)'./sqrt(cellfun(@length, xthreshold))';
h = bar(10.^model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
set(gca, 'YScale', 'Log')
% ylabel('RF gain (spike count)')
% legend('NDF 4','NDF 2', 'NDF 0');
% title('OFF')
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, 10.^model_series(:,i), 10.^(model_series(:,i)-model_error(:,i)), 10.^(model_series(:,i)+model_error(:,i)), 'k', 'linestyle', 'none');
end

