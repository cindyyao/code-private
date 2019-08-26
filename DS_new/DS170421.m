%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

datadg = load_data('/Volumes/lab/analysis/2017-04-21-0/data006-sorted/data006-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2017-04-21-0/stimuli/s06.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);
load('DS170421.mat');

% identify DS cells
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [2 3]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2017-04-21-0/data000-003-map/data000-003-map', opt);
dataflash = split_datarun(datarun, [4008 7538]);
dataflash{1}.DfParams.NDF =   [5,5,5,5,4,4,4,4,4,3,3,3,2,2,1] ; % on filter turret 
dataflash{1}.DfParams.Ftime = [1,2,4,8,2,3,4,6,8,2,4,8,2,8,2] ; % ms
dataflash{1}.DfParams.interFlashInt = [3] ; % sec

dataflash{2}.DfParams.NDF =   [5,5,5,5,4,4,4,4,4,3,3,3,2,2] ; % on filter turret 
dataflash{2}.DfParams.Ftime = [1,2,4,8,2,3,4,6,8,2,4,8,2,8] ; % ms
dataflash{2}.DfParams.interFlashInt = [3] ; % sec

dataflash{3}.DfParams.NDF =   [5,5,5,5,4,4,4,4,4,3,3,3,2,2,1] ; % on filter turret 
dataflash{3}.DfParams.Ftime = [1,2,4,8,2,3,4,6,8,2,4,8,2,8,2] ; % ms
dataflash{3}.DfParams.interFlashInt = [3] ; % sec

%% dg
% ds_id = datadg.cell_ids;
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);

delta_p = 3; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(DG{i});
rep = datadg.stimulus.repetitions;

%%
d = 1;
t = 3;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(ds_struct.U{t}, ds_struct.V{t});
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(ds_struct.U{t}, ds_struct.V{t}, x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = I;
    id_dir{i} = ds_id(idx_dir{i});
end

%%
ds_id_flash = ds_id(~flash_idx);
for ct = 1:4
    [id_dir_flash{ct}, idx_dir_flash{ct}] = intersect(ds_id_flash, id_dir{ct});
end
id_dir_flash{2}(8) = [];
idx_dir_flash{2}(8) = [];
id_dir_flash{1}(10) = [];
idx_dir_flash{1}(10) = [];
%% plot individual cell
for cc = 1:length(ds_id)
    plot_ds_raster(DG, raster_dg, cc, ds_id(cc), '', 1, 1, 0)
    pause
    close all
end

%% parameters
interFlashIntVar = 0.005; % (sec) expected precision of trigger intervals

%% find sets of flash stimuli
fig_title = {'control', 'SR', 'wash'};
for ds = 1:3
    clear trigger_set_i
    ts = 1; 
    trigger_set_i{1} = [] ;
    for a=1:length(dataflash{ds}.triggers) ; % for each trigger
        trigger_set_i{ts} = [trigger_set_i{ts},a] ; % put it in a set   
        if a<length(dataflash{ds}.triggers) ; % if its not the last trigger
            if sum(abs(dataflash{ds}.triggers(a+1)-dataflash{ds}.triggers(a)-dataflash{ds}.DfParams.interFlashInt)<interFlashIntVar) == 0 ; % next trigger is has the right interval
                ts = ts+1 ; % put it in a new set
                trigger_set_i{ts} = [] ;
            end
        end
    end

    i = 1;
    while(i<=length(trigger_set_i))
        if length(trigger_set_i{i}) < 50 || length(trigger_set_i{i}) > 70
            trigger_set_i(i) = [];
        else
            i = i+1;
        end
    end
    trigger_set_i_cd{ds} = trigger_set_i;

    %%
    ds_idx_flash = get_cell_indices(dataflash{ds}, ds_id_flash);
    bin_size = 0.02; 
    start = 0;
    XX = start+bin_size/2:bin_size:dataflash{ds}.DfParams.interFlashInt-bin_size/2;
    load('DS161208.mat', 'NdfCalibration')
    % NdfCalibration(2, :) = NdfCalibration(2, :)/10;
    Irel{ds} = (dataflash{ds}.DfParams.Ftime/1000).*NdfCalibration(2,dataflash{ds}.DfParams.NDF+1) ;
    [Irel{ds}, i] = sort(Irel{ds});

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
                dataflash{ds}.DfParams.interFlashInt, 'plot', 0);
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
    ds_flash_raster_cd{ds} = ds_flash_raster;
    ds_flash_hist_mean_cd{ds} = ds_flash_hist_mean;
    ds_flash_hist_trial_cd{ds} = ds_flash_hist_trial;
    %% get raster and psth for dark trials
    % use 400-580 second (roughly) as dark trials
    if ds == 2
        trigger = linspace(445, 495, 60);
    else
        trigger = 1:3:180;
    end

    ds_dark_raster = cell(length(ds_id_flash), 1);
    ds_dark_raster_all = cell(length(ds_id_flash), 1);
    ds_dark_hist_trial = cell(length(ds_id_flash), 1);
    ds_dark_hist_mean = cell(length(ds_id_flash), 1);
    for cc = 1:length(ds_id_flash)
        ds_dark_raster{cc} = get_raster(dataflash{ds}.spikes{ds_idx_flash(cc)}, ...
            trigger, 'stop', dataflash{ds}.DfParams.interFlashInt, ...
            'plot', 0);
        for t = 1:length(ds_dark_raster{cc})
            ds_dark_hist_trial{cc}{t} = hist(ds_dark_raster{cc}{t}, XX);
        end
        ds_dark_raster_all{cc} = sort(cell2mat(ds_dark_raster{cc}));
        ds_dark_hist_mean{cc} = hist(ds_dark_raster_all{cc}, XX)/length(trigger);
    end
    ds_dark_raster_cd{ds} = ds_dark_raster;
    ds_dark_hist_mean_cd{ds} = ds_dark_hist_mean;
    ds_dark_hist_trial_cd{ds} = ds_dark_hist_trial;
    %% 2 alternative forced choice test
    % use 400-580 second (roughly) as dark trials
    window = 1;
    bin_n = window/bin_size;

    Pc{ds} = zeros(length(ds_id_flash), length(trigger_set_i));
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
                template_flash = template_flash(1:bin_n);
                template_flash = template_flash/norm(template_flash);
                template_flash(isnan(template_flash)) = 0;
                template_dark = ds_dark_hist_sum - ds_dark_hist_trial{cc}{t};
                template_dark = template_dark(1:bin_n);
                template_dark = template_dark/norm(template_dark);
                template_dark(isnan(template_dark)) = 0;
                DV = template_flash - template_dark;
                corr_flash(t) = ds_flash_hist_trial{cc}{ts}{t}(1:bin_n) * DV(1:bin_n)';
                corr_dark(t) = ds_dark_hist_trial{cc}{t}(1:bin_n) * DV(1:bin_n)';
            end
            Pc{ds}(cc, ts) = (sum(corr_flash > corr_dark) + sum(corr_flash == corr_dark)/2)/trial_n;
%             Pc{ds}(cc, ts) = (sum(corr_flash > 0) + sum(corr_dark <= 0))/(trial_n*2);
        end
    end




    color = 'rbgk';
    figure(10)
    subplot(2,2,ds)
    for ct = 1:4
        plot(log10(Irel{ds}), Pc{ds}(idx_dir_flash{ct}, :)', 'color', color(ct))
        hold on
    end
    xlabel('log(R*/rod)')
    ylabel('probability')
    title(fig_title{ds})
    ylim([0.45 1])
    % compare across cell type
    for ct = 1:4
        Pc_dir{ds}{ct} = Pc{ds}(idx_dir_flash{ct}, :);
        Pc_dir_mean{ds}(ct, :) = mean(Pc_dir{ds}{ct}, 1);
        Pc_dir_ste{ds}(ct, :) = std(Pc_dir{ds}{ct}, [], 1)/sqrt(length(idx_dir_flash{ct}));
    end

    figure(11)
    subplot(2,2,ds)
    % Irel(5) = [];
    % Pc_dir_mean(:, 5) = [];
    % Pc_dir_ste(:, 5) = [];
    for ct = 1:4
        errorbar(log10(Irel{ds}), Pc_dir_mean{ds}(ct, :), Pc_dir_ste{ds}(ct, :), 'color', color(ct));
        hold on
    end
    legend('superior', 'anterior', 'inferior', 'posterior', 'location', 'northwest')
    xlabel('log(R*/rod)')
    ylabel('probability')
    title(fig_title{ds})
    ylim([0.45 1])

end

%% trancate Pc
Pc_truncate = Pc;
for ds = 1:3
    for cc = 1:size(Pc_truncate{ds}, 1)
        [m, i] = max(Pc_truncate{ds}(cc, :));
        Pc_truncate{ds}(cc, i:end) = ones(1, size(Pc_truncate{ds}, 2)-i+1) * m;
    end
end

color = 'rbgk';
for ds = 1:3
    figure(10)
    subplot(2,2,ds)
    for ct = 1:4
        plot(log10(Irel{ds}), Pc_truncate{ds}(idx_dir_flash{ct}, :)', 'color', color(ct))
        hold on
    end
    xlabel('log(R*/rod)')
    ylabel('probability')
    title(fig_title{ds})
    ylim([0.45 1])
    % compare across cell type
    for ct = 1:4
        Pc_trancate_dir{ds}{ct} = Pc_truncate{ds}(idx_dir_flash{ct}, :);
        Pc_trancate_dir_mean{ds}(ct, :) = mean(Pc_trancate_dir{ds}{ct}, 1);
        Pc_trancate_dir_ste{ds}(ct, :) = std(Pc_trancate_dir{ds}{ct}, [], 1)/sqrt(length(idx_dir_flash{ct}));
    end

    figure(11)
    subplot(2,2,ds)
    % Irel(5) = [];
    % Pc_dir_mean(:, 5) = [];
    % Pc_dir_ste(:, 5) = [];
    for ct = 1:4
        errorbar(log10(Irel{ds}), Pc_trancate_dir_mean{ds}(ct, :), Pc_trancate_dir_ste{ds}(ct, :), 'color', color(ct));
        hold on
    end
    legend('superior', 'anterior', 'inferior', 'posterior', 'location', 'northwest')
    xlabel('log(R*/rod)')
    ylabel('probability')
    title(fig_title{ds})
    ylim([0.45 1])
end
%% raster plot
window = 1;
a = ceil((length(trigger_set_i)+1)/2);
for cc =2:2%length(ds_id_flash);
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
            SpikeTime(SpikeTime > window) = [];
            X = [SpikeTime; SpikeTime];
            Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
            line(X, Y, 'color', 'b');
            axis([0, window, 0, trial_n]);
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
        xlim([0 window])

        b = b+3;
        ts = ts+1;
    end

%     print_close(1, [24 12], num2str(ds_id_flash(cc)))

end

% end

%% fit
threshold = 0.34;
CONDITIONS = {'control', 'SR', 'wash'};
figure
for ds = 1:3
    Irel_temp = Irel{ds};
    xdata = log10(Irel_temp)+4;
    x = linspace(min(xdata), max(xdata), 1000);

    Pc_temp = Pc_truncate{ds};
    subplot(2,2,ds)
    for ct = 1:4
        idx = idx_dir_flash{ct};
        for cc = 1:length(idx)
            ydata = Pc_temp(idx(cc), :)-0.5;
            [f, G] = fit_mm(xdata, ydata);
            fit_all{ds}{ct}{cc} = f;
            G_all{ds}{ct}{cc} = G;
            y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a);
            [~, i] = min(abs(y-threshold));
            xthreshold{ds}{ct}(cc) = x(i)-4;
            rmse{ds}{ct}(cc) = G_all{ds}{ct}{cc}.rmse;
        end
        xthreshold{ds}{ct}(xthreshold{ds}{ct} > 0.5) = [];
        ydata = Pc_dir_mean{ds}(ct, :)-0.5;
        [f, G] = fit_mm(xdata, ydata, 'Startpoints', [0.5 2 0]);
        fit_avg{ds}{ct} = f;
        G_avg{ds}{ct} = G;
        y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a);

        errorbar(log10(Irel{ds}), Pc_dir_mean{ds}(ct, :), Pc_dir_ste{ds}(ct, :), [color(ct) 'o']);
        hold on
        h(ct) = plot(x-4, y+0.5, 'color', color(ct));
    end
    xlim([-4 0.5])
    ylim([0.48 1])

    legend([h(1), h(2), h(3), h(4)], 'superior', 'anterior', 'inferior', 'posterior')
    xlabel('Log(intensity) (R*/rod)')
    ylabel('% correct')
    title(CONDITIONS{ds})
    
    color = 'brgk';
end

figure
for ds = 1:3
    for ct = 1:4
        plot(ct*ones(1, length(xthreshold{ds}{ct})), xthreshold{ds}{ct}, [color(ct) 'o'])
        hold on
        errorbar(ct+0.2, mean(xthreshold{ds}{ct}), std(xthreshold{ds}{ct})/sqrt(length(xthreshold{ds}{ct})), [color(ct) 'd']);
    end
    xlim([0.5 4.5])
    ylabel('log(R*/rod)')
    title(CONDITIONS{ds})
    xtick = {'superior'; 'anterior'; 'inferior'; 'posterior'};
    set(gca,'XTicklabel',xtick)
end

for ds = 1:3
    xthreshold_combine{ds}{1} = xthreshold{ds}{1};
    xthreshold_combine{ds}{2} = cell2mat(xthreshold{ds}(2:4));
end

% exclude outliers
for ds = 1:3
    for ct = 1:2
        notdone = 1;
        xthreshold_temp = xthreshold_combine{ds}{ct};
        while notdone
            a = length(xthreshold_temp);
            xthreshold_temp(xthreshold_temp > std(xthreshold_temp)*1.7+mean(xthreshold_temp)) = [];
            b = length(xthreshold_temp);
            if a == b
                notdone = 0;
                xthreshold_combine_clean{ds}{ct} = xthreshold_temp;
            end
        end
    end
    figure
    for ct = 1:2
        plot(ct*ones(1, length(xthreshold_combine_clean{ds}{ct})), xthreshold_combine_clean{ds}{ct}, [color(ct) 'o'])
        hold on
        errorbar(ct+0.2, mean(xthreshold_combine_clean{ds}{ct}), std(xthreshold_combine_clean{ds}{ct})/sqrt(length(xthreshold_combine_clean{ds}{ct})), [color(ct) 'd']);
    end
    xlim([0.5 2.5])
    ylabel('log(R*/rod)')
    title(CONDITIONS{ds})
    xtick = {'superior'; 'other'};
    set(gca,'XTicklabel',xtick)

end

marker = 'od';
figure
for ds = 1:2
    for ct = 1:2
        plot((ds-1)*2+ct*ones(1, length(xthreshold_combine_clean{ds}{ct})), xthreshold_combine_clean{ds}{ct}, ['k' marker(ct)], 'markersize', 7)
        hold on
    end
    errorbar((ds-1)*2+[1 2], cellfun(@mean, xthreshold_combine_clean{ds}), cellfun(@std, xthreshold_combine_clean{ds})./sqrt(cellfun(@length, xthreshold_combine_clean{ds})), 'ks-', 'markersize', 15);
end
xlim([0.5 4.5])
ylabel('log(R*/rod)') 
legend('superior', 'other')
%%
id_all = [391 2731];
stim_n = 14;
a = ceil((stim_n+1)/2);

for i = 1:2
    figure
    id = id_all(i);
    cc = find(ds_id_flash == id);
    b = 1;
    for ds = 1:3
        ds_dark_raster = ds_dark_raster_cd{ds};
        ds_flash_raster = ds_flash_raster_cd{ds};
        ts = 1;
        trial_n = length(ds_dark_raster{1});
        while ts <= length(trigger_set_i)+1
            h = subplot(3, a, b);
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
                line(X, Y, 'color', 'r');
                axis([0, window, 0, trial_n]);
                set(h, 'xticklabel',[]);
                set(h, 'yticklabel',[]);
                hold on
            end
            if b == 1
                title(num2str(ds_id_flash(cc)))
            end
            b = b+1;
            ts = ts+2;
        end
    end
end

%% create trancated Pc
Pc_dir_trancated = Pc_dir;
for i = 1:length(Pc_dir_trancated)
    for ct = 1:4
        curve_length = size(Pc_dir_trancated{i}{ct}, 2);
        for cc = 1:size(Pc_dir_trancated{i}{ct}, 1)
            [~, idx] = max(Pc_dir_trancated{i}{ct}(cc, :));
            Pc_dir_trancated{i}{ct}(cc, idx:end) = ones(1, curve_length - idx + 1) * Pc_dir_trancated{i}{ct}(cc, idx);
        end
        Pc_dir_trancated_mean{i}{ct} = mean(Pc_dir_trancated{i}{ct});
        Pc_dir_trancated_ste{i}{ct} = std(Pc_dir_trancated{i}{ct}) / sqrt(size(Pc_dir_trancated{i}{ct}, 1));
    end
end

figure
for i = 1:length(Pc_dir_trancated)
    subplot(2,2,i)
    Irel_temp = Irel{i};
    xdata = log10(Irel_temp)+4;
    x = linspace(min(xdata), max(xdata), 1000);
    for ct = 1:4
        ydata = Pc_dir_trancated_mean{i}{ct}-0.5;
        [f, G] = fit_mm(xdata, ydata, 'Startpoints', [0.5 2 0]);
        fit_avg{i}{ct} = f;
        G_avg{i}{ct} = G;
        y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a);
        errorbar(log10(Irel{i}), Pc_dir_trancated_mean{i}{ct}, Pc_dir_trancated_ste{i}{ct}, [color(ct) 'o']);
        hold on
        h(ct) = plot(x-4, y+0.5, 'color', color(ct));
    end
end

%% raster plot
window = 1;
a = ceil((max(cellfun(@length, trigger_set_i_cd))+1)/2);
for i = 2:2
    trigger_set_i = trigger_set_i_cd{i};
    for cc =2:2%length(ds_id_flash);
        ts = 1;
        b = 1;
        trial_n = length(ds_dark_raster_cd{i}{1});
        FigHandle = figure;
        set(FigHandle, 'Position', [0, 0, 1920, 1080]);
        while ts <= length(trigger_set_i)+1
            if  ts == a+1
                b = 3;
            end
            subplot(a, 8, b)
            if ts > 1
                trial_n = length(ds_flash_raster_cd{i}{cc}{ts-1});
            end
            for j = 1:trial_n
                if ts == 1
                    SpikeTime = ds_dark_raster_cd{i}{cc}{j};
                else
                    SpikeTime = ds_flash_raster_cd{i}{cc}{ts-1}{j};
                end
                SpikeTime = SpikeTime';
                SpikeTime(SpikeTime > window) = [];
                X = [SpikeTime; SpikeTime];
                Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
                line(X, Y, 'color', 'b');
                axis([0, window, 0, trial_n]);
                hold on
            end
            if b == 1
                title(num2str(ds_id_flash(cc)))
            end
            b = b+1;
            subplot(a, 8, b)
            if ts > 1
                bar(XX, ds_flash_hist_mean_cd{i}{cc}{ts-1}, 1)
            else
                bar(XX, ds_dark_hist_mean_cd{i}{cc}, 1)
            end
            xlim([0 window])

            b = b+7;
            ts = ts+1;
        end
        subplot(a, 8, [5 6 7 8 13 14 15 16 21 22 23 24 29 30 31 32])
        plot(log10(Irel{i}), Pc{i}(cc, :))
        xlim([-3.8 0.5])
        ylim([0.48 1])
    %     print_close(1, [24 12], num2str(ds_id_flash(cc)))

    end
end