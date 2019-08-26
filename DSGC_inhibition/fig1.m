%% Analyze individual datasets
[Pc_dir_0304, xthreshold_0304] = analyzeFlash160304;
[Pc_dir, xthreshold, Irel, trigger_set_i, ds_id_flash, ds_dark_raster, ds_flash_raster] = analyzeFlash160324;
[nds_Pc, nds_flash_raster, nds_flash_hist_mean, ...
    nds_dark_raster, nds_dark_hist_mean] = get_nds_pc;

for i = 1:3
    Pc_dir{i} = [Pc_dir{i}; Pc_dir_0304{i}];
end

%% trancate sensitivity curve 
for ct = 1:4
    for cc = 1:size(Pc_dir{ct}, 1)
        [m, idx] = max(Pc_dir{ct}(cc, :));
        Pc_dir{ct}(cc, idx+1:end) = ones(1, size(Pc_dir{ct}, 2)-idx) * m;
    end
end
for cc = 1:size(nds_Pc, 1)
    [m, idx] = max(nds_Pc(cc, :));
    nds_Pc(cc, idx+1:end) = ones(1, size(nds_Pc, 2)-idx) * m;
end

for ct = 1:4
    Pc_dir_mean(ct, :) = mean(Pc_dir{ct}, 1);
    Pc_dir_ste(ct, :) = std(Pc_dir{ct}, [], 1)/sqrt(size(Pc_dir{ct}, 1));
end

%% fig1A example rasters
window = 1;
id_all = [5883 7476 4295 4911];
a = ceil((length(trigger_set_i)+1)/2);
FigHandle = figure;
set(FigHandle, 'Position', [0, 0, 1920, 1080]);
cell_types = {'superior', 'anterior', 'inferior', 'posterior', 'ON sustained'};
color = 'brgk';
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
            line(X, Y, 'color', color(mod(b-1, 5) + 1));
            axis([0, window, 0, trial_n]);
            set(h, 'xticklabel',[]);
            set(h, 'yticklabel',[]);
            hold on
        end
        if b <= 5
            title(cell_types{b})
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
        line(X, Y, 'color', [.5 .5 .5]);
        axis([0, window, 0, trial_n]);
        set(h, 'xticklabel',[]);
        set(h, 'yticklabel',[]);
        hold on
    end
    if b == 5
        title('ON sustained')
    end
    b = b+5;
    ts = ts+2;
end

%% fig 1b fit
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

Pc_temp = mean(nds_Pc);
Pc_ste_temp = std(nds_Pc)/sqrt(size(nds_Pc, 1));
Irel_temp = Irel;
a = 13;
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
xlabel('R*/rod')
ylabel('Pc')
legend([h(1), h(2), h(3), h(4), h(5)], 'superior', 'anterior', 'inferior', 'posterior', 'ON alpha');

xlim([-inf 0.2])
ylim([0.48 1])

%% fig 1c response threshold

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

