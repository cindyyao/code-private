%% adaptation

load('snls130221.mat')
datarun{2}.cell_types = [];
datarun{2} = load_params(datarun{2});
display_rate = 60.35;
refresh_rate = [4 2];

cell_type = {'ON brisk transient', 'ON transient', 'OFF brisk transient', ...
    'OFF sustained'};
n = length(cell_type);
[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun{1}, ...
    datarun{2}, cell_type);

NL = cell(4, 1);

for ct = 1:4
    n = size(cell_id{ct}, 1);
    NL_temp = zeros(2, 20, n, 2);
    for cd = 1:2
        for cc = 1:n
            idx = cell_idx{ct}(cc, cd);
            gen_signals = datarun{cd}.stas.snls{idx}.gen_signal;
            gen_signals = gen_signals/std(gen_signals);
            spikes = datarun{cd}.stas.snls{idx}.spikes;

            [X, Y] = curve_from_binning(gen_signals, spikes, ...
                'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
            Y = Y*display_rate/refresh_rate(cd);
            NL_temp(1, :, cc, cd) = X;
            NL_temp(2, :, cc, cd) = Y;
        end
    end
    NL{ct} = NL_temp;
end

XX = zeros(20, 2, 4);
for ct = 1:4
    for cd = 1:2
        xx = squeeze(NL{ct}(1, :, :, cd));
        xx = xx';
        XX_temp = mean(xx);
        XX_temp(1) = max(xx(:, 1));
        XX_temp(20) = min(xx(:, 20));
        XX(:, cd, ct) = XX_temp;
    end
end

NL_Y = zeros(20, 2, 4);
for ct = 1:4
    n = size(cell_id{ct}, 1);
    NL_Y_temp = zeros(20, n, 2);
    for cd = 1:2
        for cc = 1:n
            YY_temp = interp1(NL{ct}(1, :, cc, cd), NL{ct}(2, :, cc, cd), XX(:, cd, ct));
            NL_Y_temp(:, cc, cd) = YY_temp;
        end
    end
    NL_Y_temp = squeeze(mean(NL_Y_temp, 2));
    NL_Y(:, :, ct) = NL_Y_temp;
end

        
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 800, 800]);

for ct = 1:4
    subplot(2, 2, ct)
    plot(XX(:, 1, ct), NL_Y(:, 1, ct))
    hold on 
    plot(XX(:, 2, ct), NL_Y(:, 2, ct), 'r')
    xlabel('generator signals')
    xlim([-2.5 2.5])
    ylabel('spike rate(HZ)')
    title(cell_type{ct})
    if ct == 1
        legend('NDF 4', 'NDF 0', 'location', 'northwest')
    end
    
end

    
%% circadian rod

clear all

load('snls130221.mat')
datarun11 = datarun{1};
load('snls130722-map.mat')
datarun22 = datarun2;
datarun{1} = datarun11;
datarun{2} = datarun22;

display_rate = 60.35;
refresh_rate = [4 2];

cell_type = {'ON brisk transient', 'OFF brisk transient'};
nn = length(cell_type);

NL = cell(nn, 1);

for ct = 1:nn
    NL_temp = cell(2, 1);
    for cd = 1:2
        IDX = get_cell_indices(datarun{cd}, cell_type{ct});
        n = length(IDX);
        NL_ttemp = zeros(2, 20, n);
        for cc = 1:n
            idx = IDX(cc);
            gen_signals = datarun{cd}.stas.snls{idx}.gen_signal;
            gen_signals = gen_signals/std(gen_signals);
            spikes = datarun{cd}.stas.snls{idx}.spikes;

            [X, Y] = curve_from_binning(gen_signals, spikes, ...
                'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
            Y = Y*display_rate/refresh_rate(cd);
            NL_ttemp(1, :, cc) = X;
            NL_ttemp(2, :, cc) = Y;
        end
        NL_temp{cd} = NL_ttemp;
    end
    NL{ct} = NL_temp;
end

XX = zeros(20, 2, nn);
for ct = 1:nn
    for cd = 1:2
        xx = squeeze(NL{ct}{cd}(1, :, :));
        XX_temp = mean(xx, 2);
        XX_temp(1) = max(xx(1, :));
        XX_temp(20) = min(xx(20, :));
        XX(:, cd, ct) = XX_temp;
    end
end

NL_Y = zeros(20, 2, nn);
for ct = 1:nn
    for cd = 1:2
        IDX = get_cell_indices(datarun{cd}, cell_type{ct});
        n = length(IDX);
        NL_Y_temp = zeros(20, n);
        for cc = 1:n
            YY_temp = interp1(NL{ct}{cd}(1, :, cc), NL{ct}{cd}(2, :, cc), XX(:, cd, ct));
            NL_Y_temp(:, cc) = YY_temp;
        end
        NL_Y_temp = squeeze(mean(NL_Y_temp, 2));
        NL_Y(:, cd, ct) = NL_Y_temp;
    end
    
end

        
FigHandle = figure(1);
set(FigHandle, 'Position', [100, 100, 600, 270]);

for ct = 1:nn
    subplot(1, 2, ct)
    plot(XX(:, 1, ct), NL_Y(:, 1, ct))
    hold on 
    plot(XX(:, 2, ct), NL_Y(:, 2, ct), 'r')
    xlabel('generator signals')
    xlim([-2.5 2.5])
    ylabel('spike rate(HZ)')
    title(cell_type{ct})
    if ct == 1
        legend('day', 'night', 'location', 'northwest')
    end
    
end

        
set(figure(1), 'paperpositionmode', 'auto');
print(figure(1), '-dpdf', 'avg_snl')       