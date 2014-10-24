opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);

cell_type = {'ON 1', 'ON 2', 'OFF 1', 'OFF 2', 'OFF 3', 'OFF 4'};
datarun{1} = load_data('/Analysis/xyao/2014-01-14-0/data008-map/data008-map', opt);

datarun{1} = load_java_movie(datarun{1}, '/Volumes/lab/acquisition/movie-xml/BW-15-1-0.48-11111-40x40-60.35.xml');
datarun{1} = get_sta_summaries(datarun{1}, 'all');

datarun{1} = get_snls(datarun{1}, cell_type);

%% get cell ids
n = length(cell_type);

cell_id = cell(n, 1);
cell_idx = cell(n, 1);

for i =  1:n
    for j = 1:2
        cell_id{i}{j} = get_cell_ids(datarun{j}, cell_type{i});
        cell_idx{i}{j} = get_cell_indices(datarun{j}, cell_type{i});
    end
end

%% get average snls

display_rate = 60.35;
refresh_rate = [1 1];

NL = cell(n, 1);

for ct = 1:n
    for cd = 1:2
    cell_num = length(cell_id{ct}{cd});
    NL_temp{cd} = zeros(2, 20, cell_num);    
        for cc = 1:cell_num
            idx = cell_idx{ct}{cd}(cc);
            gen_signals = datarun{cd}.stas.snls{idx}.gen_signal;
            gen_signals = gen_signals/std(gen_signals);
            spikes = datarun{cd}.stas.snls{idx}.spikes;

            [X, Y] = curve_from_binning(gen_signals, spikes, ...
                'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
            Y = Y*display_rate/refresh_rate(cd);
            NL_temp{cd}(1, :, cc) = X;
            NL_temp{cd}(2, :, cc) = Y/max(Y);
        end
        
    end
    NL{ct} = NL_temp;
end

XX = zeros(20, 2, n);
for ct = 1:n
    for cd = 1:2
        xx = squeeze(NL{ct}{cd}(1, :, :));
        XX_temp = mean(xx, 2);
        XX_temp(1) = max(xx(1, :));
        XX_temp(20) = min(xx(20, :));
        XX(:, cd, ct) = XX_temp;
    end
end

NL_Y_mean = zeros(20, 2, n);
Nonlinearity_stev = zeros(20, 2, n);

for ct = 1:n
    for cd = 1:2
    cell_num = length(cell_id{ct}{cd});
    NL_Y_temp{cd} = zeros(20, cell_num);
        for cc = 1:cell_num
            YY_temp = interp1(NL{ct}{cd}(1, :, cc), NL{ct}{cd}(2, :, cc), XX(:, cd, ct));
            NL_Y_temp{cd}(:, cc) = YY_temp;
        end
    NL_Y_mean_temp(cd, :) = squeeze(mean(NL_Y_temp{cd}, 2));
    Nonlinearity_stev_temp(cd, :) = squeeze(std(NL_Y_temp{cd}, 0, 2))/sqrt(cell_num);

    end
    NL_Y_mean(:, :, ct) = NL_Y_mean_temp';
    Nonlinearity_stev(:, :, ct) = Nonlinearity_stev_temp';
end

Nonlinearity_mean(:, :, :, 1) = XX;
Nonlinearity_mean(:, :, :, 2) = NL_Y_mean;


