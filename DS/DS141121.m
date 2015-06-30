%% load data
addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/Classification/
addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/DS' cell analysis'/

opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load data
datadg{1} = load_data('/Analysis/xyao/2014-11-21-0/data000-map/data000-map', opt);
datadg{1}.names.stimulus_path = '/Analysis/xyao/2014-11-21-0/stimuli/s00';
datadg{1} = load_stim(datadg{1}, 'user_defined_trigger_interval', 10);

% datadg{2} = load_data('/Analysis/xyao/2014-11-21-0/data002-map/data002-map', opt);
% datadg{2}.names.stimulus_path = '/Analysis/xyao/2014-11-21-0/stimuli/s02';
% datadg{2} = load_stim(datadg{2}, 'user_defined_trigger_interval', 10);

datadg{2} = load_data('/Analysis/xyao/2014-11-21-0/data004-map/data004-map', opt);
datadg{2}.names.stimulus_path = '/Analysis/xyao/2014-11-21-0/stimuli/s04';
datadg{2} = load_stim(datadg{2}, 'user_defined_trigger_interval', 10);

datadg{3} = load_data('/Analysis/xyao/2014-11-21-0/data006-map/data006-map', opt);
datadg{3}.names.stimulus_path = '/Analysis/xyao/2014-11-21-0/stimuli/s06';
datadg{3} = load_stim(datadg{3}, 'user_defined_trigger_interval', 10);

datadg{4} = load_data('/Analysis/xyao/2014-11-21-0/data008/data008', opt);
datadg{4}.names.stimulus_path = '/Analysis/xyao/2014-11-21-0/stimuli/s08';
datadg{4} = load_stim(datadg{4}, 'user_defined_trigger_interval', 10);

dataffp{1} = load_data('/Analysis/xyao/2014-11-21-0/data001-map/data001-map', opt);
dataffp{2} = load_data('/Analysis/xyao/2014-11-21-0/data003-map/data003-map', opt);
dataffp{3} = load_data('/Analysis/xyao/2014-11-21-0/data005-map/data005-map', opt);
dataffp{4} = load_data('/Analysis/xyao/2014-11-21-0/data007-map/data007-map', opt);
dataffp{5} = load_data('/Analysis/xyao/2014-11-21-0/data009-map/data009-map', opt);

%% pull out DS cells
dt = 4;

[NumSpikesCell, StimComb] = get_spikescellstim(datadg{dt},datadg{dt}.cell_ids,0);
dg_struct = dscellanalysis(NumSpikesCell, StimComb);

figure
plot(dg_struct.mag{1, 1}, dg_struct.mag{2, 1}, 'o')
hold on
[x, y] = ginput;
plot(x, y);
IN = inpolygon(dg_struct.mag{1, 1}, dg_struct.mag{2, 1}, x, y);
[~, I] = find(IN == 1);
id_init = datadg{dt}.cell_ids(I);

[C ia ib] = intersect(id_init, datadg{dt}.cell_ids);
vc = ones(length(datadg{dt}.cell_ids),1);
vc(ib) = 2; %initializing ds cells to cluster 2, everything else cluster 1

close all;
X = [];
N = [];
p = [];
X(:,1) = log(dg_struct.mag{1,1})';
X(:,2) = log(dg_struct.mag{2,1})';
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 0, 1, 0, 0, 0,0, vc);

ds_id = [];
ds_id = datadg{dt}.cell_ids(idx==2);
nonds_id = datadg{dt}.cell_ids(idx==1);

%% raster

n_dg = 4;
[raster, raster_p_sum, p_idx] = deal(cell(n_dg, 1));
for i = 1:n_dg   
    [NumSpikesCell, StimComb] = get_spikescellstim(datadg{i},ds_id,0);
    DS{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster{i} = get_ds_raster(datadg{i}, ds_id);
%     for j = 1:length(raster{i})
%         if(ds_idx(j, i))
%             raster{i}{j} = [];
%         end
%     end
end


param_p = 2; % choose which params to use to calculate prefer direction indices 
MAG_all_norm = cell(n_dg, 1);
max_r = cell(n_dg, 1);
norm_max_r = cell(n_dg, 1);

for i = 1:n_dg
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster{i}, DS{i}.angle{param_p});
    MAG_all_norm{i} = normalize_MAG(DS{i});
    rep = datadg{i}.stimulus.repetitions;
end

ll = {'NDF 4', 'NDF 2', 'NDF 1', 'NDF 0'};


%% plot cell summary

for cc = 6:length(ds_id)
    plot_ds_raster(DS, raster, cc, ds_id(cc), ll, 2, 2, 1)
end

%% ffp
n_ffp = 4;

[raster_ff, raster_ff_all] = deal(cell(n_ffp, 1));
for d = 1:n_ffp
    [raster_ff{d}, raster_ff_all{d}] = get_ffp_raster(dataffp{d}, ds_id, 3);
%     for j = 1:length(raster_ff{d})
%         if(ffp_idx(j, d))
%             raster_ff{d}{j} = [];
%             raster_ff_all{d}{j} = [];
%         end
%     end
end

for i = 1:length(ds_id) 
%     if ~isempty(raster_ff{1}{i}) || ~isempty(raster_ff{2}{i})
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 1800 800])
        for d = 1:n_ffp
        subplot(1, n_ffp, d)
        plot_ffp(raster_ff{d}, raster_ff_all{d}, i, 3)
        title([num2str(ds_id(i)) ' ' ll{d}])
        end
        
        print_close(1, [24, 12], num2str(ds_id(i)))
%     end
end

%% frequency analysis
    
duration = 8;
bin_rate = 10000;
hist_spikes = cell(n_dg, 1);
signal_length = duration*bin_rate;                
NFFT = 2^nextpow2(signal_length);
f = bin_rate/2*linspace(0,1,NFFT/2+1);
fft_spikes = cell(n_dg, 1);
[DC, F1, F2] = deal(cell(n_dg, 1));

for i = 1:n_dg
    tp = datadg{i}.stimulus.params.TEMPORAL_PERIOD;
    [DC{i}, F1{i}, F2{i}] = deal(zeros(length(ds_id), length(tp)));
    for rgc = 1:length(ds_id)
        if ~isempty(raster{i}{rgc})
        for time = 1:length(tp)
            spikes = floor(raster_p_sum{i}{rgc}{time}*bin_rate);
            spikes(spikes == 0) = 1;
            tmp_binned_spikes = zeros(1, signal_length);
            tmp_binned_spikes(spikes) = 1;
            hist_spikes{i}{rgc}(time, :) = tmp_binned_spikes;
            
            f1 = 1/tp(time)*60.35; %Hz
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
        F1{i}(rgc,:) = fund_power ./ max(DC_power);
        F2{i}(rgc,:) = sec_power ./ max(DC_power);
        DC{i}(rgc,:) = DC_power ./ max(DC_power);
        
        clear fund_power sec_power DC_power

        end
        
    end
    ratio{i} = F2{i}./F1{i};
end
%% plot 
ratio_nempty = cell(n_dg, 1);
ratio_mean = zeros(n_dg, length(datadg{1}.stimulus.params.TEMPORAL_PERIOD));
ratio_ste = zeros(n_dg, length(datadg{1}.stimulus.params.TEMPORAL_PERIOD));
for i = 1:n_dg
    ratio_nempty{i} = nan2empty(ratio{i});
    ratio_mean(i, :) = mean(ratio_nempty{i});
    ratio_ste(i, :) = std(ratio_nempty{i})/sqrt(size(ratio_nempty{i}, 1));
end

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = {'fast speed'; 'medium speed'; 'slow speed'};
model_series = [ratio_mean(1,1) ratio_mean(2,1) ratio_mean(3,1) ratio_mean(4,1) ; ratio_mean(1,2) ratio_mean(2,2) ratio_mean(3,2) ratio_mean(4,2);ratio_mean(1,3) ratio_mean(2,3) ratio_mean(3,3) ratio_mean(4,3)];   
model_error = [ratio_ste(1,1) ratio_ste(2,1) ratio_ste(3,1) ratio_ste(4,1); ratio_ste(1,2) ratio_ste(2,2) ratio_ste(3,2) ratio_ste(4,2);ratio_ste(1,3) ratio_ste(2,3) ratio_ste(3,3) ratio_ste(4,3)];
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('F2/F1 ratio')
legend('NDF 4','NDF 2', 'NDF 1', 'NDF 0', 'location', 'northeast');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

% title('F2/F1 ratio')

%% classify DSGC into subtypes (directions)
d = 4;
t = 2;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 700 700])
compass(DS{d}.U{t}, DS{d}.V{t})
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DS{d}.U{t}, DS{d}.V{t}, x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = I;
    id_dir{i} = ds_id(idx_dir{i});
end
title('DSGC')