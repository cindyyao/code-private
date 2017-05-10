% addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/DS' cell analysis'/
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);


datadg{1} = load_data('/Analysis/xyao/2016-01-15-0/data000-map/data000-map', opt);
datadg{1}.names.stimulus_path = '/Analysis/xyao/2016-01-15-0/stimuli/s00.mat';
datadg{1} = load_stim_matlab(datadg{1}, 'user_defined_trigger_interval', 10);
datadg{2} = load_data('/Analysis/xyao/2016-01-15-0/data003-map/data003-map', opt);
datadg{2}.names.stimulus_path = '/Analysis/xyao/2016-01-15-0/stimuli/s03.mat';
datadg{2} = load_stim_matlab(datadg{2}, 'user_defined_trigger_interval', 10);
datadg{3} = load_data('/Analysis/xyao/2016-01-15-0/data006-map/data006-map', opt);
datadg{3}.names.stimulus_path = '/Analysis/xyao/2016-01-15-0/stimuli/s06.mat';
datadg{3} = load_stim_matlab(datadg{3}, 'user_defined_trigger_interval', 10);
datadg{4} = load_data('/Analysis/xyao/2016-01-15-0/data009-map/data009-map', opt);
datadg{4}.names.stimulus_path = '/Analysis/xyao/2016-01-15-0/stimuli/s09.mat';
datadg{4} = load_stim_matlab(datadg{4}, 'user_defined_trigger_interval', 10);
datadg{5} = load_data('/Analysis/xyao/2016-01-15-0/data012-sorted/data012-sorted', opt);
datadg{5}.names.stimulus_path = '/Analysis/xyao/2016-01-15-0/stimuli/s12.mat';
datadg{5} = load_stim_matlab(datadg{5}, 'user_defined_trigger_interval', 10);
datadgnclean = load_data('/Analysis/xyao/2016-01-15-0/data012-sorted-noclean/data012', opt);
datadgnclean.names.stimulus_path = '/Analysis/xyao/2016-01-15-0/stimuli/s12.mat';
datadgnclean = load_stim_matlab(datadgnclean, 'user_defined_trigger_interval', 10);

datamb{1} = load_data('/Analysis/xyao/2016-01-15-0/data001-map/data001-map', opt);
datamb{1}.names.stimulus_path = '/Analysis/xyao/2016-01-15-0/stimuli/s01.mat';
datamb{1} = load_stim_matlab(datamb{1});
datamb{1}.stimulus.triggers = [datamb{1}.stimulus.triggers(1:21) datamb{1}.triggers(43) datamb{1}.stimulus.triggers(22:end)];
datamb{2} = load_data('/Analysis/xyao/2016-01-15-0/data004-map/data004-map', opt);
datamb{2}.names.stimulus_path = '/Analysis/xyao/2016-01-15-0/stimuli/s04.mat';
datamb{2} = load_stim_matlab(datamb{2});
datamb{3} = load_data('/Analysis/xyao/2016-01-15-0/data007-map/data007-map', opt);
datamb{3}.names.stimulus_path = '/Analysis/xyao/2016-01-15-0/stimuli/s07.mat';
datamb{3} = load_stim_matlab(datamb{3});
datamb{4} = load_data('/Analysis/xyao/2016-01-15-0/data010-map/data010-map', opt);
datamb{4}.names.stimulus_path = '/Analysis/xyao/2016-01-15-0/stimuli/s10.mat';
datamb{4} = load_stim_matlab(datamb{4});
datamb{5} = load_data('/Analysis/xyao/2016-01-15-0/data013-map/data013-map', opt);
datamb{5}.names.stimulus_path = '/Analysis/xyao/2016-01-15-0/stimuli/s13.mat';
datamb{5} = load_stim_matlab(datamb{5});

dataffp{1} = load_data('/Analysis/xyao/2016-01-15-0/data002-map/data002-map', opt);
dataffp{1}.triggers = dataffp{1}.triggers(2:end);
dataffp{2} = load_data('/Analysis/xyao/2016-01-15-0/data005-map/data005-map', opt);
dataffp{2}.triggers = dataffp{2}.triggers(2:end);
dataffp{3} = load_data('/Analysis/xyao/2016-01-15-0/data008-map/data008-map', opt);
dataffp{3}.triggers = dataffp{3}.triggers(2:end);
dataffp{4} = load_data('/Analysis/xyao/2016-01-15-0/data011-map/data011-map', opt);
dataffp{4}.triggers = dataffp{4}.triggers(2:end);
dataffp{5} = load_data('/Analysis/xyao/2016-01-15-0/data014-map/data014-map', opt);
dataffp{5}.triggers = dataffp{5}.triggers(2:end);


%%
n = 5;
i = 5;
[NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg{i},datadg{i}.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg{i});
params_idx = [3 4]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg{i}, ds_struct, params_idx);
title('Normalized Vector Sum Amplitude')
xlabel('High speed - 480 micron/s')
ylabel('Low speed - 240 micron/s')

%% 
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, StimComb] = get_spikescellstim(datadg{i},ds_id,0);
    DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg{i}));
    raster_dg{i} = get_ds_raster(datadg{i}, ds_id);
%     for j = 1:length(raster_dg{i})
%         if(dg_idx(j, i))
%             raster_dg{i}{j} = [];
%         end
%     end
end

delta_p = 4; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

for i = 1:n
    if ismember(i, [2 5])
        [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
    else
        [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p-3});
    end
    MAG_all_norm_dg{i} = normalize_MAG(DG{i});
    rep = datadg{i}.stimulus.repetitions;
end

ll = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};

%% plot cell summary
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [4 5]);

for cc = 1:length(ds_id)
    plot_ds_raster(DG_cut, raster_dg_cut, cc, ds_id(cc), ll, 2, 3, 1)
end

for cc = 12:12%length(ds_id)
    plot_ds_raster_one(DG_cut, raster_dg_cut, cc, ds_id(cc), 0)
end

%% plot single cell tuning curve

% use unnormalized vector sum as response
figure
v = datadg{2}.stimulus.params.SPATIAL_PERIOD./datadg{2}.stimulus.params.TEMPORAL_PERIOD*4;
subplot(1, 2, 1)
semilogx(v, exciseColumn(MAG_all_norm_dg{2}), 'b')
xlabel('micron/second')
ylabel('Response')
title(ll{2})
xlim([v(end) v(1)])

subplot(1, 2, 2)
semilogx(v, exciseColumn(MAG_all_norm_dg{5}), 'b')
xlabel('micron/second')
ylabel('Response')
title(ll{5})
xlim([v(end) v(1)])

%% classification based on speed tunning
%pca
mag_pca = MAG_all_norm_dg{5};
mag_pca = mag_pca';
[id_sub, idx_sub] = deal(cell(2, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(mag_pca);
pc1 = 1; pc2 = 2;
plot(scores(:, pc1), scores(:, pc2), 'o')
hold on
% for i = 1:2
%     [x, y] = ginput;
%     plot(x, y)
%     IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
%     [~, idx_sub{i}] = find(IN' == 1);
%     id_sub{i} = ds_id(idx_sub{i});
% end
% xlabel('1st Principal Component')
% ylabel('3rd Principal Component')
% title('NDF 0')

[x, y] = ginput;
plot(x, y)
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
title('NDF 0')

IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
[~, I] = find(IN' == 1);
id_init = ds_id(I);

[C ia ib] = intersect(id_init, ds_id);
vc = ones(length(ds_id),1);
vc(ib) = 2; %initializing ds cells to cluster 2, everything else cluster 1
close all;
X = [];
N = [];
p = [];
[idx N p] = clustering_analysis_plots(scores(:, [pc1 pc2]), 0,1, 2, 0, 1, 0, 0, 0,0, vc);
idx_sub{1} = find(idx == 1);
idx_sub{2} = find(idx == 2);
id_sub{1} = ds_id(idx_sub{1});
id_sub{2} = ds_id(idx_sub{2});

% % several cells are missing during mapping between NDF0 grating and NDF3
% % grating, simply put them into ON-OFF DSGC class here
% idx_all = 1:length(ds_id);
% idx_sub{2} = setdiff(idx_all, idx_sub{1});
% id_sub{2} = ds_id(idx_sub{2});

figure
plot(scores(idx_sub{1}, pc1), scores(idx_sub{1}, pc2), 'ro', scores(idx_sub{2}, pc1), scores(idx_sub{2}, pc2), 'bo')
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
title('NDF 0')

% 
% figure
% v = 4*datadg{2}.stimulus.params.SPATIAL_PERIOD./datadg{2}.stimulus.params.TEMPORAL_PERIOD;
% subplot(1, 2, 1)
% semilogx(v, exciseColumn(MAG_all_norm_dg{2}(:, idx_sub{1})), 'r')
% hold on
% semilogx(v, exciseColumn(MAG_all_norm_dg{2}(:, idx_sub{2})), 'b')
% xlabel('micron/second')
% ylabel('Response')
% % title(ll{2})
% xlim([v(end) v(1)])
% 
% subplot(1, 2, 2)
% figure
% semilogx(v, exciseColumn(MAG_all_norm_dg{5}(:, idx_sub{1})), 'r')
% hold on
% semilogx(v, exciseColumn(MAG_all_norm_dg{5}(:, idx_sub{2})), 'b')
% xlabel('micron/second')
% ylabel('Response')
% % title(ll{5})
% xlim([v(end) v(1)])

t = 4;
figure
compass(DG{5}.U{t}(idx_sub{1}), DG{5}.V{t}(idx_sub{1}), 'r')
hold on
compass(DG{5}.U{t}(idx_sub{2}), DG{5}.V{t}(idx_sub{2}), 'b')

%% plot average tunning curve
color = 'brk';
figure
idx_temp = [2 5];
for i = 1:2
    v = 4*datadg{idx_temp(i)}.stimulus.params.SPATIAL_PERIOD./datadg{idx_temp(i)}.stimulus.params.TEMPORAL_PERIOD;
    subplot(1, 2, i)
    for ct = 1:2
        mag_temp = exciseColumn(MAG_all_norm_dg{idx_temp(i)}(:, idx_sub{ct}));
        tuning_avg{i}(:, ct) = mean(mag_temp, 2);
        tuning_ste{i}(:, ct) = std(mag_temp, [], 2)/sqrt(size(mag_temp, 2));
        errorbar(v, tuning_avg{i}(:, ct), tuning_ste{i}(:, ct), color(ct), 'LineWidth', 2)
        hold on
    end
    set(gca, 'XScale', 'log')
    title(ll{idx_temp(i)})
    xlim([min(v) max(v)])
    xlabel('speed')
    ylabel('response')
end
legend('on-off DSGC', 'on DSGC', 'location', 'southeast')


%% compare across light level
CT = {'ON DSGC', 'ON-OFF DSGC'};
figure
v = 4*datadg{2}.stimulus.params.SPATIAL_PERIOD./datadg{2}.stimulus.params.TEMPORAL_PERIOD;
for j = 1:2
    subplot(1, 2, j)
    errorbar(v, tuning_avg{1}(:, j), tuning_ste{1}(:, j), 'b', 'LineWidth', 2)
    hold on
    errorbar(v, tuning_avg{2}(:, j), tuning_ste{2}(:, j), 'r', 'LineWidth', 2)
    title(CT{j})
    set(gca, 'XScale', 'log')
    xlim([min(v) max(v)])
    xlabel('speed')
    ylabel('response')

end
legend('NDF 3', 'NDF 0', 'location', 'northeast')

%% classify DSGC into subtypes (directions)
d = 5;
t = 3;
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

% on DS
d = 5;
t = 3;
h = figure;
dirn = 2;
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

%% get superior ds cell ei location
elec_oo = get_cell_location_ei(datadg{5}, id_dir{1});
elec_on = get_cell_location_ei(datadg{5}, id_dir_on{1});

%% DS tuning curves (drifting grating)
% all ds cells

% t = 2;
dirn = 4;
D = 5;
T = 1;

for d = 1:n
    p_direction = DG_cut{D}.angle{T}';
    xx = 0:pi/4:7*pi/4;
    xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
    xx(xx>pi) = xx(xx>pi)-2*pi;
    xx(xx<-pi) = xx(xx<-pi)+2*pi;


    subplot(2, 3, d)
    for i = 1:dirn
        for cc = 1:length(ds_id)
%             if ~dg_idx(cc, d)
            [xsort, seq] = sort(xx(cc, :));
            y_temp = DG_cut{d}.rho{T}(cc, :);
            plot(xsort, y_temp(seq), 'b')
            hold on
%             end
        end
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    title(ll{d})
    xlim([-pi pi])
end

%subtypes
clear rho_dg_mean rho_dg_ste dsi_dg_mean dsi_dg_ste
for d = 1:5
    subplot(2, 3, d)
    for i = 1:dirn
        rho_dg{d}{i} = [];
        dsi_dg{d}{i} = [];
        for cc = 1:length(idx_dir{i})
%             if ~dg_idx(idx_dir{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx_dir{i}(cc), :))>0
            [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
            y_temp = DG_cut{d}.rho{T}(idx_dir{i}(cc), :);
            plot(xsort, y_temp(seq), color(i))
            ylim([0 1])
%             pause
            hold on
            rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
            dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx_dir{i}(cc))];
%             end
        end
        rho_dg_mean{d}(i, :) = mean(rho_dg{d}{i});
        rho_dg_ste{d}(i, :) = std(rho_dg{d}{i})/sqrt(size(rho_dg{d}{i}, 1));
        dsi_dg_mean{d}(i) = mean(dsi_dg{d}{i});
        dsi_dg_ste{d}(i) = std(dsi_dg{d}{i})/sqrt(length(dsi_dg{d}{i}));
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    title(ll{d})
    xlim([-pi pi])
end
dsi_dg_mean = cell2mat(dsi_dg_mean');
dsi_dg_ste = cell2mat(dsi_dg_ste');
% plot average (cell type)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for d = 1:5
    subplot(2, 3, d)
    for i = 1:dirn
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(i));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
    title(ll{d});
end
legend(ct)

% plot average (light level)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:dirn
    subplot(2, 2, i)
    for d = 1:5
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(d));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
    title(ct{i})
end
legend(ll)
% DSI
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = ct;
model_series = [dsi_dg_mean(1,1) dsi_dg_mean(2,1) dsi_dg_mean(3,1) dsi_dg_mean(4,1) dsi_dg_mean(5,1); dsi_dg_mean(1,2) dsi_dg_mean(2,2) dsi_dg_mean(3,2) dsi_dg_mean(4,2) dsi_dg_mean(5,2);dsi_dg_mean(1,3) dsi_dg_mean(2,3) dsi_dg_mean(3,3) dsi_dg_mean(4,3) dsi_dg_mean(5,3); dsi_dg_mean(1,4) dsi_dg_mean(2,4) dsi_dg_mean(3,4) dsi_dg_mean(4,4) dsi_dg_mean(5,4)];   
model_error = [dsi_dg_ste(1,1) dsi_dg_ste(2,1) dsi_dg_ste(3,1) dsi_dg_ste(4,1) dsi_dg_ste(5,1); dsi_dg_ste(1,2) dsi_dg_ste(2,2) dsi_dg_ste(3,2) dsi_dg_ste(4,2) dsi_dg_ste(5,2);dsi_dg_ste(1,3) dsi_dg_ste(2,3) dsi_dg_ste(3,3) dsi_dg_ste(4,3) dsi_dg_ste(5,3); dsi_dg_ste(1,4) dsi_dg_ste(2,4) dsi_dg_ste(3,4) dsi_dg_ste(4,4) dsi_dg_ste(5,4)];
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('DSI')
legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

%% frequency analysis
clear ratio_avg ratio_ste ratio_dir_mean ratio_dir_ste ratio
duration = 8;
bin_rate = 10000;
hist_spikes = cell(n, 1);
signal_length = duration*bin_rate;                
NFFT = 2^nextpow2(signal_length);
f = bin_rate/2*linspace(0,1,NFFT/2+1);
fft_spikes = cell(n, 1);
[DC, F1, F2] = deal(cell(n, 1));
for i = [2 5]
%     tp = datadg{1}.stimulus.params.TEMPORAL_PERIOD;
    tp = datadg{2}.stimulus.params.TEMPORAL_PERIOD;
    [DC{i}, F1{i}, F2{i}] = deal(zeros(length(ds_id), length(tp)));
    for rgc = 1:length(ds_id)
%         if ~isempty(raster_dg_cut{i}{rgc}) && ~dg_idx(rgc, i)
%         if ~isempty(raster_dg{i}{rgc}) && ~dg_idx(rgc, i)
        for time = 1:length(tp)
%             spikes = floor(raster_p_sum_cut{i}{rgc}{time}*bin_rate);
            spikes = floor(raster_p_sum{i}{rgc}{time}*bin_rate);
            tmp_binned_spikes = zeros(1, signal_length);
            tmp_binned_spikes(spikes) = 1;
            hist_spikes{i}{rgc}(time, :) = tmp_binned_spikes;
            
            f1 = 1/tp(time); %Hz
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

%         end
        
    end
    ratio{i} = F2{i}./F1{i};
    for ct = 1:4
        ratio_dir{ct}{i} = ratio{i}(idx_dir{ct}, :);
        ratio_dir{ct}{i} = exciseRows_empty(ratio_dir{ct}{i});
        ratio_dir_mean(i, ct, :) = mean(ratio_dir{ct}{i});
        ratio_dir_ste(i, ct, :) = std(ratio_dir{ct}{i})/sqrt(size(ratio_dir{ct}{i}, 1));
    end
    for ct = 1:3
        ratio_dir_on{ct}{i} = ratio{i}(idx_dir_on_{ct}, :);
        ratio_dir_on{ct}{i} = exciseRows_empty(ratio_dir_on{ct}{i});
        ratio_dir_mean_on(i, ct, :) = mean(ratio_dir_on{ct}{i}, 1);
        ratio_dir_ste_on(i, ct, :) = std(ratio_dir_on{ct}{i}, [], 1)/sqrt(size(ratio_dir_on{ct}{i}, 1));
    end
    for ct = 1:2
        ratio_onoff{ct}{i} = ratio{i}(idx_sub_{ct}, :);
        ratio_onoff{ct}{i} = exciseRows_empty(ratio_onoff{ct}{i});
        ratio_onoff_mean(i, ct, :) = mean(ratio_onoff{ct}{i}, 1);
        ratio_onoff_ste(i, ct, :) = std(ratio_onoff{ct}{i}, [], 1)/sqrt(size(ratio_onoff{ct}{i}, 1));
    end
    ratio{i} = exciseRows_empty(ratio{i});
    ratio_mean(i, :) = mean(ratio{i});
    ratio_ste(i, :) = std(ratio{i})/sqrt(size(ratio{i}, 1));
end
for i = 1:n
    tp = datadg{1}.stimulus.params.TEMPORAL_PERIOD;
    [DC{i}, F1{i}, F2{i}] = deal(zeros(length(ds_id), length(tp)));
    for rgc = 1:length(ds_id)
        if ~isempty(raster_dg_cut{i}{rgc}) && ~dg_idx(rgc, i)
        for time = 1:length(tp)
            spikes = floor(raster_p_sum_cut{i}{rgc}{time}*bin_rate);
            tmp_binned_spikes = zeros(1, signal_length);
            tmp_binned_spikes(spikes) = 1;
            hist_spikes{i}{rgc}(time, :) = tmp_binned_spikes;
            
            f1 = 1/tp(time); %Hz
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
    for ct = 1:4
        ratio_dir{ct}{i} = ratio{i}(idx_dir{ct}, :);
        ratio_dir{ct}{i} = exciseRows_empty(ratio_dir{ct}{i});
        ratio_dir_mean(i, ct, :) = mean(ratio_dir{ct}{i});
        ratio_dir_ste(i, ct, :) = std(ratio_dir{ct}{i})/sqrt(size(ratio_dir{ct}{i}, 1));
    end
    ratio_oo{i} = ratio{i}(idx_sub{2}, :);
    ratio_oo{i} = exciseRows_empty(ratio_oo{i});
    ratio_oo_mean(i, :) = mean(ratio_oo{i});
    ratio_oo_ste(i, :) = std(ratio_oo{i})/sqrt(size(ratio_oo{i}, 1));
end

% plot average f2/f1
ct = {'posterior', 'inferior', 'anterior', 'superior'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = {'high speed'; 'low speed'};
model_series = ratio_oo_mean';
model_error = ratio_oo_ste';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('F2/F1')
legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

% plot F2/F1 tuning curve
v3 = [4800 2400 960 480 240 120 80 60 40];
figure
for i = 1:2
subplot(1, 2, i)
errorbar(v3, ratio_onoff_mean(2, i, :), ratio_onoff_ste(2, i, :), 'b')
hold on
errorbar(v3, ratio_onoff_mean(5, i, :), ratio_onoff_ste(5, i, :), 'r')
set(gca, 'XScale', 'log')
end

% F2/F1 tuning for individual direction
celltype = {'posterior', 'inferior', 'anterior', 'superior'};
figure
for ct = 1:4
    subplot(2, 2, ct)
    errorbar(v3, ratio_dir_mean(2, ct, :), ratio_dir_ste(2, ct, :), 'b')
    hold on
    errorbar(v3, ratio_dir_mean(5, ct, :), ratio_dir_ste(5, ct, :), 'r')
    set(gca, 'XScale', 'log')
    title(celltype{ct})
end

%% plot cell type specific f2/f1 
for tp = 1:2
    figure
    set(gcf, 'DefaultLineLineWidth', 1.5)
    xtick = ct;
    model_series = ratio_dir_mean(:,:,tp)';
    model_error = ratio_dir_ste(:,:,tp)';
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('F2/F1')
    legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
    if tp == 1
        title('high speed')
    else
        title('low speed')
    end
    hold on;

    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 

    groupwidth = min(0.8, numbars/(numbars+1.5));

    for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end
end


%%
duration = 573;
[raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, ~, StimComb] = get_spikescellstim_mb(datamb{i},ds_id,duration,1);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb{i} = get_mb_raster(datamb{i}, ds_id, duration);
%     for j = 1:length(raster_mb{i})
%         if(mb_idx(j, i))
%             raster_mb{i}{j} = [];
%         end
%     end
    trial_dur{i} = get_mb_trial_dur(datamb{i});
end

delta_p = 2; % choose which params to use to calculate prefer direction indices 
ll_p = 5;
MAG_all_norm_mb = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

for i = 1:n
    [raster_p_sum_mb{i}, p_idx{i}] = get_pdirection_raster(raster_mb{i}, MB{ll_p}.angle{delta_p});
    MAG_all_norm_mb{i} = normalize_MAG(MB{i});
    rep = datamb{i}.stimulus.repetitions;
end

ll = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};
ct = {'superior', 'posterior', 'inferior', 'anterior'};
close all

%% plot cell summary
Title = {'BarWidth 120, white', 'BarWidth 240, white', 'BarWidth 360, white'; 'BarWidth 120, black', 'BarWidth 240, black', 'BarWidth 360, black'};
for cc = 1:length(ds_id)
    plot_mb_raster_one(MB, raster_mb, trial_dur, cc, ds_id(cc), num2str(ds_id(cc)), 1, 1, 1)
end

%% full field pulses

n_ffp = 5;

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


%%
% load('electrode_position.mat')
% elec_n_cl = unique(ceil(ds_id/15));
% elec_n_ncl = unique(ceil(ds_id_ncl/15));
% figure
% plot_electrodes(position, 'rotation', 90, 'Alpha', false, 'label', true);
% hold on
% plot_electrodes(position(elec_n_ncl, :), 'rotation', 90, 'Alpha', false, 'elec_colors', repmat([1 0 0], length(elec_n_ncl), 1));
% % hold on
% plot_electrodes(position(elec_n_cl, :), 'rotation', 90, 'Alpha', false, 'elec_colors', repmat([0 0 1], length(elec_n_cl), 1));
% xlim([-450 450])
% ylim([-450 450])
