%% load data
addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/Classification/
addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/DS' cell analysis'/

opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load data
datadg{1} = load_data('/Analysis/xyao/2014-11-25-0/data000-map/data000-map', opt);
datadg{1}.names.stimulus_path = '/Analysis/xyao/2014-11-25-0/stimuli/s00';
datadg{1} = load_stim(datadg{1}, 'user_defined_trigger_interval', 10);

datadg{2} = load_data('/Analysis/xyao/2014-11-25-0/data003-map/data003-map', opt);
datadg{2}.names.stimulus_path = '/Analysis/xyao/2014-11-25-0/stimuli/s03';
datadg{2} = load_stim(datadg{2}, 'user_defined_trigger_interval', 10);

datadg{3} = load_data('/Analysis/xyao/2014-11-25-0/data006-map/data006-map', opt);
datadg{3}.names.stimulus_path = '/Analysis/xyao/2014-11-25-0/stimuli/s06';
datadg{3} = load_stim(datadg{3}, 'user_defined_trigger_interval', 10);

datadg{4} = load_data('/Analysis/xyao/2014-11-25-0/data008/data008', opt);
datadg{4}.names.stimulus_path = '/Analysis/xyao/2014-11-25-0/stimuli/s08';
datadg{4} = load_stim(datadg{4}, 'user_defined_trigger_interval', 10);

dataffp{1} = load_data('/Analysis/xyao/2014-11-25-0/data001-map/data001-map', opt);
dataffp{2} = load_data('/Analysis/xyao/2014-11-25-0/data004-map/data004-map', opt);
dataffp{3} = load_data('/Analysis/xyao/2014-11-25-0/data007-map/data007-map', opt);
dataffp{4} = load_data('/Analysis/xyao/2014-11-25-0/data009-map/data009-map', opt);


%%
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

ll = {'control', '300nM Strychnine', 'wash out', 'NDF 0'};


%% plot cell summary

for cc = 1:length(ds_id)
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


%% classify DSGC into subtypes (directions)
d = 1;
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
legend('ctrl','300nM strychnine', 'wash out', 'NDF 0', 'location', 'northeast');
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

%% DSI
DSI = cell(4, 1);
ds_idx_ = logical(1-ds_idx);
dsi_mean = zeros(n_dg, length(DS{1}.dsindex));
dsi_ste = zeros(n_dg, length(DS{1}.dsindex));

for i = 1:n_dg
    dsi_temp = cell2mat(DS{i}.dsindex');
    dsi_temp = dsi_temp(ds_idx_(:, i), :);
    dsi_mean(i, :) = mean(dsi_temp);
    dsi_ste(i, :) = std(dsi_temp)/sqrt(size(dsi_temp, 1));
end

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = {'fast speed'; 'medium speed'; 'slow speed'};
model_series = [dsi_mean(1,1) dsi_mean(2,1) dsi_mean(3,1) dsi_mean(4,1); dsi_mean(1,2) dsi_mean(2,2) dsi_mean(3,2) dsi_mean(4,2);dsi_mean(1,3) dsi_mean(2,3) dsi_mean(3,3) dsi_mean(4,3)];   
model_error = [dsi_ste(1,1) dsi_ste(2,1) dsi_ste(3,1) dsi_ste(4,1); dsi_ste(1,2) dsi_ste(2,2) dsi_ste(3,2) dsi_ste(4,2);dsi_ste(1,3) dsi_ste(2,3) dsi_ste(3,3) dsi_ste(4,3)];
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylim([-0.4 1.5])
ylabel('DSI')
legend('ctrl','300nM strychnine', 'wash out', 'NDF 0', 'location', 'northeast');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

title('DSI')


%% ffp psth
ctype = {'Posterior'; 'Superior'; 'Anterior'; 'Inferior'};
XX = 0.05:0.05:11.95;
figure
for i = 1:4
    idx_temp = idx_dir{i};
    for cd = 1:n_ffp
        psth{cd}{i} = [];
        subplot(4, n_ffp, n_ffp*(i-1)+cd)
        a = 0;
        for cc = 1:length(idx_temp)
            if ~isempty(raster_ff_all{cd}{idx_temp(cc)})
                h = hist(raster_ff_all{cd}{idx_temp(cc)}, XX);
                plot(XX, h/30, color(cd))
                hold on
                a = a+1;
                psth{cd}{i}(a, :) = h/30;
            end
        end
        if i == 1
            title(ll{cd})
        end
        if i == 4
            xlabel('time/sec')
        end
        if cd == 1
            ylabel(ctype{i})
        end
        xlim([0 12])    
    end
end


%% ffp of non-DS
n_ffp = 4;

[raster_ff_allcell, raster_ff_all_allcell] = deal(cell(n_ffp, 1));
for d = 1:n_ffp
    [raster_ff_allcell{d}, raster_ff_all_allcell{d}] = get_ffp_raster(dataffp{d}, datadg{4}.cell_ids, 3);
%     for j = 1:length(raster_ff_allcell{d})
%         if(ffp_idx(j, d))
%             raster_ff_allcell{d}{j} = [];
%             raster_ff_all_allcell{d}{j} = [];
%         end
%     end
end

for i = 1:length(datadg{4}.cell_ids) 
%     if ~isempty(raster_ff{1}{i}) || ~isempty(raster_ff{2}{i})
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 1800 800])
        for d = 1:n_ffp
        subplot(1, n_ffp, d)
        plot_ffp(raster_ff_allcell{d}, raster_ff_all_allcell{d}, i, 3)
        title([num2str(datadg{4}.cell_ids(i)) ' ' ll{d}])
        end
        
        print_close(1, [24, 12], num2str(datadg{4}.cell_ids(i)))
%     end
end

%% dg raster of all cells
n_dg = 4;
[raster, raster_p_sum, p_idx] = deal(cell(n_dg, 1));
for i = 1:n_dg   
    [NumSpikesCell, StimComb] = get_spikescellstim(datadg{i},datadg{4}.cell_ids,0);
    DS_all{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster_allcell{i} = get_ds_raster(datadg{i}, datadg{4}.cell_ids);
end

ll = {'control', '300nM Strychnine', 'wash out', 'NDF 0'};


%% plot cell summary

cellid = 7655;
cc = get_cell_indices(datadg{4}, cellid);
plot_ds_raster(DS_all, raster_allcell, cc, datadg{4}.cell_ids(cc), ll, 2, 2, 0)

%% individual cell ffp psth
id_ffp_eg{1} = [483 3226 3302 3305 3378 3481 3886 4202 4936 6875]; %on transient
id_ffp_eg{2} = [196 241 496 1561 1711 2326 2551 2732 3017 4623 6466]; %on sustained
id_ffp_eg{3} = [181 1637 2359 2596 2791 3228 3677 4069 4263 4396 4923]; %off transient
id_ffp_eg{4} = [3227 3271 3286 3466 4771]; %off transient small
id_ffp_eg{5} = [2371 2536 2747 3691 3693 4096 4261 4306 4531 4951 5116 5701 5896 5956 6137 6436 6856 7111]; %off sustained

XX = 0.025:0.05:11.975;
celltype = {'on transient', 'on sustained', 'off transient', 'off transient small', 'off sustained'};
for ct = 1:5
    idx_ffp_eg{ct} = get_cell_indices(datadg{4}, id_ffp_eg{ct});
    for i = 1:4
        for cc = 1:length(idx_ffp_eg{ct})
            psth{ct}(:, i, cc) = hist(raster_ff_all_allcell{i}{idx_ffp_eg{ct}(cc)}, XX);
        end
    end
end


figure
for j = 1:4
    for i = 1:5
        subplot(4, 5, (j-1)*5+i)
        plot(XX, squeeze(psth{i}(:, j, :)), 'b')
        xlim([0 12])
        if j == 1
          title(celltype{i})
        end
        if i == 1
            ylabel(ll{j})
        end
    end
end

%% off sustained response change
ct = 5;
for cd = 1:4
    for cc = 1:size(psth{ct}, 3)
        r_offs(cd, cc, 1) = 0.8*max(psth{ct}(91:120, cd, cc)) - mean(psth{ct}(61:90, cd, cc));
        r_offs(cd, cc, 2) = 0.8*max(psth{ct}(151:180, cd, cc)) - mean(psth{ct}(121:150, cd, cc));
    end
    r_offs_mean = squeeze(mean(r_offs, 2));
    r_offs_ste = squeeze(std(r_offs, 0, 2)/sqrt(size(r_offs, 2)));
    
end
    
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = {'w-g'; 'g-b'};
model_series = [r_offs_mean(1,1) r_offs_mean(2,1) r_offs_mean(3,1) r_offs_mean(4,1); r_offs_mean(1,2) r_offs_mean(2,2) r_offs_mean(3,2) r_offs_mean(4,2)];   
model_error = [r_offs_ste(1,1) r_offs_ste(2,1) r_offs_ste(3,1) r_offs_ste(4,1); r_offs_ste(1,2) r_offs_ste(2,2) r_offs_ste(3,2) r_offs_ste(4,2)];
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
% ylim([-0.4 1.5])
ylabel('off response')
legend('ctrl','300nM strychnine', 'wash out', 'NDF 0', 'location', 'northwest');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

title('off sustained')

%% off transient response change
ct = 3;
for cd = 1:4
    for cc = 1:size(psth{ct}, 3)
        r_offt(cd, cc, 1) = 0.8*max(psth{ct}(91:120, cd, cc)) - mean(psth{ct}(61:90, cd, cc));
        r_offt(cd, cc, 2) = 0.8*max(psth{ct}(151:180, cd, cc)) - mean(psth{ct}(121:150, cd, cc));
        r_offt(cd, cc, 3) = 0.8*max(psth{ct}(31:60, cd, cc)) - mean(psth{ct}(1:30, cd, cc));
        r_offt(cd, cc, 4) = 0.8*max(psth{ct}(211:240, cd, cc)) - mean(psth{ct}(181:210, cd, cc));
    end
    r_offs_mean = squeeze(mean(r_offt, 2));
    r_offs_ste = squeeze(std(r_offt, 0, 2)/sqrt(size(r_offt, 2)));
    
end
    
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = {'w-g(off)'; 'g-b(off)'; 'g-w(on)'; 'b-g(on)'};
model_series = [r_offs_mean(1,1) r_offs_mean(2,1) r_offs_mean(3,1) r_offs_mean(4,1); r_offs_mean(1,2) r_offs_mean(2,2) r_offs_mean(3,2) r_offs_mean(4,2); r_offs_mean(1,3) r_offs_mean(2,3) r_offs_mean(3,3) r_offs_mean(4,3); r_offs_mean(1,4) r_offs_mean(2,4) r_offs_mean(3,4) r_offs_mean(4,4);];   
model_error = [r_offs_ste(1,1) r_offs_ste(2,1) r_offs_ste(3,1) r_offs_ste(4,1); r_offs_ste(1,2) r_offs_ste(2,2) r_offs_ste(3,2) r_offs_ste(4,2); r_offs_ste(1,3) r_offs_ste(2,3) r_offs_ste(3,3) r_offs_ste(4,3); r_offs_ste(1,4) r_offs_ste(2,4) r_offs_ste(3,4) r_offs_ste(4,4);];
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
% ylim([-0.4 1.5])
ylabel('off response')
legend('ctrl','300nM strychnine', 'wash out', 'NDF 0', 'location', 'northwest');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

% title('off transient')

%% off transient response change
ct = 4;
for cd = 1:4
    for cc = 1:size(psth{ct}, 3)
        r_offts(cd, cc, 1) = 0.8*max(psth{ct}(91:120, cd, cc)) - mean(psth{ct}(61:90, cd, cc));
        r_offts(cd, cc, 2) = 0.8*max(psth{ct}(151:180, cd, cc)) - mean(psth{ct}(121:150, cd, cc));
        r_offts(cd, cc, 3) = 0.8*max(psth{ct}(31:60, cd, cc)) - mean(psth{ct}(1:30, cd, cc));
        r_offts(cd, cc, 4) = 0.8*max(psth{ct}(211:240, cd, cc)) - mean(psth{ct}(181:210, cd, cc));
    end
    r_offs_mean = squeeze(mean(r_offts, 2));
    r_offs_ste = squeeze(std(r_offts, 0, 2)/sqrt(size(r_offts, 2)));
    
end
    
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = {'w-g(off)'; 'g-b(off)'; 'g-w(on)'; 'b-g(on)'};
model_series = [r_offs_mean(1,1) r_offs_mean(2,1) r_offs_mean(3,1) r_offs_mean(4,1); r_offs_mean(1,2) r_offs_mean(2,2) r_offs_mean(3,2) r_offs_mean(4,2); r_offs_mean(1,3) r_offs_mean(2,3) r_offs_mean(3,3) r_offs_mean(4,3); r_offs_mean(1,4) r_offs_mean(2,4) r_offs_mean(3,4) r_offs_mean(4,4);];   
model_error = [r_offs_ste(1,1) r_offs_ste(2,1) r_offs_ste(3,1) r_offs_ste(4,1); r_offs_ste(1,2) r_offs_ste(2,2) r_offs_ste(3,2) r_offs_ste(4,2); r_offs_ste(1,3) r_offs_ste(2,3) r_offs_ste(3,3) r_offs_ste(4,3); r_offs_ste(1,4) r_offs_ste(2,4) r_offs_ste(3,4) r_offs_ste(4,4);];
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
% ylim([-0.4 1.5])
ylabel('response')
legend('ctrl','300nM strychnine', 'wash out', 'NDF 0', 'location', 'northwest');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

title('off transient small')



%% align ffp raster with triggers

for i = 1:n_ffp
    trigger_diff_temp = [diff(dataffp{i}.triggers); 0];
    trigger_diff_temp = reshape(trigger_diff_temp, 4, 25);
    trigger_diff_temp = [1.5*ones(1, 25); trigger_diff_temp(1:3, :)];
    trigger_diff_temp(2, :) = trigger_diff_temp(1, :) + trigger_diff_temp(2, :);
    trigger_diff_temp(3, :) = trigger_diff_temp(2, :) + trigger_diff_temp(3, :);
    trigger_diff_temp(4, :) = trigger_diff_temp(3, :) + trigger_diff_temp(4, :);
    trigger_diff{i} = trigger_diff_temp;
    clear trigger_diff_temp
end
    
% plot
for i = 146:146 %length(ds_id) 
%     if ~isempty(raster_ff{1}{i}) || ~isempty(raster_ff{2}{i})
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 1800 800])
        for d = 1:n_ffp
            if ~isempty(raster_ff_allcell{d}{i})
            subplot(1, n_ffp, d)
            plot_ffp(raster_ff_allcell{d}, raster_ff_all_allcell{d}, i, 3)
            hold on
            for j = 1:25
                SpikeTime = trigger_diff{d}(:, j);
                SpikeTime = SpikeTime';
                X = [SpikeTime; SpikeTime];
                Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
                line(X, Y, 'color', 'r');
                xlim([0 12])
            end
            title([num2str(datadg{4}.cell_ids(i)) ' ' ll{d}])
            end
        
%         print_close(1, [24, 12], num2str(ds_id(i)))
        end
end
         