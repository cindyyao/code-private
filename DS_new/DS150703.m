% addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/DS' cell analysis'/
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);


datadg{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data000-map/data000-map', opt);
datadg{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s00.mat';
datadg{1} = load_stim_matlab(datadg{1}, 'user_defined_trigger_interval', 10);
datadg{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data003-map/data003-map', opt);
% datadg{2} = load_data('/Analysis/xyao/2015-07-03-0/data003/data003', opt);
datadg{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s03.mat';
datadg{2} = load_stim_matlab(datadg{2}, 'user_defined_trigger_interval', 10);
datadg{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data006-map/data006-map', opt);
datadg{3}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s06.mat';
datadg{3} = load_stim_matlab(datadg{3}, 'user_defined_trigger_interval', 10);
datadg{4} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data009-map/data009-map', opt);
datadg{4}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s09.mat';
datadg{4} = load_stim_matlab(datadg{4}, 'user_defined_trigger_interval', 10);
datadg{5} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data012-map/data012-map', opt);
datadg{5}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s12.mat';
datadg{5} = load_stim_matlab(datadg{5}, 'user_defined_trigger_interval', 10);

datamb{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data001-map/data001-map', opt);
datamb{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s01.mat';
datamb{1} = load_stim_matlab(datamb{1});
datamb{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data004-map/data004-map', opt);
datamb{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s04.mat';
datamb{2} = load_stim_matlab(datamb{2});
datamb{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data007-map/data007-map', opt);
datamb{3}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s07.mat';
datamb{3} = load_stim_matlab(datamb{3});
datamb{4} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data010-map/data010-map', opt);
datamb{4}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s10.mat';
datamb{4} = load_stim_matlab(datamb{4});
datamb{5} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data014-map/data014-map', opt);
datamb{5}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/stimuli/s14.mat';
datamb{5} = load_stim_matlab(datamb{5});

dataffp{1} = load_data('/Volumes/lab/analysis/2015-07-03-0/data002-map/data002-map', opt);
dataffp{1}.triggers = dataffp{1}.triggers(2:end);
dataffp{2} = load_data('/Volumes/lab/analysis/2015-07-03-0/data005-map/data005-map', opt);
dataffp{2}.triggers = dataffp{2}.triggers(2:end);
dataffp{3} = load_data('/Volumes/lab/analysis/2015-07-03-0/data008-map/data008-map', opt);
dataffp{3}.triggers = dataffp{3}.triggers(2:end);
dataffp{4} = load_data('/Volumes/lab/analysis/2015-07-03-0/data011-map/data011-map', opt);
dataffp{4}.triggers = dataffp{4}.triggers(2:end);
dataffp{5} = load_data('/Volumes/lab/analysis/2015-07-03-0/data015-map/data015-map', opt);
dataffp{5}.triggers = dataffp{5}.triggers(2:end);

datawn = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-07-03-0/data016-map/data016-map', opt);
datawn = load_sta(datawn);
datawn = get_sta_summaries(datawn, 'all');
% datarun = load_data('/Analysis/xyao/2015-07-03-0/data016/data016', opt);
% datawn = load_data('/Volumes/lab/analysis/2015-07-03-0/data016-map/data016-map', opt);
datawn = load_ei(datawn, id_dir{4});
%%
n = 5;
i = 5;
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg{i},datadg{i}.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg{i});

params_idx = [4 5]; % which parameters to use for classification
[ds_id, nonds_id] = classify_ds(datadg{i}, ds_struct, params_idx);
title('Normalized Vector Sum Amplitude')
xlabel('High speed - 480 micron/s')
ylabel('Low speed - 240 micron/s')
%% drifting grating

% spike count
load('DS150703-1.mat')
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg{i},ds_id,0,1);
    DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg{i}));
    raster_dg{i} = get_ds_raster(datadg{i}, ds_id);
    for j = 1:length(raster_dg{i})
        if(dg_idx(j, i))
            raster_dg{i}{j} = [];
        end
    end
end

% max firing rate
load('DS150703-1.mat')
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim(datadg{i},ds_id,0,1.5);
    DG{i} = sort_direction(dscellanalysis(MaxRate, StimComb,datadg{i}));
    raster_dg{i} = get_ds_raster(datadg{i}, ds_id);
    for j = 1:length(raster_dg{i})
        if(dg_idx(j, i))
            raster_dg{i}{j} = [];
        end
    end
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
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [4 5]);

%% plot cell summary

for cc = 73:73 %length(ds_id)
    plot_ds_raster(DG_cut, raster_dg_cut, cc, ds_id(cc), ll, 2, 3, 0)
end

for cc = 73:73 %length(ds_id)
    plot_ds_raster(DG(5), raster_dg(5), cc, ds_id(cc), ll, 1, 1, 0)
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
pc1 = 1; pc2 = 3;
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

% several cells are missing during mapping between NDF0 grating and NDF3
% grating, simply put them into ON-OFF DSGC class here
% idx_all = 1:length(ds_id);
% idx_sub{2} = setdiff(idx_all, idx_sub{1});
% id_sub{2} = ds_id(idx_sub{2});
% 
% figure
% plot(scores(idx_sub{1}, pc1), scores(idx_sub{1}, pc2), 'ro', scores(idx_sub{2}, pc1), scores(idx_sub{2}, pc2), 'bo')
% xlabel('1st Principal Component')
% ylabel('3rd Principal Component')
% title('NDF 0')


figure
v = 4*datadg{2}.stimulus.params.SPATIAL_PERIOD./datadg{2}.stimulus.params.TEMPORAL_PERIOD;
subplot(1, 2, 1)
semilogx(v, exciseColumn(MAG_all_norm_dg{2}(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg{2}(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
% title(ll{2})
xlim([v(end) v(1)])

subplot(1, 2, 2)
semilogx(v, exciseColumn(MAG_all_norm_dg{5}(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg{5}(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
% title(ll{5})
xlim([v(end) v(1)])

t = 4;
figure
compass(DG{5}.U{t}(idx_sub{1}), DG{5}.V{t}(idx_sub{1}), 'r')
hold on
compass(DG{5}.U{t}(idx_sub{2}), DG{5}.V{t}(idx_sub{2}), 'b')

%% plot average tunning curve
color = 'rbk';
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
    xlabel('micron/second')
    ylabel('response')
end
legend('on DSGC', 'on-off DSGC', 'location', 'southeast')


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

%%
d = 5;
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
n = 5;
duration = 2368;
[raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, ~, StimComb] = get_spikescellstim_mb(datamb{i},ds_id,duration,1);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb{i} = get_mb_raster(datamb{i}, ds_id, duration);
    for j = 1:length(raster_mb{i})
        if(mb_idx(j, i))
            raster_mb{i}{j} = [];
        end
    end
    trial_dur{i} = get_mb_trial_dur(datamb{i}, 800, 600, 0.5);
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
ct = {'posterior', 'inferior', 'anterior', 'superior'};
close all

%% plot cell summary
Title = {'BarWidth 120, white', 'BarWidth 240, white', 'BarWidth 360, white'; 'BarWidth 120, black', 'BarWidth 240, black', 'BarWidth 360, black'};
for cc = 1:length(ds_id)
    plot_mb_raster_one(MB, raster_mb, trial_dur, cc, ds_id(cc), Title, 2, 3, 1)
end


%% Direction tuning (moving bar)
% all ds cells
color = 'brgkc';

% t = 2;
% ct_on = {'Inferior', 'anterior', 'superior'};
% idx_dir = idx_dir_on;
% id_dir = id_dir_on;
dirn = 4;
D = 5;
T = 2;

for d = 1:n 
    p_direction = MB{D}.angle{T}';
    xx = 0:pi/4:7*pi/4;
    xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
    xx(xx>pi) = xx(xx>pi)-2*pi;
    xx(xx<-pi) = xx(xx<-pi)+2*pi;


    subplot(2, 3, d)
    for i = 1:dirn
        for cc = 1:length(ds_id)
            if ~mb_idx(cc, d)
            [xsort, seq] = sort(xx(cc, :));
            y_temp = MB{d}.rho{T}(cc, :);
            plot(xsort, y_temp(seq), 'b')
            hold on
            end
        end
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    title(ll{d})
    xlim([-pi pi])
end

%subtypes
clear rho_mb_mean rho_mb_ste dsi_mb_mean dsi_mb_ste rho_mb dsi_mb
for t = 1:2
    for bw = 1:3
        for cl = 1:2
            p_direction = MB{D}.angle{t,bw,cl}';
            xx = 0:pi/4:7*pi/4;
            xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
            xx(xx>pi) = xx(xx>pi)-2*pi;
            xx(xx<-pi) = xx(xx<-pi)+2*pi;

            figure
            for d = 1:5
                subplot(2, 3, d)
                for i = 1:dirn
                    rho_mb{d}{i}{t,bw,cl} = [];
                    RHO_mb{d}{i}{t,bw,cl} = [];
                    dsi_mb{d}{i}{t,bw,cl} = [];
                    for cc = 1:length(idx_dir{i})
                        if ~mb_idx(idx_dir{i}(cc), d) && sum(MB{d}.rho{t, bw,cl}(idx_dir{i}(cc), :))>0
                        [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
                        y_temp = MB{d}.rho{t,bw,cl}(idx_dir{i}(cc), :);
                        Y_temp = MB{d}.RHO{t,bw,cl}(idx_dir{i}(cc), :);
                        plot(xsort, y_temp(seq), color(i))
                        ylim([0 1])
            %             pause
                        hold on
                        rho_mb{d}{i}{t,bw,cl} = [rho_mb{d}{i}{t,bw,cl}; y_temp(seq)];
                        RHO_mb{d}{i}{t,bw,cl} = [RHO_mb{d}{i}{t,bw,cl}; Y_temp(seq)];
                        dsi_mb{d}{i}{t,bw,cl} = [dsi_mb{d}{i}{t,bw,cl}; MB{d}.dsindex{t,bw,cl}(idx_dir{i}(cc))];
                        end
                    end
                    rho_mb_mean{t,bw,cl}{d}(i, :) = mean(rho_mb{d}{i}{t,bw,cl},1);
                    rho_mb_ste{t,bw,cl}{d}(i, :) = std(rho_mb{d}{i}{t,bw,cl},[],1)/sqrt(size(rho_mb{d}{i}{t,bw,cl}, 1));
                    RHO_mb_mean{t,bw,cl}{d}(i, :) = mean(RHO_mb{d}{i}{t,bw,cl},1);
                    RHO_mb_ste{t,bw,cl}{d}(i, :) = std(RHO_mb{d}{i}{t,bw,cl},[],1)/sqrt(size(RHO_mb{d}{i}{t,bw,cl}, 1));
                    dsi_mb_mean{t,bw,cl}{d}(i) = mean(dsi_mb{d}{i}{t,bw,cl});
                    dsi_mb_ste{t,bw,cl}{d}(i) = std(dsi_mb{d}{i}{t,bw,cl})/sqrt(length(dsi_mb{d}{i}{t,bw,cl}));
                end
                xlabel('direction (rad)')
                ylabel('normalized response')
                title(ll{d})
                xlim([-pi pi])
            end
            dsi_mb_mean_all{t,bw,cl} = cell2mat(dsi_mb_mean{t,bw,cl}');
            dsi_mb_ste_all{t,bw,cl} = cell2mat(dsi_mb_ste{t,bw,cl}');            
        end
    end
end


% plot average (cell type)
for t = 1:2
    for bw = 1:3
        for cl = 1:2
            h = figure;
            set(h, 'Position', [1 1 1520,1080])
            for d = 1:5
                subplot(2, 3, d)
                for i = 1:dirn
                    errorbar(xsort/pi*180, rho_mb_mean{t,bw,cl}{d}(i, :), rho_mb_ste{t,bw,cl}{d}(i, :), color(i));
                    hold on
                end
                xlabel('degrees')
                ylabel('normalized mean spike number')
                title(ll{d});
            end
            legend(ct)
        end
    end
end

% plot average (light level)
for t = 1:2
    for bw = 1:3
        for cl = 1:2
            h = figure;
            set(h, 'Position', [1,1,980,1080])            
            for i = 1:dirn
                subplot(2, 2, i)
                for d = 1:5
                    errorbar(xsort/pi*180, rho_mb_mean{t,bw,cl}{d}(i, :), rho_mb_ste{t,bw,cl}{d}(i, :), color(d));
                    hold on
                end
                ylim([0 1])
                xlabel('degrees')
                ylabel('normalized average response')
                title(ct{i})
            end
            legend(ll)
        end
    end
end
% DSI

for t = 1:2
    for bw = 1:3
        for cl = 1:2
            clear dsi_mb_mean dsi_mb_ste
            dsi_mb_mean = dsi_mb_mean_all{t,bw,cl};
            dsi_mb_ste = dsi_mb_ste_all{t,bw,cl};
            FigHandle = figure;
            set(FigHandle, 'Position', [100, 100, 1000, 500]);

            xtick = ct;
            model_series = [dsi_mb_mean(1,1) dsi_mb_mean(2,1) dsi_mb_mean(3,1) dsi_mb_mean(4,1) dsi_mb_mean(5,1); dsi_mb_mean(1,2) dsi_mb_mean(2,2) dsi_mb_mean(3,2) dsi_mb_mean(4,2) dsi_mb_mean(5,2);dsi_mb_mean(1,3) dsi_mb_mean(2,3) dsi_mb_mean(3,3) dsi_mb_mean(4,3) dsi_mb_mean(5,3); dsi_mb_mean(1,4) dsi_mb_mean(2,4) dsi_mb_mean(3,4) dsi_mb_mean(4,4) dsi_mb_mean(5,4)];   
            model_error = [dsi_mb_ste(1,1) dsi_mb_ste(2,1) dsi_mb_ste(3,1) dsi_mb_ste(4,1) dsi_mb_ste(5,1); dsi_mb_ste(1,2) dsi_mb_ste(2,2) dsi_mb_ste(3,2) dsi_mb_ste(4,2) dsi_mb_ste(5,2);dsi_mb_ste(1,3) dsi_mb_ste(2,3) dsi_mb_ste(3,3) dsi_mb_ste(4,3) dsi_mb_ste(5,3); dsi_mb_ste(1,4) dsi_mb_ste(2,4) dsi_mb_ste(3,4) dsi_mb_ste(4,4) dsi_mb_ste(5,4)];
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
        end
    end
end

%% fitting
xdata = MB{1}.theta{1}(1, :);
for i = 1:4
    for j = 1:5
        ydata = RHO_mb_mean{1,2,1}{j}(i, :);
        [f, g] = fit_cos(xdata, ydata);
        Ymax{j}(i) = f.ymax;
        Phi{j}(i) = f.phi;
        Alpha{j}(i) = f.alpha;
        B{j}(i) = f.b;

        yfit = f.ymax.*(0.5+0.5*cos(xdata+f.phi)).^f.alpha+f.b;
        figure(1)
        plot(xdata, ydata, 'b')
        hold on
        plot(xdata, yfit, 'r')
    end
end

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
            if ~dg_idx(cc, d)
            [xsort, seq] = sort(xx(cc, :));
            y_temp = DG_cut{d}.rho{T}(cc, :);
            plot(xsort, y_temp(seq), 'b')
            hold on
            end
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
            if ~dg_idx(idx_dir{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx_dir{i}(cc), :))>0
            [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
            y_temp = DG_cut{d}.rho{T}(idx_dir{i}(cc), :);
            plot(xsort, y_temp(seq), color(i))
            ylim([0 1])
%             pause
            hold on
            rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
            dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx_dir{i}(cc))];
            end
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

%% on-off vs on
%subtypes
clear rho_dg_mean rho_dg_ste dsi_dg_mean dsi_dg_ste
for d = 1:5
    subplot(2, 3, d)
    for i = 1:2
        rho_dg{d}{i} = [];
        dsi_dg{d}{i} = [];
        for cc = 1:length(idx_sub{i})
            if ~dg_idx(idx_sub{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx_sub{i}(cc), :))>0
            [xsort, seq] = sort(xx(idx_sub{i}(cc), :));
            y_temp = DG_cut{d}.rho{T}(idx_sub{i}(cc), :);
            plot(xsort, y_temp(seq), color(i))
            ylim([0 1])
%             pause
            hold on
            rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
            dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx_sub{i}(cc))];
            end
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
    for i = 1:2
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(i));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
    title(ll{d});
end
legend('on DSGC', 'on-off DSGC')

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
% delete cell 5239, looks weird, check back later
idx_dir_on_ = idx_dir_on;
idx_dir_on_{3}(5) = [];
idx_sub_ = idx_sub;
idx_sub_{1}(13) = [];

for i = [2 5]
%     tp = datadg{1}.stimulus.params.TEMPORAL_PERIOD;
    tp = datadg{2}.stimulus.params.TEMPORAL_PERIOD;
    [DC{i}, F1{i}, F2{i}] = deal(zeros(length(ds_id), length(tp)));
    for rgc = 1:length(ds_id)
%         if ~isempty(raster_dg_cut{i}{rgc}) && ~dg_idx(rgc, i)
        if ~isempty(raster_dg{i}{rgc}) && ~dg_idx(rgc, i)
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

        end
        
    end
    ratio{i} = F2{i}./F1{i};
%     for ct = 1:4
%         ratio_dir{ct}{i} = ratio{i}(idx_dir{ct}, :);
%         ratio_dir{ct}{i} = exciseRows_empty(ratio_dir{ct}{i});
%         ratio_dir_mean(i, ct, :) = mean(ratio_dir{ct}{i});
%         ratio_dir_ste(i, ct, :) = std(ratio_dir{ct}{i})/sqrt(size(ratio_dir{ct}{i}, 1));
%     end
%     for ct = 1:3
%         ratio_dir_on{ct}{i} = ratio{i}(idx_dir_on_{ct}, :);
%         ratio_dir_on{ct}{i} = exciseRows_empty(ratio_dir_on{ct}{i});
%         ratio_dir_mean_on(i, ct, :) = mean(ratio_dir_on{ct}{i}, 1);
%         ratio_dir_ste_on(i, ct, :) = std(ratio_dir_on{ct}{i}, [], 1)/sqrt(size(ratio_dir_on{ct}{i}, 1));
%     end
%     for ct = 1:2
%         ratio_onoff{ct}{i} = ratio{i}(idx_sub_{ct}, :);
%         ratio_onoff{ct}{i} = exciseRows_empty(ratio_onoff{ct}{i});
%         ratio_onoff_mean(i, ct, :) = mean(ratio_onoff{ct}{i}, 1);
%         ratio_onoff_ste(i, ct, :) = std(ratio_onoff{ct}{i}, [], 1)/sqrt(size(ratio_onoff{ct}{i}, 1));
%     end
%     ratio{i} = exciseRows_empty(ratio{i});
%     ratio_mean(i, :) = mean(ratio{i});
%     ratio_ste(i, :) = std(ratio{i})/sqrt(size(ratio{i}, 1));
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
%     for ct = 1:4
%         ratio_dir{ct}{i} = ratio{i}(idx_dir{ct}, :);
%         ratio_dir{ct}{i} = exciseRows_empty(ratio_dir{ct}{i});
%         ratio_dir_mean(i, ct, :) = mean(ratio_dir{ct}{i});
%         ratio_dir_ste(i, ct, :) = std(ratio_dir{ct}{i})/sqrt(size(ratio_dir{ct}{i}, 1));
%     end
%     ratio_oo{i} = ratio{i}(idx_sub{2}, :);
%     ratio_oo{i} = exciseRows_empty(ratio_oo{i});
%     ratio_oo_mean(i, :) = mean(ratio_oo{i});
%     ratio_oo_ste(i, :) = std(ratio_oo{i})/sqrt(size(ratio_oo{i}, 1));
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


%% frequency doubling of moving bar response
color = 'brgkm';
% get spike histgram of preferred direction
bin_size = 0.2;
for i = 1:n
    for t = 1:2
        for bw = 1:3
            for cl = 1:2
                dur = max(trial_dur{i}(t,bw));
                XX = bin_size/2:bin_size:dur;
                hist_mb_p{i}{bw,t,cl} = zeros(length(ds_id), length(XX));
                for cc = 1:length(ds_id)
                    if ~mb_idx(cc, i) && ~isempty(raster_p_sum_mb{i}{cc})
                        hist_mb_p{i}{bw,t,cl}(cc, :) = hist(raster_p_sum_mb{i}{cc}{bw,t,cl}, XX);
                    end
                end
            end
        end
    end
end

%plot
% x = [4 4 4 4]; y = [4 4 5 5];
x = [2 3 2]; y = [2 3 3];
for bw = 3:3
    for t = 1:1
        for cl = 2:2
            for ct = 1:3
                figure
                set(gcf, 'Position', [1,1,1080,1080])
                for cc = 1:length(idx_dir_on{ct})
                    subplot(x(ct), y(ct), cc)
            %         subplot(3, 4, cc-3)
                    dur = max(trial_dur{1}(t,bw));
                    XX = bin_size/2:bin_size:dur;
            %         for i = 4:4
                    for i = 1:n
                        plot(XX, 1.25*hist_mb_p{i}{bw,t,cl}(idx_dir_on{ct}(cc), :), color(i))
                        hold on
                    end
                    title(num2str(ds_id(idx_dir_on{ct}(cc))))
                    xlim([0 dur])
                end
                legend('NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0')
            end
        end
    end
end

x = [3 3 5 4]; y = [4 4 5 4];
for bw = 1:1
    for t = 1:1
        for cl = 1:1
            for ct = 1:4
                figure
                set(gcf, 'Position', [1,1,1080,1080])
                for cc = 1:length(idx_dir{ct})
                    subplot(x(ct), y(ct), cc)
            %         subplot(3, 4, cc-3)
                    dur = max(trial_dur{1}(t,bw));
                    XX = bin_size/2:bin_size:dur;
            %         for i = 4:4
                    for i = 1:n
                        plot(XX, 1.25*hist_mb_p{i}{bw,t,cl}(idx_dir{ct}(cc), :), color(i))
                        hold on
                    end
                    title(num2str(ds_id(idx_dir{ct}(cc))))
                    xlim([0 dur])
                end
                legend('NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0')
            end
        end
    end
end

    
%%c F2/F1 ratio
clear ratio_ste ratio_onoff ratio_avg max_onoff
refresh_rate = 60.3578;
range = 1/3;

for t = 1:1
    for bw = 1:1
        for cl = 1:1
            delta = datamb{1}.stimulus.params.DELTA(t);
            bar_width = datamb{1}.stimulus.params.BAR_WIDTH(bw);
            bar_time = bar_width/delta/refresh_rate;
            for ct = 1:4
                for cc = 1:length(idx_dir{ct})
                    figure(1)
                    dur = max(trial_dur{1}(t,bw));
                    XX = bin_size/2:bin_size:dur;
                    for i = 1:n
                        plot(XX, hist_mb_p{i}{bw,t,cl}(idx_dir{ct}(cc), :), color(i))
                        hold on
                    end
                    title(num2str(ds_id(idx_dir{ct}(cc))))
                    [x, ~] = ginput;
                    on_time{ct}(cc, t,bw) = x;
                    left_on = max(round((on_time{ct}(cc, t,bw)-bar_time*range)/bin_size), 1);
                    right_on = min(round((on_time{ct}(cc, t,bw)+bar_time*range)/bin_size), size(hist_mb_p{1}{bw,t,cl}, 2));
                    left_off = max(round((on_time{ct}(cc, t,bw)+bar_time*(1-range))/bin_size), 1);
                    right_off = min(round((on_time{ct}(cc, t,bw)+bar_time*(1+range))/bin_size), size(hist_mb_p{1}{bw,t,cl}, 2));
                    for i = 1:n
                        max_onoff{i}{ct}{bw,t,cl}(cc, 1) = max(hist_mb_p{i}{bw,t,cl}(idx_dir{ct}(cc), left_on:right_on));
                        max_onoff{i}{ct}{bw,t,cl}(cc, 2) = max(hist_mb_p{i}{bw,t,cl}(idx_dir{ct}(cc), left_off:right_off));
                    end
                    hold off
                end
                legend('NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0')
            end

            for i = 1:5
                for ct = 1:4
                    for cc = 1:size(max_onoff{i}{ct}{bw,t,cl}, 1)
                        if cl == 1
                            ratio_onoff{ct}{bw,t,cl}{i}(cc) = max_onoff{i}{ct}{bw,t,cl}(cc, 2)/max_onoff{i}{ct}{bw,t,cl}(cc, 1);
                        elseif cl == 2
                            ratio_onoff{ct}{bw,t,cl}{i}(cc) = max_onoff{i}{ct}{bw,t,cl}(cc, 1)/max_onoff{i}{ct}{bw,t,cl}(cc, 2);
                        end
                    end
                    index = logical(isnan(ratio_onoff{ct}{bw,t,cl}{i}) + isinf(ratio_onoff{ct}{bw,t,cl}{i}));
                    ratio_onoff{ct}{bw,t,cl}{i}(index) = [];
                    % delete cell 4231
%                     if i == 4
%                         ratio_onoff{1}{1}{4}(5) = [];
%                     end
                    % % % % 
                    ratio_avg{ct,bw,t,cl}(i) = mean(ratio_onoff{ct}{bw,t,cl}{i});
                    ratio_ste{ct,bw,t,cl}(i) = std(ratio_onoff{ct}{bw,t,cl}{i})/sqrt(length(ratio_onoff{ct}{bw,t,cl}{i}));
                end
            end
        end
    end
end


%%
delta = datamb{1}.stimulus.params.DELTA;
bar_width = datamb{1}.stimulus.params.BAR_WIDTH;
color_ = {'white bar', 'black bar'};

for t = 1:2
    for bw = 1:3
        for cl = 1:2
            figure
            set(gcf, 'DefaultLineLineWidth', 1.5)
            xtick = {'posterior', 'inferior', 'anterior', 'superior'};
            model_series = cell2mat(ratio_avg(:,bw,t,cl));
            model_error = cell2mat(ratio_ste(:,bw,t,cl));
            h = bar(model_series);
            set(h,'BarWidth',1); % The bars will now touch each other

            set(gca,'XTicklabel',xtick)
            ylabel('OFF/ON')
            legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
            hold on;
            title(['bw:' num2str(bar_width(bw)) '  speed:' num2str(delta(t)) '  ' color_{cl}])

            numgroups = size(model_series, 1); 
            numbars = size(model_series, 2); 

            groupwidth = min(0.8, numbars/(numbars+1.5));

            for i = 1:numbars
            % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
            x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
            errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
            end
        end
    end
end

%% full field pulses

n_ffp = 5;
bin_size = 0.1;
XX = bin_size/2:bin_size:6-bin_size/2;

[raster_ff, raster_ff_all] = deal(cell(n_ffp, 1));
for d = 1:n_ffp
    [raster_ff{d}, raster_ff_all{d}] = get_ffp_raster(dataffp{d}, ds_id, 3);
    for j = 1:length(raster_ff{d})
        if(ffp_idx(j, d))
            raster_ff{d}{j} = [];
            raster_ff_all{d}{j} = [];
        end
        step_raster = get_raster(raster_ff_all{d}{j}, [0 6], 'plot', 0);
        hist_on = hist(step_raster{1}, XX);
        hist_off = hist(step_raster{2}, XX);
        ratio_ffp{d}(j) = max(hist_off)/max(hist_on);
    end
    ratio_temp = ratio_ffp{d};
    ratio_temp(isnan(ratio_temp)) = [];
    ratio_ffp_avg(d) = mean(ratio_temp);
    ratio_ffp_ste(d) = std(ratio_temp)/sqrt(length(ratio_temp));
end

figure
errorbar([0:4], ratio_ffp_avg, ratio_ffp_ste)
title('ffp off on ratio')
xlabel('light level')
ylabel('off on ratio')
xlim([-0.5 4.5])

for i = 1:1%length(ds_id) 
%     if ~isempty(raster_ff{1}{i}) || ~isempty(raster_ff{2}{i})
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 1800 800])
        for d = 1:n_ffp
        subplot(1, n_ffp, d)
        plot_ffp(raster_ff{d}, raster_ff_all{d}, i, 3)
%         title([num2str(ds_id(i)) ' ' ll{d}])
        end
        
%         print_close(1, [24, 12], num2str(ds_id(i)))
%     end
end

%% ffp psth
repeats = 25;
dur = 3;
dt = 0.1;
XX = dt/2:dt:4*dur-dt/2;
for i = 1:n
    for cc = 1:length(ds_id)
        psth_ffp{i}(cc, :) = hist(raster_ff_all{i}{cc}, XX)/dt/repeats;
        psth_norm_ffp{i}(cc, :) = psth_ffp{i}(cc, :)/max(psth_ffp{i}(cc, :));
    end
end

% on cell
for i = 1:n
    for ct = 1:3
        temp_psth = psth_norm_ffp{i}(idx_dir_on{ct}, :);
        psth_norm_ffp_avg{i}(ct,:) = mean(temp_psth(~ffp_idx(idx_dir_on{ct}, i), :), 1);
        psth_norm_ffp_ste{i}(ct,:) = std(temp_psth(~ffp_idx(idx_dir_on{ct}, i), :), [], 1)/sqrt(sum(~ffp_idx(idx_dir_on{ct})));
        temp_psth = psth_ffp{i}(idx_dir_on{ct}, :);
        psth_ffp_avg{i}(ct,:) = mean(temp_psth(~ffp_idx(idx_dir_on{ct}, i), :), 1);
        psth_ffp_ste{i}(ct,:) = std(temp_psth(~ffp_idx(idx_dir_on{ct}, i), :), [], 1)/sqrt(sum(~ffp_idx(idx_dir_on{ct})));
    end
    temp_psth = psth_norm_ffp{i}(cell2mat(idx_dir_on), :);
    psth_norm_ffp_avgall{i} = mean(temp_psth(~ffp_idx(cell2mat(idx_dir_on), i), :), 1);
    psth_norm_ffp_steall{i} = std(temp_psth(~ffp_idx(cell2mat(idx_dir_on), i), :), [],1)/sqrt(sum(~ffp_idx(cell2mat(idx_dir_on))));
    temp_psth = psth_ffp{i}(cell2mat(idx_dir_on), :);
    psth_ffp_avgall{i} = mean(temp_psth(~ffp_idx(cell2mat(idx_dir_on), i), :), 1);
    psth_ffp_steall{i} = std(temp_psth(~ffp_idx(cell2mat(idx_dir_on), i), :), [],1)/sqrt(sum(~ffp_idx(cell2mat(idx_dir_on))));
end

% unnormalized ffp psth
for i = 1:n
cell_type = {'inferior', 'anterior', 'superior'};
figure
set(gcf, 'Position', [1,1,880,1080])
subplot(5,1,1)
imagesc([0.5 1 0.5 0])
set(gca, 'xtick', [], 'ytick', [])
colormap gray
title('stimulus')
for ct = 1:3
    for cc = 1:length(id_dir_on{ct})
        subplot(5,1,5)
        plot(XX, psth_ffp{i}(idx_dir_on{ct}(cc), :), 'color', [0.5 0.9 1])
        hold on
        subplot(5,1,ct+1)
        plot(XX, psth_ffp{i}(idx_dir_on{ct}(cc), :), 'color', [0.5 0.9 1])
        hold on
    end
    plot(XX, psth_ffp_avg{i}(ct,:), 'b');
    title(cell_type{ct})
end
subplot(5,1,5)
plot(XX, psth_ffp_avgall{i}, 'b')
title('all')
xlabel('second')
ylabel('spike/second')

% normalized ffp psth
figure
set(gcf, 'Position', [1,1,880,1080])
subplot(5,1,1)
imagesc([0.5 1 0.5 0])
set(gca, 'xtick', [], 'ytick', [])
colormap gray
title('stimulus')
for ct = 1:3
    for cc = 1:length(id_dir_on{ct})
        subplot(5,1,5)
        plot(XX, psth_norm_ffp{i}(idx_dir_on{ct}(cc), :), 'color', [0.5 0.9 1])
        hold on
        subplot(5,1,ct+1)
        plot(XX, psth_norm_ffp{i}(idx_dir_on{ct}(cc), :), 'color', [0.5 0.9 1])
        hold on
    end
    plot(XX, psth_norm_ffp_avg{i}(ct,:), 'b');
    title(cell_type{ct})
end
subplot(5,1,5)
plot(XX, psth_norm_ffp_avgall{i}, 'b')
title('all')
xlabel('second')
ylabel('normalized firing rate')    
end

% on-off cell
for i = 1:n
    for ct = 1:4
        temp_psth = psth_norm_ffp{i}(idx_dir{ct}, :);
        psth_norm_ffp_avg{i}(ct,:) = mean(temp_psth(~ffp_idx(idx_dir{ct}, i), :), 1);
        psth_norm_ffp_ste{i}(ct,:) = std(temp_psth(~ffp_idx(idx_dir{ct}, i), :), [], 1)/sqrt(sum(~ffp_idx(idx_dir{ct})));
        temp_psth = psth_ffp{i}(idx_dir{ct}, :);
        psth_ffp_avg{i}(ct,:) = mean(temp_psth(~ffp_idx(idx_dir{ct}, i), :), 1);
        psth_ffp_ste{i}(ct,:) = std(temp_psth(~ffp_idx(idx_dir{ct}, i), :), [], 1)/sqrt(sum(~ffp_idx(idx_dir{ct})));
    end
    temp_psth = psth_norm_ffp{i}(cell2mat(idx_dir), :);
    psth_norm_ffp_avgall{i} = mean(temp_psth(~ffp_idx(cell2mat(idx_dir), i), :), 1);
    psth_norm_ffp_steall{i} = std(temp_psth(~ffp_idx(cell2mat(idx_dir), i), :), [],1)/sqrt(sum(~ffp_idx(cell2mat(idx_dir))));
    temp_psth = psth_ffp{i}(cell2mat(idx_dir), :);
    psth_ffp_avgall{i} = mean(temp_psth(~ffp_idx(cell2mat(idx_dir), i), :), 1);
    psth_ffp_steall{i} = std(temp_psth(~ffp_idx(cell2mat(idx_dir), i), :), [],1)/sqrt(sum(~ffp_idx(cell2mat(idx_dir))));
end

% unnormalized ffp psth
for i = 1:n
cell_type = {'posterior', 'inferior', 'anterior', 'superior'};
figure
set(gcf, 'Position', [1,1,880,1080])
subplot(6,1,1)
imagesc([0.5 1 0.5 0])
set(gca, 'xtick', [], 'ytick', [])
colormap gray
title('stimulus')
for ct = 1:4
    for cc = 1:length(id_dir{ct})
        if(~ffp_idx(idx_dir{ct}(cc), i))
        subplot(6,1,6)
        plot(XX, psth_ffp{i}(idx_dir{ct}(cc), :), 'color', [0.5 0.9 1])
        hold on
        subplot(6,1,ct+1)
        plot(XX, psth_ffp{i}(idx_dir{ct}(cc), :), 'color', [0.5 0.9 1])
        hold on
        end
    end
    plot(XX, psth_ffp_avg{i}(ct,:), 'b');
    title(cell_type{ct})
end
subplot(6,1,6)
plot(XX, psth_ffp_avgall{i}, 'b')
title('all')
xlabel('second')
ylabel('spike/second')

% normalized ffp psth
figure
set(gcf, 'Position', [1,1,880,1080])
subplot(6,1,1)
imagesc([0.5 1 0.5 0])
set(gca, 'xtick', [], 'ytick', [])
colormap gray
title('stimulus')
for ct = 1:4
    for cc = 1:length(id_dir{ct})
        if(~ffp_idx(idx_dir{ct}(cc), i))
        subplot(6,1,6)
        plot(XX, psth_norm_ffp{i}(idx_dir{ct}(cc), :), 'color', [0.5 0.9 1])
        hold on
        subplot(6,1,ct+1)
        plot(XX, psth_norm_ffp{i}(idx_dir{ct}(cc), :), 'color', [0.5 0.9 1])
        hold on
        end
    end
    plot(XX, psth_norm_ffp_avg{i}(ct,:), 'b');
    title(cell_type{ct})
end
subplot(6,1,6)
plot(XX, psth_norm_ffp_avgall{i}, 'b')
title('all')
xlabel('second')
ylabel('normalized firing rate')    
end    

%% max firing rate speed tuning

max_rate_prefer = cell(5,1);
max_rate_null = cell(5,1);
max_rate_prefer_sub = cell(5, 1);
max_rate_null_sub = cell(5, 1);
n_idx = ceil(mod(p_idx{5}+4-0.1,8));
for i = [2 5]
    clear max_rate_prefer_temp max_rate_null_temp
    for tp = 1:length(DG{2}.RHO)
        for cc = 1:length(ds_id)
            if ~dg_idx(cc, i)
                max_rate_prefer_temp{cc}(tp) = DG{i}.RHO{tp}(cc, p_idx{5}(cc));
                max_rate_null_temp{cc}(tp) = DG{i}.RHO{tp}(cc, n_idx(cc));
            end
        end
    end
    max_rate_prefer{i} = max_rate_prefer_temp';
    max_rate_null{i} = max_rate_null_temp';
    for sub = 1:2
        max_rate_prefer_sub{i}{sub} = max_rate_prefer{i}(idx_sub{sub});
        max_rate_null_sub{i}{sub} = max_rate_null{i}(idx_sub{sub});
    end
end

% unnormalized individual tuning curve
color = 'br';
for i = [2 5]
    figure
    subplot(1, 2, 1)
    for sub = 1:2
        tuning_p = cell2mat(max_rate_prefer_sub{i}{sub});
        semilogx(repmat(v', 1, size(tuning_p, 1)), tuning_p', 'color', color(sub))
        hold on
    end
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')
    ylabel('Hz')
    subplot(1, 2, 2)
    for sub = 1:2
        tuning_n = cell2mat(max_rate_null_sub{i}{sub});
        semilogx(repmat(v', 1, size(tuning_n, 1)), tuning_n', 'color', color(sub))
        hold on
    end 
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')
    ylabel('Hz')

end

% unnormalized average tuning curve
for i = [2 5]
    figure
    subplot(1, 2, 1)
    for sub = 1:2
        tuning_p = cell2mat(max_rate_prefer_sub{i}{sub});
        errorbar(v, mean(tuning_p),std(tuning_p)/sqrt(size(tuning_p, 1)), 'color', color(sub))
        hold on
    end
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')
    ylabel('Hz')
    legend('ON', 'ON-OFF', 'Location', 'northwest')
    title('preferred direction')
    subplot(1, 2, 2)
    for sub = 1:2
        tuning_n = cell2mat(max_rate_null_sub{i}{sub});
        errorbar(v, mean(tuning_n), std(tuning_n)/sqrt(size(tuning_n, 1)), 'color', color(sub))
        hold on
    end
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')
    ylabel('Hz')
    legend('ON', 'ON-OFF', 'Location', 'northwest')
    title('null direction')
end

% normalized individual tuning curve
color = 'br';
for i = [2 5]
    figure
    subplot(1, 2, 1)
    for sub = 1:2
        tuning_p = cell2mat(max_rate_prefer_sub{i}{sub});
        tuning_p_n = tuning_p./repmat(max(tuning_p, [], 2), 1,9);
        tuning_p_n = nan2empty(tuning_p_n);
        semilogx(repmat(v', 1, size(tuning_p_n, 1)), tuning_p_n', 'color', color(sub))
        hold on
    end
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')
    subplot(1, 2, 2)
    for sub = 1:2
        tuning_n = cell2mat(max_rate_null_sub{i}{sub});
        tuning_n_n = tuning_n./repmat(max(tuning_n, [], 2), 1,9);
        tuning_n_n = nan2empty(tuning_n_n);
        semilogx(repmat(v', 1, size(tuning_n_n, 1)), tuning_n_n', 'color', color(sub))
        hold on
    end 
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')

end

% normalized average tuning curve
for i = [2 5]
    figure
    subplot(1, 2, 1)
    for sub = 1:2
        tuning_p = cell2mat(max_rate_prefer_sub{i}{sub});
        tuning_p_n = tuning_p./repmat(max(tuning_p, [], 2), 1,9);
        tuning_p_n = nan2empty(tuning_p_n);
        errorbar(v, mean(tuning_p_n),std(tuning_p_n)/sqrt(size(tuning_p_n, 1)), 'color', color(sub))
        hold on
    end
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')
    legend('ON', 'ON-OFF', 'Location', 'northwest')
    title('preferred direction')
    subplot(1, 2, 2)
    for sub = 1:2
        tuning_n = cell2mat(max_rate_null_sub{i}{sub});
        tuning_n_n = tuning_n./repmat(max(tuning_n, [], 2), 1,9);
        tuning_n_n = nan2empty(tuning_n_n);
        errorbar(v, mean(tuning_n_n), std(tuning_n_n)/sqrt(size(tuning_n_n, 1)), 'color', color(sub))
        hold on
    end
    set(gca,'xscale','log')
    xlim([min(v) max(v)])
    xlabel('micron/second')
    legend('ON', 'ON-OFF', 'Location', 'northwest')
    title('null direction')
end
%% save ei movie
save_ei_movie(datadg{5}, 873, '/Users/xyao/Desktop/cell873', 'frames', [10:20])

%% cross correlation
duration = 3600;
bin_size = 0.0005;
max_lag = 40;
xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
for c1 = 1:length(id_dir{4})
    FigHandle = figure;
    set(FigHandle, 'Position', [1 1 2000 2000])
    for c2 = 1:length(id_dir{4})
        if c1 ~= c2
            id1 = id_dir{4}(c1);
            id2 = id_dir{4}(c2);
            idx1 = get_cell_indices(datawn, id1);
            idx2 = get_cell_indices(datawn, id2);
            spikes1 = datawn.spikes{idx1};
            spikes1_TF= ceil(spikes1/bin_size);
            spikes1 = zeros(duration/bin_size, 1);
            spikes1(spikes1_TF) = 1;
            
            spikes2 = datawn.spikes{idx2};
            spikes2_TF= ceil(spikes2/bin_size);
            spikes2 = zeros(duration/bin_size, 1);
            spikes2(spikes2_TF) = 1;
            
            A = xcov(spikes1, spikes2, max_lag, 'coeff');
            subplot(4, 4, c2)
            bar(xx, A)
            title([num2str(id1) '  ' num2str(id2)])
            xlim([-0.02 0.02])

%             pause
        end
    end
    print_close(1, [24 12], num2str(c1));
end


%% neighboring pairs
pos = datawn.ei.position;
mode = 'neg';
neighbors = [];
for cc1 = 1:length(id_dir{4})
    for cc2 = cc1+1:length(id_dir{4})
        id1 = id_dir{4}(cc1);
        idx1 = get_cell_indices(datawn, id1);
        ei1 = datawn.ei.eis{idx1};
        com1 = ei_com_xy(ei1, pos, 30*3, mode);
        id2 = id_dir{4}(cc2);
        idx2 = get_cell_indices(datawn, id2);
        ei2 = datawn.ei.eis{idx2};
        com2 = ei_com_xy(ei2, pos, 30*3, mode);
        if pdist([com1;com2]) < 150
            neighbors = [neighbors; id1 id2];
        end
    end
end

coms = [];
for cc = 1:length(id_dir{4})
    id = id_dir{4}(cc);
    idx = get_cell_indices(datawn, id);
    ei = datawn.ei.eis{idx};
    com = ei_com_xy(ei, pos, 30*3, mode);
    coms = [coms; com];
end

figure
for cc = 1:length(id_dir{4})
    plot(coms(cc, 1), coms(cc, 2),'ko')
    hold on
%     text(coms(cc, 1)+5, coms(cc, 2)+5, num2str(id_dir{1}(cc)), 'FontSize', 10)
    
end
for cp = 1:size(corr_cells)
    idx1 = find(id_dir{4} == corr_cells(cp, 1));
    idx2 = find(id_dir{4} == corr_cells(cp, 2));
    plot([coms(idx1, 1), coms(idx2, 1)], [coms(idx1, 2), coms(idx2, 2)], 'k');
end
%% fit direction tuning curves
T = 1;
for dir = 1:4
    figure
    for i = 1:5
        CC = 1;
        for cc = 1:length(id_dir{dir})
            xdata = DG_cut{i}.theta{1}(idx_dir{dir}(cc), :);
            ydata = circshift(DG_cut{i}.rho{T}(idx_dir{dir}(cc), :), -2, 2);
            if ~(dg_idx(idx_dir{dir}(cc), i)) && sum(ydata) ~= 0
                [f, g] = fit_cos(xdata, ydata);
                Ymax{dir}{i}(CC) = f.ymax;
                Phi{dir}{i}(CC) = f.phi;
                Alpha{dir}{i}(CC) = f.alpha;
                Width{dir}{i}(CC) = acos(2*0.5^(1/f.alpha) - 1)/pi*360;
                B{dir}{i}(CC) = f.b;
                CC = CC + 1;
                
                xfit = linspace(0, 2*pi, 100);
                yfit = f.ymax.*(0.5+0.5*cos(xfit+f.phi)).^f.alpha+f.b;
                subplot(4, 6, cc)
                plot(xdata, ydata, 'b')
                hold on
                plot(xfit, yfit, 'r')
                ylim([0 1.1])
                width = acos(2 * (0.5.^(1/f.alpha) - 0.5))/pi*360;
                title(['width = ', num2str(width)])


            end
        end
        Ymax_mean(dir, i) = robust_mean(Ymax{dir}{i});
        Phi_mean(dir, i) = robust_mean(Phi{dir}{i});
        Alpha_mean(dir, i) = robust_mean(Alpha{dir}{i});
        B_mean(dir, i) = robust_mean(B{dir}{i});
    end
end


for dir = 1:4
    figure

    for i = 1:5
        subplot(2, 3, i)
        hist(Width{dir}{i})
        title(num2str(mean(Width{dir}{i})))
    end
end

for i = 1:5
    for dir = 1:4
        WidthMean(dir, i) = mean(Width{dir}{i});
        WidthSte(dir, i) = std(Width{dir}{i})/sqrt(length(Width{dir}{i}));
    end
end

marker = 'xosd';
figure
for dir = 1:4
    errorbar(0:4, WidthMean(dir, :), WidthSte(dir, :), 'Color', 'k', 'Marker', marker(dir), 'MarkerSize', 10)
    hold on
end
legend('superior', 'anterior', 'inferior', 'posterior')
ylim([0 300])
xlabel('light level')
ylabel('tuning width (degree)')
xlim([-0.5 4.5])

%% plot sta
figure
plot_rf(datawn, 6857);
figure
plot_rf(datawn, 7759);
figure
plot_rf(datawn, 5956);

%% SNR
ctr = 1;
color = 'brgk';
for ll = 1:5;
    rep = datamb{ll}.stimulus.repetitions;
    for ct = 1:4
        [SpikeN_temp{ct}, SpikeN_mean_temp{ct}, SpikeN_std_temp{ct}] = deal(cell(length(id_dir{ct}), 1));
        for cc = 1:length(id_dir{ct})
            raster_temp = raster_mb{ll}{idx_dir{ct}(cc)};
            if ~isempty(raster_temp)
                SpikeN_temp{ct}{cc} = squeeze(cellfun('length',raster_temp(1,ctr,:,:)));
                SpikeN_mean_temp{ct}{cc} = mean(SpikeN_temp{ct}{cc},2);
                SpikeN_std_temp{ct}{cc} = std(SpikeN_temp{ct}{cc},[],2);
            end
        end
    end
    SpikeN{ll} = SpikeN_temp;
    SpikeN_mean{ll} = SpikeN_mean_temp;
    SpikeN_std{ll} = SpikeN_std_temp;
end

clear SpikeN_mean_all SpikeN_std_all
for ll = 1:5
    SpikeN_mean_all(:, ll) = cat(1,SpikeN_mean{ll}{:});
    SpikeN_std_all(:, ll) = cat(1,SpikeN_std{ll}{:});
end
SpikeN_mean_all(any(cellfun(@isempty, SpikeN_mean_all), 2), :) = [];
SpikeN_std_all(any(cellfun(@isempty, SpikeN_std_all), 2), :) = [];

[rmax, i] = cellfun(@max, SpikeN_mean_all);
for ll = 1:5
    for cc = 1:size(SpikeN_mean_all, 1)
        stdmax(cc, ll) = SpikeN_std_all{cc, ll}(i(cc, ll));
    end
end

SNRmax = rmax./stdmax;
SNRmax = exciseRows_empty(SNRmax);

figure
plot(SNRmax(:, 1), SNRmax(:, 5), 'ko')
hold on
plot([0 15], [0 15], 'k--')
errorbar(mean(SNRmax(:, 1)), mean(SNRmax(:, 5)), std(SNRmax(:, 5)/sqrt(size(SNRmax, 1))), 'ro')
herrorbar(mean(SNRmax(:, 1)), mean(SNRmax(:, 5)), std(SNRmax(:, 1)/sqrt(size(SNRmax, 1))), 'ro')

% xlim([0 8])
% ylim([0 8])
xlabel('low light level SNR')
ylabel('high light level SNR')


%% subtract bgnd firing
% estimate background activity
n = 5;
interval = 1;
raster_interval_dg = deal(cell(n, 1));
for i = 1:n    
    raster_interval_dg{i} = get_ds_interval_raster(datadg{i}, ds_id, interval);
    for j = 1:length(raster_interval_dg{i})
        if iscell(raster_interval_dg{i}{j})
            bgfr{i}(j) = mean(cellfun(@length, raster_interval_dg{i}{j}))/interval;
        end
    end
    for dir = 1:4
        bgfr_ct{i}{dir} = [];
        for cc = 1:length(idx_dir{dir})
            if ~dg_idx(idx_dir{dir}(cc), i) && sum(DG_cut{i}.rho{1}(idx_dir{dir}(cc), :))>0
                bgfr_ct{i}{dir} = [bgfr_ct{i}{dir} bgfr{i}(idx_dir{dir}(cc))];
            end
        end
        bgfr_ct_mean(i, dir) = mean(bgfr_ct{i}{dir});
        bgfr_ct_ste(i, dir) = std(bgfr_ct{i}{dir})/sqrt(length(bgfr_ct{i}{dir}));
    end
end

celltype = {'posterior', 'inferior', 'anterior', 'superior'};marker = 'xosd';
figure
for i = 1:4
    errorbar([0:4], bgfr_ct_mean(:, i), bgfr_ct_ste(:, i), 'Color', 'k', 'Marker', marker(i), 'MarkerSize', 10);
    hold on
end
xlim([-0.5 4.5])
xlabel('R*/rod/s')
ylabel('background firing (Hz)')
legend(celltype)

n = 5;
% spike count
DG_bgnd = cell(n, 1);
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim(datadg{i},ds_id,0, 0.05);
    DG_bgnd{i} = sort_direction(dscellanalysis_bgnd_subtract(NumSpikesCell, StimComb, bgfr{i}));
end
[DG_bgnd_cut, ~, ~] = cut_dg(DG_bgnd, raster_dg, raster_p_sum, 2, [4 5]);

