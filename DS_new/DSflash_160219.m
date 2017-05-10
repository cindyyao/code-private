cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);


datadg = load_data('/Volumes/lab/analysis/2016-02-19-0/data004/data004', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-02-19-0/stimuli/s04.mat';
datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);

% dataflash = load_data('/Volumes/lab/analysis/2016-02-19-0/data001/data001', opt);
dataflash = load_data('/Volumes/lab/analysis/2016-02-19-0/data000-001-map/data000-001-map', opt);
dataflash.DfParams.NDF =   [5,5,5,4,4,4,4,4,3,3,3,2,2,1,0] ; % on filter turret 
dataflash.DfParams.Ftime = [2,4,8,2,3,4,6,8,2,4,8,2,4,2,2] ; % ms
dataflash.DfParams.interFlashInt = [3] ; % sec
dataflash.triggers(160:161) = [];

[NumSpikesCell, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [4 5]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

load('DS160219.mat')
ds_id_flash = ds_id(flash_idx);
ds_id_flash(26) = 6522;
ds_idx_flash = get_cell_indices(dataflash, ds_id_flash);

%% 
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0, 1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);
%     for j = 1:length(raster_dg{i})
%         if(dg_idx(j, i))
%             raster_dg{i}{j} = [];
%         end
%     end

delta_p = 4; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(DG{i});
rep = datadg.stimulus.repetitions;

%%
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [4 5]);

for cc = 1:length(ds_id)
    plot_ds_raster(DG_cut, raster_dg_cut, cc, ds_id(cc), '', 1, 1, 1)
end

% for cc = 2:2 %length(ds_id)
%     plot_ds_raster(DG, raster_dg, cc, ds_id(cc), '', 1, 1, 0)
% end

%% classification based on speed tunning
L = 1;
% mag_pca = MAG_all_norm_dg{L};
% mag_pca = mag_pca';
% [id_sub, idx_sub] = deal(cell(2, 1));
% 
% FigHandle = figure;
% set(FigHandle, 'Position', [1 1 380 400])
% 
% [~,scores,~,~] = princomp(mag_pca);
% pc1 = 1; pc2 = 2;
% plot(scores(:, pc1), scores(:, pc2), 'o')
% hold on
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

% automatic classification
% pca
mag_pca = MAG_all_norm_dg{L};
mag_pca = mag_pca';
mag_pca = nan2empty(mag_pca);

[id_sub, idx_sub] = deal(cell(2, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(mag_pca);
pc1 = 1; pc2 = 2;
plot(scores(:, pc1), scores(:, pc2), 'o')
hold on
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

figure
plot(scores(idx_sub{1}, pc1), scores(idx_sub{1}, pc2), 'ro', scores(idx_sub{2}, pc1), scores(idx_sub{2}, pc2), 'bo')
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
title('NDF 0')
% 
% several cells are missing during mapping between NDF0 grating and NDF3
% grating, simply put them into ON-OFF DSGC class here
id_sub{2} = setdiff(ds_id, id_sub{1});
[~, idx_sub{1}] = intersect(ds_id, id_sub{1});
[~, idx_sub{2}] = intersect(ds_id, id_sub{2});



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

t = 3;
figure
compass(DG{1}.U{t}(idx_sub{1}), DG{1}.V{t}(idx_sub{1}), 'r')
hold on
compass(DG{1}.U{t}(idx_sub{2}), DG{1}.V{t}(idx_sub{2}), 'b')

%%
d = 1;
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
id_dir_flash{2} = unique([id_dir_flash{2} 6522]);
idx_dir_flash{2} = unique([idx_dir_flash{2}; 29]);
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

trigger = 400:3:580;
ds_dark_raster = cell(length(ds_id_flash), 1);
ds_dark_raster_all = cell(length(ds_id_flash), 1);
ds_dark_hist_trial = cell(length(ds_id_flash), 1);
ds_dark_hist_mean = cell(length(ds_id_flash), 1);
for cc = 1:length(ds_id_flash)
    ds_dark_raster{cc} = get_raster(dataflash.spikes{ds_idx_flash(cc)}, ...
        trigger, 'stop', dataflash.DfParams.interFlashInt, ...
        'plot', 0);
    for t = 1:length(ds_dark_raster{cc})
        ds_dark_hist_trial{cc}{t} = hist(ds_dark_raster{cc}{t}, XX);
    end
    ds_dark_raster_all{cc} = sort(cell2mat(ds_dark_raster{cc}));
    ds_dark_hist_mean{cc} = hist(ds_dark_raster_all{cc}, XX)/length(trigger);
end

%% 2 alternative forced choice test
% use 400-580 second (roughly) as dark trials
window = 2;
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


color = 'rbgk';
figure
for ct = 1:4
    h_temp = plot(log10(Irel), Pc(idx_dir_flash{ct}, :)', 'color', color(ct));
    h(ct) = h_temp(1);
    hold on
end
for i = [5 17 21 28]
    h_temp= plot(log10(Irel), Pc(i, :)', 'c');
    h(5) = h_temp(1);
end
legend([h(1), h(2), h(3), h(4), h(5)], 'superior', 'anterior', 'inferior', 'posterior', 'ON')
% legend(h(1),'superior', h(2), 'anterior', h(3), 'inferior', h(4), 'posterior', h(5), 'ON')


% compare across cell type
for ct = 1:4
    Pc_dir{ct} = Pc(idx_dir_flash{ct}, :);
    Pc_dir_mean(ct, :) = mean(Pc_dir{ct});
    Pc_dir_ste(ct, :) = std(Pc_dir{ct})/sqrt(length(idx_dir_flash{ct}));
end

figure
for ct = 1:4
    errorbar(log10(Irel), Pc_dir_mean(ct, :), Pc_dir_ste(ct, :), 'color', color(ct));
    hold on
end
legend('superior', 'A(P)', 'inferior', 'P(A)')
xlabel('log(R*/rod)')
ylabel('probability')
%% intensity-response curve
window = 1;
for cc = 1:length(ds_id_flash)
    for ts = 1:length(trigger_set_i)
%         response(cc, ts) = sum(ds_flash_hist_mean{cc}{ts}(1:window/bin_size)) - sum(ds_flash_hist_mean{cc}{ts}(end-window/bin_size+1:end));
        response(cc, ts) = sum(ds_flash_hist_mean{cc}{ts}(1:window/bin_size)) - sum(ds_dark_hist_mean{cc}(1:window/bin_size));
    end
end
 
response_norm = response./repmat(max(response, [], 2), 1, length(trigger_set_i));

figure
for ct = 1:4
    plot(log10(Irel), response_norm(idx_dir{ct}, :), 'color', color(ct));
    hold on
end
%% raster plot

a = ceil((length(trigger_set_i)+1)/2);
for cc =[17 21 28]%length(ds_id_flash);
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
