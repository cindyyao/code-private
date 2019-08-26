cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

datadg = load_data('/Volumes/lab/analysis/2017-06-15-0/data017-sorted/data017-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2017-06-15-0/stimuli/s17.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [1 2]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

datarun = load_data('/Volumes/lab/analysis/2017-06-15-0/datadg-map/datadg-map', opt);
time_points = [480 1130 1780];
dataDG(1:4) = split_datarun(datarun, time_points);
dataDG{1}.names.stimulus_path = '/Volumes/lab/analysis/2017-06-15-0/stimuli/s02.txt';
dataDG{1} = load_stim(dataDG{1}, 'user_defined_trigger_interval', 10);

dataDG{2}.names.stimulus_path = '/Volumes/lab/analysis/2017-06-15-0/stimuli/s05.txt';
dataDG{2} = load_stim(dataDG{2}, 'user_defined_trigger_interval', 10);
dataDG{3}.names.stimulus_path = '/Volumes/lab/analysis/2017-06-15-0/stimuli/s09.txt';
dataDG{3} = load_stim(dataDG{3}, 'user_defined_trigger_interval', 10);
dataDG{4}.names.stimulus_path = '/Volumes/lab/analysis/2017-06-15-0/stimuli/s13.txt';
dataDG{4} = load_stim(dataDG{4}, 'user_defined_trigger_interval', 10);

dataDG{5} = datadg;

datarun = load_data('/Volumes/lab/analysis/2017-06-15-0/dataflash-map/dataflash-map', opt);
time_points = [2044 4342];
dataflash(1:3) = split_datarun(datarun, time_points);

dataflash{1}.DfParams.NDF =   [4,4,4,3,3,3,3,2,2,2] ; % on filter turret 
dataflash{1}.DfParams.Ftime = [2,4,8,2,4,6,8,2,4,8] ; % ms
dataflash{1}.DfParams.interFlashInt = [3] ; % sec

dataflash{2}.DfParams.NDF =   [4,4,4,3,3,3,3,3,2,2,2] ; % on filter turret 
dataflash{2}.DfParams.Ftime = [2,4,8,2,3,4,6,8,2,8,4] ; % ms
dataflash{2}.DfParams.interFlashInt = [3] ; % sec

dataflash{3}.DfParams.NDF =   [4,4,4,3,3,3,3,2,2,2,1,1] ; % on filter turret 
dataflash{3}.DfParams.Ftime = [2,4,8,2,4,6,8,2,4,8,2,8] ; % ms
dataflash{3}.DfParams.interFlashInt = [3] ; % sec

%%
d = 5;
t = 2;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(DG{d}.U{t}, DG{d}.V{t});
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}, DG{d}.V{t}, x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = I;
    id_dir{i} = ds_id(idx_dir{i});
end

%% DG
n = 5;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, ~, StimComb] = get_spikescellstim(dataDG{i},ds_id,0,1);
    DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,dataDG{i}));
    raster_dg{i} = get_ds_raster(dataDG{i}, ds_id);
    for j = 1:length(raster_dg{i})
        if(dg_idx(j, i))
            raster_dg{i}{j} = [];
        end
    end
end

delta_p = 2; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{1}.angle{delta_p});
    [raster_n_sum{i}, n_idx{i}] = get_ndirection_raster(raster_dg{i}, DG{1}.angle{delta_p});
    MAG_all_norm_dg{i} = normalize_MAG(DG{i});
    rep = dataDG{i}.stimulus.repetitions;
end

ll = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};

%% plot cell summary
conditions = {'NDF 4 ctr', 'NDF 4 SR', 'NDF 4 SR+TPMPA', 'NDF 4 wash', 'NDF 0'};
for dir = 1:4
    for cc = 1:length(id_dir{dir})
        plot_ds_raster(DG, raster_dg, idx_dir{dir}(cc), id_dir{dir}(cc), conditions, 2, 3, 1)
    end
end


%%
ds_id_flash = ds_id(~flash_idx);
for ct = 1:4
    [id_dir_flash{ct}, idx_dir_flash{ct}] = intersect(ds_id_flash, id_dir{ct});
end

%% parameters
interFlashIntVar = 0.005; % (sec) expected precision of trigger intervals

%% find sets of flash stimuli
% trigger = {[1:3:300], [601:3:900], [601:3:900], [1:3:300]};
for ds = 1:3
clear trigger_set_i
ts = 1; 
trigger_set_i{1} = [] ;
for a=1:length(dataflash{ds}.triggers) ; % for each trigger
    trigger_set_i{ts} = [trigger_set_i{ts},a] ; % put it in a set   
    if a<length(dataflash{ds}.triggers) ; % if its not the last trigger
        if sum(abs(dataflash{ds}.triggers(a+1)-dataflash{ds}.triggers(a)-dataflash{ds}.DfParams.interFlashInt)<interFlashIntVar)>0 ; % next trigger is has the right interval
            ts = ts ; % keep it in the set
        else
            ts = ts+1 ; % put it in a new set
            trigger_set_i{ts} = [] ;
        end
    end
end

i = 1;
while(i<=length(trigger_set_i))
    if length(trigger_set_i{i}) < 50
        trigger_set_i(i) = [];
    else
        i = i+1;
    end
end


%%
% interval = dataflash{ds}.DfParams.interFlashInt;
interval = 3;
ds_idx_flash = get_cell_indices(dataflash{ds}, ds_id_flash);
bin_size = 0.01; 
start = 0;
XX = start+bin_size/2:bin_size:dataflash{ds}.DfParams.interFlashInt-bin_size/2;
load('DS161208.mat', 'NdfCalibration')
% NdfCalibration(2, :) = NdfCalibration(2, :)/10;
Irel = (dataflash{ds}.DfParams.Ftime/1000).*NdfCalibration(2,dataflash{ds}.DfParams.NDF+1) ;
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
        raster_1cell{ts} = get_raster(dataflash{ds}.spikes{ds_idx_flash(cc)}, ...
            dataflash{ds}.triggers(trigger_set_i{ts}), 'start', start, 'stop', ...
            interval, 'plot', 0);
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

trigger = 1:3:200;
ds_dark_raster = cell(length(ds_id_flash), 1);
ds_dark_raster_all = cell(length(ds_id_flash), 1);
ds_dark_hist_trial = cell(length(ds_id_flash), 1);
ds_dark_hist_mean = cell(length(ds_id_flash), 1);
for cc = 1:length(ds_id_flash)
    ds_dark_raster{cc} = get_raster(dataflash{ds}.spikes{ds_idx_flash(cc)}, ...
        trigger, 'stop', interval, ...
        'plot', 0);
    for t = 1:length(ds_dark_raster{cc})
        ds_dark_hist_trial{cc}{t} = hist(ds_dark_raster{cc}{t}, XX);
    end
    ds_dark_raster_all{cc} = sort(cell2mat(ds_dark_raster{cc}));
    ds_dark_hist_mean{cc} = hist(ds_dark_raster_all{cc}, XX)/length(trigger);
end
DS_dark_raster_all(:, ds) = ds_dark_raster_all;

%% 2 alternative forced choice test
% use 400-580 second (roughly) as dark trials
clear Pc Pc_dir Pc_dir_mean Pc_dir_ste
Irel = (dataflash{ds}.DfParams.Ftime/1000).*NdfCalibration(2,dataflash{ds}.DfParams.NDF+1) ;
[Irel, i] = sort(Irel);

window = 1;
bin_n = window/bin_size;

Pc = zeros(length(ds_id_flash), length(trigger_set_i));
for cc = 1:length(ds_id_flash)
    for ts = 1:length(trigger_set_i)
%         trial_n = length(ds_flash_raster{1}{ts});
        trial_n = 50;
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
figure(1)
subplot(2,2,ds)
for ct = 1:4
    plot(log10(Irel), Pc(idx_dir_flash{ct}, :)', 'color', color(ct))
    hold on
end
xlabel('log(R*/rod)')
ylabel('probability')
title('wash')
% compare across cell type
for ct = 1:4
    Pc_dir{ct} = Pc(idx_dir_flash{ct}, :);
    Pc_dir_mean(ct, :) = mean(Pc_dir{ct}, 1);
    Pc_dir_ste(ct, :) = std(Pc_dir{ct}, [], 1)/sqrt(length(idx_dir_flash{ct}));
end

figure(2)
subplot(2,2,ds)

for ct = 1:4
    errorbar(log10(Irel), Pc_dir_mean(ct, :), Pc_dir_ste(ct, :), 'color', color(ct));
    hold on
end
legend('superior', 'Anterior', 'inferior', 'posterior')
xlabel('log(R*/rod)')
ylabel('probability')
title('wash')





%% raster plot

a = ceil((length(trigger_set_i)+1)/2);
for cc =1:length(ds_id_flash);
    ts = 1;
    b = 1;
    trial_n = length(ds_dark_raster{1});
    FigHandle = figure(10);
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

    print_close(10, [24 12], [num2str(ds_id_flash(cc)) '_' num2str(ds)])

end

end


for ds = 1:3
clear trigger_set_i
ts = 1; 
trigger_set_i{1} = [] ;
for a=1:length(dataflash{ds}.triggers) ; % for each trigger
    trigger_set_i{ts} = [trigger_set_i{ts},a] ; % put it in a set   
    if a<length(dataflash{ds}.triggers) ; % if its not the last trigger
        if sum(abs(dataflash{ds}.triggers(a+1)-dataflash{ds}.triggers(a)-dataflash{ds}.DfParams.interFlashInt)<interFlashIntVar)>0 ; % next trigger is has the right interval
            ts = ts; % keep it in the set
        else
            ts = ts+1; % put it in a new set
            trigger_set_i{ts} = [] ;
        end
    end
end

i = 1;
while(i<=length(trigger_set_i))
    if length(trigger_set_i{i}) < 50
        trigger_set_i(i) = [];
    else
        i = i+1;
    end
end

%%
% interval = dataflash{ds}.DfParams.interFlashInt;
interval = 3;
ds_idx_flash = get_cell_indices(dataflash{ds}, ds_id_flash);
bin_size = 0.01; 
start = 0;
XX = start+bin_size/2:bin_size:dataflash{ds}.DfParams.interFlashInt-bin_size/2;
load('DS161208.mat', 'NdfCalibration')
% NdfCalibration(2, :) = NdfCalibration(2, :)/10;
Irel = (dataflash{ds}.DfParams.Ftime/1000).*NdfCalibration(2,dataflash{ds}.DfParams.NDF+1) ;
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
        raster_1cell{ts} = get_raster(dataflash{ds}.spikes{ds_idx_flash(cc)}, ...
            dataflash{ds}.triggers(trigger_set_i{ts}), 'start', start, 'stop', ...
            interval, 'plot', 0);
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

d = 200;
ds_dark_fr = zeros(length(ds_id_flash), 1);
for cc = 1:length(ds_id_flash)
    ds_dark_fr(cc) = sum(dataflash{ds}.spikes{ds_idx_flash(cc)} < d)/d;
end
ds_dark_fr_all(:, ds) = ds_dark_fr;
%%
trigger = 1:3:180;
ds_dark_raster = cell(length(ds_id_flash), 1);
ds_dark_raster_all = cell(length(ds_id_flash), 1);
ds_dark_hist_trial = cell(length(ds_id_flash), 1);
ds_dark_hist_mean = cell(length(ds_id_flash), 1);
for cc = 1:length(ds_id_flash)
    ds_dark_raster{cc} = get_raster(dataflash{ds}.spikes{ds_idx_flash(cc)}, ...
        trigger, 'stop', interval, ...
        'plot', 0);
    for t = 1:length(ds_dark_raster{cc})
        ds_dark_hist_trial{cc}{t} = hist(ds_dark_raster{cc}{t}, XX);
    end
    ds_dark_raster_all{cc} = sort(cell2mat(ds_dark_raster{cc}));
    ds_dark_hist_mean{cc} = hist(ds_dark_raster_all{cc}, XX)/length(trigger);
end
DS_dark_raster_all(:, ds) = ds_dark_raster_all;
%% sensitivity
Irel = (dataflash{ds}.DfParams.Ftime/1000).*NdfCalibration(2,dataflash{ds}.DfParams.NDF+1) ;
[Irel, i] = sort(Irel);
trial_n = 50;

window = 0.01;
bin_n = window/bin_size;
sensitivity{ds} = zeros(length(ds_id_flash), length(trigger_set_i));
% sensitivity{ds} = zeros(length(ds_id_flash), 11);
for cc = 1:length(ds_id_flash)
    for ts = 1:length(trigger_set_i)
%         sensitivity{ds}(cc, ts) = max(conv(ds_flash_hist_mean{cc}{ts}, ones(bin_n, 1), 'valid')) - ds_dark_fr(cc) * window;
        sensitivity{ds}(cc, ts) = max(conv(ds_flash_hist_mean{cc}{ts}, ones(bin_n, 1), 'valid')) - max(conv(ds_dark_hist_mean{cc}, ones(bin_n, 1), 'valid'));
        sensitivity_norm{ds} = sensitivity{ds}./repmat(max(sensitivity{ds}, [], 2), 1, length(trigger_set_i));
%         sensitivity_norm{ds} = sensitivity{ds}./repmat(max(sensitivity{ds}, [], 2), 1, 11);
    end
end


color = 'rbgk';
figure
for ct = 1:4
    plot(log10(Irel), sensitivity_norm{ds}(idx_dir_flash{ct}, :)', 'color', color(ct))
    hold on
end
xlabel('log(R*/rod)')
ylabel('response')

% compare across cell type
for ct = 1:4
    sensitivity_dir{ds}{ct} = sensitivity_norm{ds}(idx_dir_flash{ct}, :);
    sensitivity_dir_mean{ds}(ct, :) = mean(sensitivity_dir{ds}{ct});
    sensitivity_dir_ste{ds}(ct, :) = std(sensitivity_dir{ds}{ct})/sqrt(length(idx_dir_flash{ct}));
end

figure

for ct = 1:4
    errorbar(log10(Irel), sensitivity_dir_mean{ds}(ct, :), sensitivity_dir_ste{ds}(ct, :), 'color', color(ct));
    hold on
end
legend('superior', 'Anterior', 'inferior', 'posterior')
xlabel('log(R*/rod)')
ylabel('response')

end
