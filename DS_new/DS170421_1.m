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

datarun = load_data('/Volumes/lab/analysis/2017-04-21-0/data000-003-map/data000-003-map', opt);
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

%% find trigger sets
clear trigger_set_i
interFlashIntVar = 0.005; % (sec) expected precision of trigger intervals
for ds = 1:3
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
    if length(trigger_set_i{i}) < 60 || length(trigger_set_i{i}) > 70
        trigger_set_i(i) = [];
    else
        i = i+1;
    end
end

%%
ds_idx_flash = get_cell_indices(dataflash{ds}, ds_id_flash);
bin_size = 0.02; 
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

%% raster plot

a = ceil((length(trigger_set_i)+1)/2);
for cc =1:length(ds_id_flash);
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

    print_close(1, [24 12], [num2str(ds_id_flash(cc)) '_' num2str(ds)])

end
end