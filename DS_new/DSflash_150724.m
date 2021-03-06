opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);


datadg = load_data('/Analysis/jcafaro/2015-07-24-0/data009/data009', opt);
datadg.names.stimulus_path = '/Analysis/jcafaro/2015-07-24-0/stimuli/s009';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10.001);

dataflash = load_data('/Analysis/jcafaro/2015-07-24-0/data002-map/data002-map', opt);
dataflash.DfParams.NDF =   [5,5,5,4,4,4,2,2,2,3,0,1] ; % on filter turret 
dataflash.DfParams.Ftime = [4,8,2,2,4,8,2,4,8,2,2,2] ; % ms
dataflash.DfParams.interFlashInt = [3] ; % sec

[NumSpikesCell, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);
params_idx = [1 2]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);
title('Normalized Vector Sum Amplitude')
xlabel('High speed - 480 micron/s')
ylabel('Low speed - 240 micron/s')

ds_id_flash = intersect(ds_id, dataflash.cell_ids);
ds_id_flash([1 27 28 32]) = [];
ds_idx_flash = get_cell_indices(dataflash, ds_id_flash);

%% classify DSGC into subtypes (directions)

[NumSpikesCell, StimComb] = get_spikescellstim(datadg,ds_id_flash,0);
DG = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
raster_dg = get_ds_raster(datadg, ds_id_flash);

t = 1;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(DG.U{t}, DG.V{t});
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG.U{t}, DG.V{t}, x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = I;
    id_dir{i} = ds_id_flash(idx_dir{i});
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
load /Volumes/lab/Experiments/Calibration/NdfCalibration
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
% use 180-300 second (roughly) as dark trials

trigger = 180:3:320;
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
% use 180-300 second (roughly) as dark trials
window = 2;
bin_n = window/bin_size;

Pc = zeros(length(ds_id_flash), length(trigger_set_i));
for cc = 1:length(ds_id_flash)
    for ts = 1:length(trigger_set_i)
        corr_flash = zeros(length(ds_flash_raster{1}{ts}), 1);
        corr_dark = zeros(length(ds_flash_raster{1}{ts}), 1);
        temp = ds_dark_hist_trial{cc}(1:length(ds_flash_raster{1}{ts}));
        ds_dark_hist_sum = sum(cell2mat(temp'));
        ds_flash_hist_sum = sum(cell2mat(ds_flash_hist_trial{cc}{ts}'));
        for t = 1:length(ds_flash_raster{1}{ts})
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
        Pc(cc, ts) = (sum(corr_flash > 0) + sum(corr_dark <= 0))/(length(ds_flash_raster{1}{ts})*2);
    end
end


color = 'rbgk';
figure
for ct = 1:4
    plot(log10(Irel), Pc(idx_dir{ct}, :)', 'color', color(ct))
    hold on
end

% compare across cell type
for ct = 1:4
    Pc_dir{ct} = Pc(idx_dir{ct}, :);
    Pc_dir_mean(ct, :) = mean(Pc_dir{ct});
    Pc_dir_ste(ct, :) = std(Pc_dir{ct})/sqrt(length(idx_dir{ct}));
end

figure
for ct = 1:4
    errorbar(log10(Irel), Pc_dir_mean(ct, :), Pc_dir_ste(ct, :), 'color', color(ct));
    hold on
end
legend('superior', 'A(P)', 'inferior', 'P(A)')
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

    print_close(1, [24 12], num2str(ds_id_flash(cc)))

end

