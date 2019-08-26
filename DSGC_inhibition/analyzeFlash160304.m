function [Pc_dir, xthreshold] = analyzeFlash160304
cd /Volumes/dusom_fieldlab/All_Staff/lab/Development/matlab/private/xyao/matlab/code/DS_new
path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2016-03-04-0/';
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

dataflash = load_data(strcat(path, 'data000-map/data000-map'), opt);
dataflash.DfParams.NDF =   [5,5,5,4,4,4,4,4,3,3,3,2,2,1,0] ; % on filter turret 
dataflash.DfParams.Ftime = [2,4,8,2,3,4,6,8,2,4,8,2,4,2,2] ; % ms
dataflash.DfParams.interFlashInt = [3] ; % sec

load('DS160304.mat', 'ds_id', 'id_dir', 'idx_dir', 'id_dir_on', 'idx_dir_on', 'flash_idx', 'NdfCalibration')
ds_id_flash = ds_id(flash_idx);
ds_idx_flash = get_cell_indices(dataflash, ds_id_flash);

datawn = load_data(strcat(path, 'data002-map/data002-map'), opt);

for ct = 1:4
    [id_dir_flash{ct}, idx_dir_flash{ct}] = intersect(ds_id_flash, id_dir{ct});
end
for ct = 1:3
    [id_dir_on_flash{ct}, idx_dir_on_flash{ct}] = intersect(ds_id_flash, id_dir_on{ct});
end

%% find sets of flash stimuli

interFlashIntVar = 0.005; % (sec) expected precision of trigger intervals

ts=1 ; 
trigger_set_i{1} = [] ;
for a=1:length(dataflash.triggers) ; % for each trigger
    trigger_set_i{ts} = [trigger_set_i{ts},a] ; % put it in a set   
    if a<length(dataflash.triggers) && ... % if its not the last trigger
            abs(dataflash.triggers(a+1)-dataflash.triggers(a)-dataflash.DfParams.interFlashInt) > interFlashIntVar ; 
        ts = ts + 1 ; % keep it in the set
        trigger_set_i{ts} = [] ;
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

trigger = 150:3:330;
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
window = 1; % second
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
            template_flash = template_flash(1:bin_n);
            template_flash = template_flash/norm(template_flash);
            template_flash(isnan(template_flash)) = 0;
            template_dark = ds_dark_hist_sum - ds_dark_hist_trial{cc}{t};
            template_dark = template_dark(1:bin_n);
            template_dark = template_dark/norm(template_dark);
            template_dark(isnan(template_dark)) = 0;
            DV = template_flash - template_dark;
            corr_flash(t) = ds_flash_hist_trial{cc}{ts}{t}(1:bin_n) * DV(1:bin_n)';
            corr_dark(t) = ds_dark_hist_trial{cc}{t}(1:bin_n) * DV(1:bin_n)';
        end
        Pc(cc, ts) = (sum(corr_flash > corr_dark) + sum(corr_flash == corr_dark)/2)/trial_n;
%         Pc(cc, ts) = (sum(corr_flash > 0) + sum(corr_dark <= 0))/(trial_n*2);
    end
end


color = 'rbgk';
% figure
% for ct = 1:3
%     plot(log10(Irel), Pc(idx_dir_flash{ct}, :)', 'color', color(ct))
%     hold on
% end
% xlabel('log(R*/rod)')
% ylabel('probability')

% compare across cell type
for ct = 1:3
    Pc_dir{ct} = Pc(idx_dir_flash{ct}, :);
    Pc_dir_mean(ct, :) = mean(Pc_dir{ct}, 1);
    Pc_dir_ste(ct, :) = std(Pc_dir{ct}, [], 1)/sqrt(length(idx_dir_flash{ct}));
end

% figure
% for ct = 1:3
%     errorbar(log10(Irel), Pc_dir_mean(ct, :), Pc_dir_ste(ct, :), 'color', color(ct));
%     hold on
% end
% legend('superior', 'Anterior', 'inferior')
% xlabel('log(R*/rod)')
% ylabel('probability')

%% fit
Pc_temp = Pc;
Irel_temp = Irel;

for ct = 1:4;
    idx = idx_dir_flash{ct};
    for cc = 1:length(idx)
        ydata = Pc_temp(idx(cc),:)-0.5;
        xdata = log10(Irel_temp)+4;
        [f, G] = fit_mm(xdata, ydata, 'Upper', [0.5, 10, max(xdata)]);
        fit_all{ct}{cc} = f;
        G_all{ct}{cc} = G;
    end
end

xdata = log10(Irel_temp)+4;
x = linspace(min(xdata), max(xdata), 1000);

threshold = 0.34;
for ct = 1:3
    idx = idx_dir_flash{ct};
    for cc = 1:length(idx)
        f = fit_all{ct}{cc};
        y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a);
        [~, i] = min(abs(y-threshold));
%         xthreshold{ct}(cc) = 10^(x(i)-4);
        xthreshold{ct}(cc) = x(i)-4;
        rmse{ct}(cc) = G_all{ct}{cc}.rmse;
    end
end

% save('DS160304-1.mat', 'Pc_dir_0304', '-append')
