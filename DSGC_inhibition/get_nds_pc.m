function [nds_Pc, nds_flash_raster, nds_flash_hist_mean, nds_dark_raster, nds_dark_hist_mean] = get_nds_pc

cd /Users/xyao/matlab/code-private/DS_new/
path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2016-03-04-0/';
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

dataflash = load_data(strcat(path, 'data000-map-002/data000-map-002'), opt);
dataflash.DfParams.NDF =   [5,5,5,4,4,4,4,4,3,3,3,2,2,1,0] ; % on filter turret 
dataflash.DfParams.Ftime = [2,4,8,2,3,4,6,8,2,4,8,2,4,2,2] ; % ms
dataflash.DfParams.interFlashInt = [3] ; % sec

datawn = load_data(strcat(path, 'data002/data002'), opt);
load('DS160324.mat', 'NdfCalibration')
%% get raster and psth for stimulus trials

id_nds_flash{1} = intersect(get_cell_ids(datawn,'ON transient'), dataflash.cell_ids);
idx_nds_flash{1} = get_cell_indices(dataflash, id_nds_flash{1});
id_nds_flash{2} = intersect(get_cell_ids(datawn,'ON sustained'), dataflash.cell_ids);
idx_nds_flash{2} = get_cell_indices(dataflash, id_nds_flash{2});
id_nds_flash{3} = intersect(get_cell_ids(datawn,'OFF transient'), dataflash.cell_ids);
idx_nds_flash{3} = get_cell_indices(dataflash, id_nds_flash{3});
id_nds_flash{4} = intersect(get_cell_ids(datawn,'OFF sustained'), dataflash.cell_ids);
idx_nds_flash{4} = get_cell_indices(dataflash, id_nds_flash{4});

nds_i = 2;

% trigger_set_i{3} = trigger_set_i{3}(1:30);
% trigger_set_i{4} = trigger_set_i{4}(30:end);
% trigger_set_i{5} = trigger_set_i{5}(1:40);
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

nds_flash_raster = cell(length(id_nds_flash{nds_i}), 1);
nds_flash_hist_trial = cell(length(id_nds_flash{nds_i}), 1);
nds_flash_hist_mean = cell(length(id_nds_flash{nds_i}), 1);
for cc = 1:length(id_nds_flash{nds_i})
    raster_1cell = cell(length(trigger_set_i), 1);
    raster_1cell_all = cell(length(trigger_set_i), 1);
    hist_1cell = cell(length(trigger_set_i), 1);
    hist_1cell_all = cell(length(trigger_set_i), 1);
    for ts = 1:length(trigger_set_i)
        raster_1cell{ts} = get_raster(dataflash.spikes{idx_nds_flash{nds_i}(cc)}, ...
            dataflash.triggers(trigger_set_i{ts}), 'start', start, 'stop', ...
            dataflash.DfParams.interFlashInt, 'plot', 0);
        for t = 1:length(raster_1cell{ts})
            hist_1cell{ts}{t} = hist(raster_1cell{ts}{t}, XX);
        end
        raster_1cell_all{ts} = sort(cell2mat(raster_1cell{ts}));
        hist_1cell_all{ts} = hist(raster_1cell_all{ts}, XX)/length(trigger_set_i{ts});
    end
    nds_flash_raster{cc} = raster_1cell(i, :);
    nds_flash_hist_trial{cc} = hist_1cell(i, :);
    nds_flash_hist_mean{cc} = hist_1cell_all(i, :);
end

%% get raster and psth for dark trials
% use 400-580 second (roughly) as dark trials

trigger = 150:3:330;
nds_dark_raster = cell(length(id_nds_flash{nds_i}), 1);
nds_dark_raster_all = cell(length(id_nds_flash{nds_i}), 1);
nds_dark_hist_trial = cell(length(id_nds_flash{nds_i}), 1);
nds_dark_hist_mean = cell(length(id_nds_flash{nds_i}), 1);
for cc = 1:length(id_nds_flash{nds_i})
    nds_dark_raster{cc} = get_raster(dataflash.spikes{idx_nds_flash{nds_i}(cc)}, ...
        trigger, 'stop', dataflash.DfParams.interFlashInt, ...
        'plot', 0);
    for t = 1:length(nds_dark_raster{cc})
        nds_dark_hist_trial{cc}{t} = hist(nds_dark_raster{cc}{t}, XX);
    end
    nds_dark_raster_all{cc} = sort(cell2mat(nds_dark_raster{cc}));
    nds_dark_hist_mean{cc} = hist(nds_dark_raster_all{cc}, XX)/length(trigger);
end

%% 2 alternative forced choice test
% use 400-580 second (roughly) as dark trials
window = 2;
bin_n = window/bin_size;

nds_Pc = zeros(length(id_nds_flash{nds_i}), length(trigger_set_i));
for cc = 1:length(id_nds_flash{nds_i})
    for ts = 1:length(trigger_set_i)
        trial_n = min(length(nds_flash_raster{1}{ts}), 60);
%         trial_n = 60;
        corr_flash = zeros(trial_n, 1);
        corr_dark = zeros(trial_n, 1);
        temp = nds_dark_hist_trial{cc}(1:trial_n);
        nds_dark_hist_sum = sum(cell2mat(temp'));
        nds_flash_hist_sum = sum(cell2mat(nds_flash_hist_trial{cc}{ts}'));
        for t = 1:trial_n
            template_flash = nds_flash_hist_sum - nds_flash_hist_trial{cc}{ts}{t};
            template_flash = template_flash/norm(template_flash);
            template_flash(isnan(template_flash)) = 0;
            template_dark = nds_dark_hist_sum - nds_dark_hist_trial{cc}{t};
            template_dark = template_dark/norm(template_dark);
            template_dark(isnan(template_dark)) = 0;
            DV = template_flash - template_dark;
            corr_flash(t) = nds_flash_hist_trial{cc}{ts}{t}(1:bin_n) * DV(1:bin_n)';
            corr_dark(t) = nds_dark_hist_trial{cc}{t}(1:bin_n) * DV(1:bin_n)';
        end
%         nds_Pc(cc, ts) = (sum(corr_flash > 0) + sum(corr_dark <= 0))/(trial_n*2);
        nds_Pc(cc, ts) = (sum(corr_flash > corr_dark) + sum(corr_flash == corr_dark)/2)/trial_n;
    end
end


% color = 'rbgk';
% figure
% for i = 1:size(nds_Pc, 1)
%     plot(log10(Irel), nds_Pc(i, :), 'color', 'r')
%     hold on
% %     pause
% end
% xlabel('log(R*/rod)')
% ylabel('probability')
% 
% figure
% errorbar(log10(Irel), mean(nds_Pc), std(nds_Pc)/sqrt(size(nds_Pc, 1)));

