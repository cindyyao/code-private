for ds = 1:2
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
    if length(trigger_set_i{i}) < 60
        trigger_set_i(i) = [];
    else
        i = i+1;
    end
end

%%
% interval = dataflash{ds}.DfParams.interFlashInt;
interval = 3;
ds_idx_flash = get_cell_indices(dataflash{ds}, ds_id_flash);
bin_size = 0.001; 
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

d = 300;
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
trial_n = 60;

window = 0.002;
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


%%
for drug = 1:2
    ds_dark_fr_dir{1}(:, drug) = ds_dark_fr_all(idx_dir_flash{1}, drug);
    ds_dark_fr_dir{2}(:, drug) = ds_dark_fr_all(cell2mat(idx_dir_flash(2:4)'), drug);
    for ct = 1:2
        ds_dark_fr_dir_mean(ct, drug) = mean(ds_dark_fr_dir{ct}(:, drug));
        ds_dark_fr_dir_ste(ct, drug) = std(ds_dark_fr_dir{ct}(:, drug))/sqrt(size(ds_dark_fr_dir{ct}, 1));
    end
end

figure
errorbar(ds_dark_fr_dir_mean(2, 1), ds_dark_fr_dir_mean(1, 1), ds_dark_fr_dir_ste(1, 1), 'bo')
hold on
errorbar(ds_dark_fr_dir_mean(2, 2), ds_dark_fr_dir_mean(1, 2), ds_dark_fr_dir_ste(1, 2), 'ro')
herrorbar(ds_dark_fr_dir_mean(2, 1), ds_dark_fr_dir_mean(1, 1), ds_dark_fr_dir_ste(2, 1), 'bo')
herrorbar(ds_dark_fr_dir_mean(2, 2), ds_dark_fr_dir_mean(1, 2), ds_dark_fr_dir_ste(2, 2), 'ro')
legend('control', 'SR')
plot([0 6], [0 6], 'k--')
xlabel('others (Hz)')
ylabel('superior (Hz)')
%%

drug = 1;
for dir = 1:4
    for cc = 1:length(id_dir_mb{dir})
        null_temp{dir}(cc) = (length(raster_n_sum{drug}{idx_dir_mb{dir}(cc)}{ctr})/trial_dur - bgnd_firing{dir}(drug, cc))/datamb{drug}.stimulus.repetitions;
    end
end
null{1} = null_temp{1};
null{2} = cell2mat(null_temp(2:4));


drug = 1;
ctr = 5;
for dir = 1:4
    for cc = 1:length(id_dir_mb{dir})
        prefer_temp{dir}(cc) = (length(raster_p_sum{drug}{idx_dir_mb{dir}(cc)}{ctr})/trial_dur - bgnd_firing{dir}(drug, cc))/datamb{drug}.stimulus.repetitions;
    end
end
prefer{1} = prefer_temp{1};
prefer{2} = cell2mat(prefer_temp(2:4));