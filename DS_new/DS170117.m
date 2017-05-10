%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load drifting grating data
datadg = load_data('/Volumes/lab/analysis/2017-01-17-0/data006-sorted/data006-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2017-01-17-0/stimuli/s06.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

% identify DS cells
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [1 3]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

% data000:3794s
datarun = load_data('/Volumes/lab/analysis/2017-01-17-0/data000-001-map/data000-001-map', opt);
dataflash = split_datarun(datarun, 3794);
dataflash{1}.DfParams.NDF =   [5,5,4,4,4,3,3,3,3,3,2,2,2,1,1] ; % on filter turret 
dataflash{1}.DfParams.Ftime = [2,8,2,4,8,2,3,4,6,8,2,4,8,2,8] ; % ms
dataflash{1}.DfParams.interFlashInt = [3] ; % sec

dataflash{2}.DfParams.NDF =   [5,5,5,4,4,4,3,3,3,3,3,2,2,2,1,1] ; % on filter turret 
dataflash{2}.DfParams.Ftime = [2,4,8,2,4,8,2,3,4,6,8,2,4,8,2,8] ; % ms
dataflash{2}.DfParams.interFlashInt = [3] ; % sec

datarun = load_data('/Volumes/lab/analysis/2017-01-17-0/data000-001-map-all/data000-001-map-all', opt);
DataFlashAll = split_datarun(datarun, 3794);
DataFlashAll{1}.DfParams.NDF =   [5,5,4,4,4,3,3,3,3,3,2,2,2,1,1] ; % on filter turret 
DataFlashAll{1}.DfParams.Ftime = [2,8,2,4,8,2,3,4,6,8,2,4,8,2,8] ; % ms
DataFlashAll{1}.DfParams.interFlashInt = [3] ; % sec

DataFlashAll{2}.DfParams.NDF =   [5,5,5,4,4,4,3,3,3,3,3,2,2,2,1,1] ; % on filter turret 
DataFlashAll{2}.DfParams.Ftime = [2,4,8,2,4,8,2,3,4,6,8,2,4,8,2,8] ; % ms
DataFlashAll{2}.DfParams.interFlashInt = [3] ; % sec

% moving flashing squares
% data002:2400s
datarun = load_data('/Volumes/lab/analysis/2017-01-17-0/data002-003-map/data002-003-map', opt);
datafs = split_datarun(datarun, 2400);
datafs{1}.names.stimulus_path = '/Volumes/lab/analysis/2017-01-17-0/stimuli/s02.mat';
datafs{1} = load_stim_mfs(datafs{1});
datafs{2}.names.stimulus_path = '/Volumes/lab/analysis/2017-01-17-0/stimuli/s03.mat';
datafs{2} = load_stim_mfs(datafs{2});

load('DS170117.mat')
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

%% parameters
interFlashIntVar = 0.005; % (sec) expected precision of trigger intervals

%% find sets of flash stimuli
clear trigger_set_i
ds = 1;
ts = 1 ; 
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
    if length(trigger_set_i{i}) < 60
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
NdfCalibration(2, :) = NdfCalibration(2, :)/10;
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

trigger = 1:3:300;
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

% exclude 2 bad data points
Irel(7:8) = [];
Pc(:, 7:8) = [];


color = 'rbgk';
figure
for ct = 1:4
    plot(log10(Irel), Pc(idx_dir_flash{ct}, :)', 'color', color(ct))
    hold on
end
xlabel('log(R*/rod)')
ylabel('probability')
title('control')
% compare across cell type
for ct = 1:4
    Pc_dir{ct} = Pc(idx_dir_flash{ct}, :);
    Pc_dir_mean(ct, :) = mean(Pc_dir{ct}, 1);
    Pc_dir_ste(ct, :) = std(Pc_dir{ct}, [], 1)/sqrt(length(idx_dir_flash{ct}));
end

figure
for ct = 1:4
    errorbar(log10(Irel), Pc_dir_mean(ct, :), Pc_dir_ste(ct, :), 'color', color(ct));
    hold on
end
legend('superior', 'Anterior', 'inferior', 'posterior')
xlabel('log(R*/rod)')
ylabel('probability')
title('control')
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

%% find sets of flash stimuli
clear trigger_set_i
ds = 2;
ts = 1 ; 
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
    if length(trigger_set_i{i}) < 60
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
NdfCalibration(2, :) = NdfCalibration(2, :)/10;
Irel = (dataflash{ds}.DfParams.Ftime/1000).*NdfCalibration(2,dataflash{ds}.DfParams.NDF+1) ;
Irel(2) = [];
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

trigger = 1:3:250;
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

%% 2 alternative forced choice test
% use 400-580 second (roughly) as dark trials
clear Pc Pc_dir Pc_dir_mean Pc_dir_ste

window = 1;
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
    plot(log10(Irel), Pc(idx_dir_flash{ct}, :)', 'color', color(ct))
    hold on
end
xlabel('log(R*/rod)')
ylabel('probability')
title('PTX')

% compare across cell type
for ct = 1:4
    Pc_dir{ct} = Pc(idx_dir_flash{ct}, :);
    Pc_dir_mean(ct, :) = mean(Pc_dir{ct}, 1);
    Pc_dir_ste(ct, :) = std(Pc_dir{ct}, [], 1)/sqrt(length(idx_dir_flash{ct}));
end

figure
for ct = 1:4
    errorbar(log10(Irel), Pc_dir_mean(ct, :), Pc_dir_ste(ct, :), 'color', color(ct));
    hold on
end
legend('superior', 'Anterior', 'inferior', 'posterior')
xlabel('log(R*/rod)')
ylabel('probability')
title('PTX')
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

%% plot rfs
for ds = 1:2
    cell_type = {'superior', 'anterior', 'inferior', 'posterior'};
    field_width = 17; field_height = 17;
    subregion = 0;
    stop = 0.5; %second
    %     fs_raster{i} = get_fs_raster(datafs{i}, ds_id, 'stop', 0.5);
    fs_raster = get_fs_raster(datafs{ds}, ds_id);
    for cc = 1:length(ds_id)
        if fs_idx(cc)
            fs_raster{cc} = [];
        end
    end
    fs_spike = fs_get_spike(fs_raster);
    [rf_all{ds}, rf_std] = get_fs_rf(fs_spike, field_width, field_height,subregion);
end

%%
for cc = 1:length(ds_id)
    figure(1)
    set(gcf, 'Position', [1 1 900 600])
    id = ds_id(cc);
    if ~isempty(rf_all{1}{cc})
%             rf = padarray(rf_all{i}{cc},[7,7]);
        rf = rf_all{2}{cc};

        subplot(2,3,1)
        imagesc(sum(rf,3))
        colormap gray
        axis image
        axis off
        title('control')

        subplot(2,3,2)
        imagesc(rf(:,:,1))
        colormap gray
        axis image
        axis off
        title('on')

        subplot(2,3,3)
        imagesc(rf(:,:,2))
        colormap gray
        axis image
        axis off
        title('off')
        
        rf = rf_all{1}{cc};

        subplot(2,3,4)
        imagesc(sum(rf,3))
        colormap gray
        axis image
        axis off
        title('PTX')

        subplot(2,3,5)
        imagesc(rf(:,:,1))
        colormap gray
        axis image
        axis off
        title('on')

        subplot(2,3,6)
        imagesc(rf(:,:,2))
        colormap gray
        axis image
        axis off
        title('off')


    end
    print_close(1,[15 10],num2str(id))
end

%% fit RF with Gaussian                
PixelArea = (30*4)^2/10^6;
% fit and compute rf area
for ds = 1:2
    for dir = 1:4
        clear rf_area_temp
        for onoff = 1:2
            rf_area_temp = [];
            for cc = 1:length(id_dir{dir})
                if ~fs_idx(idx_dir{dir}(cc))
                    data = rf_all{ds}{idx_dir{dir}(cc)}(:, :, onoff);
    %                     data = rf_wt{dir}{cc}{onoff};
                    if sum(sum(data > mean(data(:))+5*std(data(:))))>0
    %                         figure(100)
    %                         imagesc(data)
    %                         colormap gray
    %                         pause

                        params = fit_2d_gaussian(data);
    %                     Gaussian_params{ll}{dir}{cc}{onoff} = params;
                        rf_area_temp = [rf_area_temp params.xd * params.yd * pi * PixelArea];
                    end
                end
            end
            rf_area{ds}{dir}{onoff} = rf_area_temp;
        end
    end
end

% exclude outliers
stdn = 100;
for ds = 1:2
    for dir = 1:4
        for onoff = 1:2
            notdone = 1;
            rf_area_temp = rf_area{ds}{dir}{onoff};
            while notdone
                a = length(rf_area_temp);
                rf_area_temp(rf_area_temp > std(rf_area_temp)*stdn + mean(rf_area_temp)) = [];
                b = length(rf_area_temp);
                if a == b
                    notdone = 0;
                    rf_area_clean{ds}{dir}{onoff} = rf_area_temp;
                end
            end
            rf_area_clean_mean{ds}{onoff}(dir) = mean(rf_area_clean{ds}{dir}{onoff});
            rf_area_clean_ste{ds}{onoff}(dir) = std(rf_area_clean{ds}{dir}{onoff})/sqrt(length(rf_area_clean{ds}{dir}{onoff}));
        end
    end
end

% plot 
color = 'brgkc';
figure
for onoff = 1:2
    subplot(1,2,onoff)
    for dir = 1:4
        ll = 1;
        n = length(rf_area_clean{ll}{dir}{onoff});
        h{ll} = plot((dir-1)*5+ll*ones(n,1), rf_area_clean{ll}{dir}{onoff}, [color(ll) 'o']);
        hold on
        ll = 2;
        n = length(rf_area_clean{ll}{dir}{onoff});
        h{ll} = plot((dir-1)*5+ll*ones(n,1), rf_area_clean{ll}{dir}{onoff}, [color(ll) 'o']);
    end
%     set(gca, 'yscale', 'log')
    legend([h{1}(1), h{2}(1)], 'PTX', 'control')
    if onoff == 1
        title('ON')
    else
        title('OFF')
    end
    ylabel('RF area (mm^2)')
    set(gca, 'xtick', [])
%         ylim([0 0.3])

end

%% dim flash responses of all recorded GCs
bin_size = 0.02; 
start = 0;
ds= 1;
XX = start+bin_size/2:bin_size:DataFlashAll{ds}.DfParams.interFlashInt-bin_size/2;
load('DS161208.mat', 'NdfCalibration')
NdfCalibration(2, :) = NdfCalibration(2, :)/10;
Irel = (DataFlashAll{ds}.DfParams.Ftime/1000).*NdfCalibration(2,DataFlashAll{ds}.DfParams.NDF+1) ;
[Irel, i] = sort(Irel);

flash_raster = cell(length(DataFlashAll{ds}.cell_ids), 1);
flash_hist_trial = cell(length(DataFlashAll{ds}.cell_ids), 1);
flash_hist_mean = cell(length(DataFlashAll{ds}.cell_ids), 1);
for cc = 1:length(DataFlashAll{ds}.cell_ids)
    raster_1cell = cell(length(trigger_set_i), 1);
    raster_1cell_all = cell(length(trigger_set_i), 1);
    hist_1cell = cell(length(trigger_set_i), 1);
    hist_1cell_all = cell(length(trigger_set_i), 1);
    for ts = 1:length(trigger_set_i)
        raster_1cell{ts} = get_raster(DataFlashAll{ds}.spikes{cc}, ...
            DataFlashAll{ds}.triggers(trigger_set_i{ts}), 'start', start, 'stop', ...
            DataFlashAll{ds}.DfParams.interFlashInt, 'plot', 0);
        for t = 1:length(raster_1cell{ts})
            hist_1cell{ts}{t} = hist(raster_1cell{ts}{t}, XX);
        end
        raster_1cell_all{ts} = sort(cell2mat(raster_1cell{ts}));
        hist_1cell_all{ts} = hist(raster_1cell_all{ts}, XX)/length(trigger_set_i{ts});
    end
    flash_raster{cc} = raster_1cell(i, :);
    flash_hist_trial{cc} = hist_1cell(i, :);
    flash_hist_mean{cc} = hist_1cell_all(i, :);
end

%% get raster and psth for dark trials
% use 400-580 second (roughly) as dark trials

trigger = 1:3:300;
dark_raster = cell(length(DataFlashAll{ds}.cell_ids), 1);
dark_raster_all = cell(length(DataFlashAll{ds}.cell_ids), 1);
dark_hist_trial = cell(length(DataFlashAll{ds}.cell_ids), 1);
dark_hist_mean = cell(length(DataFlashAll{ds}.cell_ids), 1);
for cc = 1:length(DataFlashAll{ds}.cell_ids)
    dark_raster{cc} = get_raster(DataFlashAll{ds}.spikes{cc}, ...
        trigger, 'stop', DataFlashAll{ds}.DfParams.interFlashInt, ...
        'plot', 0);
    for t = 1:length(dark_raster{cc})
        dark_hist_trial{cc}{t} = hist(dark_raster{cc}{t}, XX);
    end
    dark_raster_all{cc} = sort(cell2mat(dark_raster{cc}));
    dark_hist_mean{cc} = hist(dark_raster_all{cc}, XX)/length(trigger);
end

%% 2 alternative forced choice test
% use 400-580 second (roughly) as dark trials
clear Pc Pc_dir Pc_dir_mean Pc_dir_ste
Irel = (DataFlashAll{ds}.DfParams.Ftime/1000).*NdfCalibration(2,DataFlashAll{ds}.DfParams.NDF+1) ;
[Irel, i] = sort(Irel);

window = 1;
bin_n = window/bin_size;

Pc = zeros(length(ds_id_flash), length(trigger_set_i));
for cc = 1:length(ds_id_flash)
    for ts = 1:length(trigger_set_i)
%         trial_n = length(ds_flash_raster{1}{ts});
        trial_n = 60;
        corr_flash = zeros(trial_n, 1);
        corr_dark = zeros(trial_n, 1);
        temp = dark_hist_trial{cc}(1:trial_n);
        ds_dark_hist_sum = sum(cell2mat(temp'));
        ds_flash_hist_sum = sum(cell2mat(flash_hist_trial{cc}{ts}'));
        for t = 1:trial_n
            template_flash = ds_flash_hist_sum - flash_hist_trial{cc}{ts}{t};
            template_flash = template_flash/norm(template_flash);
            template_flash(isnan(template_flash)) = 0;
            template_dark = ds_dark_hist_sum - dark_hist_trial{cc}{t};
            template_dark = template_dark/norm(template_dark);
            template_dark(isnan(template_dark)) = 0;
            DV = template_flash - template_dark;
            corr_flash(t) = flash_hist_trial{cc}{ts}{t}(1:bin_n) * DV(1:bin_n)';
            corr_dark(t) = dark_hist_trial{cc}{t}(1:bin_n) * DV(1:bin_n)';
        end
        Pc(cc, ts) = (sum(corr_flash > 0) + sum(corr_dark <= 0))/(trial_n*2);
    end
end

% exclude 2 bad data points
Irel(7:8) = [];
Pc(:, 7:8) = [];


color = 'rbgk';
figure
for ct = 1:4
    plot(log10(Irel), Pc(idx_dir_flash{ct}, :)', 'color', color(ct))
    hold on
end
xlabel('log(R*/rod)')
ylabel('probability')
title('control')
% compare across cell type
for ct = 1:4
    Pc_dir{ct} = Pc(idx_dir_flash{ct}, :);
    Pc_dir_mean(ct, :) = mean(Pc_dir{ct}, 1);
    Pc_dir_ste(ct, :) = std(Pc_dir{ct}, [], 1)/sqrt(length(idx_dir_flash{ct}));
end

figure
for ct = 1:4
    errorbar(log10(Irel), Pc_dir_mean(ct, :), Pc_dir_ste(ct, :), 'color', color(ct));
    hold on
end
legend('superior', 'Anterior', 'inferior', 'posterior')
xlabel('log(R*/rod)')
ylabel('probability')
title('control')
%% raster plot

a = ceil((length(trigger_set_i)+1)/2);
for cc =1:length(DataFlashAll{ds}.cell_ids);
    ts = 1;
    b = 1;
    trial_n = length(dark_raster{1});
    FigHandle = figure;
    set(FigHandle, 'Position', [0, 0, 1920, 1080]);
    while ts <= length(trigger_set_i)+1
        if  ts == a+1
            b = 3;
        end
        subplot(a, 4, b)
        if ts > 1
            trial_n = length(flash_raster{cc}{ts-1});
        end
        for j = 1:trial_n
            if ts == 1
                SpikeTime = dark_raster{cc}{j};
            else
                SpikeTime = flash_raster{cc}{ts-1}{j};
            end
            SpikeTime = SpikeTime';
            X = [SpikeTime; SpikeTime];
            Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
            line(X, Y, 'color', 'b');
            axis([0, 3, 0, trial_n]);
            hold on
        end
        if b == 1
            title(num2str(DataFlashAll{ds}.cell_ids(cc)))
        end
        b = b+1;
        subplot(a, 4, b)
        if ts > 1
            bar(XX, flash_hist_mean{cc}{ts-1}, 1)
        else
            bar(XX, dark_hist_mean{cc}, 1)
        end
        xlim([0 3])

        b = b+3;
        ts = ts+1;
    end

    print_close(1, [24 12], [num2str(DataFlashAll{ds}.cell_ids(cc)) '_ctrl'])

end
