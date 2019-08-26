%% data001 NDF 3

wn_datapath = '/Analysis/xyao/2013-11-07-0/data000/data000';
dg_datapath = '/Analysis/xyao/2013-11-07-0/data001/data001';

opts = struct('load_neurons', true, 'load_ei', true, 'load_params', true);

datarun{1} = load_data(wn_datapath, opts);
datarun{2} = load_data(dg_datapath, opts);

% load stimulus for drifting gratings data
datarun{2}.names.stimulus_path = '/Analysis/xyao/2013-11-07-0/stimuli/s01';
datarun{2} = load_stim(datarun{2});

% define some parameters of the stimulus
mic_contrast = 0.12;
trial_duration = 8; % sec
start_time =0; % start at the beginning of each trial
bin_rate = 10000; % how to bin the data in units of Hz.

%% data005 NDF 0

wn_datapath = '/Analysis/xyao/2013-11-07-0/data002/data002';
dg_datapath = '/Analysis/xyao/2013-11-07-0/data005/data005';

opts = struct('load_neurons', true, 'load_ei', true, 'load_params', true);

datarun{1} = load_data(wn_datapath, opts);
datarun{2} = load_data(dg_datapath, opts);

% load stimulus for drifting gratings data
datarun{2}.names.stimulus_path = '/Analysis/xyao/2013-11-07-0/stimuli/s05';
datarun{2} = load_stim(datarun{2});

% define some parameters of the stimulus
mic_contrast = 0.12;
trial_duration = 8; % sec
start_time =0; % start at the beginning of each trial
bin_rate = 10000; % how to bin the data in units of Hz.

%%
% CALCULATIONS BEGIN HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Map cells

cell_type = {'ON brisk transient', 'ON transient', 'ON sustained', 'OFF brisk transient', ...
    'OFF transient', 'OFF sustained', 'OFF slow'};



n = 7;

master_cell_indices = get_cell_indices(datarun{1}, cell_type{n});


[cell_list, failed_cells] = map_ei(datarun{1}, datarun{2}, 'master_cell_type', cell_type{n},...
            'verbose', true, 'corr_threshold', 0.95);

mapped_cntr = 0;
clear slave_IDs
for rgc = 1:length(master_cell_indices)
    if ~isempty(cell_list{master_cell_indices(rgc)})
      mapped_cntr = mapped_cntr + 1;
      slave_IDs(mapped_cntr) = cell_list{master_cell_indices(rgc)};
    end
end 
   
% analyze grating response for the successfully mapped cells
slave_cell_indices = get_cell_indices(datarun{2}, slave_IDs)

% organize the triggers to index to the correct spike sequences
spat_periods = datarun{2}.stimulus.params.SPATIAL_PERIOD;
num_trials = length(datarun{2}.stimulus.trials);
condition_triggers = zeros(datarun{2}.stimulus.repetitions, length(spat_periods));

% get the trigger indices that correspond to the above conditions
for trl = 1:num_trials
    tmp_trial = datarun{2}.stimulus.trials(trl);
    if tmp_trial.RGB(1) == mic_contrast
        temp_spat_per = tmp_trial.SPATIAL_PERIOD;
        spat_per_index = find(datarun{2}.stimulus.params.SPATIAL_PERIOD == temp_spat_per);
        rep = ceil(trl./ length(datarun{2}.stimulus.combinations));
        condition_triggers(rep, spat_per_index) = trl;
    else
        continue
    end
end

% get tuning curves for each rgc

tuning_curves = zeros(length(slave_cell_indices), length(datarun{2}.stimulus.params.SPATIAL_PERIOD));

for rgc = 1:length(slave_cell_indices)
    
    % initialize some variables and keep track of stuff
    temp_cell_index = slave_cell_indices(rgc);
    tmp_spikes = datarun{2}.spikes{temp_cell_index};
    sig_length = trial_duration*bin_rate;
    hist_spikes = zeros(length(spat_periods),sig_length);

    % get spike rasters and sum them over trials
    for spp = 1:length(spat_periods)
        spp_raster = get_raster(tmp_spikes, datarun{2}.stimulus.triggers(condition_triggers(:,spp)), 'stop', trial_duration, 'start', start_time,'plot', false);
        for rep = 1:length(spp_raster)
            tmp_times = spp_raster{rep};
            if isempty(tmp_times)
                break
            end
            tmp_times = floor(tmp_times * bin_rate);
            if tmp_times(1) == 0
                tmp_times = tmp_times(2:end);
            end
            tmp_binned_spikes = zeros(1,trial_duration*bin_rate);
            tmp_binned_spikes(tmp_times) = 1;
            hist_spikes(spp, :) = hist_spikes(spp,:) + tmp_binned_spikes;
        end            
    end

    % calculates fourier transform and gets power at fundamental frequency (f1)
    NFFT = 2^nextpow2(sig_length);
    f = bin_rate/2*linspace(0,1,NFFT/2+1);
    f1 = 2; %Hz
    f_diff = f - 2;
    [~,f1_index] = min(abs(f_diff));
    clear tmp_fft fft_spikes
    for spp = 1:length(spat_periods)
        tmp_fft = fft(hist_spikes(spp,:), NFFT)./ sig_length;
        fft_spikes(spp,:) = 2*abs(tmp_fft(1:NFFT/2+1));
        fund_power(spp) = sum(fft_spikes(spp,f1_index:f1_index+2));
    end

    % stores info for this cell into the matrix tuning curves
    tuning_curves(rgc,:) = fund_power ./ max(fund_power);
end

% plot tuning curve portraits
side_square = ceil(sqrt(length(slave_cell_indices)));
figure
for rgc = 1:length(slave_cell_indices)
    subplot(side_square, side_square, rgc)
    semilogx(800./datarun{2}.stimulus.params.SPATIAL_PERIOD, tuning_curves(rgc,:), '-o')
    axis([0.3 300 0 1]) 
end

%% fit these tuning curves

sp_freq = 800./datarun{2}.stimulus.params.SPATIAL_PERIOD;

% HACK Should be removed eventually
sp_indices = [1 3 5 6 7 8 9 10 11 12];

cell_numb = length(slave_cell_indices)+1;
if cell_numb <= ceil(sqrt(cell_numb))*floor(sqrt(cell_numb))
    dimx = floor(sqrt(cell_numb));
    dimy = ceil(sqrt(cell_numb));
else
    dimx = ceil(sqrt(cell_numb));
    dimy = dimx;
end


tuning_fits = zeros(length(slave_cell_indices), length(sp_freq(sp_indices)));
fit_params = cell(length(slave_cell_indices),1);
surround_strengths = zeros(length(slave_cell_indices), 1);
figure
for rgc = 1:length(slave_cell_indices)

    [temp_tuning_fit, temp_fit_params] = fit_spatial_tuning(sp_freq(sp_indices), tuning_curves(rgc,sp_indices), 'verbose', false);
    temp_fit_params
    tuning_fits(rgc,:) = temp_tuning_fit;
    fit_params{rgc} = temp_fit_params;
    
    % get volume of center and surround
    surround_volume = 2 * pi * temp_fit_params.surround_gain * temp_fit_params.surround_radius^2;
    center_volume = 2 * pi * temp_fit_params.center_gain * temp_fit_params.center_radius^2;
    surround_strengths(rgc) = surround_volume / center_volume;

        id = datarun{2}.cell_ids(slave_cell_indices(rgc));
        %side_square = ceil(sqrt(length(slave_cell_indices)));
        subplot(dimx, dimy, rgc)
        semilogx(sp_freq, tuning_curves(rgc,:), 'ko', sp_freq(sp_indices), temp_tuning_fit, 'r-')
        title([num2str(id),'  ' ,num2str(surround_strengths(rgc))])
        axis([0.3 300 0 1]) 
    
end

subplot(dimx, dimy, rgc+1)
hist(surround_strengths, [0:0.1:1.5])
axis([0 1.5 0 5])
title(cell_type{n})






%% Map cells

cell_type = {'ON brisk transient', 'ON transient', 'ON sustained', 'OFF brisk transient', ...
    'OFF transient', 'OFF sustained', 'OFF slow'};

surround_strengths = zeros(1, length(cell_type));

sp_indices = cell(5, 1);
sp_indices{1} = [1 2 3 4 5 6 7 8 9 10 11 12];
sp_indices{2} = [1 2 3 4 5 6 7 8 9 10 11 12];
sp_indices{3} = [1 2 3 4 5 6 7 8 9 10 11 12];
sp_indices{4} = [1 2 5 6 8 9 10 11 12];
sp_indices{5} = [1 3 5 6 7 8 9 10 11 12];
sp_indices{6} = [1 3 5 7 8 9 10 11 12];
sp_indices{7} = [1 2 3 4 5 6 7 8 9 10 11 12];

% sp_indices{1} = [1 2 4 5 6 9 10];
% sp_indices{2} = [1 3 5 8 9 10];
% sp_indices{3} = [1 3 5 7 8 9 10];
% sp_indices{4} = [1 3 5 6 7 8 9 10];
% sp_indices{5} = [1 2 3 4 5 6 7 8 9 10];


figure

for i = 1:length(cell_type)
master_cell_indices = get_cell_indices(datarun{1}, cell_type(i));

[cell_list, failed_cells] = map_ei(datarun{1}, datarun{2}, 'master_cell_type', cell_type(i),...
            'verbose', true, 'corr_threshold', 0.95);

mapped_cntr = 0;
clear slave_IDs
for rgc = 1:length(master_cell_indices)
    if ~isempty(cell_list{master_cell_indices(rgc)})
      mapped_cntr = mapped_cntr + 1;
      slave_IDs(mapped_cntr) = cell_list{master_cell_indices(rgc)};
    end
end 
   
% analyze grating response for the successfully mapped cells
slave_cell_indices = get_cell_indices(datarun{2}, slave_IDs)

% organize the triggers to index to the correct spike sequences
spat_periods = datarun{2}.stimulus.params.SPATIAL_PERIOD;
num_trials = length(datarun{2}.stimulus.trials);
condition_triggers = zeros(datarun{2}.stimulus.repetitions, length(spat_periods));

% get the trigger indices that correspond to the above conditions
for trl = 1:num_trials
    tmp_trial = datarun{2}.stimulus.trials(trl);
    if tmp_trial.RGB(1) == mic_contrast
        temp_spat_per = tmp_trial.SPATIAL_PERIOD;
        spat_per_index = find(datarun{2}.stimulus.params.SPATIAL_PERIOD == temp_spat_per);
        rep = ceil(trl./ length(datarun{2}.stimulus.combinations));
        condition_triggers(rep, spat_per_index) = trl;
    else
        continue
    end
end

% get tuning curves for each rgc

tuning_curves = zeros(length(slave_cell_indices), length(datarun{2}.stimulus.params.SPATIAL_PERIOD));

for rgc = 1:length(slave_cell_indices)
    
    % initialize some variables and keep track of stuff
    temp_cell_index = slave_cell_indices(rgc);
    tmp_spikes = datarun{2}.spikes{temp_cell_index};
    sig_length = trial_duration*bin_rate;
    hist_spikes = zeros(length(spat_periods),sig_length);

    % get spike rasters and sum them over trials
    for spp = 1:length(spat_periods)
        spp_raster = get_raster(tmp_spikes, datarun{2}.stimulus.triggers(condition_triggers(:,spp)), 'stop', trial_duration, 'start', start_time,'plot', false);
        for rep = 1:length(spp_raster)
            tmp_times = spp_raster{rep};
            if isempty(tmp_times)
                break
            end
            tmp_times = floor(tmp_times * bin_rate);
            if tmp_times(1) == 0
                tmp_times = tmp_times(2:end);
            end
            tmp_binned_spikes = zeros(1,trial_duration*bin_rate);
            tmp_binned_spikes(tmp_times) = 1;
            hist_spikes(spp, :) = hist_spikes(spp,:) + tmp_binned_spikes;
        end            
    end

    % calculates fourier transform and gets power at fundamental frequency (f1)
    NFFT = 2^nextpow2(sig_length);
    f = bin_rate/2*linspace(0,1,NFFT/2+1);
    f1 = 2; %Hz
    f_diff = f - 2;
    [~,f1_index] = min(abs(f_diff));
    clear tmp_fft fft_spikes
    for spp = 1:length(spat_periods)
        tmp_fft = fft(hist_spikes(spp,:), NFFT)./ sig_length;
        fft_spikes(spp,:) = 2*abs(tmp_fft(1:NFFT/2+1));
        fund_power(spp) = sum(fft_spikes(spp,f1_index:f1_index+2));
    end

    % stores info for this cell into the matrix tuning curves
    tuning_curves(rgc,:) = fund_power ./ max(fund_power);
end



% calculate mean tuning curve

tuning_curves_all = sum(tuning_curves);
tuning_curves = tuning_curves_all/max(tuning_curves_all);


% fit this mean tuning curve

sp_freq = 800./datarun{2}.stimulus.params.SPATIAL_PERIOD;


subplot(2, ceil(length(cell_type)/2), i)

[tuning_fits, fit_params] = fit_spatial_tuning(sp_freq(sp_indices{i}), tuning_curves(sp_indices{i}), 'verbose', false);
  
    
% get volume of center and surround
surround_volume = 2 * pi * fit_params.surround_gain * fit_params.surround_radius^2;
center_volume = 2 * pi * fit_params.center_gain * fit_params.center_radius^2;
surround_strengths(i) = surround_volume / center_volume;

semilogx(sp_freq, tuning_curves, 'ko', sp_freq(sp_indices{i}), tuning_fits, 'r-')
title(cell_type(i))
axis([0.3 300 0 1]) 
    
end

%%
figure
x = [1:5]; y = [ss1; ss2]'; 
bar(x,y)
set(gca,'xticklabel',cell_type)
legend('with melatonin', 'without melatonin', 'location', 'northwest')
ylabel('surround center volume ratio')
title('contrast 0.24')



%% scotopic level

load('dg131107.mat');
load('dg130701.mat', 'surround_ndf0_12')
load('dg130722.mat', 'surround_ndf4_24')

% compare 'on brisk transient, on transient, off brisk transient'
% idx for WT: [1 2 3]
% idx for VGAT: [1 2 4]
cell_type = {'on brisk transient', 'on transient', 'off brisk transient'};

    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 1300, 650]);
    
    s_wt = surround_ndf4_24;
    s_vgat = surround_strengths_rod{2};
    
    xtick = cell_type;
    model_series = [s_wt.surround_mean(1) s_vgat.surround_mean(1); s_wt.surround_mean(2) s_vgat.surround_mean(2); s_wt.surround_mean(3) s_vgat.surround_mean(4)];   
    model_error = [s_wt.surround_stev(1) s_vgat.surround_stev(1); s_wt.surround_stev(2) s_vgat.surround_stev(2); s_wt.surround_stev(3) s_vgat.surround_stev(4)];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('Surround Center Volume Ratio')
    legend('WT','VGAT-/-', 'location', 'northeast');
    hold on;
 
    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 

    groupwidth = min(0.8, numbars/(numbars+1.5));
    
    for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end
    title('scotopic level')


% compare 'on brisk transient, on transient, off brisk transient', 'OFF transient', 'OFF sustained', 'OFF slow'
% idx for WT: [2 1 4 7 8 9]
% idx for VGAT: [1 2 4 5 6 7]

cell_type = {'on brisk transient', 'on transient', 'off brisk transient', 'off transient', 'off sustained', 'off slow'};

    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 1300, 650]);
    
    clear s_wt s_vgat
    s_wt = surround_ndf0_12;
    s_vgat = surround_strengths_cone{1};
    
    xtick = cell_type;
    model_series = [s_wt.surround_mean(2) s_vgat.surround_mean(1); s_wt.surround_mean(1) s_vgat.surround_mean(2); s_wt.surround_mean(4) s_vgat.surround_mean(4); s_wt.surround_mean(7) s_vgat.surround_mean(5);s_wt.surround_mean(8) s_vgat.surround_mean(6);s_wt.surround_mean(9) s_vgat.surround_mean(7)];   
    model_error = [s_wt.surround_stev(2) s_vgat.surround_stev(1); s_wt.surround_stev(1) s_vgat.surround_stev(2); s_wt.surround_stev(4) s_vgat.surround_stev(4); s_wt.surround_stev(7) s_vgat.surround_stev(5);s_wt.surround_stev(8) s_vgat.surround_stev(6);s_wt.surround_stev(9) s_vgat.surround_stev(7)];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('Surround Center Volume Ratio')
    legend('WT','VGAT-/-', 'location', 'northeast');
    hold on;
 
    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 

    groupwidth = min(0.8, numbars/(numbars+1.5));
    
    for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end
    title('photopic level')



% compare 'on brisk transient, on transient, on sustained, off brisk transient', 'OFF transient', 'OFF sustained'
% idx for ndf3: [1 2 3 4 5 6]
% idx for ndf0: [1 2 3 4 5 6]

cell_type = {'on brisk transient', 'on transient', 'on sustained', 'off brisk transient', 'off transient', 'off sustained'};

    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 1300, 650]);
    
    s_ndf3 = surround_strengths_rod{2};
    s_ndf0 = surround_strengths_cone{2};
    
    xtick = cell_type;
    model_series = [s_ndf3.surround_mean(1) s_ndf0.surround_mean(1); s_ndf3.surround_mean(2) s_ndf0.surround_mean(2); s_ndf3.surround_mean(3) s_ndf0.surround_mean(3); s_ndf3.surround_mean(4) s_ndf0.surround_mean(4);s_ndf3.surround_mean(5) s_ndf0.surround_mean(5);s_ndf3.surround_mean(6) s_ndf0.surround_mean(6)];   
    model_error = [s_ndf3.surround_stev(1) s_ndf0.surround_stev(1); s_ndf3.surround_stev(2) s_ndf0.surround_stev(2); s_ndf3.surround_stev(3) s_ndf0.surround_stev(3); s_ndf3.surround_stev(4) s_ndf0.surround_stev(4);s_ndf3.surround_stev(5) s_ndf0.surround_stev(5);s_ndf3.surround_stev(6) s_ndf0.surround_stev(6)];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('Surround Center Volume Ratio')
    legend('NDF 3','NDF 0', 'location', 'northeast');
    hold on;
 
    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 

    groupwidth = min(0.8, numbars/(numbars+1.5));
    
    for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end
    title('VGAT-/-')



