%% data 004 NDF 4

wn_datapath = '/Volumes/lab/Analysis/2013-07-01-0/data008/data008';
dg_datapath = '/Analysis/xyao/2013-07-01-0/data004/data004';

opts = struct('load_neurons', true, 'load_ei', true, 'load_params', true);

datarun{1} = load_data(wn_datapath, opts);
datarun{2} = load_data(dg_datapath, opts);

% load stimulus for drifting gratings data
datarun{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2013-07-01-0/stimuli/s04';
datarun{2} = load_stim(datarun{2});

% define some parameters of the stimulus
mic_contrast = 0.36;
trial_duration = 8; % sec
start_time =0; % start at the beginning of each trial
bin_rate = 10000; % how to bin the data in units of Hz.

%% data 007 NDF 0

wn_datapath = '/Volumes/lab/Analysis/2013-07-01-0/data008/data008';
dg_datapath = '/Analysis/xyao/2013-07-01-0/data007/data007';

opts = struct('load_neurons', true, 'load_ei', true, 'load_params', true);

datarun{1} = load_data(wn_datapath, opts);
datarun{2} = load_data(dg_datapath, opts);

% load stimulus for drifting gratings data
datarun{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2013-07-01-0/stimuli/s07';
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

% cell_type = {'on transient', 'on brisk transient', 'off brisk transient', ...
%     'off brisk transient large', 'off sustained'};

cell_type = {'on transient', 'on brisk transient', 'on slow transient', 'off brisk transient', 'off brisk transient2'...
    'off brisk transient large', 'off transient', 'off sustained', 'off slow'};

n = 8; % 1:length(cell_type)

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

% fit these tuning curves

sp_freq = 800./datarun{2}.stimulus.params.SPATIAL_PERIOD;


tuning_fits = zeros(length(slave_cell_indices), length(sp_freq));
fit_params = cell(length(slave_cell_indices),1);
surround_strengths = zeros(length(slave_cell_indices), 1);
figure
for rgc = 1:length(slave_cell_indices)

    [temp_tuning_fit, temp_fit_params] = fit_spatial_tuning(sp_freq, tuning_curves(rgc,:), 'verbose', false);
    tuning_fits(rgc,:) = temp_tuning_fit;
    fit_params{rgc} = temp_fit_params;
    
    % get volume of center and surround
    surround_volume = 2 * pi * temp_fit_params.surround_gain * temp_fit_params.surround_radius^2;
    center_volume = 2 * pi * temp_fit_params.center_gain * temp_fit_params.center_radius^2;
    surround_strengths(rgc) = surround_volume / center_volume;

        height = 2;
        width = 2;
        id = datarun{2}.cell_ids(slave_cell_indices(rgc));
        %side_square = ceil(sqrt(length(slave_cell_indices)));
        subplot(height, width, rgc)
        semilogx(sp_freq, tuning_curves(rgc,:), 'ko', sp_freq, temp_tuning_fit, 'r-')
        title([num2str(id),'  ' ,num2str(surround_strengths(rgc))])
        axis([0.3 300 0 1]) 
    
end

subplot(height, width, rgc+1)
hist(surround_strengths, [0:0.1:1.5])
axis([0 1.5 0 5])

% surround{n} = surround_strengths;




%% Map cells

cell_type = {'on transient', 'on brisk transient', 'off brisk transient', ...
    'off brisk transient large', 'off sustained'};

surround_strengths = zeros(1, length(cell_type));
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

[tuning_fits, fit_params] = fit_spatial_tuning(sp_freq, tuning_curves, 'verbose', false);
  
    
% get volume of center and surround
surround_volume = 2 * pi * fit_params.surround_gain * fit_params.surround_radius^2;
center_volume = 2 * pi * fit_params.center_gain * fit_params.center_radius^2;
surround_strengths(i) = surround_volume / center_volume;

semilogx(sp_freq, tuning_curves, 'ko', sp_freq, tuning_fits, 'r-')
title(cell_type(i))
axis([0.3 300 0 1]) 
    
end

%%
figure
x = [1:5]; y = [surround_strengths{3}.ndf4; surround_strengths{3}.ndf0; surroud_strengths_nd4_temp; surroud_strengths_nd0_temp]'; 
bar(x,y)
set(gca,'xticklabel',cell_type)
legend('NDF 4 night', 'NDF 0 night', 'NDF 4 day', 'NDF 0 day', 'location', 'northwest')
ylabel('surround center volume ratio')
title('contrast 0.24')