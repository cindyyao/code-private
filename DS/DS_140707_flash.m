opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
datarun1 = load_data('/Analysis/xyao/2014-07-07-0/data000-map/data000-map', opt);

% find the start time of stimulus 
trigger_diff = diff(datarun1.triggers); 
start_I = [1; find(trigger_diff >3.1)+1]; % after turning on the output, the duration of each trial is 3sec.
start_time = datarun1.triggers(start_I); % convert the index of triggers to time

% delete the wrong triggers
start_time([2 14]) = []; 
start_I([2 14]) = [];


% calculate stimuli strengths
load('LED_calibration.mat')
load('DS140707.mat', 'stimuli_setting')


for j = 1:size(stimuli_setting, 1)
    I1 = find(Intensity == stimuli_setting(j, 1));
    I2 = find(Duration == stimuli_setting(j, 2));
    stimuli(j) = LED_calibration(I1, I2)*NDF_factor(stimuli_setting(j, 3));
end
[stimuli, seq] = sort(stimuli);
delta = diff(stimuli);
delta(delta ~= 0) = 1;
delta = [1 delta];
k = 0;
for j = 1:length(stimuli)
    if delta(j) == 1
        k = k+1;
        stimuli_nrpt{k} = [];
        stimuli_nrpt{k} = [stimuli_nrpt{k} j];
    else
        stimuli_nrpt{k} = [stimuli_nrpt{k} j];
    end
end
start_time = start_time(seq);
start_I = start_I(seq);
        


% get cell ids
flash_id = datarun1.cell_ids;

% get raster & psth

begin_time = 0:3:117;

bin_size = 0.0125;
tau = 4*bin_size;
tt = -3*tau:bin_size:3*tau;
filter = exp(-tt.^2/(2*tau^2));
circ = (length(filter)-1)/2;


for cn = 1:length(flash_id)
    idx = get_cell_indices(datarun1, flash_id(cn));
    spike = datarun1.spikes{idx};
    for st = 1:length(stimuli)
        if i == 3 && st < 17
            trigger = datarun1.triggers(start_I(st):start_I(st)+79);
        else
            trigger = datarun1.triggers(start_I(st):start_I(st)+39);
        end
        raster{cn}{st} = get_raster(spike, trigger, 'plot', false, 'stop', 3);
        raster_all{cn}{st} = sort(cell2mat(raster{cn}{st})); %spikes from all trials are added together
        for t = 1:length(raster{cn}{st})
            clear psth_temp
            psth_temp = hist(raster{cn}{st}{t}, bin_size/2:bin_size:3-bin_size/2);
            psth_temp = [psth_temp(end-circ+1:end) psth_temp psth_temp(1:circ)];
            psth_flash{cn}{st}{t} = conv(psth_temp, filter, 'valid');
        end
    end
%             spike_number_mean_temp = mean(cellfun(@length, raster{i}{cn}));
    for j = 1:length(stimuli_nrpt)
%             spikes_number_mean{i}{cn}(j) = mean(spike_number_mean_temp(stimuli_nrpt{i}{j}));
        raster_all_nrpt{cn}{j} = cell2mat(raster_all{cn}(stimuli_nrpt{j})');
    end
    raster_dark{cn} = get_raster(spike, begin_time, 'plot', false, 'stop', 3);
    for t = 1:length(raster_dark{cn})
        clear psth_temp
        psth_temp = hist(raster_dark{cn}{t}, bin_size/2:bin_size:3-bin_size/2);
        psth_temp = [psth_temp(end-circ+1:end) psth_temp psth_temp(1:circ)];
        psth_dark{cn}{t} = conv(psth_temp, filter, 'valid');
    end

end

%% ideal observer

for cn = 1:length(flash_id)
    for st = 1:length(stimuli)
        for t = 1:length(raster{cn}{st})
            V = sum(cell2mat(psth_flash{cn}{st}')) - psth_flash{cn}{st}{t};
            template1 = V/norm(V);
            V = sum(cell2mat(psth_dark{cn}(1:length(raster{cn}{st}))')) - psth_dark{cn}{t};
            template0 = V/norm(V);
            template = template1 - template0;
            corr_flash{cn}{st}(t) = psth_flash{cn}{st}{t}*template';
            corr_dark{cn}{st}(t) = psth_dark{cn}{t}*template';
        end
        % two alternative forced choice method 1
%             Pc_temp(st) = (sum(corr_flash{cn}{st}>0) + sum(corr_flash{cn}{st}==0)/2 + ...
%                 sum(corr_dark{cn}{st}<0) + sum(corr_dark{cn}{st}==0)/2)/ ...
%                 (length(corr_flash{cn}{st})+length(corr_dark{cn}{st}));
        % two alternative forced choice method 2
%             Pc_temp(st) = (sum(corr_flash{cn}{st}>0) + sum(corr_flash{cn}{st}==0)/2)/ ...
%                 length(corr_flash{cn}{st});
        % two alternative forced choice (0 spike --> darkness)
          Pc_temp(st) = (sum(corr_flash{cn}{st}>0) + ...
            sum(corr_dark{cn}{st}<0) + sum(corr_dark{cn}{st}==0))/ ...
            (length(corr_flash{cn}{st})+length(corr_dark{cn}{st}));
%             % two interval forced choice
%             Pc_temp(st) = sum((corr_flash{cn}{st}(1:40) - corr_dark{cn}{st}) > 0)/40;

    end
    for st = 1:length(stimuli_nrpt)
        Pc{cn}(st) = mean(Pc_temp(stimuli_nrpt{st}));
    end
end
    

