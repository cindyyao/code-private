function datarun = load_stim_matlab(datarun, varargin)

% datarun = load_stim_matlab(datarun)
% xyao
% 2014-10-11

p = inputParser;
p.addParamValue('user_defined_trigger_interval', [], @isnumeric); 
p.addParamValue('user_defined_trigger_set', [], @isnumeric);
p.addParamValue('user_defined_trigger_interval_error', 0.1, @isnumeric);

p.parse(varargin{:});
params = p.Results;

datarun.stimulus = [];
% check if stimulus file path is defined
if isfield(datarun.names, 'stimulus_path')
    stim_path = datarun.names.stimulus_path;
else
    fprintf('\t stimulus path not recognized. Please define stimulus path and try again. \n');
    return
end

load(stim_path)
datarun.stimulus.repetitions = stim_out.repeats;
datarun.stimulus.trial_list = stim_out.trial_list;

if ~isempty(p.Results.user_defined_trigger_set);
    triggers = datarun.triggers(p.Results.user_defined_trigger_set);
elseif ~isempty(p.Results.user_defined_trigger_interval)
    % define new stimulus trigger interval
    stim_interval = p.Results.user_defined_trigger_interval; % units are seconds
    % find the triggers that are closest to occuring subsequent to these intervals
    triggers = datarun.triggers;
    a = 1;
    while a<length(triggers)
        interval = triggers(a+1) - triggers(a);
        if mod(interval, stim_interval)<params.user_defined_trigger_interval_error || mod(interval, stim_interval)>stim_interval-params.user_defined_trigger_interval_error
            a = a+1;
        else
            triggers(a+1) = [];
        end
    end
else
    triggers = datarun.triggers(2:2:end);
end

datarun.stimulus.triggers = triggers';
trial_n = length(stim_out.trial_list);

if (strcmp(stim_out.type, 'MG') || strcmp(stim_out.type, 'CG'))
    datarun.stimulus.params.SPATIAL_PERIOD = stim_out.spatial_period;
    datarun.stimulus.params.TEMPORAL_PERIOD = stim_out.temporal_period;
    datarun.stimulus.params.DIRECTION = stim_out.direction;
    datarun.stimulus.params.BACK_RGB = stim_out.back_rgb;
    datarun.stimulus.params.RGB = stim_out.rgb;
    trials(1:trial_n) = struct;
    for i = 1:trial_n
        trials(i).SPATIAL_PERIOD = stim_out.trials(i).spatial_period;
        trials(i).TEMPORAL_PERIOD = stim_out.trials(i).temporal_period;
        trials(i).DIRECTION = stim_out.trials(i).direction;
        trials(i).BACK_RGB = stim_out.trials(i).back_rgb;
        trials(i).RGB = stim_out.trials(i).rgb;
    end
    datarun.stimulus.trials = trials;
    datarun.stimulus.combinations = trials(1:trial_n/stim_out.repeats);
elseif (strcmp(stim_out.type, 'MB'))
    if isfield(stim_out,'x_start')
        datarun.stimulus.x_start = stim_out.x_start;
        datarun.stimulus.x_end = stim_out.x_end;
        datarun.stimulus.y_start = stim_out.y_start;
        datarun.stimulus.y_end = stim_out.y_end;
    end
    datarun.stimulus.params.BAR_WIDTH = stim_out.bar_width;
    datarun.stimulus.params.DELTA = stim_out.delta;
    datarun.stimulus.params.DIRECTION = stim_out.direction;
    datarun.stimulus.params.BACK_RGB = stim_out.back_rgb;
    datarun.stimulus.params.RGB = stim_out.rgb;
    trials(1:trial_n) = struct;
    for i = 1:trial_n
        trials(i).BAR_WIDTH = stim_out.trials(i).bar_width;
        trials(i).DELTA = stim_out.trials(i).delta;
        trials(i).DIRECTION = stim_out.trials(i).direction;
        trials(i).BACK_RGB = stim_out.trials(i).back_rgb;
        trials(i).RGB = stim_out.trials(i).rgb;
    end
    datarun.stimulus.trials = trials;
    datarun.stimulus.combinations = trials(1:trial_n/stim_out.repeats);
% elseif (strcmp(stim_out.type, 'MS'))
%     datarun.stimulus.x_vertices = stim_out.x_vertices_all;
%     datarun.stimulus.y_vertices = stim_out.y_vertices_all;
% else
%     fprintf('\t stimulus type not recognized. \n');
%     return
end

if length(datarun.stimulus.trial_list) ~= length(datarun.stimulus.triggers)
    disp(sprintf('!!!!!!!!!!!!\n load_stim: ERROR %d triggers - %d stimuli\n!!!!!!!!!!!!',length(datarun.stimulus.triggers), length(datarun.stimulus.trial_list)));
    trepeats=floor(length(datarun.stimulus.triggers)/length(datarun.stimulus.combinations));
    disp(sprintf('\n asume run ended early  delete last repetions: %d rep -> %d rep\n',datarun.stimulus.repetitions, trepeats)); 
    datarun.stimulus.repetitions=trepeats;
    datarun.stimulus.trials=datarun.stimulus.trials(1:datarun.stimulus.repetitions*length(datarun.stimulus.combinations));
    datarun.stimulus.trial_list=datarun.stimulus.trial_list(1:datarun.stimulus.repetitions*length(datarun.stimulus.combinations));
    datarun.stimulus.triggers=datarun.stimulus.triggers(1:datarun.stimulus.repetitions*length(datarun.stimulus.combinations));

end
