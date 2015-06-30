function [datarun] = delete_last_repeat(datarun)

trial_n = length(datarun.stimulus.combinations)*(datarun.stimulus.repetitions - 1);
cutout_t = datarun.stimulus.triggers(trial_n + 1);
datarun.stimulus.trial_list = datarun.stimulus.trial_list(1:trial_n);
datarun.stimulus.triggers = datarun.stimulus.triggers(1:trial_n);
datarun.stimulus.trials = datarun.stimulus.trials(1:trial_n);
datarun.stimulus.repetitions = datarun.stimulus.repetitions - 1;

for cc = 1:length(datarun.spikes)
    datarun.spikes{cc} = datarun.spikes{cc}(datarun.spikes{cc} < cutout_t);
end

