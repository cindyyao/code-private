function datarun = load_mfs_161019(datarun)

datarun.stimulus = [];
% check if stimulus file path is defined
if isfield(datarun.names, 'stimulus_path')
    stim_path = datarun.names.stimulus_path;
else
    fprintf('\t stimulus path not recognized. Please define stimulus path and try again. \n');
    return
end

current_path = pwd;
cd('/Users/xyao/Field-lab/Photons/Classes')
load(stim_path)
cd(current_path)

datarun.stimulus.repetitions = stim_out.repeats - 1;
datarun.stimulus.trial_list = stim_out.trial_list(stim_out.trial_num+1:end);

start_i = stim_out.trial_num*(stim_out.repeats-1)*2;
datarun.triggers = datarun.triggers(end-start_i+1:end);
triggers = datarun.triggers(1:2:end);

datarun.stimulus.triggers = triggers';
trial_n = length(stim_out.trial_list);