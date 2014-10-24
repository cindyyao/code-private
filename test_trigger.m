%% load data
opt = struct('load_params', 1,'load_neurons', 1);

% load data
datarun_t{1} = load_data('/Analysis/xyao/trigger_test_2/data000/data000', opt);


datarun_t{2} = load_data('/Analysis/xyao/trigger_test_2/data001/data001', opt);
datarun_t{2}.names.stimulus_path = '/RSM/test_mb.mat';
datarun_t{2} = load_stim_matlab(datarun_t{2});



datarun_t{3} = load_data('/Analysis/xyao/trigger_test_2/data002/data002', opt);
datarun_t{3}.names.stimulus_path = '/Analysis/xyao/trigger_test_2/stimuli/s02.mat';
datarun_t{3} = load_stim_matlab(datarun_t{3});

datarun_t{4} = load_data('/Analysis/xyao/trigger_test_2/data003/data003', opt);
datarun_t{5} = load_data('/Analysis/xyao/trigger_test_2/data004/data004', opt);
datarun_t{6} = load_data('/Analysis/xyao/trigger_test_2/data005/data005', opt);