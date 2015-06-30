opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

data_trigger{1} = load_data('/Analysis/xyao/2014-10-14-0/data000-mapn/data000-map', opt);
data_trigger{1}.names.stimulus_path = '/Analysis/xyao/2014-10-14-0/stimuli/s00.mat';
data_trigger{1} = load_stim_matlab(data_trigger{1});

data_trigger{2} = load_data('/Analysis/xyao/2014-10-14-0/data003/data003', opt);
data_trigger{2}.names.stimulus_path = '/Analysis/xyao/2014-10-14-0/stimuli/s03.mat';
data_trigger{2} = load_stim_matlab(data_trigger{2});

data_trigger{3} = load_data('/Analysis/xyao/2014-10-28-0/data003/data003', opt);
data_trigger{3}.names.stimulus_path = '/Analysis/xyao/2014-10-28-0/stimuli/s03.mat';
data_trigger{3} = load_stim_matlab(data_trigger{3}, 'user_defined_trigger_interval', 10);

data_trigger{4} = load_data('/Analysis/xyao/2014-10-28-0/data000-map/data000-map', opt);
data_trigger{4}.names.stimulus_path = '/Analysis/xyao/2014-10-28-0/stimuli/s00.mat';
data_trigger{4} = load_stim_matlab(data_trigger{4}, 'user_defined_trigger_interval', 10);

figure
for i = 1:4
   trigger_interval{i} = diff(data_trigger{i}.stimulus.triggers)*60.35;
   subplot(4, 1, i)
   plot(trigger_interval{i})
   xlabel('trial#')
   ylabel('trigger interval(frames)')
   ylim([min(trigger_interval{i}) max(trigger_interval{i})])
   if i == 1
       title('Drifting Grating')
   end
end

data_ffp{1} = load_data('/Analysis/xyao/2014-10-14-0/data002-mapn/data002-map', opt);
data_ffp{1}.triggers = data_ffp{1}.triggers(2:end);
data_ffp{2} = load_data('/Analysis/xyao/2014-10-14-0/data006-mapn/data006-map', opt);
data_ffp{2}.triggers = data_ffp{2}.triggers(2:end);
data_ffp{3} = load_data('/Analysis/xyao/2014-10-28-0/data002/data002', opt);
data_ffp{3}.triggers = data_ffp{3}.triggers(2:end);
data_ffp{4} = load_data('/Analysis/xyao/2014-10-28-0/data005/data005', opt);
data_ffp{4}.triggers = data_ffp{4}.triggers(2:end);

figure
for i = 1:4
   trigger_interval{i} = diff(data_ffp{i}.triggers)*60.35;
   subplot(4, 1, i)
   plot(trigger_interval{i})
   xlabel('trial#')
   ylabel('trigger interval(frames)')
   ylim([min(trigger_interval{i}) max(trigger_interval{i})])
   if i == 1
       title('Full Field Pulses')
   end
end

