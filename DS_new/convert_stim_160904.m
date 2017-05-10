function datarun = convert_stim_160904(datarun)
stim.params.DIRECTION = datarun.stimulus.params.direction;
stim.params.DELTA = datarun.stimulus.params.delta;
stim.params.BAR_WIDTH = datarun.stimulus.params.bar_width;
stim.params.BACK_RGB = {[0.2 0.2 0.2]};
for i = 1:length(datarun.stimulus.params.rgb)
    stim.params.RGB{i} = datarun.stimulus.params.rgb(i)*ones(1,3);
end
for i = 1:length(datarun.stimulus.trials)
    stim.trials(i).DELTA = datarun.stimulus.trials(i).delta;
    stim.trials(i).BAR_WIDTH = datarun.stimulus.trials(i).bar_width;
    stim.trials(i).DIRECTION = datarun.stimulus.trials(i).direction;
    stim.trials(i).BACK_RGB = [0.2 0.2 0.2];
    stim.trials(i).RGB = datarun.stimulus.trials(i).rgb*ones(1,3);
end

for i = 1:length(datarun.stimulus.combinations)
    stim.combinations(i).BAR_WIDTH = datarun.stimulus.combinations(i).bar_width;
    stim.combinations(i).DELTA = datarun.stimulus.combinations(i).delta;
    stim.combinations(i).DIRECTION = datarun.stimulus.combinations(i).direction;
    stim.combinations(i).BACK_RGB = [0.2 0.2 0.2];
    stim.combinations(i).RGB = datarun.stimulus.combinations(i).rgb*ones(1,3);
end
datarun.stimulus.params = stim.params;
datarun.stimulus.trials = stim.trials;
datarun.stimulus.combinations = stim.combinations;
end