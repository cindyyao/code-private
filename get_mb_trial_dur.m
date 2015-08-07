function trial_dur = get_mb_trial_dur(datarun)

% trial_dur: M delta x N bar_width

display_width = 800; 
display_height = 600;
refresh_rate = 60.35;
delta = datarun.stimulus.params.DELTA;
bar_width = datarun.stimulus.params.BAR_WIDTH;

for dt = 1:length(delta)
    for bw = 1:length(bar_width)
        trial_dur(dt, bw) = (sqrt(display_width^2+display_height^2)+bar_width(bw))/delta(dt)/refresh_rate;
    end
end
end