function trial_dur = get_mb_trial_dur(datarun, width, height, delay)

% trial_dur: M delta x N bar_width

if isfield(datarun.stimulus,'x_start')
    display_width = datarun.stimulus.x_end - datarun.stimulus.x_start;
    display_height = datarun.stimulus.y_end - datarun.stimulus.y_start;
else
    display_width = width;
    display_height = height;
end

refresh_rate = 60.35;
delta = datarun.stimulus.params.DELTA;
bar_width = datarun.stimulus.params.BAR_WIDTH;

for dt = 1:length(delta)
    for bw = 1:length(bar_width)
        trial_dur(dt, bw) = (sqrt(display_width^2+display_height^2)+bar_width(bw))/delta(dt)/refresh_rate + delay;
    end
end
end