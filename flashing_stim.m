%% Rectangular Flashing Pulses
fprintf('\n\n<strong> Rectangular Pulses: any sequence. </strong>\n');
clear parameters stimulus;

parameters.class = 'FP';
parameters.frames = 60;
parameters.delay_frames = 60;
parameters.back_rgb = repmat([1 1 1]*0.5, 4, 1);
parameters.x_start = 100;  parameters.x_end = 700;
parameters.y_start = 0;   parameters.y_end = 600;

stixel_size = 15;
subregions = [2 2];
maps = generate_flashing_map(parameters, stixel_size, subregions);

rgb = [1 1 1]*0.5;
for i=1:size(rgb,1)
    stimulus = make_stimulus(parameters, 'rgb', repmat(rgb(i, :), 4, 1), def_params);
    display_stimulus(stimulus);
end