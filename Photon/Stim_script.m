% cd('/Users/alexth/test4/Photons/Utils/mex_functions/')
% mex Draw_Random_Frame_opt.c

%% Initialization

cd('/Users/alexth/test4')

path2save = '/Users/alexth/test4/Photons/saved_stim/2015-08-25-1';
screen_number = 2;
def_params = initialize_display('OLED', screen_number);

% set gamma
scale = [0.9936    1.0018    0.9958];
power = [2.8349    2.8767    2.8196];
offset = [0.0006   -0.0021    0.0004];
set_gamma_from_fit_params(scale, power, offset);

%% Gamma calibration

steps = 17;
run_gamma_flashes('linear', steps, [0 0 1]); % FIRST PARAMETER will reset gamma to linear!

% linear measurements (INSERT HERE)
r = [0.01e-3 0.01e-3 0.06e-3 0.35e-3 1.10e-3 2.40e-3 4.34e-3 6.95e-3 10.30e-3 13.90e-3 18.27e-3 2.33e-2 2.95e-2 3.72e-2 4.63e-2 5.68e-2 6.84e-2];
g = [0e-3 0e-3 0.01e-2 0.03e-2 0.11e-2 0.27e-2 0.52e-2 0.87e-2 1.33e-2 1.82e-2 2.42e-2 3.13e-2 4.01e-2 5.08e-2 6.36e-2 7.77e-2 9.33e-2];
b = [0.01e-3 0.01e-3 0.09e-3 0.51e-3 1.58e-3 3.42e-3 6.12e-3 9.71e-3 14.21e-3 19.05e-3 2.5e-2 3.23e-2 4.10e-2 5.15e-2 6.41e-2 7.83e-2 9.40e-2];
w = [0.01e-3 0.02e-3 0.18e-3 1.11e-3 3.66e-3 8.15e-3 14.77e-3 2.37e-2 3.5e-2 4.75e-2 6.48e-2 8.53e-2 11.06e-2 14.08e-2 17.46e-2 0.212 0.255];

[scale, power, offset] = fit_gamma([r',g',b']);

% control
run_gamma_flashes([scale, power, offset], steps, [1 0 0]);

% corrected measurements (INSERT HERE, then proceed below)
r = [];
g = [];
b = [];
w = [];

% plot
x = linspace(0,255,steps);
figure
subplot(1,2,1); plot(x, r, 'r', x, g, 'g', x, b, 'b');
subplot(1,2,2); plot(x, w, 'k-o', x, r+g+b, 'm-*'); legend('measured', 'sum')

%% Set background

% white
mglClearScreen([1 1 1]);
mglFlush

% gray
mglClearScreen([1 1 1]*0.5);
mglFlush

% black
mglClearScreen([1 1 1]*0);
mglFlush


%% Focus Squares

fprintf('\n\n<strong> Focus squares. </strong>\n');
clear parameters stimulus;

parameters.class = 'FS';
stimulus = make_stimulus(parameters,def_params);
display_stimulus(stimulus);


%% Rectangular Flashing Pulses
fprintf('\n\n<strong> Rectangular Pulses: any sequence. </strong>\n');
clear parameters stimulus;

parameters.class = 'FP';
parameters.frames = 30;
parameters.delay_frames = 30;
parameters.back_rgb = [1 1 1]*0.5;
parameters.x_start = 120;  parameters.x_end = 620;
parameters.y_start = 50;   parameters.y_end = 400;

rgb = [1 0 1; 1 1 1; 0 1 0; -1 -1 -1]*0.5;
for i=1:size(rgb,1)
    stimulus = make_stimulus(parameters, 'rgb', rgb(i,:), def_params);
    save_parameters(stimulus, path2save, 'data000');
    display_stimulus(stimulus);
end


%% Cone-Isolating Pulse
fprintf('\n\n<strong> Cone-Isolating Pulse </strong>\n');
clear parameters stimulus;

parameters.class = 'PL';
parameters.back_rgb = [1 1 1]*0.5;
parameters.map_file_name = '/Users/alexth/test4/Photons/Maps/map_data034/map_data034.txt';

s_params = read_stim_lisp_output_hack('/Users/alexth/test4/Photons/Maps/s36_test'); % read S file

for i=2:size(s_params,2)
    trial_params = combine_parameters(parameters, s_params{1}, s_params{i});
    stimulus{i-1} = make_stimulus(trial_params, def_params);
    save_parameters(stimulus{i-1}, path2save, 'data000');
end

time_stamps = cell(1,length(stimulus));

for i=1:length(stimulus)
    time_stamps{i} = display_stimulus(stimulus{i});
end

for i=1:length(stimulus)
    save_time_stamps(time_stamps{i}, path2save, 'data000');
end

%% Moving bar

fprintf('\n\n<strong> Moving bar. </strong>\n');
clear parameters stimulus;

parameters.class = 'MB';
parameters.back_rgb = [1 1 1]*0.5;
parameters.rgb = -[1, 1, 1]*0.48;
parameters.bar_width = 30;
parameters.direction = 45;
parameters.delta = 2;  % pixels per frame
parameters.x_start = 100;  parameters.x_end = 300;
parameters.y_start = 100;   parameters.y_end = 350;
parameters.frames = 200;
parameters.delay_frames = 30;
parameters.repeats = 4;

direction = [0 30 45];% 60 90 120 135 150 180 210 225 240 270 300 315 330];
for i = 1:length(direction)
    stimulus{i} = make_stimulus(parameters,'direction', direction(i), def_params);    
    save_parameters(stimulus{i}, path2save, 'data000');
end

time_stamps = cell(1,length(stimulus));

for i=1:length(stimulus)
    time_stamps{i} = display_stimulus(stimulus{i});
end

for i=1:length(stimulus)
    save_time_stamps(time_stamps(i), path2save, 'data000'); % as one file!
end



%% Moving Grating
% gray
mglClearScreen([1 1 1]*0.5);
mglFlush


fprintf('\n\n<strong> Moving Grating. </strong>\n');
clear parameters stimulus

parameters.class = 'MG';
parameters.spatial_modulation = 'sine'; % sine or square
parameters.rgb = [1 1 1]*0.48;
parameters.back_rgb = [1 1 1]*0.5;
parameters.frames = 1; % presentation of each grating, frames
% parameters.x_start = 100;  parameters.x_end = 540;
% parameters.y_start = 100;   parameters.y_end = 380;
% parameters.x_start = 51;  parameters.x_end = 690;
% parameters.y_start = 51;   parameters.y_end = 530;
parameters.x_start = 1;  parameters.x_end = 640;
parameters.y_start = 1;   parameters.y_end = 480;
parameters.spatial_phase = 0; % pixels - offset from sin(0)
parameters.temporal_period = 180;  % frames (how long it takes to drift one period)
parameters.spatial_period = 160; % pixels
parameters.orientation = 180; % pixels

stimulus = make_stimulus(parameters, def_params);
display_stimulus(stimulus, 'erase', 0);



% s_params = read_stim_lisp_output_hack('/Users/alexth/test4/Photons/Maps/gratings/s03'); % read S file
% for i=2:size(s_params,2)
%     trial_params = combine_parameters(parameters, s_params{1}, s_params{i});
%     stimulus{i-1} = make_stimulus(trial_params, def_params);
%     save_parameters(stimulus{i-1}, path2save, 'data000');
% end

% 
% orientation = [0 30 45 60 90 120 135 150 180 210 225 240 270 300 315 330];
% for i = 1:length(orientation)
%     stimulus{i} = make_stimulus(parameters,'orientation', orientation(i), def_params);    
%     save_parameters(stimulus{i}, path2save, 'data000');
% end
% 
% time_stamps = cell(1,length(stimulus));
% 
% for i=1:length(stimulus)
%     time_stamps{i} = display_stimulus(stimulus{i});
% end
% 
% for i=1:length(stimulus)
%     save_time_stamps(time_stamps(i), path2save, 'data000'); % as one file!
% end

%% Counterphase Grating

fprintf('\n\n<strong> Counterphase Grating. </strong>\n');
clear parameters stimulus

parameters.class = 'CG';
parameters.spatial_modulation = 'sine'; % sine or square
parameters.temporal_modulation = 'sine'; % sine or square

parameters.rgb = [1 1 1]*0.48;
parameters.back_rgb = [1 1 1]*0.5;
parameters.frames = 120; % presentation of each grating, frames
parameters.x_start = 1;  parameters.x_end = 640;
parameters.y_start = 1;   parameters.y_end = 480;
parameters.spatial_phase = 0; % pixels - offset
parameters.temporal_period = 480;  % frames (how long it takes to return to initial phase)
parameters.spatial_period = 160; % pixels
parameters.orientation = 180; % pixels

stimulus = make_stimulus(parameters, def_params);
display_stimulus(stimulus, 'erase', 0);
% 
% % S file: copy logic from cone pulses
% orientation = [0 30 45];% 60 90 120 135 150 180 210 225 240 270 300 315 330];
% for i = 1:length(orientation)
%     stimulus{i} = make_stimulus(parameters,'orientation', orientation(i), def_params);    
%     save_parameters(stimulus{i}, path2save, 'data000');
% end
% 
% time_stamps = cell(1,length(stimulus));
% 
% for i=1:length(stimulus)
%     time_stamps{i} = display_stimulus(stimulus{i});
% end
% 
% for i=1:length(stimulus)
%     save_time_stamps(time_stamps(i), path2save, 'data000'); % as one file!
% end



%% Random Noise

%cd('/Users/alexth/test4/Photons/Utils/mex_functions/')
%mex Draw_Random_Frame_opt.c

fprintf('\n\n<strong> Random Noise </strong>\n');
clear parameters stimulus

parameters.class = 'RN';
parameters.back_rgb = [1 1 1]*0.5;
parameters.rgb = [1 1 1]*0.48;
parameters.interval = 1;
parameters.seed = 11111;
parameters.independent = 0;
parameters.binary = 1;
parameters.probability = 1;
parameters.delay_frames = 0;
% 
% parameters.x_start = 101;  parameters.x_end = 500;
% parameters.y_start = 101;   parameters.y_end = 500;
% parameters.stixel_width = 40;   parameters.stixel_height = 40;
% parameters.field_width = 10;  parameters.field_height = 10;

% parameters.x_start = 1;  parameters.x_end = 640;
% parameters.y_start = 81;   parameters.y_end = 400;
% parameters.stixel_width = 80;   parameters.stixel_height = 80;
% parameters.field_width = 8;  parameters.field_height = 4;

parameters.x_start = 1;  parameters.x_end = 320;
parameters.y_start = 1;   parameters.y_end = 320;
parameters.stixel_width = 10;   parameters.stixel_height = 10;
parameters.field_width = 32;  parameters.field_height = 32;

% parameters.x_start = 1;  parameters.x_end = 600;
% parameters.y_start = 1;   parameters.y_end = 600;
% parameters.stixel_width = 1;   parameters.stixel_height = 1;
% parameters.field_width = 600;  parameters.field_height = 600;

% parameters.map_file_name = '/Users/alexth/test4/Photons/Maps/2011-12-13-2_f04_vorcones/map-0000.txt';

parameters.frames = 10;%ceil(1*60*60);  % min * refresh rate (ceil it?) * 60(sec in min) - duration of each repetition!

parameters.jitter = 0;

stimulus = make_stimulus(parameters, def_params);
% save_parameters(stimulus, path2save, 'data000');
time_stamps = cell(1,10);
% for i=1:10
    time_stamps{i} = display_stimulus(stimulus, 'erase',0);
% end
% this is to test timing
% figure
% for i=1:10
%     subplot(3,4,i)
%     a = diff(time_stamps{i})*1000;
%     plot(a)
%     xlabel('flush')
%     ylabel('time per flush, ms')
% end

%% Raw Movie

fprintf('\n\n<strong> Raw Movie </strong>\n');
clear parameters stimulus;

parameters.class = 'RM';
parameters.back_rgb = [1 1 1]*0.5;
parameters.x_start = 200; % x_end and y_end wil depend on movie size (and stixel size)!
parameters.y_start = 200;
parameters.stixel_width = 1;   parameters.stixel_height = 1;
parameters.frames = ceil(0.05*60*60);  % min * refresh rate (ceil it?) * 60(sec in min) - duration of each repetition!
parameters.start_frame = 1; % >0
parameters.interval = 1;
parameters.flip = 1;  % 1 = normal; 2 = vertical flip; 3 = horizontal flip; 4 = vertical + horizontal flip
parameters.reverse = 0;   % 1 = backward (reverse), 0 = forward
parameters.movie_name = '/Users/alexth/test4/Photons/Movies/test_5_A.rawMovie';

stimulus = make_stimulus(parameters, def_params);
save_parameters(stimulus, path2save, 'data000');
time_stamps = display_stimulus(stimulus, 'trigger_interval', 100, 'wait_key',1, 'erase', 0);


%% photographic mapping

fprintf('\n\n<strong> Random Noise : Binary. </strong>\n');
clear parameters stimulus;

parameters.class = 'RN';
parameters.rgb = [1 1 1]*0.48;
parameters.independent = 0;
parameters.seed = 11111;
parameters.frames = 1;
parameters.x_start = 101;  parameters.x_end = 420;
parameters.y_start = 101;   parameters.y_end = 420;
% large 
parameters.stixel_width = 32;   parameters.stixel_height = 32;
parameters.field_width = 10;  parameters.field_height = 10;

% medium 
parameters.stixel_width = 10;   parameters.stixel_height = 10;
parameters.field_width = 32;  parameters.field_height = 32;

% small 
parameters.stixel_width = 1;   parameters.stixel_height = 1;
parameters.field_width = 320;  parameters.field_height = 320;

stimulus = make_stimulus(parameters, def_params);
display_stimulus(stimulus, 'erase',0);


%%
Stop_Photons

make_normal_table(8)

