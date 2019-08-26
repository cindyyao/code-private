%% Initialization

my_path = '/Users/stimulus/Desktop/Photons-duke';

addpath(genpath(my_path))
cd(my_path)

path2save = [my_path, '/saved_stim'];
screen_number = 2;
def_params = initialize_display('OLED', screen_number);
% mglMoveWindow(1, 700)

% real refresh rate
% mglTestRefresh(2)

% set gamma CRT nov 2015
% scale = [0.9998    1.0072    1.0019];
% power = [2.7807    2.8437    2.7429];
% offset = [-0.0017   -0.0043   -0.0002];
% set_gamma_from_fit_params(scale, power, offset);

% set gamma OLED nov 2015
scale = [1.1156    1.0919    1.0921];
power = [1.1799    1.2878    1.2614];
offset = [-0.1228   -0.0961   -0.0955];
set_gamma_from_fit_params(scale, power, offset);
%% Run after canceling a stimulus

mglClearScreen(0.5);
mglFlush
mglClearScreen(0.5);
mglFlush


%% Focus Squares

fprintf('\n\n<strong> Focus squares. </strong>\n');
clear parameters stimulus;
 
parameters.class = 'FS';
stimulus = make_stimulus(parameters,def_params);
display_stimulus(stimulus); 

%% Random Noise. Daq triggers every 1s.

fprintf('\n\n<strong> Random Noise </strong>\n');
clear parameters stimulus

parameters.class = 'RN';
parameters.back_rgb = [1 1 1]*0.5;
parameters.rgb = [1 1 1]*0.48;
parameters.seed = 11111;
parameters.binary = 1;
parameters.probability = 1;
parameters.jitter = 0;
parameters.delay_frames = 0;

%%%%%%%%%%%%%% OLED %%%%%%%%%%%%%% 
parameters.x_start = 10;  parameters.x_end = 790;
parameters.y_start = 0;   parameters.y_end = 600;

%%%%%%%%%%%%%% CRT %%%%%%%%%%%%%% 
% parameters.x_start = 1;  parameters.x_end = 640;
% parameters.y_start = 1;   parameters.y_end = 480;

parameters.independent = 0;
parameters.interval = 6;
parameters.stixel_width = 30;
parameters.frames = 60*3700;

parameters.stixel_height = parameters.stixel_width;
parameters.field_width = (parameters.x_end-parameters.x_start+1)/parameters.stixel_width;  
parameters.field_height = (parameters.y_end-parameters.y_start+1)/parameters.stixel_height;

% For Voronoi, set stixel_height and stixel_width to 1 and pass a map path
% parameters.map_file_name = [my_path, '/Maps/2011-12-13-2_f04_vorcones/map-0000.txt'];

stimulus = make_stimulus(parameters, def_params);

time_stamps = display_stimulus(stimulus, 'wait_trigger', 1);

%% Noise repeats. Daq triggers every 1s.

fprintf('\n\n<strong> Random Noise Repeats </strong>\n');
clear parameters stimulus time_stamps

num_rep=10;
parameters.class = 'RN';
parameters.back_rgb = [1 1 1]*0.5;
parameters.rgb = [1 1 1]*0.48;
parameters.seed = 11111;
parameters.binary = 1;
parameters.probability = 1;
parameters.jitter = 0;
parameters.delay_frames = 10;

%%%%%%%%%%%%%% OLED %%%%%%%%%%%%%% 
parameters.x_start = 10;  parameters.x_end = 790;
parameters.y_start = 0;   parameters.y_end = 600;

%%%%%%%%%%%%%% CRT %%%%%%%%%%%%%% 
% parameters.x_start = 1;  parameters.x_end = 640;
% parameters.y_start = 1;   parameters.y_end = 480;

parameters.independent = 0;
parameters.interval = 6;
parameters.stixel_width = 30;
parameters.frames = 60*10;

parameters.stixel_height = parameters.stixel_width;
parameters.field_width = (parameters.x_end-parameters.x_start+1)/parameters.stixel_width;  
parameters.field_height = (parameters.y_end-parameters.y_start+1)/parameters.stixel_height;

% For Voronoi, set stixel_height and stixel_width to 1 and pass a map path
% parameters.map_file_name = [my_path, '/Maps/2011-12-13-2_f04_vorcones/map-0000.txt'];

stimulus = make_stimulus(parameters, def_params);
for i=1:num_rep
    time_stamps{i} = display_stimulus(stimulus, 'wait_trigger',1);
end

%% Moving Grating S File write. Saved as MG.txt in /Users/stimulus/Desktop/Photons-duke/saved_stim/2016-03-28-0
% 
% fprintf('\n\n<strong> Moving Grating. </strong>\n');
% clear parameters stimulus
% 
% parameters.class = 'MG';
% parameters.spatial_modulation = 'square'; % sine or square
% parameters.rgb = [1 1 1]*0.48;
% parameters.back_rgb = [1 1 1]*0.5;
% parameters.frames = 8*60; % presentation of each grating, frames
% parameters.x_start = 0;  parameters.x_end = 790;
% parameters.y_start = 0;   parameters.y_end = 600;
% % parameters.direction = 45;
% parameters.delay=0;
% 
% variable_parameters = randomize_parameters('direction', [0 30 60 90 120 150 180 210 240 270 300 330], 'temporal_period', [120 240], 'spatial_period', [240], 'nrepeats',4);
% path2file = write_s_file(parameters, variable_parameters);
% s_params = read_s_file(path2file);

% see second option example in "S File read"
% for i=2:size(s_params,2)
%     trial_params = combine_parameters(s_params{1}, s_params{i});
%     stimulus{i-1} = make_stimulus(trial_params, def_params);
%     display_stimulus(stimulus{i-1},'wait_trigger',1);
% end 

%%%%%%%%%% clean up %%%%%%%%%% 
% for i=1:length(stimulus)
%     for j=1:stimulus{i}.temporal_period
%         mglDeleteTexture(stimulus{i}.texture{j});
%     end
% end

%% Moving Grating S File read. Daq triggers every 10s.

fprintf('\n\n<strong> Moving Grating. </strong>\n');
clear parameters stimulus

parameters.class = 'MG';
parameters.spatial_modulation = 'sine'; % sine or square
parameters.back_rgb = [1 1 1]*0.5;

s_params = read_s_file([my_path, '/saved_stim/2016-03-28-0/mg.txt']);

%%%%%%%%%% OPTION 1: make/display in loop %%%%%%%%%%%%%%%
% for i=2:size(s_params,2)
%     trial_params = combine_parameters(parameters, s_params{1}, s_params{i});
%     stimulus{i-1} = make_stimulus(trial_params, def_params);
%     display_stimulus(stimulus{i-1});
% end

%%%%%%%%%% OPTION 2: pre-make all, display all %%%%%%%%%% 
for i=2:size(s_params,2)
    trial_params = combine_parameters(parameters, s_params{1}, s_params{i});
    stimulus{i-1} = make_stimulus(trial_params, def_params);
end
for i=1:length(stimulus)
    display_stimulus(stimulus{i},'wait_trigger',0);
end

%%%%%%%%%% clean up %%%%%%%%%% 
% for i=1:length(stimulus)
%     for j=1:stimulus{i}.temporal_period
%         mglDeleteTexture(stimulus{i}.texture{j});
%     end
% end
