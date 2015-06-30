
datarun = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data004/data004');
% load the path of all the data files in data004??
datarun = load_neurons(datarun);
% neurons file contain spike time information
% triggers: times that stimuli being deliverd to the retina
% duration: in second
% sampling_rate: in HZ
% spikes: each cell contains an array that represent spike time of a neuron
% cell_ids: same with that in "Vision"
% channels: electrode id of the neurons???
datarun = load_sta(datarun);
% sta files contain spike trigged average information
% stas:
%       depth: number of frames that used to compute the sta
%       stas: every cell in it contains 4-D data with size 40*40(image size)*1(BW)*15(depth) 
%             
datarun = load_params(datarun);
% vision: 
%         sta_fits: each cell contains RF parameters of 1 neuron
%         cell_types: cell types defined in "Vision"


datarun = get_sta_summaries(datarun, 'all');


plot_rf(datarun, 61, 'foa', 2, 'com', true, 'autozoom', true)

%  temp_image = zeros(mov_height, mov_width);
%         for fm = 1:length(frame_id)
%             temp_image = temp_image + squeeze(mov(:,:,frame_id(fm)-2));
%             imagesc(temp_image); colormap gray
%             drawnow
%         end
%         imagesc(temp_image); colormap gray
%         %drawnow        
%   

datarun = load_data('/Analysis/xyao/2013-04-09-0/data000/data000');
datarun = load_neurons(datarun);

datarun.names.stimulus_path = '/lab/Experiments/Array/Analysis/2013-04-09-0/stimuli/s00';
datarun = load_stim(datarun, 'user_defined_trigger_set', [1:2:480], 'find_trigger', 0);

plot_rf_summaries
% Plot collection of RFs in a variety of ways

sta_fit_function
% This low level function receives a list of params (the order of which is
% important) and makes a spatial-temporal-chromatic STA based on those
% parameters.


load_txt_cell_types


temp_marks_sta15 = significant_stixels(sta15, 'thresh', 5, 'time', 'max');
 

time_course_from_sta

% plot more than one curves with the colors you want
p.addParamValue('colors', 'krbg')
for tt = 1:length(cell_types)

    % GATHER THE NEEDED DATA
    
    % get the cell type
    cell_type = cell_types(tt);
    % get the color
    color = colors(mod(tt - 1,length(colors))+1);
end

%% 

screen_size = [22 12];
% name = {'type_1', 'type_2'};
for i =1:4
set(figure(i), 'paperpositionmode', 'auto');
% set(figure(1), 'PaperPosition', [-0.5 -0.25 22 10]);
set(gcf, 'PaperUnits', 'inch');
set(figure(i), 'PaperSize', screen_size);
print(figure(i), '-dpdf', num2str(i))
end

%%

inpolygon
grpstats
clearvars -except

%% 
i = 2;

% set(figure(i), 'paperpositionmode', 'auto');
% set(figure(i), 'PaperSize', 'auto');
print(figure(i), '-dpdf', 'vector_sum_plot')
