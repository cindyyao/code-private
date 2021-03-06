function [vision_movie] = compute_vision_movie_(movie_xml_file, triggers, monitorFrequency)

% gets white noise movie from vision
% see COMPUTE_STA  

% import java classes
import('edu.ucsc.neurobiology.vision.util.*');
import('edu.ucsc.neurobiology.vision.matlab.*');

% grab sampling frequency
sampling_frequency = VisionParams.SAMPLES_PER_MILLISECOND * 1000;

% convert to integers
triggers = int32(triggers * sampling_frequency);

% compute the movieed
if nargin < 3
    vision_movie = Matlab.computeMovie(movie_xml_file, triggers);
else
    vision_movie = Matlab.computeMovie(movie_xml_file, triggers, monitorFrequency);    
end
