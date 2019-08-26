opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);
path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-11-0/';
datarun = load_data(strcat(path, 'data004/data004'), opt);

% path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-05-22-0/';
% datatemp = load_data(strcat(path, 'data000/data000'), opt);

FIRST_TRIGGER = 0.95;
DELTA_TRIGGER = 1.66;
DURATION = 1800;
SAMPLING_RATE = 20000;
triggers = [FIRST_TRIGGER];
while triggers(end) + DELTA_TRIGGER < DURATION
    triggers(length(triggers) + 1) = triggers(end) + DELTA_TRIGGER;
end

movie_path = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-15-2-0.48-11111-53x40-60.35.xml';
vision_movie = compute_vision_movie_(movie_path, triggers);
sta_depth = 15;

% spikes = javaArray('java.lang.Integer', size(datarun.spikes, 1), 100000);
% 
% % spikes = datarun.spikes;
% for i = 1:length(spikes)
%     for j = 1:length(datarun.spikes{i})
%         spikes(i,j) = java.lang.Integer(datarun.spikes{i}(j) * SAMPLING_RATE);
%     end
% end

spikes = datarun.spikes;
for i = 1:length(spikes)
    spikes{i} = int32(spikes{i} * SAMPLING_RATE);
end
stas = compute_stas(spikes, int32(FIRST_TRIGGER*20000), vision_movie, int32(sta_depth));


%% Kiersten's code
%stim
FIRST_TRIGGER = 0.95;
DELTA_TRIGGER = 1.66;
DURATION = 1800;
SAMPLING_RATE = 20000;
movie_path = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-15-2-0.48-11111-53x40-60.35.xml';

triggers = [FIRST_TRIGGER];
while triggers(end) + DELTA_TRIGGER < DURATION
    triggers(length(triggers) + 1) = triggers(end) + DELTA_TRIGGER;
end
triggers_fine = [FIRST_TRIGGER];
while triggers_fine(end) + DELTA_TRIGGER/1000 < DURATION
    triggers_fine(length(triggers_fine) + 1) = triggers_fine(end) + DELTA_TRIGGER/1000;
end
[mov,~,~,dur,refresh] = get_movie(movie_path,triggers,54200);
movie_white_noise = squeeze(mov(:,:,1,:));
slen = 54100;
mov1 = movie_white_noise(:,:,1:slen);
nky = size(movie_white_noise,1);
nkx = size(movie_white_noise,2);
Stim0 = reshape(mov1,nky*nkx,[])'; %nky = size(movie_white_noise,1) nkx = size(movie_white_noise,2)
Stim = bsxfun(@minus,Stim0,mean(Stim0)); %subtract mean

%spikes
spikes = datarun.spikes{1};
nframes = 54200;
upsamplefactor = 10;
% [binned_spikes_coarse,binned_spikes_fine,dtStim,~] = bin_spikes( spikes,dataruns.(wn_runs_names{j}).triggers,nframes,upsamplefactor); %just bin up spikes according to frame rate
binned_spikes_fine = histc(spikes, triggers_fine);
rlen = upsamplefactor*slen; 
sps = binned_spikes_fine(1:rlen); 
%rebin coarsely 
sps_coarse = sum(reshape(sps,[],slen),1)';


%get sta
nkt = 30; %number of frames to go back in time
sta = simpleSTC(Stim,sps_coarse,nkt);

%% Suva's code
% data000 - working example

N_SPIKES_STA = 10000;


datarunpath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-11-0/data004/data004';
moviechunksfolder = 'path_to_xml_file_folder' ;
interval = 2;

% STA parameters
headerCapacity = int32(10000);
width = int32(53);
height = int32(40);
staDepth = int32(15);
staOffset = int32(0);
stixelwidth = 15;
stixelheight = 15;
refreshtime = 3;
stafilepath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-11-0/data004/data004.sta';
staFile = edu.ucsc.neurobiology.vision.io.STAFile(stafilepath, headerCapacity, width, height, staDepth, staOffset, stixelwidth, stixelwidth, refreshtime);

% Load the data run
datarun = load_data(datarunpath);
datarun = load_neurons(datarun);
samples_to_frames = [datarun.names.rrs_prefix '.stf'];    

% RRS works in samples, but it is more natural to work in samples for us.
% Let's convert the datarun to samples.
datarun = convert_datarun_times_to_samples(datarun);

% Get the time of the image refreshes from the ttls
t_frames = time_imrefresh_from_ttls(datarun.triggers);

% Use that to cache a map of samples to frame indices. 
% This only needs to be calculated once per dataset and it's slow, so 
% comment out the following two lines if you need to run the sta 
% calculation more than once.
map_samples_to_frames(1:length(t_frames), t_frames, datarun.duration, samples_to_frames);


% Progress bar
fprintf('Calculating STAs.\n');
fprintf([repmat('.',1,80) '\n']);
ndots = 0;

% Get the STAs
ncells = length(datarun.cell_ids);
for k = 1:ncells
    % Update progress bar
    if (k/ncells)*80 > ndots
        fprintf('.');
        ndots = ndots + 1;
    end
    
    % Get spike times and trim
    st = datarun.spikes{k};
    if length(st) > N_SPIKES_STA
        st = st(1:N_SPIKES_STA);
    end
    
    % Get cell ID
    cellid = datarun.cell_ids(k);
    
    % Calculate STA frame indices
    staind = compute_ta_ind(st, samples_to_frames);
    
    % Remap frame indices to movie indices
    % staindremapped = remap_event_indices(staind, ...
    %     1:length(t_frames), ...
    %     reshape(repmat(1:(length(t_frames)/interval), interval, 1), length(t_frames), 1).');
    % For WN remapping is simply dividing by interval
    staindremapped = ceil(staind/interval);
    
    % Calculate STA from movie indices
    [sta, e_sta] = compute_ta_from_ind(staindremapped, moviechunksfolder);

    % Convert the cell array STA to a Vision STA
    vsta = cell_array_to_vision_sta(sta, e_sta, refreshtime);
    
    % Add the STA to the STA file
    staFile.addSTA(cellid, vsta)
end
fprintf('\nSTA calculation done.\n');
staFile.close()


