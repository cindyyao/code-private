%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);
datadg = load_data('/Analysis/xyao/2012-10-31-0/data002-map/data002-map', opt);
datadg.names.stimulus_path = '/Analysis/xyao/2012-10-31-0/stimuli/s02';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);
datagwn = load_data('/Analysis/xyao/2012-10-31-0/data004/data004', opt);
datagwn = load_sta(datagwn);
datagwn = get_rf_coms(datagwn, datagwn.cell_ids);
%% stimulus distribution

movie_path = '/Volumes/lab/acquisition/movie-xml/BW-gaussian-8-4-0.16-11111.xml';
display_frame_rate = 120.0841;
movie_refresh = 4;
movie_frame_rate = display_frame_rate / movie_refresh;
movie_frames = floor(movie_frame_rate * datagwn.duration);
mov = get_movie(movie_path, datagwn.triggers, movie_frames);
mov = squeeze(mov(:,:,1,:));

%% get spike-triggered ensemble
mov_height = 40;
mov_width = 80;
depth = 15;
movie_frame_duration = 1/movie_frame_rate;

cell_id = 361;
cell_idx = get_cell_indices(datagwn, cell_id);
spikes = datagwn.spikes{cell_idx};
spikes = spikes(spikes>depth*movie_frame_duration); % make sure can get enough frames preceding the spike
spike_n = length(spikes);
ste = zeros(mov_height, mov_width, depth+1, spike_n); 
frame_id = ceil(spikes/movie_frame_duration);
for i = 1:spike_n
    ste(:,:,:,i) = mov(:,:,frame_id(i)-depth:frame_id(i));
end
    
sta = mean(ste, 4);  
com = round(datagwn.stas.rf_coms{cell_idx});

bin_size = 0.01;
XX = bin_size/2:bin_size:1-bin_size/2;
h_s = hist(squeeze(mov(com(2),com(2),:)),XX);
h_s_n = h_s/sum(h_s);
for i = 1:depth+1
    h_ste(:, i) = hist(squeeze(ste(com(2),com(1),i,:)), XX);
end
h_ste_n = h_ste./repmat(sum(h_ste), length(XX),1);

%%

% distribution
figure
for i = 1:depth+1
    plot(XX, h_s_n, 'b',XX, h_ste_n(:, i), 'r')
    pause
end

s_n_cum = cumsum(h_s_n);
s_ste_cum = cumsum(h_ste_n);

% cumulative distribution
figure
for i = 1:16
    plot(XX, s_n_cum, 'b', XX, s_ste_cum(:,i), 'r')
    pause
end

% diff of distribution
h_diff = h_ste_n - repmat(h_s_n', 1, 16); 
figure
for i = 1:16
    plot(XX, h_diff(:,i))
    pause
end

% diff of cumulative distrubution
cum_diff = s_ste_cum - repmat(s_n_cum', 1, 16);
figure
for i = 1:16
    plot(XX, cum_diff(:,i), XX, zeros(100,1));
    pause
end

sta_(:,:,1,:) = sta;
sta = repmat(sta_, [1,1,3,1]);

plot_sta_(sta)


    

