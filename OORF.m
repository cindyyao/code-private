%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1, 'load_sta', 1);

datarun = load_data('/Volumes/lab/analysis/2016-06-17-0/data000-map/data000-map', opt);

movie_path = '/Volumes/lab/acquisition/movie-xml/BW-10-2-0.48-11111-60x60-60.35.xml';
display_frame_rate = 60.35;
movie_refresh = 2;
movie_frame_rate = display_frame_rate / movie_refresh;
movie_frames = length(datarun.triggers)*100/movie_refresh;
mov = get_movie(movie_path, datarun.triggers, movie_frames);
mov = mov(:, :, 1, :);
mov = squeeze(mov);
mov_height = 60;
mov_width = 60;
depth = 15;

frame_T = 1/movie_frame_rate;

% STA = cell(cell_numb, 2);

frame_time = zeros(1, (length(datarun.triggers) - 1)*100/movie_refresh);
for j = 1:length(datarun.triggers)-1
    frame_time(100/movie_refresh*(j-1)+1:100/movie_refresh*j+1) = linspace(datarun.triggers(j), datarun.triggers(j+1), 100/movie_refresh+1);
end


cell_id = 81;
idx = get_cell_indices(datarun, cell_id);

spike_time = datarun.spikes{idx};
spike_time = spike_time(spike_time>frame_time(depth)& spike_time<=datarun.triggers(end)); 
% make sure can get enough frames preceding the spike
        
frame_id = zeros(1, length(spike_time));
for n = 1:length(spike_time)
    t = frame_time(frame_time-spike_time(n)<0);
    frame_id(n) = find(frame_time == t(end));
end
                
ste = zeros(mov_height, mov_width, length(spike_time), depth);
for m = 1:depth
    ste(:, :, :, m) = mov(:, :, frame_id-m+1);
end

sta = mean(ste, 3);
sta = squeeze(sta);

figure
for t = 1:15
    imagesc(sta(:,:,t))
    t
    pause
end

figure
plot(squeeze(sta(37,41,:)))

center_ste = squeeze(ste(37,41,:,:));
[coeff, score] = pca(center_ste);
figure
plot(score(:, 1), score(:, 2), 'o')