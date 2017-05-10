clear all
clc
close all

datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2013-02-21-0/data004/data004');
datarun = load_neurons(datarun);
interval_time = .05;

% cell_numb = 1;
cell_numb = length(datarun.cell_ids);
for i = 1:cell_numb
    spike_time = datarun.spikes{i};
    train_size = length(datarun.spikes{i});
    interval = spike_time(2:train_size) - spike_time(1:train_size-1);
    interval = [interval; 0];    
    spike2 = interval < interval_time;
    % pick out the first spike of spike pairs whose interval time is < 50ms
    spike2_time = spike_time.*spike2;
    spike2_time(spike2_time == 0) = [];
    datarun.spikes{i,2} = spike2_time;
end

movie_path = '/Volumes/lab/acquisition/movie-xml/BW-15-2-0.48-11111-40x40-60.35.xml';
display_frame_rate = 60.35;
movie_refresh = 2;
movie_frame_rate = display_frame_rate / movie_refresh;
movie_frames = floor(movie_frame_rate * datarun.duration);
mov = get_movie(movie_path, datarun.triggers, movie_frames);
mov = mov(:, :, 1, :);
mov = squeeze(mov);
mov_height = 40;
mov_width = 40;
depth = 15;

frame_T = 1/movie_frame_rate;

% STA = cell(cell_numb, 2);

frame_time = zeros(1, (length(datarun.triggers) - 1)*50);
for j = 1:length(datarun.triggers)-1
    frame_time(50*(j-1)+1:50*j+1) = linspace(datarun.triggers(j), datarun.triggers(j+1), 51);
end


for k = 1:2
    for i = 1:5 %cell_numb
        spike_time = datarun.spikes{i, k};
        spike_time = spike_time(spike_time>frame_time(depth)& spike_time<=datarun.triggers(end)); 
        % make sure can get enough frames preceding the spike
        
        frame_id = zeros(1, length(spike_time));
        for n = 1:length(spike_time)
            t = frame_time(frame_time-spike_time(n)<0);
            frame_id(n) = find(frame_time == t(end));
        end
                
        sts = zeros(mov_height, mov_width, length(spike_time), depth);
        for m = 1:depth
            sts(:, :, :, m) = mov(:, :, frame_id-m+1);
        end
        
        sta = mean(sts, 3);
        sta = squeeze(sta);
        STA{i, k} = sta;        
               
    end
    
end
save('STA2.mat', 'STA')

figure;
cell_n = 25;
% time course difference: 3, 25
% amplitude difference: 1, 2, 4, 6
% 5
sta = STA{cell_n, 1};
[~, b] = max(sta(1:end));
z = ceil(b/1600);
y = ceil((b-(z-1)*1600)/40);
x = b-(z-1)*1600-40*(y-1);
a = 15:-1:1;
plot(a, squeeze(sta(x, y, :)), 'r');
hold on

sta = STA{cell_n, 1};
for h = 1:depth;
    imshow(sta(:, :, h));
    pause
end
[~, b] = max(sta(1:end));
z = ceil(b/1600);
y = ceil((b-(z-1)*1600)/40);
x = b-(z-1)*1600-40*(y-1);
a = 15:-1:1;
plot(a, squeeze(sta(x, y, :)), 'b');
legend('STA', '2STA', 'location', 'NorthWest')
title(['cell number = ' num2str(cell_n)])

sta_diff = STA{cell_n, 1} - STA{cell_n, 2};
for h = 1:depth
    colormap gray
    imagesc(sta(:, :, h));
    pause
end