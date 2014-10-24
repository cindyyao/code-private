clear all
clc
close all

datarun = load_data('/lab/Experiments/Array/Analysis/2013-02-21-0/data004/data004');
datarun = load_neurons(datarun);
interval_time = .05;
cell_numb = 1;
% cell_numb = length(datarun.cell_ids);
for i = 1:cell_numb
    spike_time = datarun.spikes{i};
    train_size = length(datarun.spikes{i});
    interval = spike_time(2:train_size) - spike_time(1:train_size-1);
    interval = [interval; 0];
    spike2 = interval > interval_time;
    % pick out the first spike of spike pairs whose interval time is < 50ms
    spike2_time = spike_time.*spike2;
    spike2_time(spike2_time == 0) = [];
    datarun.spikes{i,2} = spike2_time;
end

% plot_sta(datarun, 33);

movie_path = '/lab/acquisition/movie-xml/BW-15-2-0.48-11111-40x40-60.35.xml';
display_frame_rate = 60;
movie_refresh = 2;
movie_frame_rate = display_frame_rate / movie_refresh;
movie_frames = movie_frame_rate * datarun.duration;
mov = get_movie(movie_path, datarun.triggers, movie_frames);
mov = mov(:,:,1,:);



mov_height = 40;
mov_width = 40;
depth = 15;

% rep = 30;
% A = 1:frame_numb;
% mov_fine = zeros(mov_height, mov_width, 1, frame_numb*rep);
% for i = 1:rep
%     mov_fine(:, :, 1, rep*(A-1)+i) = mov(:, :, 1, :);
% end


    
movie_frame_duration = 1/movie_frame_rate;

STA = zeros(mov_height, mov_width, depth+1, cell_numb, 2);
rep = 30;
for k = 1:2
    for i = 1:cell_numb
        spikes = datarun.spikes{i, k};
        spikes = spikes(spikes>depth*movie_frame_duration); % make sure can get enough frames preceding the spike
        spike_numb = size(spikes);
        sts = zeros(mov_height, mov_width, depth+1, spike_numb);
        frame_id = ceil(spikes/movie_frame_duration);
        for j = 1:depth+1
            frame = frame_id-j+1;
            sts(:, :, j, :) = mov(:, :, 1, frame);
        end
        
        sts_fine = zeros(mov_height, mov_width, rep*(depth+1), spike_numb);
        A = 1:depth+1;
        for m = 1:rep
            sts_fine(:, :, rep*(A-1)+m, :) = sts(:, :, :, :);
        end
        
        frame_id = ceil((spikes/movie_frame_duration - fix(spikes/movie_frame_duration))*rep);
        sts_fine1 = zeros(mov_height, mov_width, rep*depth, spike_numb);
        for n = 1:spike_numb
            sts_fine1(:, :, :, n) = sts_fine(:, :, frame_id:frame_id+rep*depth, :);
        end
                
        stsm_fine1 = mean(sts_fine1, 4);
        
        stsm = zeros(mov_height, mov_width, depth, spike_numb);
        for h = 1:depth
            stsm(:, :, h, :) = mean(stsm_fine1(:, :, rep*(h-1)+1:rep*h, :), 3);
        end
                
        STA(:, :, :, i, k) = stsm;
    end
    
end

STA_diff = STA(:, :, :, :, 1) - STA(:, :, :, :, 2);
colormap(gray);
imshow(STA(:, :, 3, 5, 1));

    

