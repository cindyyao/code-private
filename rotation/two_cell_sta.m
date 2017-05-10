datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2013-02-21-0/data004/data004');
datarun = load_neurons(datarun);
interval = 0.01;
cell_numb = length(datarun.cell_ids);

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

frame_time = zeros(1, (length(datarun.triggers) - 1)*50);
for j = 1:length(datarun.triggers)-1
    frame_time(50*(j-1)+1:50*j+1) = linspace(datarun.triggers(j), datarun.triggers(j+1), 51);
end

two_cell_spike = cell(cell_numb, cell_numb);
two_cell_sta = cell(cell_numb, cell_numb);

load cell_pair.mat
for k = 1:1
    cell_1 = find(datarun.cell_ids == cell_pair(k, 1));
    cell_2 = find(datarun.cell_ids == cell_pair(k, 2));
    % on brisk-transient cell pairs:
    % [33 7666] [5476 5583] [661 1036] [1681 2118] [1066 1186]

    % on transient cell pairs:
    % [2252 2357] [3364 3557] [7292 7426] [5342 5869] [3047 6722]

    % off brisk-transient cell pairs:
    % [2206 2463] [2881 3241] [1951 2012] [4156 4546] [4531 5161]

    % off transient cell pairs:
    % [976 1351] [167 558] [7517 7728] [5915 6095] [5371 6886]

    % off sustained cell pairs:
    % [541 993] [106 7231] [3917 4111]* [5057 5312] [4891 5357]
    spike_1 = datarun.spikes{cell_1, 1};
    spike_2 = datarun.spikes{cell_2, 1};
    train_size_1 = length(spike_1);
    
    T = [];
    for i = 1:train_size_1
        p = find(spike_2<spike_1(i)+interval & spike_2>spike_1(i)-interval);
        if isempty(p) == 0;
           for j = 1: length(p);
               t = min(spike_1(i), spike_2(p(j)));
               T = [T; t];               
           end
                              
        end
    end
    
    T = sort(T);
    delta = T(2:end) - T(1:end-1);
    delta(delta ~= 0) = 1;
    delta = [1; delta];
    T = T.*delta;
    T(T == 0) = [];
    
    two_cell_spike{cell_1, cell_2} = T;
    two_cell_spike{cell_2, cell_1} = T;

    spike_time = two_cell_spike{cell_1, cell_2};
    spike_time = spike_time(spike_time>frame_time(depth)& spike_time<=datarun.triggers(end)); 
    % make sure can get enough frames preceding the spike
        
    frame_id = zeros(1, length(spike_time));
    for n = 1:length(spike_time)
        t = frame_time(frame_time-spike_time(n)<0);
        frame_id(n) = find(frame_time == t(end));
    end


    sts = zeros(mov_height, mov_width, 1, length(spike_time), depth);
    for m = 1:depth
        sts(:, :, 1, :, 16-m) = mov(:, :, frame_id-m+1);
        fprintf(num2str(m))
    end
        
    sta = mean(sts, 4);
    sta = sta- mean(sta(:));
    sta = reshape(sta, mov_height, mov_width, 1, depth);
    two_cell_sta{cell_1, cell_2} = sta;  
    two_cell_sta{cell_2, cell_1} = sta; 
    
%     save('two_cell_temp.mat', 'two_cell_spike', 'two_cell_sta')
    
    fprintf(num2str(k))
end


figure;
sta1 = STA{cell_1, 1};
[~, b] = max(sta1(1:end));
z = ceil(b/1600);
y = ceil((b-(z-1)*1600)/40);
x = b-(z-1)*1600-40*(y-1);
a = 15:-1:1;
plot(a, squeeze(sta1(x, y, :)), 'r');
hold on

sta2 = STA{cell_2, 1};
[~, b] = max(sta2(1:end));
z = ceil(b/1600);
y = ceil((b-(z-1)*1600)/40);
x = b-(z-1)*1600-40*(y-1);
a = 15:-1:1;
plot(a, squeeze(sta2(x, y, :)), 'g');

sta3 = two_cell_sta{cell_1, cell_2};
[~, b] = max(sta3(1:end));
z = ceil(b/1600);
y = ceil((b-(z-1)*1600)/40);
x = b-(z-1)*1600-40*(y-1);
a = 15:-1:1;
plot(a, squeeze(sta3(x, y, :)), 'b');

legend('cell 1 STA',two_cell_sta 'cell 2 STA', 'two cell STA', 'location', 'NorthWest')
title(['cell number = ' num2str(cell_1) ' and ' num2str(cell_2)])
figure;
for h = 1:depth
    colormap gray
    imagesc(sta(:, :, 1,  h));
    pause
end

% fit_para = cell(25, 3);
% sta_fit = cell(25, 3);


for i = 1:25    
    cell_1 = find(datarun.cell_ids == cell_pair(i, 1));
    cell_2 = find(datarun.cell_ids == cell_pair(i, 2));
    sta = two_cell_sta{cell_1, cell_2};
    sta_temp = fit_sta(sta, 'fit_scale_one', false, 'fit_scale_two', false, 'fit_tau_one', false, 'fit_tau_two', false);
    sta_fit{i, 5} = sta_temp;
    fit_para_temp = zeros(1, 5);
    fit_para_temp(1) = sta_temp.center_point_x;
    fit_para_temp(2) = sta_temp.center_point_y;
    fit_para_temp(3) = sta_temp.center_sd_x;
    fit_para_temp(4) = sta_temp.center_sd_y;
    fit_para_temp(5) = sta_temp.center_rotation_angle;
    fit_para{i, 5} = fit_para_temp;
end

for i = 1:25
    for j = 1:2
        cell_n = find(datarun.cell_ids == cell_pair(i, j));
        sta = STA{cell_n, 1};
%         sta_temp = fit_sta(sta, 'fit_scale_one', false, 'fit_scale_two', false, 'fit_tau_one', false, 'fit_tau_two', false);
        sta_temp = fit_sta(sta, 'fit_surround', true, 'fit_surround_sd_scale', true, 'fit_surround_amp_scale', true, 'verbose', true);
        sta_fit{i, j} = sta_temp;
        fit_para_temp = zeros(1, 5);
        fit_para_temp(1) = sta_temp.center_point_x;
        fit_para_temp(2) = sta_temp.center_point_y;
        fit_para_temp(3) = sta_temp.center_sd_x;
        fit_para_temp(4) = sta_temp.center_sd_y;
        fit_para_temp(5) = sta_temp.center_rotation_angle;
        fit_para{i, j} = fit_para_temp;
    end
end

figure;
pair_numb = 21;
[X1, Y1] = drawEllipse(fit_para{pair_numb, 1});
[X2, Y2] = drawEllipse(fit_para{pair_numb, 2});
[X3, Y3] = drawEllipse(fit_para{pair_numb, 5});

plot(X1, Y1, 'r');
hold on
plot(X2, Y2, 'g');
plot(X3, Y3, 'b');
legend('cell 1', 'cell 2', 'two-cell')
title(['cell number = ' num2str(cell_pair(pair_numb, 1)) ...
    ' and ' num2str(cell_pair(pair_numb, 2))])

center_error = zeros(25, 1);
for i = 1:25
    center_distance = norm(fit_para{i, 1}(1:2) - fit_para{i, 2}(1:2));
    midpoint = norm((fit_para{i, 1}(1:2) + fit_para{i, 2}(1:2))/2 - fit_para{i, 5}(1:2));
    error = midpoint/center_distance;
    center_error(i) = error;
end

figure;
hist(center_error, 20);