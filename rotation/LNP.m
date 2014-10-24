%% stimulus: binary BW 1800 s

movie_rate = 30;
T = 1800; % in sec i.e. 30min
movie_size = movie_rate*T;
movie_height = 40;
movie_width = 40;
depth = 15;
movie = rand(movie_height, movie_width, movie_size);
movie(movie<.5) = 0;
movie(movie>=.5) = 1;

%% generate linear filter

x = 1:movie_width;
y = 1:movie_height;
[XX, YY] = meshgrid(x, y);
A1 = 1;
A2 = 0.2;
sigma_c = 1;
sigma_s = 2;
x1 = 15;
y1 = 15;
sf1 = A1*exp(-((XX-x1).^2+(YY-y1).^2)/(2*sigma_c^2)) - A2*exp(-((XX-x1).^2+(YY-y1).^2)/(2*sigma_s^2));
Ds1 = repmat(sf1, [1, 1, 15]);

x2 = 17;
y2 = 15;
sf2 = A1*exp(-((XX-x2).^2+(YY-y2).^2)/(2*sigma_c^2)) - A2*exp(-((XX-x2).^2+(YY-y2).^2)/(2*sigma_s^2));
Ds2 = repmat(sf2, [1, 1, 15]);

figure;
mesh(sf1);
hold on
mesh(sf2);

p1 = 35;
p2 = 15;
tau1 = 2/30;
tau2 = 4/30;
n = 6;
dt = 1/30;
t = 0:dt:0.5;
t = t(1:depth);

tf = zeros(1, 1, depth);
tf(1, 1, :) = p1*(t/tau1).^n.*exp(-n*(t/tau1-1)) - p2*(t/tau2).^n.*exp(-n*(t/tau2-1));
Dt = repmat(tf, [movie_width movie_height 1]);

figure
plot(squeeze(tf))

D1 = Ds1.*Dt;
D2 = Ds2.*Dt;

%% generate response
D1r = zeros(movie_width, movie_height, depth);
for i = 1:depth
    D1r(:, :, i) = D1(:, :, depth+1-i);
end

n = movie_size-depth+1;
r1 = zeros(1, n);
for i = 1:n
    r_temp = movie(:, :, i:i+depth-1).*D1r;
    r1(i) = sum(r_temp(:));
end
r1(r1<0) = 0; % add non-linearity



D2r = zeros(movie_width, movie_height, depth);
for i = 1:15
    D2r(:, :, i) = D2(:, :, depth+1-i);
end

r2 = zeros(1, n);
for i = 1:n
    r_temp = movie(:, :, i:i+depth-1).*D2r;
    r2(i) = sum(r_temp(:));
end
r2(r2<0) = 0; % add non-linearity


%% generate spike trains

dt = 1/900;
N = T/dt;
r1r = repmat(r1, 30, 1);
r1r = reshape(r1r, 1, n*30);
spike1 = r1r*dt>rand(1, n*30);

r2r = repmat(r2, 30, 1);
r2r = reshape(r2r, 1, n*30);
spike2 = r2r*dt>rand(1, n*30);

frame_id1 = ceil(find(spike1 == 1)/30)+depth-1;
frame_id2 = ceil(find(spike2 == 1)/30)+depth-1;

spike_time1 = find(spike1 == 1)*dt;
spike_time2 = find(spike2 == 1)*dt;


%% calculate one-cell-sta

sts1 = zeros(movie_width, movie_height, length(frame_id1), depth);
for m = 1:depth
    sts1(:, :, :, m) = movie(:, :, frame_id1-15+m);
end
sta1 = mean(sts1, 3);
sta1 = sta1-mean(sta1(:));

sts2 = zeros(movie_width, movie_height, length(frame_id2), depth);
for m = 1:depth
    sts2(:, :, :, m) = movie(:, :, frame_id2-15+m);
end
sta2 = mean(sts2, 3);
sta2 = sta2-mean(sta2(:));

% figure;
% for h = 1:depth
%     colormap gray
%     imagesc(sta1(:, :, 1, h));
%     pause
% end


%% calculate two-cell-sta
interval = 0.03;
train_size = length(spike_time1);
T = [];
for i = 1:train_size
    p = find(spike_time2<spike_time1(i)+interval & spike_time2>spike_time1(i)-interval);
    if isempty(p) == 0;
       for j = 1: length(p);
           t = min(spike_time1(i), spike_time2(p(j)));
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

frame_id = ceil(T*movie_rate)+depth-1;
sts = zeros(movie_width, movie_height, length(frame_id), depth);
for m = 1:depth
    sts(:, :, :, m) = movie(:, :, frame_id-15+m);
end
sta = mean(sts, 3);
sta = sta-mean(sta(:));

% figure;
% for h = 1:depth
%     colormap gray
%     imagesc(sta(:, :, 1, h));
%     pause
% end

%% fit sta
opt = struct('fit_surround_amp_scale', true, 'fit_surround_sd_scale', true, 'verbose', true);
sta_temp = fit_sta(sta1, opt);
fit_para1 = zeros(1, 5);
fit_para1(1) = sta_temp.center_point_x;
fit_para1(2) = sta_temp.center_point_y;
fit_para1(3) = sta_temp.center_sd_x;
fit_para1(4) = sta_temp.center_sd_y;
fit_para1(5) = sta_temp.center_rotation_angle;

sta_temp = fit_sta(sta2, opt);
fit_para2 = zeros(1, 5);
fit_para2(1) = sta_temp.center_point_x;
fit_para2(2) = sta_temp.center_point_y;
fit_para2(3) = sta_temp.center_sd_x;
fit_para2(4) = sta_temp.center_sd_y;
fit_para2(5) = sta_temp.center_rotation_angle;

sta_temp = fit_sta(sta, opt);
fit_para = zeros(1, 5);
fit_para(1) = sta_temp.center_point_x;
fit_para(2) = sta_temp.center_point_y;
fit_para(3) = sta_temp.center_sd_x;
fit_para(4) = sta_temp.center_sd_y;
fit_para(5) = sta_temp.center_rotation_angle;

[X1, Y1] = drawEllipse(fit_para1);
[X2, Y2] = drawEllipse(fit_para2);
[X3, Y3] = drawEllipse(fit_para);

figure;
plot(X1, Y1, 'r');
hold on
plot(X2, Y2, 'g');
plot(X3, Y3, 'b');
legend('cell 1', 'cell 2', 'two-cell')
axis([13 19 13 19]);
