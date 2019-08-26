%% single cell model
%% basc function of conductance
t = 0:2000;
coeff = [600 100 200 10 1 100 100];
samplerate = 1000;
% speed = 850; % um/s
% bar_size = 300; %um
% t_lag = round(bar_size/speed*samplerate);
% input conductance
% Gexc1 = 1*EPSCfit(coeff, t);
% Gexc = Gexc1 + circshift(Gexc1, t_lag, 2);
Gexc = 5*EPSCfit(coeff, t);
figure
plot(t, Gexc, 'LineWidth',2)
xlabel('time (ms)'); ylabel('conductance (nS)')


Ginh_min = 1*Gexc;
Ginh_max = 3*Gexc;
Gleak = 4; %
V_depolarize = 0;
alpha = 1.2 ;
ndir = 9;
theta = linspace(-pi, pi, ndir);
Inh_tuning = (0.5 + 0.5*cos(theta-pi)).^alpha;

%% spike tuning curve
A = 1;
% V_trace = LIFgj([Gexc; Gexc],[Ginh; Ginh],Gleak,repmat(Ggap, 2, 2),0, samplerate,[]);
% V_trace = LIFgj(Gexc,Ginh,Gleak,Ggap,0, samplerate,[]);
V_trace = cell(length(Inh_tuning), length(A));
spike_counts = zeros(length(Inh_tuning), length(A));
for i = 1:length(A)
    for trial = 1:length(Inh_tuning)
        Ginh = A(i)*(Ginh_min + (Ginh_max - Ginh_min)*Inh_tuning(trial));
        Ginh_p(trial, i)= max(Ginh);
        Gexc_p(trial, i)= max(Gexc);
        V_trace{trial, i} = LIFgj(Gexc,Ginh,Gleak,0,0, samplerate,[]);
        [spike_times, spike_count] = get_spike_time(V_trace{trial, i}, V_depolarize, samplerate);
        spike_counts(trial, i) = spike_count';
    end
end
figure
subplot(2,1,1)
plot(theta/pi*180, Gexc_p, 'LineWidth',2)
hold on
plot(theta/pi*180, Ginh_p, 'LineWidth',2)
xlim([-180, 180])
xlabel('direction')
ylabel('peak conductance (nS)')
legend('Exc', 'Inh')
subplot(2,1,2)
plot(theta/pi*180, spike_counts, 'LineWidth',2, 'color', 'k')
% plot(theta/pi*180, spike_counts/max(spike_counts(:)), 'LineWidth',2)
xlim([-180, 180])
xlabel('direction')
ylabel('spike count')

%% Voltage traces of prefer and null direction
figure
for i = 1:size(V_trace{5}, 1)
    plot(t,V_trace{5}(i, :), 'Linewidth', 2, 'color', 'k')
    hold on
end
for i = 1:size(V_trace{1}, 1)
    plot(t,V_trace{1}(i, :), 'Linewidth', 2, 'color', [1 1 1]*0.6)
    hold on
end
legend('pref', 'null')
xlabel('time (ms)')
ylabel('Voltage (mV)')
xlim([0 1000])
ylim([-70 60])


%% effect of Ginh on spike tuning curve
% gain change
A = [1 .8 .6 .4 .3];
for i = 1:length(A)
    for trial = 1:length(Inh_tuning)
        Ginh = A(i)*(Ginh_min + (Ginh_max - Ginh_min)*Inh_tuning(trial));
        Ginh_p(trial, i)= max(Ginh);
        Gexc_p(trial, i)= max(Gexc);
        V_trace{trial, i} = LIFgj(Gexc,Ginh,Gleak,0,0, samplerate,[]);
        [spike_times, spike_count] = get_spike_time(V_trace{trial, i}, V_depolarize, samplerate);
        spike_counts(trial, i) = spike_count';
    end
end

figure
subplot(1,3,1)
plot(theta/pi*180, Ginh_p, 'LineWidth',2)
xlim([-180, 180])
xlabel('direction')
ylabel('peak conductance (nS)')
legend(['A = ' num2str(A(1))], ['A = ' num2str(A(2))], ['A = ' num2str(A(3))], ...
    ['A = ' num2str(A(4))], ['A = ' num2str(A(5))])
subplot(1,3,2)
plot(theta/pi*180, spike_counts, 'LineWidth',2)
% plot(theta/pi*180, spike_counts/max(spike_counts(:)), 'LineWidth',2)
xlim([-180, 180])
xlabel('direction')
ylabel('spike count')

subplot(1,3,3)
% plot(theta/pi*180, spike_counts, 'LineWidth',2)
plot(theta/pi*180, spike_counts./repmat(max(spike_counts), length(theta), 1), 'LineWidth',2)
xlim([-180, 180])
xlabel('direction')
ylabel('normalized spike count')


% DC change
A = [1 .8 .6 .4 .3];
for i = 1:length(A)
    for trial = 1:length(Inh_tuning)
        Ginh = A(i)*(Ginh_min) + (Ginh_max - Ginh_min)*Inh_tuning(trial);
        Ginh_p(trial, i)= max(Ginh);
        Gexc_p(trial, i)= max(Gexc);
        V_trace{trial, i} = LIFgj(Gexc,Ginh,Gleak,0,0, samplerate,[]);
        [spike_times, spike_count] = get_spike_time(V_trace{trial, i}, V_depolarize, samplerate);
        spike_counts(trial, i) = spike_count';
    end
end

figure
subplot(1,3,1)
plot(theta/pi*180, Ginh_p, 'LineWidth',2)
xlim([-180, 180])
xlabel('direction')
ylabel('peak conductance (nS)')
legend(['A = ' num2str(A(1))], ['A = ' num2str(A(2))], ['A = ' num2str(A(3))], ...
    ['A = ' num2str(A(4))], ['A = ' num2str(A(5))])
subplot(1,3,2)
plot(theta/pi*180, spike_counts, 'LineWidth',2)
% plot(theta/pi*180, spike_counts/max(spike_counts(:)), 'LineWidth',2)
xlim([-180, 180])
xlabel('direction')
ylabel('spike count')

subplot(1,3,3)
% plot(theta/pi*180, spike_counts, 'LineWidth',2)
plot(theta/pi*180, spike_counts./repmat(max(spike_counts), length(theta), 1), 'LineWidth',2)
xlim([-180, 180])
xlabel('direction')
ylabel('normalized spike count')

%% Network model
% 2 neuron model
Ginh = 1*(Ginh_min) + (Ginh_max - Ginh_min)*Inh_tuning(5);
Gexc_cells = [Gexc; zeros(1, length(t))];
Ginh_cells = 0.3*[Ginh; zeros(1, length(t))];

figure
subplot(1,2,1)
Ggap = [0 1;1 0]*1;
V_trace = LIFgj(Gexc_cells,Ginh_cells,Gleak,Ggap,0, samplerate,[]);
plot(t, V_trace, 'Linewidth', 2)
xlim([200 1000])
ylim([-70 60])
legend('cell 1', 'cell 2')
xlabel('time (ms)')
ylabel('voltage (mV)')
title('weak Ggap')

subplot(1,2,2)
Ggap = [0 1;1 0]*3;
V_trace = LIFgj(Gexc_cells,Ginh_cells,Gleak,Ggap,0, samplerate,[]);
plot(t, V_trace, 'Linewidth', 2)
xlim([200 1000])
ylim([-70 60])
legend('cell 1', 'cell 2')
xlabel('time (ms)')
ylabel('voltage (mV)')
title('strong Ggap')

%% multi neuron model
clear V_trace spike_counts
load('LIFmodel.mat')
speed = 850; % um/s
Ggap = [0 1 0 0 0 0 0 1;
        1 0 1 1 0 0 0 0;
        0 1 0 0 1 1 0 1;
        0 1 0 0 1 0 0 0;
        0 0 1 1 0 1 0 0;
        0 0 1 0 1 0 1 0;
        0 0 0 0 0 1 0 1;
        1 0 1 0 0 0 1 0]*1;
% Ggap = zeros(8,8);
    
cell_locations_1 = cell_locations - repmat(cell_locations(1, :), size(cell_locations, 1), 1);
theta = linspace(-pi, pi, ndir);
directions = [cos(theta); sin(theta)];
time_off = cell_locations_1*directions/speed;
time_off = round((time_off - repmat(min(time_off), size(time_off, 1), 1))*samplerate);
A = 0.3;
G = [0 0.01];
for g = 1:length(G)
    for i = 1:length(A)
        for trial = 1:length(theta)
            Ginh = A(i)*(Ginh_min + (Ginh_max - Ginh_min)*Inh_tuning(trial));
%             Ginh = Ginh_min + A(i)*(Ginh_max - Ginh_min)*Inh_tuning(trial);
            Ginh_p(trial, i)= max(Ginh);
            Gexc_p(trial, i)= max(Gexc);
            gexc = [];
            ginh = [];
            for cc = 1:length(cell_locations_1)
                gexc(cc, :) = circshift(Gexc, time_off(cc, trial), 2);
                ginh(cc, :) = circshift(Ginh, time_off(cc, trial), 2);
            end
            V_trace{trial, i} = LIFgj(gexc,ginh,Gleak,Ggap*G(g),0, samplerate,[]);
            [spike_times, spike_count] = get_spike_time(V_trace{trial, i}, V_depolarize, samplerate);
            spike_counts(g, i, trial, :) = spike_count;
        end
    end
end

figure
for g = 1:length(G)
    spike_unnorm = squeeze(spike_counts(g, 1, :, :));
    spike_norm = spike_unnorm./repmat(max(spike_unnorm), size(spike_unnorm, 1), 1);
    plot(theta/pi*180, mean(spike_norm, 2), 'Linewidth', 2)
    hold on
end
ylim([0 1])
xlabel('direction')
ylabel('normalized spike count')
legend('uncoupled', 'coupled')
figure
for g = 1:length(G)
    plot(theta/pi*180, squeeze(spike_counts(g, 1, :, 3))/max(squeeze(spike_counts(g, 1, :, 3))), 'Linewidth', 2)
    hold on
end
ylim([0 1])
xlabel('direction')
ylabel('normalized spike count')
legend('uncoupled', 'coupled')

%% change E/I ratio
A = 1;
V_trace = cell(length(Inh_tuning), length(A));
spike_counts = zeros(length(Inh_tuning), length(A));
for i = 1:length(A)
    for trial = 1:length(Inh_tuning)
        ginh = A(i)*(Ginh_min + (Ginh_max - Ginh_min)*Inh_tuning(trial));
        gexc = A(i)*Gexc;
        Ginh_p(trial, i)= max(ginh);
        Gexc_p(trial, i)= max(gexc);
        V_trace{trial, i} = LIFgj(gexc,ginh,Gleak,0,0, samplerate,[]);
        [spike_times, spike_count] = get_spike_time(V_trace{trial, i}, V_depolarize, samplerate);
        spike_counts(trial, i) = spike_count';
    end
end
figure
subplot(2,1,1)
plot(theta/pi*180, Gexc_p, 'LineWidth',2)
hold on
plot(theta/pi*180, Ginh_p, 'LineWidth',2)
xlim([-180, 180])
xlabel('direction')
ylabel('peak conductance (nS)')
legend('Exc', 'Inh')
subplot(2,1,2)
plot(theta/pi*180, spike_counts/max(spike_counts), 'LineWidth',2, 'color', 'k')
% plot(theta/pi*180, spike_counts/max(spike_counts(:)), 'LineWidth',2)
xlim([-180, 180])
xlabel('direction')
ylabel('spike count')

%% add noise to the model

clear V_trace spike_counts
load('LIFmodel.mat')
speed = 850; % um/s
Ggap = [0 1 0 0 0 0 0 1;
        1 0 1 1 0 0 0 0;
        0 1 0 0 1 1 0 1;
        0 1 0 0 1 0 0 0;
        0 0 1 1 0 1 0 0;
        0 0 1 0 1 0 1 0;
        0 0 0 0 0 1 0 1;
        1 0 1 0 0 0 1 0]*1;
% Ggap = zeros(8,8);
figure
for cc = 1:length(cell_locations)
    plot(cell_locations(cc, 1), cell_locations(cc, 2),'ko')
    hold on
    text(cell_locations(cc, 1)+5, cell_locations(cc, 2)+5, num2str(cc), 'FontSize', 10)
    
end

cell_locations_1 = cell_locations - repmat(cell_locations(1, :), size(cell_locations, 1), 1);
theta = linspace(-pi, pi, ndir);
directions = [cos(theta); sin(theta)];
time_off = cell_locations_1*directions/speed;
time_off = round((time_off - repmat(min(time_off), size(time_off, 1), 1))*samplerate);
A = 0.5;
G = [0 1];
for repeat = 1:30
for g = 1:2
    for i = 1:length(A)
        for trial = 1:length(theta)
            Ginh = A(i)*(Ginh_min + (Ginh_max - Ginh_min)*Inh_tuning(trial));
            Ginh_p(trial, i)= max(Ginh);
            Gexc_p(trial, i)= max(Gexc);
            gexc = [];
            ginh = [];
            for cc = 1:length(cell_locations_1)
                gexc_temp = circshift(Gexc, time_off(cc, trial), 2);
                gexc(cc, :) = 0.02*sqrt(gexc_temp).*randn(size(gexc_temp))+gexc_temp+0.01*randn(size(gexc_temp));
                ginh_temp = circshift(Ginh, time_off(cc, trial), 2);
                ginh(cc, :) = 0.02*sqrt(ginh_temp).*randn(size(ginh_temp))+ginh_temp+0.01*randn(size(ginh_temp));
            end
            V_trace{trial, i} = LIFgj(gexc,ginh,Gleak,Ggap*G(g),0, samplerate,[]);
            [spike_times, spike_count] = get_spike_time(V_trace{trial, i}, V_depolarize, samplerate);
            spike_counts(repeat, g, i, trial, :) = spike_count;
        end
    end
end
end

spike_counts_mean = squeeze(mean(spike_counts));
% figure
% for g = 1:length(G)
%     plot(squeeze(spike_counts_mean(g, :, 3))/max(squeeze(spike_counts_mean(g, :, 3))))
%     hold on
% end
% ylim([0 1])
% legend('uncoupled', 'coupled')
figure
for g = 1:length(G)
    plot(mean(squeeze(spike_counts_mean(g, :, :)), 2)/max(mean(squeeze(spike_counts_mean(g, :, :)), 2)))
    hold on
end
ylim([0 1])
legend('uncoupled', 'coupled')

%% generate larger network
n = 10; % layers of cells
spacing = 120; % um
x = [1:n]'*spacing;
cell_coordinates = [];
y = spacing/2*sqrt(3)*ones(n,1);
for i = 1:n
    if mod(i, 2) == 1
        xtemp = x;
    else
        xtemp = x + spacing/2;
    end
    ytemp = i*y;
    cell_coordinates = [cell_coordinates;[xtemp ytemp]];
end
jitter = 20*randn(size(cell_coordinates));
% jitter = 0;
cell_coordinates = cell_coordinates + jitter;
figure(1)
plot(cell_coordinates(:, 1), cell_coordinates(:, 2), 'o')
hold on

Ggap = zeros(n^2);
for i = 1:n^2-1
    for j = i+1:n^2
        if norm(cell_coordinates(i, :) - cell_coordinates(j, :)) < spacing + 30
            Ggap(i,j) = 1;
            Ggap(j,i) = 1;
            plot(cell_coordinates([i,j],1), cell_coordinates([i,j],2))
        end
    end
end
%% mosaic positions
spacing = 150; % um
DSdensity = 120 ; % /mm^2 (from rabbit Vaney 1994 ~20-50)
exclusion_mean = 150 ; % dendrite length (um)/(um/pix)
exclusion_sigma = 20 ; 
grid_size = [550 550] ; % subtracting exclusion mean to keep cells of edges
num_cells = ceil((prod(grid_size)/1000^2)*DSdensity) ; % pix size of grid*(um/pix)/(um/mm)*density/mm^2
cell_coordinates = make_serial_exclusion_mosaic(grid_size, num_cells, exclusion_mean, exclusion_sigma) ; %
figure
plot(cell_coordinates(:, 1), cell_coordinates(:, 2), 'o')
hold on
Ggap = zeros(num_cells);
for i = 1:num_cells-1
    for j = i+1:num_cells
        if norm(cell_coordinates(i, :) - cell_coordinates(j, :)) < spacing
            Ggap(i,j) = 1;
            Ggap(j,i) = 1;
            plot(cell_coordinates([i,j],1), cell_coordinates([i,j],2))
        end
    end
end

%% multi neuron model
clear V_trace spike_counts
speed = 850; % um/s
% Ggap = zeros(8,8);
    
theta = linspace(-pi, pi, ndir);
directions = [cos(theta); sin(theta)];
time_off = cell_coordinates*directions/speed;
time_off = round((time_off - repmat(min(time_off), size(time_off, 1), 1))*samplerate);
A = 0.47; % 0.66
G = [0 1 1.5];
for g = 1:length(G)
    for i = 1:length(A)
        for trial = 1:length(theta)
            Ginh = A(i)*(Ginh_min + (Ginh_max - Ginh_min)*Inh_tuning(trial));
%             Ginh = Ginh_min + A(i)*(Ginh_max - Ginh_min)*Inh_tuning(trial);
            Ginh_p(trial, i)= max(Ginh);
            Gexc_p(trial, i)= max(Gexc);
            gexc = [];
            ginh = [];
            for cc = 1:length(cell_coordinates)
                gexc(cc, :) = circshift(Gexc, time_off(cc, trial), 2);
                ginh(cc, :) = circshift(Ginh, time_off(cc, trial), 2);
            end
            V_trace{trial, i} = LIFgj(gexc,ginh,Gleak,Ggap*G(g),0, samplerate,[]);
            [spike_times, spike_count] = get_spike_time(V_trace{trial, i}, V_depolarize, samplerate);
            spike_counts(g, i, trial, :) = spike_count;
        end
    end
end
% 
% figure
% for g = 1:length(G)
%     spike_unnorm = squeeze(spike_counts(g, 1, :, :));
%     spike_norm = spike_unnorm./repmat(max(spike_unnorm), size(spike_unnorm, 1), 1);
%     plot(theta/pi*180, mean(spike_norm, 2), 'Linewidth', 2)
%     hold on
% end
% ylim([0 1])
% xlabel('direction')
% ylabel('normalized spike count')
% legend('uncoupled', 'coupled')
lower = 100;
upper = 550-lower;
aa = cell_coordinates(:, 1) > lower & cell_coordinates(:, 1) < upper ...
     & cell_coordinates(:, 2) > lower & cell_coordinates(:, 2) < upper;


figure
for g = 1:length(G)
    spike_unnorm = squeeze(spike_counts(g, 1, :, :));
%     spike_unnorm = squeeze(spike_counts(g, 1, :, aa));
    spike_unnorm_mean = mean(spike_unnorm, 2);
    spike_norm = spike_unnorm_mean/max(spike_unnorm_mean);
%     plot(theta/pi*180, spike_norm, 'Linewidth', 2)
    plot(theta/pi*180, spike_unnorm_mean, 'Linewidth', 2)
    hold on
end
xlabel('direction')
ylabel('normalized spike count')
legend('uncoupled', 'coupled')

figure
for g = 1:length(G)
%     spike_unnorm = squeeze(spike_counts(g, 1, :, :));
    spike_unnorm = squeeze(spike_counts(g, 1, :, aa));
    spike_unnorm_mean = mean(spike_unnorm, 2);
    spike_norm = spike_unnorm_mean/max(spike_unnorm_mean);
%     plot(theta/pi*180, spike_norm, 'Linewidth', 2)
    plot(theta/pi*180, spike_unnorm_mean, 'Linewidth', 2)
    hold on
end
xlabel('direction')
ylabel('normalized spike count')
legend('uncoupled', 'coupled')


figure
for g = 1:length(G)
    spike_unnorm = squeeze(spike_counts(g, 1, :, aa));
    spike_unnorm_mean = mean(spike_unnorm, 2);
    spike_norm = spike_unnorm_mean/max(spike_unnorm_mean);
    plot(theta/pi*180, spike_norm, 'Linewidth', 2)
%     plot(theta/pi*180, spike_unnorm_mean, 'Linewidth', 2)
    hold on
end
ylim([0 1])
xlabel('direction')
ylabel('normalized spike count')
legend('uncoupled', 'coupled')

% figure
% for g = 1:length(G)
%     plot(theta/pi*180, squeeze(spike_counts(g, 1, :, 3))/max(squeeze(spike_counts(g, 1, :, 3))), 'Linewidth', 2)
%     hold on
% end
% ylim([0 1])
% xlabel('direction')
% ylabel('normalized spike count')
% legend('uncoupled', 'coupled')
