%% basc function of conductance
t = 1:30;
coeff = [0.822 1 4 3 1 2 2];
samplerate = 1000;
frame_rate = 50;
Gexc = 1*EPSCfit(coeff, t);
figure
plot(t, Gexc, 'LineWidth',2)
xlabel('time (ms)'); ylabel('conductance (nS)')

T = 1:1000;
figure
st = repmat(random('poiss', 1, 1, length(T)*frame_rate/1000), samplerate/frame_rate, 1);
st = st(:);
a =conv(st,Gexc', 'same');
plot(a, 'LineWidth',2)
xlabel('time (ms)'); ylabel('conductance (nS)')

Gleak = 4; %
V_depolarize = 0;

%% mosaic positions
spacing = 80; % um
DSdensity = 20 ; % /mm^2 (from rabbit Vaney 1994 ~20-50)
exclusion_mean = 150/4 ; % dendrite length (um)/(um/pix)
exclusion_sigma = 20/4 ; 
grid_size = [300 300]-(exclusion_mean*2) ; % subtracting exclusion mean to keep cells of edges
num_cells = ceil((prod(grid_size)*4^2/1000^2)*DSdensity) ; % pix size of grid*(um/pix)/(um/mm)*density/mm^2
cell_coordinates = make_serial_exclusion_mosaic(grid_size, num_cells, exclusion_mean, exclusion_sigma) ; %
figure
plot(cell_coordinates(:, 1), cell_coordinates(:, 2), 'o')
hold on
Ggap = zeros(num_cells);
neighbors = [];
for i = 1:num_cells-1
    for j = i+1:num_cells
        if norm(cell_coordinates(i, :) - cell_coordinates(j, :)) < spacing
            Ggap(i,j) = 1;
            Ggap(j,i) = 1;
            neighbors = [neighbors; i j];
            plot(cell_coordinates([i,j],1), cell_coordinates([i,j],2))
        end
    end
end

%% Voltage traces
G = 0.1:0.02:0.3;
T = 1:1000*120;
poiss_mean = 2.5;
clear gexc ginh
for cc = 1:length(cell_coordinates)
    st = repmat(random('poiss', poiss_mean, 1, length(T)*frame_rate/1000), samplerate/frame_rate, 1);
    st = st(:);
    gexc(cc, :) = conv(st,Gexc', 'same');
end
ginh = zeros(size(gexc));

for g = 1:length(G)
    V_trace = LIFgj(gexc,ginh,Gleak,Ggap*G(g),0, samplerate,[]);
    [spike_times{g}, ~] = get_spike_time(V_trace, V_depolarize, samplerate);
    g
end

%%
duration = 120;
bin_size = 0.001;
max_lag = 10;
A_all = []; A_all_coeff = [];

for cp = 1:size(neighbors, 1)
%     xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
%     FigHandle = figure;
%     set(FigHandle, 'Position', [1 1 900 200])
    for g = 1:length(G)
    %     id1 = corr_cells(cp, 1);
    %     id2 = corr_cells(cp, 2);


        spikes1 = spike_times{g}{neighbors(cp, 1)};
        spikes1_TF= ceil(spikes1/bin_size);
        spikes1 = zeros(duration/bin_size, 1);
        spikes1(spikes1_TF) = 1;

        spikes2 = spike_times{g}{neighbors(cp, 2)};
        spikes2_TF= ceil(spikes2/bin_size);
        spikes2 = zeros(duration/bin_size, 1);
        spikes2(spikes2_TF) = 1;

        A = xcorr(spikes1, spikes2, max_lag)/duration;
        A_all(cp, :, g) = A;
        
        A_coeff = xcorr(spikes1, spikes2, max_lag, 'coeff');
        A_all_coeff(cp, :, g) = A_coeff;
%     % %     A = xcorr(spikes1, spikes2, max_lag);
%     %     [h, filteredA] = find_smallest_h(A);
%     %     [bootstat,bootsam] = bootstrp(N,@find_smallest_h_hist,rude(round(filteredA), 1:max_lag*2+1), max_lag);
%     %     p = sum(bootstat > h)/N;
%         subplot(1, length(G), g)
%     %     if p(c1, c2) < 0.05 && ratio(c1, c2) > 2 && maxi(c1, c2) > 0.75*max_lag && maxi(c1, c2) < 1.25*max_lag+1 && maxv(c1, c2) > 10
%            bar(xx, A, 'k')
%     %     else
%     %        bar(xx, A, 'r')
%     %     endclose
%         xlim([-0.01 0.01])
%         ylim([0 inf])
        

    end
%     print_close(1, [12,4], [num2str(id1) '  ' num2str(id2)])
end
t2p = squeeze(max(A_all, [], 2)) - squeeze(A_all(:, 11, :));

figure
subplot(1,2,1)
for i = 1:size(t2p,1)
    plot(G, t2p(i, :))
    hold on
end
xlabel('Ggap')
ylabel('Hz')
title('correlated firing rate')

t2p_mean = mean(t2p);
t2p_ste = std(t2p)/sqrt(size(t2p,1));
subplot(1,2,2)
errorbar(G, t2p_mean, t2p_ste);
xlabel('Ggap')
ylabel('Hz')
title('correlated firing rate')
xlim([0 0.6])



t2p_coeff = squeeze(max(A_all_coeff, [], 2)) - squeeze(A_all_coeff(:, 11, :));

figure
subplot(1,2,1)
for i = 1:size(t2p_coeff,1)
    plot(G, t2p_coeff(i, :))
    hold on
end
xlabel('Ggap')
ylabel('coeff')
title('correlation coefficient')

t2p_coeff_mean = mean(t2p_coeff);
t2p_coeff_ste = std(t2p_coeff)/sqrt(size(t2p_coeff,1));
subplot(1,2,2)
errorbar(G, t2p_coeff_mean, t2p_coeff_ste);
xlabel('Ggap')
ylabel('coeff')
title('correlation coefficient')



%% G = 0.32
%% Voltage traces
G = 0.15;
T = 1:1000*120;
poiss_mean = 2.5:0.05:3;
clear gexc ginh

for fr = 1:length(poiss_mean)
    for cc = 1:length(cell_coordinates)
        st = repmat(random('poiss', poiss_mean(fr), 1, length(T)*frame_rate/1000), samplerate/frame_rate, 1);
        st = st(:);
        gexc(cc, :) = conv(st,Gexc', 'same');
    end
    ginh = zeros(size(gexc));

    V_trace = LIFgj(gexc,ginh,Gleak,Ggap*G,0, samplerate,[]);
    [spike_times{fr}, ~] = get_spike_time(V_trace, V_depolarize, samplerate);
    fr
end

for fr = 1:length(poiss_mean)
    fr_mean(fr) = mean(cellfun(@length, spike_times{fr}))/120;
end
%%
duration = 120;
bin_size = 0.001;
max_lag = 10;
A_all = []; A_all_coeff = [];

for cp = 1:size(neighbors, 1)
%     xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
%     FigHandle = figure;
%     set(FigHandle, 'Position', [1 1 900 200])
    for g = 1:length(poiss_mean)
    %     id1 = corr_cells(cp, 1);
    %     id2 = corr_cells(cp, 2);


        spikes1 = spike_times{g}{neighbors(cp, 1)};
        spikes1_TF= ceil(spikes1/bin_size);
        spikes1 = zeros(duration/bin_size, 1);
        spikes1(spikes1_TF) = 1;

        spikes2 = spike_times{g}{neighbors(cp, 2)};
        spikes2_TF= ceil(spikes2/bin_size);
        spikes2 = zeros(duration/bin_size, 1);
        spikes2(spikes2_TF) = 1;

        A = xcorr(spikes1, spikes2, max_lag)/duration;
        A_all(cp, :, g) = A;
        
        A_coeff = xcorr(spikes1, spikes2, max_lag, 'coeff');
        A_all_coeff(cp, :, g) = A_coeff;
%     % %     A = xcorr(spikes1, spikes2, max_lag);
%     %     [h, filteredA] = find_smallest_h(A);
%     %     [bootstat,bootsam] = bootstrp(N,@find_smallest_h_hist,rude(round(filteredA), 1:max_lag*2+1), max_lag);
%     %     p = sum(bootstat > h)/N;
%         subplot(1, length(G), g)
%     %     if p(c1, c2) < 0.05 && ratio(c1, c2) > 2 && maxi(c1, c2) > 0.75*max_lag && maxi(c1, c2) < 1.25*max_lag+1 && maxv(c1, c2) > 10
%            bar(xx, A, 'k')
%     %     else
%     %        bar(xx, A, 'r')
%     %     endclose
%         xlim([-0.01 0.01])
%         ylim([0 inf])
        

    end
%     print_close(1, [12,4], [num2str(id1) '  ' num2str(id2)])
end
t2p = squeeze(max(A_all, [], 2)) - squeeze(A_all(:, 11, :));

figure
subplot(1,2,1)
for i = 1:size(t2p,1)
    plot(fr_mean, t2p(i, :))
    hold on
end
xlabel('firing rate')
ylabel('Hz')
title('correlated firing rate')

t2p_mean = mean(t2p);
t2p_ste = std(t2p)/sqrt(size(t2p,1));
subplot(1,2,2)
errorbar(fr_mean, t2p_mean, t2p_ste);
xlabel('firing rate')
ylabel('Hz')
title('correlated firing rate')
% xlim([0 0.6])


t2p_coeff = squeeze(max(A_all_coeff, [], 2)) - squeeze(A_all_coeff(:, 11, :));

figure
subplot(1,2,1)
for i = 1:size(t2p_coeff,1)
    plot(fr_mean, t2p_coeff(i, :))
    hold on
end
xlabel('firing rate')
ylabel('coeff')
title('correlation coefficient')

t2p_coeff_mean = mean(t2p_coeff);
t2p_coeff_ste = std(t2p_coeff)/sqrt(size(t2p_coeff,1));
subplot(1,2,2)
errorbar(fr_mean, t2p_coeff_mean, t2p_coeff_ste);
xlabel('firing rate')
ylabel('coeff')
title('correlation coefficient')
% xlim([0 0.6])