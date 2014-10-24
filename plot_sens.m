% draw raster and psth
cell_id = 2206;
cell_n = find(datarun.cell_ids == cell_id);
raster = raster_all_cells{cell_n, 1};
N_total = [];

for i = 1:light_level
    t = spike_added{cell_n}{i};
    X = 0.025:0.05:trial_duration-0.025; % 50ms/bin
    N = hist(t, X);
    N_total = [N_total N];
end

figure;
for i = 1:light_level
    subplot(light_level, 2, 2*(i-1)+1);
    for j = 1:repeat
        SpikeTime = raster{i}{j};
        SpikeTime = SpikeTime';
        X = [SpikeTime; SpikeTime];
        Y = [ones(1, length(SpikeTime))*(j-0.9); ...
            ones(1, length(SpikeTime))*j];
        line(X, Y, 'color', 'b');
        xlim([0 trial_duration])
        hold on
    end
    ylabel([num2str(params(i))])
    if i == 1
        title(['cell id = ' num2str(cell_id)]);
    end
    
    subplot(light_level, 2, 2*i);        
    t = spike_added{cell_n}{i};
    X = 0.025:0.05:trial_duration-0.025; % 50ms/bin
    hist(t, X);     
    xlim([0 trial_duration])  
    
    M = max(N);
    if isempty(M) == 0;
        ylim([0,max(N_total)]);
    end

end



% plot sensitivity curve

% off GC

cell_id = 2206;
cell_n = datarun.cell_ids == cell_id;
sens= zeros(light_level, 1);
for j = 1:light_level
    t = spike_added{cell_n}{j};
    X = 0.025:0.05:trial_duration-0.025; % 50ms/bin
    N = hist(t, X);
    sens(j) = mean(N) - min(N);
end
m_sens = max(sens(:));
sens = sens/m_sens;

figure;
plot(params, sens);
title(['cell id = ' num2str(cell_id)]);
xlabel('light level')



% on GC
cell_id = 2942;
cell_n = find(datarun.cell_ids == cell_id);
sens= zeros(light_level, 1);
for j = 1:light_level
    t = spike_added{cell_n}{j};
    X = 0.025:0.05:trial_duration-0.025; % 50ms/bin
    N = hist(t, X);
    sens(j) = max(N) - mean(N);
end
m_sens = max(sens(:));
sens = sens/m_sens;

figure;
plot(params, sens);
title(['cell id = ' num2str(cell_id)]);
xlabel('light level')
