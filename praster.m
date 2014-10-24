cell_n = 4;
raster = raster_all_cells{cell_n, 1};
for i = 1:9
    subplot(9, 2, 2*(i-1)+1);
    for j = 1:24
        SpikeTime = raster{25*(i-1)+j, 1};
        SpikeTime = SpikeTime';
        X = [SpikeTime; SpikeTime];
        Y = [ones(1, length(SpikeTime))*(j-0.9); ...
            ones(1, length(SpikeTime))*j];
        line(X, Y, 'color', 'b');
        axis([0,12,0,25])
        hold on
    end
    subplot(9, 2, 2*i);
    t = spike_added{cell_n, i};
    N = hist(t, 100);
    hist(t, 100);     
    xlim([0 12])  
    
    M = max(N);
    if isempty(M) == 0;
        ylim([0,M+10]);
    end

end