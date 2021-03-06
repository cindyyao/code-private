function [] = plot_ds_raster_one(DS, raster, cell_idx, cell_id, print_out)
% generate subplot indices for DS raster plot with multiple dataruns. 
% So far assume that there are 8 directions.
% DS: DS structure with all datasets.
% raster: rasters of all datasets.
% print: 1:print as pdf and close
%        0:don't print

[idx, xx, yy] = subplot_idx(1, 1);
ii = find(~cellfun(@isempty,raster{1}),1); % get the idx of 1st non-empty cell in raster
tpn = size(raster{1}{ii}, 2);
color = 'kbrgc';
tt = DS{1}.theta{1}(1, :);
for time = 1:tpn
    FigHandle = figure;
    set(FigHandle, 'Position', [1 1 800 800])
%     set(FigHandle, 'Position', [1 1 1400 600])
%     set(FigHandle, 'Position', [1 1 1620 1080])
%     set(FigHandle, 'Position', [1 1 1080 1080])
    repeats = 0;
    for j = 1:length(DS)
        if ~isempty(raster{j}{cell_idx}) && time <= length(DS{j}.rho)
            h = subplot(xx, yy, idx(1)); polar(tt, DS{j}.rho{time}(cell_idx, :), color(j));
            polar_theta_off(h)
            hold on
            for i = 2:9
                subplot(xx, yy, idx(i)); 
                plot_raster(squeeze(raster{j}{cell_idx}(1, time, i-1, :)), 0, 8, 'color', color(j), 'first_trial', repeats+1)
                if i == 4
                    title(num2str(cell_id))
                end 
                if mod(idx(i), yy) == 1
                    ylabel('trial number')
                end
                if idx(i) > yy*(xx-1)
                    xlabel('time (s)')
                end
            end
            repeats = repeats + size(raster{j}{ii}, 5);
        end
    end
     
    if print_out
        name = [num2str(cell_id) '_' num2str(time) '_dg'];
        screen_size = [24 12];
        set(figure(1), 'paperpositionmode', 'auto');
        set(gcf, 'PaperUnits', 'inch');
        set(figure(1), 'PaperSize', screen_size);
        print(figure(1), '-dpdf', name)
        close
    end
end

        

