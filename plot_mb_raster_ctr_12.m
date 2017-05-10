function [] = plot_mb_raster_ctr_12(MB, raster, trial_dur, cell_idx, cell_id, Title, x, y, print_out)
% generate subplot indices for DS raster plot with multiple dataruns. 
% So far assume that there are 12 directions.
% DS: DS structure with all datasets.
% raster: rasters of all datasets.
% print: 1:print as pdf and close
%        0:don't print

[idx, xx, yy] = subplot_idx_12(x, y);
ii = find(~cellfun(@isempty,raster{1}),1); % get the idx of 1st non-empty cell in raster
tpn = size(raster{1}{ii}, 2);
tt = MB{1}.theta{1}(1, :);
for time = 1:tpn
    FigHandle = figure;
%     set(FigHandle, 'Position', [1 1 1400 600])
    set(FigHandle, 'Position', [1 1 1520 1080])
%     set(FigHandle, 'Position', [1 1 1080 1080])
    for j = 1:length(MB)
        for ctr = 1:size(MB{j}.rho, 3)
%             FigHandle = figure;
%             %     set(FigHandle, 'Position', [1 1 1400 600])
%             set(FigHandle, 'Position', [1 1 1520 1080])
%             %     set(FigHandle, 'Position', [1 1 1080 1080])

            if ~isempty(raster{j}{cell_idx}) && time <= length(MB{j}.rho)
                jj = (j-1)*size(MB{j}.rho, 3) + ctr;
%                 jj = 1;
%                 jj = j;
                h = subplot(xx, yy, idx(jj, 13:16)) ; polar(tt([1:end 1]), MB{j}.rho{1,time,ctr}(cell_idx, ([1:end 1])));
                polar_theta_off(h)
                for i = 1:12
                    h = subplot(xx, yy, idx(jj, i)); plot_raster(squeeze(raster{j}{cell_idx}(1, time, i, ctr,:)), 0, trial_dur{j}(time))
                    set(h, 'xticklabel', [])
                    set(h, 'yticklabel', [])
                end
            end
        end
    end
     
    if print_out
        name = [num2str(cell_id) '_' num2str(time) '_' Title];
        screen_size = [24 12];
        set(figure(1), 'paperpositionmode', 'auto');
        set(gcf, 'PaperUnits', 'inch');
        set(figure(1), 'PaperSize', screen_size);
        print(figure(1), '-dpdf', name)
        close
    end
end

        

