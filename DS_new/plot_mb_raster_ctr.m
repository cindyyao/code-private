function [] = plot_mb_raster_ctr(MB, raster, trial_dur, cell_idx, cell_id, Title, x, y, print_out)
% generate subplot indices for DS raster plot with multiple dataruns. 
% So far assume that there are 8 directions.
% DS: DS structure with all datasets.
% raster: rasters of all datasets.
% print: 1:print as pdf and close
%        0:don't print

[idx, xx, yy] = subplot_idx(x, y);
ii = find(~cellfun(@isempty,raster{1}),1); % get the idx of 1st non-empty cell in raster
tpn = size(raster{1}{ii}, 2);
tt = MB{1}.theta{1}(1, :);
for time = 1:tpn
    FigHandle = figure;
%     set(FigHandle, 'Position', [1 1 1400 600])
%     set(FigHandle, 'Position', [1 1 2000 1080])
%     set(FigHandle, 'Position', [1 1 1080 1080])
    set(FigHandle, 'Position', [1 1 2000 1080])
    for j = 1:length(MB)
        for ctr = 1:size(MB{j}.rho, 3)
            if ~isempty(raster{j}{cell_idx}) && time <= length(MB{j}.rho)
                jj = (j-1)*size(MB{j}.rho, 3) + ctr;
                h = subplot(xx, yy, idx(jj, 1)) ; polar(tt, MB{j}.rho{1,time,ctr}(cell_idx, :));
                polar_theta_off(h)
                for i = 2:9
                    h = subplot(xx, yy, idx(jj, i)); plot_raster(squeeze(raster{j}{cell_idx}(1, time, i-1, ctr,:)), 0, trial_dur{j}(time))
%                     if idx(jj, i) == 2
%                         title('high contrast')
%                     elseif idx(jj, i) == 6
%                         title('medium contrast')
%                     elseif idx(jj, i) == 10
%                         title('low contrast')
%                     end 
%                     if idx(jj, i) == 12
%                         ylabel('control')
%                     elseif idx(jj, i) == 56
%                         ylabel('drug')
%                     elseif idx(jj, i) == 100 
%                         ylabel('wash')
%                     end
%                     if idx(jj, i) > yy*(xx-1)
%                         xlabel('time (s)')
%                     else
                        set(h, 'xticklabel', [])
%                     end
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

        

