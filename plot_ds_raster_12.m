function [] = plot_ds_raster_12(DS, raster, cell_idx, cell_id, Title, print_out)
% generate subplot indices for DS raster plot with multiple dataruns. 
% So far assume that there are 8 directions.
% DS: DS structure with all datasets.
% raster: rasters of all datasets.
% print: 1:print as pdf and close
%        0:don't print

idx = [8 4 3 2 1 5 9 13 14 15 16 12];
ii = find(~cellfun(@isempty,raster{1}),1); % get the idx of 1st non-empty cell in raster
tpn = size(raster{1}{ii}, 2);
tt = DS{1}.theta{1}(1, :);
for time = 1:tpn
    FigHandle = figure;
     set(FigHandle, 'Position', [1 1 400 400])
%     set(FigHandle, 'Position', [1 1 1400 600])
%     set(FigHandle, 'Position', [1 1 1620 1080])
%     set(FigHandle, 'Position', [1 1 900 900])
    for j = 1:length(DS)
        if ~isempty(raster{j}{cell_idx}) && time <= length(DS{j}.rho)
            h = subplot(4, 4, [6 7 10 11]); 
            u_temp = DS{j}.U{end-time+1}(cell_idx);
            v_temp = DS{j}.V{end-time+1}(cell_idx);
            alim = max(sqrt(u_temp^2+v_temp^2), 3);
            P = polar(0, alim);
            set(P, 'Visible', 'off')
            hold on
            compass(DS{j}.U{end-time+1}(cell_idx), DS{j}.V{end-time+1}(cell_idx), 'r');
            polar(tt, DS{j}.rho{end-time+1}(cell_idx, :), 'b');
            polar_theta_off(h)
            for i = 1:12
                subplot(4, 4, idx(j, i)); plot_raster(squeeze(raster{j}{cell_idx}(1, end-time+1, i, :)), 0, 8)
            end
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

        

