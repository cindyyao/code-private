opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);

% load data
datarun{1} = load_data('/Analysis/xyao/2014-01-14-0/data009-map/data009-map', opt);
datarun{1}.names.stimulus_path = '/Analysis/xyao/2014-01-14-0/stimuli/s09';
datarun{1} = load_stim(datarun{1}, 'user_defined_trigger_interval', 10);

datarun{2} = load_data('/Analysis/xyao/2014-01-14-0/data008-map/data008-map', opt);


[NumSpikesCell, StimComb] = get_spikescellstim(datarun{1},datarun{1}.cell_ids,0);
[mag MAG dsindex magmax magave angle rho RHO theta num U V] = dscellanalysis(NumSpikesCell, StimComb);

% pull out DS cells

figure
plot(mag{1, 1}, mag{2, 1}, 'o')
hold on
[x, y] = ginput;
plot(x, y);
xlabel('TP 30')
ylabel('TP 60')
title('Vector sum plot')

IN = inpolygon(mag{1, 1}, mag{2, 1}, x, y);
[~, I] = find(IN == 1);
id = datarun{1}.cell_ids(I);

for j = 1:2
    [theta_seq, I_seq] = sort(theta{j}(1, :));
    r = rho{j};
    R = RHO{j};
    rho_seq = r(:, I_seq);
    RHO_seq = R(:, I_seq);
    rho{j} = rho_seq;
    RHO{j} = RHO_seq;
end

tt = datarun{1}.stimulus.params.DIRECTION*pi/180;
tp = datarun{1}.stimulus.params.TEMPORAL_PERIOD;
sp = datarun{1}.stimulus.params.SPATIAL_PERIOD;

raster = get_ds_raster(datarun{1}, datarun{1}.cell_ids);

DS_type = {'DS on1', 'DS on2', 'DS off1', 'DS off3'};
nDS_type = {'ON 1', 'ON 2', 'OFF 1', 'OFF 2', 'OFF 3'};

for i = 1:length(nDS_type)
    nDS_id{i} = intersect(get_cell_ids(datarun{2}, nDS_type{i}), datarun{1}.cell_ids);
end

x = 3;
y = 6;
% name = cell(length(id)*length(tp), 1);
for ct = 4:length(nDS_type)
    for cn = 1:length(nDS_id{ct})
        idx = get_cell_indices(datarun{1}, nDS_id{ct}(cn));
        FigHandle = figure;
        set(FigHandle, 'Position', get(0, 'ScreenSize'))
        if ~isempty(raster{idx})
        subplot(x, y, 8); polar(tt, rho{1}(idx, :))
        subplot(x, y, 9); plot_raster(squeeze(raster{idx}(1, 1, 1, :)), 0, 8)
        subplot(x, y, 3); plot_raster(squeeze(raster{idx}(1, 1, 2, :)), 0, 8)
        subplot(x, y, 2); plot_raster(squeeze(raster{idx}(1, 1, 3, :)), 0, 8); title('TP 32')
        subplot(x, y, 1); plot_raster(squeeze(raster{idx}(1, 1, 4, :)), 0, 8); title([nDS_type{ct} ' ' num2str(nDS_id{ct}(cn))])
        subplot(x, y, 7); plot_raster(squeeze(raster{idx}(1, 1, 5, :)), 0, 8)
        subplot(x, y, 13); plot_raster(squeeze(raster{idx}(1, 1, 6, :)), 0, 8)
        subplot(x, y, 14); plot_raster(squeeze(raster{idx}(1, 1, 7, :)), 0, 8)
        subplot(x, y, 15); plot_raster(squeeze(raster{idx}(1, 1, 8, :)), 0, 8)
        
        subplot(x, y, 11); polar(tt, rho{2}(idx, :))
        subplot(x, y, 12); plot_raster(squeeze(raster{idx}(1, 2, 1, :)), 0, 8)
        subplot(x, y, 6); plot_raster(squeeze(raster{idx}(1, 2, 2, :)), 0, 8)
        subplot(x, y, 5); plot_raster(squeeze(raster{idx}(1, 2, 3, :)), 0, 8); title('TP 128')
        subplot(x, y, 4); plot_raster(squeeze(raster{idx}(1, 2, 4, :)), 0, 8)
        subplot(x, y, 10); plot_raster(squeeze(raster{idx}(1, 2, 5, :)), 0, 8)
        subplot(x, y, 16); plot_raster(squeeze(raster{idx}(1, 2, 6, :)), 0, 8)
        subplot(x, y, 17); plot_raster(squeeze(raster{idx}(1, 2, 7, :)), 0, 8)
        subplot(x, y, 18); plot_raster(squeeze(raster{idx}(1, 2, 8, :)), 0, 8)
        end
        
screen_size = [24 12];
set(figure(1), 'paperpositionmode', 'auto');
set(gcf, 'PaperUnits', 'inch');
set(figure(1), 'PaperSize', screen_size);
print(figure(1), '-dpdf', [nDS_type{ct} ' ' num2str(nDS_id{ct}(cn))])
close

%         name{length(tp)*(cc-1)+t} = [num2str(id(cc)) '_' num2str(tp(t))];
    end
end

