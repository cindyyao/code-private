function raster_all = get_ds_interval_raster(datarun, cell_id, interval)
% raster_all = get_ds_raster(datarun, cell_id)
% datarun: drifting grating dataset
% raster_all: S cells, M spatial period x N temporal period x R direction x P repeats
% for each cell
%
% xyao
% 2013-12-16

triggers = datarun.stimulus.triggers;
% index = grp_stim(datarun);

raster_all = cell(length(cell_id), 1);
for rgc = 1:length(cell_id)
    if ismember(cell_id(rgc), datarun.cell_ids)
    idx = get_cell_indices(datarun, cell_id(rgc));
    spike = datarun.spikes{idx};
    raster = get_raster(spike, triggers-interval, 'stop', interval,'plot', false);
    raster_all{rgc} = raster;
    end
end

    