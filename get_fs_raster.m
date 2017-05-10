function raster_all = get_fs_raster(datarun, cell_id, varargin)

p = inputParser;
p.addRequired('datarun', @isstruct);
p.addRequired('cell_id', @isnumeric);

p.addParameter('stop', [], @isnumeric);
p.parse(datarun, cell_id, varargin{:});
datarun = p.Results.datarun;
cell_id = p.Results.cell_id;
stop = p.Results.stop;

% trigger = datarun.triggers(2:end);
trigger = datarun.triggers;
index = grp_stim_fs(datarun);

raster_all = cell(length(cell_id), 1);
for i = 1:length(cell_id)
    if ~isempty(intersect(datarun.cell_ids, cell_id(i)))
        idx = get_cell_indices(datarun, cell_id(i));
        raster = get_raster(datarun.spikes{idx}, trigger, 'plot', false, 'stop', stop);
        raster_all{i} = raster(index);
    end
end