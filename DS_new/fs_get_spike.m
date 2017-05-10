function spikes = fs_get_spike(raster)

spikes = cell(length(raster),1);
for i = 1:length(raster)
    if ~isempty(raster{i})
        spikes{i} = cellfun('length', raster{i});
    end
end
end