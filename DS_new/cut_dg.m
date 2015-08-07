function [DG_less, raster_dg_less, raster_p_sum_less] = cut_dg(DG, raster_dg, raster_p_sum, speed_number, idx)

for i = 1:length(DG)
    if length(DG{i}.mag) > speed_number
        DG{i}.mag = DG{i}.mag(idx);
        DG{i}.MAG = DG{i}.MAG(idx);
        DG{i}.dsindex = DG{i}.dsindex(idx);
        DG{i}.angel = DG{i}.angle(idx);
        DG{i}.rho = DG{i}.rho(idx);
        DG{i}.RHO = DG{i}.RHO(idx);
        DG{i}.theta = DG{i}.theta(idx);
        DG{i}.num_t = DG{i}.num_t(idx);
        DG{i}.U = DG{i}.U(idx);
        DG{i}.V = DG{i}.V(idx);
        
        for cc = 1:length(raster_dg{i})
            if ~isempty(raster_dg{i}{cc})
                raster_dg{i}{cc} = raster_dg{i}{cc}(:, idx, :, :);
            end
            if ~isempty(raster_p_sum{i}{cc})
                raster_p_sum{i}{cc} = raster_p_sum{i}{cc}(idx);
            end
        end
        
    end
end
DG_less = DG;
raster_dg_less = raster_dg;
raster_p_sum_less = raster_p_sum;
end
        
        



