function [raster_sum, p_idx, raster_sum_all] = get_pdirection_raster(raster, pdirection)

idx = find(~cellfun(@isempty,raster),1); % get the idx of 1st non-empty cell in raster
spn = size(raster{idx}, 1);
tpn = size(raster{idx}, 2);
dirn = size(raster{idx}, 3);
cln = size(raster{idx}, 4);
pdirection(pdirection<0) = pdirection(pdirection<0) + 2*pi;
p_idx = round(pdirection/(2*pi/dirn))+1;
p_idx(p_idx == dirn+1) = 1;

raster_sum = cell(length(raster), 1);
raster_sum_all = cell(length(raster), 1);
for cc = 1:length(raster)
    if ~isempty(raster{cc})
        for sp = 1:spn
            for tp = 1:tpn
                for cl = 1:cln
                    raster_sum{cc}{sp, tp, cl} = sort(cell2mat(squeeze(raster{cc}(sp, tp, p_idx(cc),cl, :))));
                    raster_sum_all{cc}(sp, tp, cl, :) = squeeze(raster{cc}(sp, tp, p_idx(cc),cl, :));
                end
            end
        end
    end
end

end
    

