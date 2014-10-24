function [corr_corrected] = shift_predictor_dg(datarun, cell1_id, cell2_id, bin_size, maxlag, dir)

raster1_t = get_ds_raster(datarun, cell1_id);
raster1_t = raster1_t{1};
raster2_t = get_ds_raster(datarun, cell2_id);
raster2_t = raster2_t{1};
edges = 0:bin_size:length(datarun.stimulus.triggers)*10;
max_lag = round(maxlag/bin_size); 



if dir == 0
    rsize = size(raster1_t);
    combn = prod(rsize(1:3));
    raster1 = [];
    raster2 = [];
    for i = 1:datarun.stimulus.repetitions
        raster1 = [raster1; reshape(raster1_t(:, :, :, i), combn, 1)];
        raster2 = [raster2; reshape(raster2_t(:, :, :, i), combn, 1)];
    end
else
    rsize = size(raster1_t);
    combn = prod(rsize(1:2));
    raster1 = [];
    raster2 = [];
    for i = 1:datarun.stimulus.repetitions
        raster1 = [raster1; reshape(raster1_t(:, :, dir, i), combn, 1)];
        raster2 = [raster2; reshape(raster2_t(:, :, dir, i), combn, 1)];
    end
end

for j = 1:length(raster1)
    raster1_temp{j} = raster1{j}+(j-1)*8;
    raster2_temp{j} = raster2{j}+(j-1)*8;
end
raster1_all = cell2mat(raster1_temp');
raster2_all = cell2mat(raster2_temp');
raster1_hc = histc(raster1_all, edges);
raster2_hc = histc(raster2_all, edges);

% calculate shift predictor
shift_index = [combn+1:length(raster1) 1:combn];
raster2_shift_temp = raster2;
for i = 2:datarun.stimulus.repetitions
    raster2_shift_temp = raster2_shift_temp(shift_index);
    for j = 1:length(raster1)
        raster2_shift_temp2{j} = raster2_shift_temp{j}+(j-1)*8;
    end
    raster2_shift(i-1, :) = cell2mat(raster2_shift_temp2');
    raster2_shift_hc(i-1, :) = histc(raster2_shift(i-1, :), edges);
    shift_predictor(i-1, :) = xcorr(raster1_hc, raster2_shift_hc(i-1, :), max_lag);
end

shift_predictor_mean = mean(shift_predictor);
corr = xcorr(raster1_hc, raster2_hc, max_lag);
corr_corrected = corr' - shift_predictor_mean;