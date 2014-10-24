function [corr_corrected] = shift_predictor_ffp(datarun, cell1_id, cell2_id, bin_size)

cell1_idx = get_cell_indices(datarun, cell1_id);
cell2_idx = get_cell_indices(datarun, cell2_id);
spikes1 = datarun.spikes{cell1_idx};
spikes2 = datarun.spikes{cell2_idx};
edges = 0:bin_size:datarun.duration;
maxlag = round(0.1/bin_size); 
spikes1_hc = histc(spikes1, edges);
spikes2_hc = histc(spikes2, edges); 

n = length(datarun.triggers)/4;
trigger = zeros(1, n);
for i = 1:n
    trigger(i) = datarun.triggers(4*(i-1)+1);
end
trigger = trigger - trigger(1);

for i = 1:n-1
    spike2_temp = spikes2 - trigger(i+1);
    spike2_temp(spike2_temp < 0) = spike2_temp(spike2_temp < 0) + datarun.duration;
    spike2_temp = sort(spike2_temp);
    spike_shifted(i, :) = spike2_temp;
    spikes2_temp_hc = histc(spike2_temp, edges); 
    shift_predictor(i, :) = xcorr(spikes1_hc, spikes2_temp_hc, maxlag);
end

shift_predictor_mean = mean(shift_predictor);
corr = xcorr(spikes1_hc, spikes2_hc, maxlag);
corr_corrected = corr' - shift_predictor_mean;
    