function [spike_time, spike_count] = get_spike_time(V_trace, V_depolarize, samplerate)

spike_time = cell(size(V_trace, 1), 1);
spike_count = zeros(size(V_trace, 1), 1);
for cc = 1:size(V_trace, 1)
    spike_time{cc} = (find(V_trace(cc, :) == V_depolarize) - 1)/samplerate;
    spike_count(cc) = length(spike_time{cc});
end

end

