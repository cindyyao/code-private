spike1_T = datarun.spikes{49};
spike1_TF= ceil(spike1_T/0.01);
spike1 = zeros(spike1_TF(end), 1);
spike1(spike1_TF) = 1;

spike2_T = datarun.spikes{74};
spike2_TF= ceil(spike2_T/0.01);
spike2 = zeros(spike2_TF(end), 1);
spike2(spike2_TF) = 1;

A = xcorr(spike1, spike2, 20);
plot(A)