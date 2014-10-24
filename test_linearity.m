function test_linearity(obj)

% not working properly, need to figure out the meaning of the three gamma
% parameters for each gun.

repetition = size(obj.optometer_raw_data, 1)/51;
optometer_readout = obj.optometer_raw_data(:, 4);
optometer_readout = reshape(optometer_readout, 17, 3, repetition);
if repetition > 1
    optometer_readout_avg = mean(optometer_readout, 3);
end
light_intensity = linspace(0, 1, 17);
params = obj.red_params;
fit{1} = params(1) + params(2)*light_intensity.^params(3);
params = obj.green_params;
fit{2} = params(1) + params(2)*light_intensity.^params(3);
params = obj.blue_params;
fit{3} = params(1) + params(2)*light_intensity.^params(3);

for i = 1:3
    figure(i)
    plot(light_intensity, optometer_readout_avg(:, 1), 'k')
    hold on
    plot(light_intensity, fit{i}, 'r')
    legend('raw data', 'fit')
end


