function [theta_e_mean, theta_e_ste] = run_poiss_MLdecoder_megf10(phi, alpha2, r_base, ymax, repeat_n)
%% training
THETA = linspace(0, 2*pi, 100); % direction
% repeat_n = 100;

% Output = zeros(length(THETA), repeat_n, 4);
for dir = 1:length(THETA)
    theta = THETA(dir);
    gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2;
    modelOutput(dir,:) = ymax.*gain + r_base;
end


%% test performance

% direction decoding
clear theta_e theta_e_mean theta_e_ste

for dir = 1:length(THETA)
    for repeat = 1:repeat_n
        theta = THETA(dir);
        gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2;
        output = ymax.*gain + r_base;
        for i = 1:length(output)
            output(i) = random('poiss', output(i));
%             output(i) = random('bino', round(2*output(i)), 0.9);
%             output(i) = round(max(random('norm', output(i), sqrt(1.34*output(i)^0.66)), 0));
        end

        output_r = repmat(output, length(THETA), 1);
        p = exp(-(output_r - modelOutput).^2./(2*modelOutput.^2))./sqrt(2*pi*modelOutput.^2);
%         p = modelOutput.^output_r.*exp(-modelOutput)./factorial(output_r);
        p = prod(p, 2);
        [~,i] = max(p);
        error = abs(THETA(i)-theta)/pi*180;
        if error > 180
            error = 360-error;
        end
        theta_e(repeat,dir) = error;
    end
end
theta_e_mean = mean(theta_e);
theta_e_ste = std(theta_e)/sqrt(repeat_n);

% ERROR = mean(theta_e_mean);
% ERROR = theta_e_mean;
% figure
% errorbar(THETA, theta_e_mean, theta_e_ste, 'color', 'k'); 
% xlim([-0.5 4.5])
% xlabel('Log(luminance)(R*/rod/s)')
% ylabel('mean error (degree)')

end