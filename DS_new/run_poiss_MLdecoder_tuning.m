function [ERROR, theta_e] = run_poiss_MLdecoder_tuning(phi, alpha2, r_base, ymax)
%% training
THETA = linspace(0, 2*pi, 100); % direction
% repeat_n = 100;

% Output = zeros(length(THETA), repeat_n, 4);
for dir = 1:length(THETA)
    theta = THETA(dir);
    gain = ymax.*(0.5 + 0.5 * cos(theta + phi)).^alpha2;
    modelOutput(dir,:) = gain + r_base;
end

%% test performance

% direction decoding
clear theta_e theta_e_mean theta_e_ste

repeat_n = 1000;
for dir = 1:length(THETA)
    for repeat = 1:repeat_n
        theta = THETA(dir);
        gain = ymax.*(0.5 + 0.5 * cos(theta + phi)).^alpha2;
        output = gain + r_base;
        for i = 1:length(output)
            output(i) = random('poiss', output(i));
        end

        output_r = repmat(output, length(THETA), 1);
%         p = exp(-(output_r - Output_mean).^2./(2*Output_std.^2))./sqrt(2*pi*Output_std.^2);
        p = modelOutput.^output_r.*exp(-modelOutput)./factorial(output_r);
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

ERROR = mean(theta_e_mean);

end