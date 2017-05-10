function [theta_e] = run_MLdecoder(c50, alpha1, sigma1, phi, alpha2, r_base, sigma2, ymax, ctr)
%% training
THETA = linspace(0, 2*pi, 100); % direction
repeat_n = 100;

Output = zeros(length(THETA), repeat_n, 4);
for dir = 1:length(THETA)
    for repeat = 1:repeat_n
        theta = THETA(dir);
        input = ymax.*(ctr.^alpha1./(ctr.^alpha1 + c50.^alpha1) + sigma1.*random('norm', 0, 1, [1, 4]));
%         gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2+r_base;
%         output = max(max(input.*gain, 0) + sigma2.*random('norm', 0, 1, [1, 4]), 0);

        gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2;
        output = max(max(input.*gain, 0) + r_base + sigma2.*random('norm', 0, 1, [1, 4]), 0);
        Output(dir, repeat, :) = floor(output);
    end
end

Output_mean = squeeze(mean(Output, 2));
Output_std = squeeze(std(Output, [], 2));
Output_var = squeeze(var(Output, [], 2));
Output_ste = Output_std/sqrt(repeat_n);


figure
for cc = 1:4
    errorbar(THETA, Output_mean(:, cc), Output_ste(:, cc))
    hold on
end
xlim([min(THETA)-0.1 max(THETA)+0.1])
xlabel('direction')
ylabel('spike #')
figure
scatter(Output_mean(:), Output_var(:))
xlabel('mean')
ylabel('variance')

%% test performance

% direction decoding
clear theta_e theta_e_mean theta_e_ste

repeat_n = 100;
for dir = 1:length(THETA)
    for repeat = 1:repeat_n
        theta = THETA(dir);
        input = ymax.*(ctr.^alpha1./(ctr.^alpha1 + c50.^alpha1) + sigma1.*random('norm', 0, 1, [1, 4]));
        gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2;
        output = max(max(input.*gain, 0) + r_base + sigma2.*random('norm', 0, 1, [1, 4]), 0);
%         gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2+r_base;
%         output = max(max(input.*gain, 0) + sigma2.*random('norm', 0, 1, [1, 4]), 0);

        output_r = repmat(output, length(THETA), 1);
        p = exp(-(output_r - Output_mean).^2./(2*Output_std.^2))./sqrt(2*pi*Output_std.^2);
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

mean(theta_e_mean)
figure
errorbar(THETA, theta_e_mean, theta_e_ste, 'color', 'k'); 
% xlim([-0.5 4.5])
% xlabel('Log(luminance)(R*/rod/s)')
% ylabel('mean error (degree)')

end