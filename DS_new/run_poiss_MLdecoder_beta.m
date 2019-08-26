function [ERROR, theta_e_mean] = run_poiss_MLdecoder_beta(c50, alpha1, phi, alpha2, r_base, ymax, c, repeat_n, beta)
%% training
x = log10(c);
THETA = linspace(0, 2*pi, 17); % direction
THETA = THETA(1:end-1);
% repeat_n = 100;

% Output = zeros(length(THETA), repeat_n, 4);
for dir = 1:length(THETA)
    theta = THETA(dir);
    input = ymax.* x.^alpha1./(x.^alpha1 + c50.^alpha1);
    gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2.*beta+(1-beta);
    modelOutput(dir,:) = input.*gain + r_base;
%     for repeat = 1:repeat_n
%         for i = 1:size(modelOutput,2)
%             output(i) = random('poiss', modelOutput(dir, i));
%         end
%         Output(dir, repeat, :) = floor(output);
%     end
end

% Output_mean = squeeze(mean(Output, 2));
% Output_std = squeeze(std(Output, [], 2));
% Output_var = squeeze(var(Output, [], 2));
% Output_ste = Output_std/sqrt(repeat_n);


% figure
% for cc = 1:4
%     plot(THETA, modelOutput(:, cc))
%     hold on
% end
% xlim([min(THETA)-0.1 max(THETA)+0.1])
% xlabel('direction')
% ylabel('spike #')
% figure
% scatter(Output_mean(:), Output_var(:))
% xlabel('mean')
% ylabel('variance')

%% test performance

% direction decoding
clear theta_e theta_e_mean theta_e_ste

for repeat = 1:repeat_n
    for dir = 1:length(THETA)
        theta = THETA(dir);
        input = ymax.* x.^alpha1./(x.^alpha1 + c50.^alpha1);
        gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2.*beta+(1-beta);
        output = input.*gain + r_base;
        for i = 1:length(output)
            output(i) = random('poiss', output(i));
%             output(i) = random('bino', round(2*output(i)), 0.9);
        end

        output_r = repmat(output, length(THETA), 1);
%         p = exp(-(output_r - Output_mean).^2./(2*Output_std.^2))./sqrt(2*pi*Output_std.^2);
        p = modelOutput.^output_r.*exp(-modelOutput)./factorial(output_r);
        p = prod(p, 2);
        p = p + randn(size(p))*10^-16;
        [~,i] = max(p);
        prediction(repeat,dir) = THETA(i);
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
% figure
% errorbar(THETA, theta_e_mean, theta_e_ste, 'color', 'k'); 
% xlim([-0.5 4.5])
% xlabel('Log(luminance)(R*/rod/s)')
% ylabel('mean error (degree)')

end