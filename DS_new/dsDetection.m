%%
c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
Alpha2 = [0.8 0.8 0.8 0.8; 0.8 2 2 2; 2 2 2 2]; % tuning width 
r_base = [1 1 1 1]*0.2; % background firing
ymax = [1 1 1 1]*5;

contrast = [5 10 20 40 80 150 300];
repeat_n = 1000;

for width = 1:size(Alpha2, 1)
    alpha2 = Alpha2(width, :);
    for ctr = 1:length(contrast)
        c = contrast(ctr);
        x(1) = 0;
        x(2) = log10(c);
        THETA = linspace(0, 2*pi, 50); % direction
%         for dir = 1:length(THETA)
%             theta = THETA(dir);
%             input = ymax.* x(1).^alpha1./(x(1).^alpha1 + c50.^alpha1);
%             gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2;
%             modelOutput(dir,:) = input.*gain + r_base;
%         end

        correct = zeros(repeat_n, length(THETA));

        for dir = 1:length(THETA)
            for repeat = 1:repeat_n
                for j = 1:2
                    theta = THETA(dir);
                    input = ymax.* x(j).^alpha1./(x(j).^alpha1 + c50.^alpha1);
                    gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2;
                    output = input.*gain + r_base;
                    for i = 1:length(output)
                        output(i) = random('poiss', output(i));
                    end

%                     output_r = repmat(output, length(THETA), 1);
            %         p = exp(-(output_r - Output_mean).^2./(2*Output_std.^2))./sqrt(2*pi*Output_std.^2);
                    p = r_base.^output.*exp(-r_base)./factorial(output);
                    p = prod(p, 2);
%                     p = p/sum(p);
%                     l(j) = max(p);

                    l(j) = p;
                end
                if l(1) > l(2)
                    correct(repeat, dir) = 1;
                elseif l(1) == l(2)
                    correct(repeat, dir) = 0.5;
%                 else
%                     pause
                end
            end
        end
        Precision(width, ctr) = sum(correct(:))/repeat_n/length(THETA);
        ctr
    end
end
figure
for i = 1:3
plot(log10(contrast), Precision(i, :))
hold on
end
legend('all broad', 'one broad', 'all narrow')
xlabel('log(contrast)')
ylabel('Percentage correct')
ylim([0.5 1])


%% likelihood ratio

c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
Alpha2 = [0.8 0.8 0.8 0.8; 0.8 0.8 0.8 2; 0.8 0.8 2 2; 0.8 2 2 2; 2 2 2 2]; % tuning width 
r_base = [1 1 1 1]*0.2; % background firing
ymax = [1 1 1 1]*5;

contrast = [5 10 20 40 80 150 300];
repeat_n = 1000;

for width = 1:size(Alpha2, 1)
    alpha2 = Alpha2(width, :);
    for ctr = 1:length(contrast)
        c = contrast(ctr);
        x(1) = 0;
        x(2) = log10(c);
        THETA = linspace(0, 2*pi, 50); % direction
        for dir = 1:length(THETA)
            theta = THETA(dir);
            input = ymax.* x(2).^alpha1./(x(2).^alpha1 + c50.^alpha1);
            gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2;
            modelOutput(dir,:) = input.*gain + r_base;
        end

        correct = zeros(repeat_n, length(THETA));

        for dir = 1:length(THETA)
            for repeat = 1:repeat_n
                for j = 1:2
                    theta = THETA(dir);
                    input = ymax.* x(j).^alpha1./(x(j).^alpha1 + c50.^alpha1);
                    gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2;
                    output = input.*gain + r_base;
                    for i = 1:length(output)
                        output(i) = random('poiss', output(i));
                    end

                    output_r = repmat(output, length(THETA), 1);
            %         p = exp(-(output_r - Output_mean).^2./(2*Output_std.^2))./sqrt(2*pi*Output_std.^2);
                    p = r_base.^output.*exp(-r_base)./factorial(output);
                    p_null = prod(p, 2);
                    p = modelOutput.^output_r.*exp(-modelOutput)./factorial(output_r);
                    p_stim = max(prod(p, 2));
                    l(j) = p_stim/p_null;
%                     p = p/sum(p);
%                     l(j) = max(p);
                end
                if l(1) < l(2)
                    correct(repeat, dir) = 1;
                elseif l(1) == l(2)
                    correct(repeat, dir) = 0.5;
%                 else
%                     pause
                end
            end
        end
        Precision_temp(width, ctr) = sum(correct(:))/repeat_n/length(THETA);
        ctr
    end
end
figure
for i = 1:5
plot(log10(contrast), Precision(i, :))
hold on
end
% legend('all broad', 'one broad', 'all narrow')
xlabel('log(contrast)')
ylabel('Percentage correct')
ylim([0.5 1])


