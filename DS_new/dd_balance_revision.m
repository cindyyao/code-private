%% discrimination and detection task (precision likelihood ratio)

c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
a = 1.8; b = 0.8;
Alpha2 = [a a a a ; b a a a; b b a a; b b b a; b b b b]; % tuning width 
r_base = [1 1 1 1]*0.02; % background firing
ymax = [1 1 1 1]*5;

contrast = [40 50];
repeat_n = 1000;

for simrepeat = 1:10
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
%                             output(i) = random('bino', round(2*output(i)), 0.9);
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
            Precision_ratio_cone{simrepeat}(width, ctr) = sum(correct(:))/repeat_n/length(THETA);
            disp(['repeat ' num2str(simrepeat) ' detection task: tuning width ' num2str(width) ' contrast ' num2str(ctr) ' finished.'])
        end
    end

    %
    for i = 1:5
        alpha2 = Alpha2(i, :);
        for ii = 1:length(contrast)
            % run_MLdecoder
            [error_cone{simrepeat}(i, ii), Theta_e{simrepeat}{i, ii}] = run_poiss_MLdecoder(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ii),repeat_n);
            disp(['repeat ' num2str(simrepeat) ' discrimination task: tuning width ' num2str(i) ' contrast ' num2str(ii) ' finished.'])
        end
    end

    save('dsmodel_temp5.mat', 'Precision_ratio_cone', 'error_cone')
    % dsmodel_temp1.mat: contrast = [40 50 60], r_base = 0.1, a = 2; b = 0.8;
    % dsmodel_temp2.mat: contrast = [40 50 60], r_base = 0.02, a = 2; b = 0.8; 
    % dsmodel_temp3.mat: contrast = [40 50 60], r_base = 1, a = 2; b = 0.8;
    % dsmodel_temp4.mat: contrast = [40 50], r_base = 0.02, a = 1.8; b = 0.8;
    % dsmodel_temp5.mat: contrast = [40 50], r_base = 0.02, a = 1.9; b = 0.8;
end
%%

temp(1, 1, :) = error_cone;
error_cone_all = mean(cell2mat(temp), 3);

temp(1, 1, :) = Precision_ratio_cone;
Precision_cone_all = mean(cell2mat(temp), 3);

error_all = error_cone_all;
precision_all = Precision_cone_all - 0.5;

figure
c1 = 1;
subplot(2, 1, 1)
% yyaxis left
plot(0:4, (90-error_all(:, c1))/max(90-error_all(:, c1)), 'ro-');
hold on
xx = linspace(-1, 5, 100);
yy = (error_all(1, c1)-error_all(2, c1))/max(90-error_all(:, c1))*xx + (90 - error_all(1, c1))/max(90 - error_all(1, c1));
plot(xx, yy, 'r--')
% ylim([0 90])

% yyaxis right
% precision_all = precision_all - 0.5;
plot(0:4, precision_all(:, c1)/max(precision_all(:, c1)), 'bo-')
xx = linspace(-1, 5, 100);
yy = (precision_all(2, c1)-precision_all(1, c1))/max(precision_all(:, c1))*xx + precision_all(1, c1)/max(precision_all(:, c1));
plot(xx, yy, 'b--')
% ylim([0.5 1])
xlim([-0.5 4.5])
