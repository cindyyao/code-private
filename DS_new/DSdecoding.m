load('dsmodel.mat')
% ERROR_rod/ERROR_cone: 7 contrasts X 20 tuning width
% error_rod/error_cone: (5 models X 7 contrasts) X 10 repeats
% Precision_ratio_rod/Precision_ratio_cone: (5 models X 7 contrasts) X 10 repeats
%% optimal tuning width
clear error theta_e
ctr = [5 10 20 40 80 150 300]; % contrast
theta = unique(sort([linspace(10,50,4)/180*pi linspace(50,80,10)/180*pi linspace(80,170,8)/180*pi]));
Alpha2 = log(0.5)./log(0.5+0.5*cos(theta));


for i = 1:length(ctr)
    for j = 1:length(Alpha2)
        % NDF 4
        % CRF 
        c50 = [1.8 1.8 1.8 1.8];
        alpha1 = [7 7 7 7];
        % tuning curve
        phi = [0 1/2 1 3/2]*pi; % preferred direction
        alpha2 = [1 1 1 1]*Alpha2(j); % tuning width 
        r_base = [1 1 1 1]; % background firing
        ymax = [1 1 1 1]*25; % ymax = 5 for rod light level
        % run_MLdecoder
        [error(i,j), theta_e{i,j}] = run_poiss_MLdecoder(c50,alpha1,phi,alpha2,r_base,ymax,ctr(i));
        [i,j]
    end
end 


%% discrimination and detection task (precision likelihood ratio)

c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
a = 2; b = 1.2;
Alpha2 = [a a a a ; b a a a; b b a a; b b b a; b b b b]; % tuning width 
r_base = [1 1 1 1]*0.2; % background firing
ymax = [1 1 1 1]*5;

contrast = [5 10 20 40 80 150 300];
repeat_n = 1000;

for simrepeat = 1:1
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
            Precision_ratio_cone{simrepeat}(width, ctr) = sum(correct(:))/repeat_n/length(THETA);
            disp(['repeat ' num2str(simrepeat) ' detection task: tuning width ' num2str(width) ' contrast ' num2str(ctr) ' finished.'])
        end
    end

    %
    for i = 1:5
        alpha2 = Alpha2(i, :);
        for ii = 1:length(contrast)
            % run_MLdecoder
            [error_cone{simrepeat}(i, ii), Theta_e{simrepeat}{i, ii}] = run_poiss_MLdecoder(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ii));
            disp(['repeat ' num2str(simrepeat) ' discrimination task: tuning width ' num2str(i) ' contrast ' num2str(ii) ' finished.'])
        end
    end

    save('temp.mat', 'Precision_ratio_cone', 'error_cone')
end

%% figure
%% optimal tuning width
figure
subplot(1, 2, 1)
ctr = [5 10 20 40 80 150 300];
ERROR = ERROR_rod;
for i = 1:length(ctr)
    plot(theta/pi*180*2, ERROR(i, :), 'color', [1 1 1]/(i+1));
    hold on
end
[~, optimal] = min(ERROR');
optimalWidth = zeros(length(ctr),1);
for i = 1:length(ctr)
    optimalWidth(i) = ERROR(i,optimal(i));
end
plot(theta(optimal)/pi*180*2, optimalWidth, 'r*')
xlabel('tuning width / degree')
ylabel('prediction error / degree')
ylim([0 100])
xlim([0 400])


subplot(1, 2, 2)
ERROR = ERROR_cone;
for i = 1:length(ctr)
    plot(theta/pi*180*2, ERROR(i, :), 'color', [1 1 1]/(i+1));
    hold on
end
[~, optimal] = min(ERROR');
optimalWidth = zeros(length(ctr),1);
for i = 1:length(ctr)
    optimalWidth(i) = ERROR(i,optimal(i));
end
plot(theta(optimal)/pi*180*2, optimalWidth, 'r*')
xlabel('tuning width / degree')
ylabel('prediction error / degree')
ylim([0 100])
xlim([0 400])

%% task performance
contrast = [5 10 20 40 80 150 300];
rep = 1;
figure
subplot(1, 2, 1)
for i = 1:5
    plot(contrast, error_rod{rep}(i, :))
    hold on 
end
set(gca, 'xscale', 'log')
xlim([4 400])
ylim([0 90])
xlabel('contrast')
ylabel('error(degree)')
legend('0 broad', '1 broad', '2 broad', '3 broad', '4 broad')

subplot(1, 2, 2)
for i = 1:5
    plot(contrast, error_cone{rep}(i, :))
    hold on 
end
set(gca, 'xscale', 'log')
xlim([4 400])
ylim([0 90])
xlabel('contrast')
ylabel('error(degree)')
legend('0 broad', '1 broad', '2 broad', '3 broad', '4 broad')

figure
subplot(1, 2, 1)
for i = 1:5
    plot(contrast, Precision_ratio_rod{rep}(i, :))
    hold on 
end
set(gca, 'xscale', 'log')
xlim([4 400])
ylim([0.5 1])
xlabel('contrast')
ylabel('2AFC performance')
legend('0 broad', '1 broad', '2 broad', '3 broad', '4 broad')

subplot(1, 2, 2)
for i = 1:5
    plot(contrast, Precision_ratio_cone{rep}(i, :))
    hold on 
end
set(gca, 'xscale', 'log')
xlim([4 400])
ylim([0.5 1])
xlabel('contrast')
ylabel('2AFC performance')
legend('0 broad', '1 broad', '2 broad', '3 broad', '4 broad')


%% estimation of ymax from recorded data
load('/Users/xyao/matlab/data/DS160506.mat', 'fittingDs')
alpha = 7; c50 = 1.8; x = log10(80);
f = x^alpha/(x^alpha + c50^alpha);
index = [1 5];
for i = 1:2
    ymax_temp = [];
    for dir = 1:4
        ymax_temp = [ymax_temp fittingDs.ymax{dir}{index(i)}];
    end
    ymax_data{i} = ymax_temp/f;
end
   
figure
scatter(ymax_data{1}, ymax_data{2}, 'MarkerEdgeColor', [1 1 1]*0.7);
hold on
plot([0 100], [0 100], 'k--')
errorbar(mean(ymax_data{1}), mean(ymax_data{2}), std(ymax_data{2})/sqrt(length(ymax_data{2})), 'ro')
herrorbar(mean(ymax_data{1}), mean(ymax_data{2}), std(ymax_data{1})/sqrt(length(ymax_data{1})), 'ro')
xlim([0 55])
ylim([0 55])
xlabel('ymax at rod light level')
ylabel('ymax at cone light level')

%% heatmap
alpha = 7; c50 = 1.8; phi = 0; alpha2 = 2; r_base = 0;
contrast = 10.^linspace(log10(1), log10(1000), 100);
direction = linspace(-pi, pi, 100);
[xx, yy] = meshgrid(direction, log10(contrast));
zz = yy.^alpha./(yy.^alpha + c50.^alpha).*(0.5 + 0.5 * cos(xx + phi)).^alpha2 + r_base;
figure
s = surf(xx, 10.^yy, zz);
s.EdgeColor = 'none';
set(gca, 'yscale', 'log')
xlim([-pi pi])
ylim([min(contrast) max(contrast)])
xlabel('direction (rad)')
ylabel('contrast (%)')
caxis([0 1])
colorbar

%% discrimination model
direction = linspace(0, 2*pi, 100); phi = pi; r_base = 0.2; alpha2 = 2;
figure
for dir = 1:4
    y = 10 * (0.5 + 0.5 * cos(direction - phi)).^alpha2 + r_base;
    plot(direction, y);
    hold on
    direction = direction + pi/2;
    phi = phi + pi/2;
end
%% percentage improvement

% load('dsDetection.mat')
contrast = [5 10 20 40 80 150 300];
error = error_cone;
Precision_ratio = Precision_ratio_cone;

rep = 1;
figure
for i = 1:5
    plot(contrast, error{rep}(i, :))
    hold on 
end
set(gca, 'xscale', 'log')
xlim([4 400])
xlabel('contrast')
ylabel('error(degree)')
legend('0 broad', '1 broad', '2 broad', '3 broad', '4 broad')

figure
for i = 1:5
    plot(contrast, Precision_ratio{rep}(i, :))
    hold on 
end
set(gca, 'xscale', 'log')
xlim([4 400])
xlabel('contrast')
ylabel('2AFC performance')
legend('0 broad', '1 broad', '2 broad', '3 broad', '4 broad')


for rep = 1:1
    error_total_area = sum(error{rep}(5, :) - error{rep}(1, :));
    detection_total_area = sum(Precision_ratio{rep}(5, :) - Precision_ratio{rep}(1, :));

    for i = 1:4
        error_improv_rate(rep, i) = sum(error{rep}(i + 1, :) - error{rep}(i, :))/error_total_area;
        detection_improv_rate(rep, i) = sum(Precision_ratio{rep}(i + 1, :) - Precision_ratio{rep}(i, :))/detection_total_area;
    end
end

ratio = detection_improv_rate./error_improv_rate;

figure
errorbar(1:4, mean(error_improv_rate), std(error_improv_rate)/sqrt(10))
xlabel('# of broad tuning curve')
ylabel('% worse in discrimination task (beta)')
xlim([0 5])

figure
errorbar(1:4, mean(detection_improv_rate), std(detection_improv_rate)/sqrt(10))
xlabel('# of broad tuning curve')
ylabel('% better in detection task (alpha)')
xlim([0 5])


figure
errorbar(1:4, mean(ratio), std(ratio)/sqrt(10))
hold on
plot([0 5], [1 1], 'k--')
% ylim([0.5 1.7])
xlabel('# of broad tuning curve')
ylabel('improvemnet index (alpha/beta)')
xlim([0.5 4.5])

%%
for i = 1:10
    error_rod_40(i, :) = error_rod{i}(:, 4);
    Precision_rod_40(i, :) = Precision_ratio_rod{i}(:, 4);
end

figure
errorbar(mean(error_rod_40), mean(Precision_rod_40), std(Precision_rod_40)/sqrt(10), 'k');
hold on
herrorbar(mean(error_rod_40), mean(Precision_rod_40), std(error_rod_40)/sqrt(10), 'k');
ylim([0.86 0.96])
xlim([40 50])
xlabel('prediction error (deg)')
ylabel('correct rate')
title('n = 3')