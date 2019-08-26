clear SpikeN_mean_all SpikeN_std_all
for ll = 1:3
    SpikeN_mean_all(:, ll) = cat(1,SpikeN_mean{ll}{:});
    SpikeN_std_all(:, ll) = cat(1,SpikeN_std{ll}{:});
end
SpikeN_mean_all(any(cellfun(@isempty, SpikeN_mean_all), 2), :) = [];
% SpikeN_mean_all = cell2mat(SpikeN_mean_all);
SpikeN_std_all(any(cellfun(@isempty, SpikeN_std_all), 2), :) = [];
% SpikeN_std_all = cell2mat(SpikeN_std_all);
% SNR = SpikeN_mean_all./SpikeN_std_all;
% SNR = exciseRows_empty(SNR);


[rmax, i] = cellfun(@max, SpikeN_mean_all);
for ll = 1:3
    for cc = 1:size(SpikeN_mean_all, 1)
        stdmax(cc, ll) = SpikeN_std_all{cc, ll}(i(cc, ll));
    end
end

SNRmax = rmax./stdmax;
SNRmax = exciseRows_empty(SNRmax);
figure
plot(SNRmax(:, 1), SNRmax(:, 3), 'ko')
hold on
plot([0 9], [0 9], 'k--')
errorbar(mean(SNRmax(:, 1)), mean(SNRmax(:, 3)), std(SNRmax(:, 3)/sqrt(size(SNRmax, 1))), 'ro')
herrorbar(mean(SNRmax(:, 1)), mean(SNRmax(:, 3)), std(SNRmax(:, 1)/sqrt(size(SNRmax, 1))), 'ro')

xlim([0 9])
ylim([0 9])
xlabel('low light level SNR')
ylabel('high light level SNR')

%%
load('dsmodel-1.mat')
load('dsmodel.mat', 'ERROR_cone')
load('dsmodel170504.mat', 'Precision_ratio_cone')
ERROR_cone_mean = ERROR_cone';
Precision_cone_mean = Precision_ratio_cone;

ymax = 5; alpha1 = 7; c50 = 1.8;
ctr = log10([10 20 40 80 150 300]);
input1 = ymax.* ctr.^alpha1./(ctr.^alpha1 + c50.^alpha1);

ymax = 25; alpha1 = 7; c50 = 1.8;
ctr = log10([5 10 20 40 80 150 300]);
input2 = ymax.* ctr.^alpha1./(ctr.^alpha1 + c50.^alpha1);

[input, i] = sort([input1 input2]);
ERROR = [ERROR_rod_mean ERROR_cone'];
ERROR = ERROR(:, i);
ERROR = ERROR(:, 4:end);

Precision_ratio = [Precision_rod_mean Precision_ratio_cone];
Precision_ratio = Precision_ratio(:, i);
Precision_ratio = Precision_ratio(:, 4:end);

theta = unique(sort([linspace(10,50,4)/180*pi linspace(50,80,10)/180*pi linspace(80,170,8)/180*pi]));
[~, optimal] = min(ERROR);
for i = 1:size(ERROR,2)
    optimalWidth(i) = ERROR(optimal(i), i);
    optimaPrecision(i) = Precision_ratio(optimal(i), i);
end

figure
c1 = 3;
subplot(1, 2, 1)
[ax, h1, h2] = plotyy(theta/pi*360, ERROR(:, c1), theta/pi*360, Precision_ratio(:, c1));
hold on
plot(theta(optimal(c1))/pi*180*2, optimalWidth(c1), 'r*')
set(ax(1),'ylim',[0 90])
set(ax(2),'ylim',[0.5 1])
set(ax(1),'xlim',[0 360])
set(ax(2),'xlim',[0 360])

c2 = 9;
subplot(1, 2, 2)
[ax, h1, h2] = plotyy(theta/pi*360, ERROR(:, c2), theta/pi*360, Precision_ratio(:, c2));
hold on
plot(theta(optimal(c2))/pi*180*2, optimalWidth(c2), 'r*')
set(ax(1),'ylim',[0 90])
set(ax(2),'ylim',[0.5 1])
set(ax(1),'xlim',[0 360])
set(ax(2),'xlim',[0 360])



figure
c1 = 3;
subplot(1, 2, 1)
yyaxis left
plot(theta/pi*360, ERROR(:, c1), 'k');
hold on
plot(theta(optimal(c1))/pi*180*2, optimalWidth(c1), 'r*')
ylim([0 90])
yyaxis right
plot(theta/pi*360, Precision_ratio(:, c1), 'k');
ylim([0.5 1])
xlim([0 360])

c1 = 9;
subplot(1, 2, 2)
yyaxis left
plot(theta/pi*360, ERROR(:, c2), 'k');
hold on
plot(theta(optimal(c2))/pi*180*2, optimalWidth(c2), 'r*')
ylim([0 90])
yyaxis right
plot(theta/pi*360, Precision_ratio(:, c2), 'k');
ylim([0.5 1])
xlim([0 360])

load('DS150603.mat', 'Width')
width_temp = cat(1,Width{:});
for i = 1:5
    width_temp{2, i} = cell2mat(width_temp(2:4, i)');
end
width = width_temp(1:2, :);
width_mean = cellfun(@mean, width);
width_ste = cellfun(@std, width)./sqrt(cellfun(@length, width));


SNRmax_mean = mean(SNRmax);
SNRmax_ste = std(SNRmax)/sqrt(size(SNRmax, 1));
SNRmax_mean(2) = [];
SNRmax_ste(2) = [];
figure
plot(sqrt(input(4:end)), theta(optimal)/pi*180*2, 'k')
hold on
plot(sqrt(input(4:end)), ones(10, 1)*360, 'k')
errorbar(reshape(repmat(SNRmax_mean, 2, 1), 4, 1), [width_mean(:, 1); width_mean(:, 5)], [width_ste(:, 1); width_ste(:, 5)], 'ko')
herrorbar(reshape(repmat(SNRmax_mean, 2, 1), 4, 1), [width_mean(:, 1); width_mean(:, 5)], reshape(repmat(SNRmax_ste, 2, 1), 4, 1), 'ko')

xlabel('SNR')
ylabel('optimal tuning width')
ylim([90 380])

load dsmodel.mat
temp(1, 1, :) = error_rod;
error_rod_all = mean(cell2mat(temp), 3);
temp(1, 1, :) = error_cone;
error_cone_all = mean(cell2mat(temp), 3);
temp(1, 1, :) = Precision_ratio_rod;
Precision_rod_all = mean(cell2mat(temp), 3);
temp(1, 1, :) = Precision_ratio_cone;
Precision_cone_all = mean(cell2mat(temp), 3);


ymax = 5; alpha1 = 7; c50 = 1.8;
ctr = log10([5 10 20 40 80 150 300]);
input1 = ymax.* ctr.^alpha1./(ctr.^alpha1 + c50.^alpha1);

ymax = 30; alpha1 = 7; c50 = 1.8;
ctr = log10([5 10 20 40 80 150 300]);
input2 = ymax.* ctr.^alpha1./(ctr.^alpha1 + c50.^alpha1);

[input, i] = sort([input1 input2]);
error_all = [error_rod_all error_cone_all];
error_all = error_all(:, i);
precision_all = [Precision_rod_all Precision_cone_all];
precision_all = precision_all(:, i);

%%
h = figure;
c1 = 6;
subplot(2, 1, 1)
yyaxis left
plot(0:4, 90-error_all(:, c1), 'bo-');
hold on
xx = linspace(-1, 5, 100);
yy = (error_all(1, c1)-error_all(2, c1))*xx - error_all(1, c1) + 90;
plot(xx, yy, 'b--')
ylim([0 90])

yyaxis right
plot(0:4, precision_all(:, c1), 'ro-')
xx = linspace(-1, 5, 100);
yy = (precision_all(2, c1)-precision_all(1, c1))*xx + precision_all(1, c1);
plot(xx, yy, 'r--')
ylim([0.5 1])
xlim([-0.5 4.5])

c1 = 11;
subplot(2, 1, 2)
yyaxis left
plot(0:4, 90-error_all(:, c1), 'bo-');
hold on
xx = linspace(-1, 5, 100);
yy = (error_all(1, c1)-error_all(2, c1))*xx - error_all(1, c1) + 90;
plot(xx, yy, 'b--')
ylim([0 90])

yyaxis right
plot(0:4, precision_all(:, c1), 'ro-')
xx = linspace(-1, 5, 100);
yy = (precision_all(2, c1)-precision_all(1, c1))*xx + precision_all(1, c1);
plot(xx, yy, 'r--')
ylim([0.5 1])
xlim([-0.5 4.5])

%%
h = figure;
c1 = 6;
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

c1 = 11;
subplot(2, 1, 2)
% yyaxis left
plot(0:4, (90-error_all(:, c1))/max(90-error_all(:, c1)), 'ro-');
hold on
xx = linspace(-1, 5, 100);
yy = (error_all(1, c1)-error_all(2, c1))/max(90-error_all(:, c1))*xx + (90 - error_all(1, c1))/max(90 - error_all(1, c1));
plot(xx, yy, 'r--')
% ylim([0 90])

% yyaxis right
plot(0:4, precision_all(:, c1)/max(precision_all(:, c1)), 'bo-')
xx = linspace(-1, 5, 100);
yy = (precision_all(2, c1)-precision_all(1, c1))/max(precision_all(:, c1))*xx + precision_all(1, c1)/max(precision_all(:, c1));
plot(xx, yy, 'b--')
% ylim([0.5 1])
xlim([-0.5 4.5])

%% 

thresh = size(Precision_ratio, 1) - sum(Precision_ratio>0.99) + 1;
% theta = [theta 2*pi];
figure
plot(sqrt(input(4:end)), theta(optimal)/pi*180*2, 'k')
hold on
plot(sqrt(input(4:end)), ones(10, 1)*360, 'k')
plot(sqrt(input(4:end)), theta(thresh)/pi*180, 'k')

errorbar(reshape(repmat(SNRmax_mean, 2, 1), 4, 1), [width_mean(:, 1); width_mean(:, 5)], [width_ste(:, 1); width_ste(:, 5)], 'ko')
herrorbar(reshape(repmat(SNRmax_mean, 2, 1), 4, 1), [width_mean(:, 1); width_mean(:, 5)], reshape(repmat(SNRmax_ste, 2, 1), 4, 1), 'ko')

xlabel('SNR')
ylabel('optimal tuning width')
% ylim([90 380])


%%
contrast = {'5%', '10%', '20%', '40%', '80%', '150%', '300%'};
h = figure;
for c1 = 1:7
    subplot(2, 7, c1)
    % yyaxis left
    plot(0:4, (90-error_rod_all(:, c1))/max(90-error_rod_all(:, c1)), 'bo-');
    hold on
    xx = linspace(-1, 5, 100);
    yy = (error_rod_all(1, c1)-error_rod_all(2, c1))/max(90-error_rod_all(:, c1))*xx + (90 - error_rod_all(1, c1))/max(90 - error_rod_all(:, c1));
    plot(xx, yy, 'b--')
    % ylim([0 90])
    title(contrast{c1})
    % yyaxis right
    plot(0:4, Precision_rod_all(:, c1)/max(Precision_rod_all(:, c1)), 'ro-')
    xx = linspace(-1, 5, 100);
    yy = (Precision_rod_all(2, c1)-Precision_rod_all(1, c1))/max(Precision_rod_all(:, c1))*xx + Precision_rod_all(1, c1)/max(Precision_rod_all(:, c1));
    plot(xx, yy, 'r--')
    % ylim([0.5 1])
    xlim([-0.5 4.5])

    subplot(2, 7, c1+7)
    % yyaxis left
    plot(0:4, (90-error_cone_all(:, c1))/max(90-error_cone_all(:, c1)), 'bo-');
    hold on
    xx = linspace(-1, 5, 100);
    yy = (error_cone_all(1, c1)-error_cone_all(2, c1))/max(90-error_cone_all(:, c1))*xx + (90 - error_cone_all(1, c1))/max(90 - error_cone_all(:, c1));
    plot(xx, yy, 'b--')
    % ylim([0 90])

    % yyaxis right
    plot(0:4, Precision_cone_all(:, c1)/max(Precision_cone_all(:, c1)), 'ro-')
    xx = linspace(-1, 5, 100);
    yy = (Precision_cone_all(2, c1)-Precision_cone_all(1, c1))/max(Precision_cone_all(:, c1))*xx + Precision_cone_all(1, c1)/max(Precision_cone_all(:, c1));
    plot(xx, yy, 'r--')
    % ylim([0.5 1])
    xlim([-0.5 4.5])
end