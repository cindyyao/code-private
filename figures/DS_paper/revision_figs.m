load('DS161017.mat', 'response_max_norm')
drug = 1;
response_max_norm = response_max_norm{drug};
response_max_norm{2} = [response_max_norm{2}; response_max_norm{3}; response_max_norm{4}];
response_max_norm(3:4) = [];
ctr = 5; % contrast = 80%
for i = 1:2
    response_mean_ndf0(i, :) = mean(response_max_norm{i}(:, :, ctr));
    response_ste_ndf0(i, :) = std(response_max_norm{i}(:, :, ctr))/sqrt(size(response_max_norm{i}, 1));
end

load('DS170629.mat', 'response_max_norm_mean', 'response_max_norm_ste')
response_max_norm_mean = response_max_norm_mean{drug};
response_max_norm_ste = response_max_norm_ste{drug};
ctr = 5; % contrast = 80%
for i = 1:2
    response_mean_ndf3(i, :) = response_max_norm_mean{i}{ctr};
    response_ste_ndf3(i, :) = response_max_norm_ste{i}{ctr};
end

theta = linspace(-pi, pi, 9);
theta = theta(2:9)/pi*180;
figure
for i = 1:2
    subplot(1, 2, i)
    errorbar(theta, response_mean_ndf3(i, :), response_ste_ndf3(i, :), 'k')
    hold on
    errorbar(theta, response_mean_ndf0(i, :), response_ste_ndf0(i, :), 'color', [1 1 1]*.5)
    legend('NDF 3', 'NDF 0')
    xlabel('direction')
    ylabel('response')
    if i == 1
        title('superior')
    else
        title('others')
    end
    xlim([-200 200])
end

figure
subplot(1,2,1)
for i = 1:2
    errorbar(theta, response_mean_ndf3(i, :), response_ste_ndf3(i, :))
    hold on
end
subplot(1,2,2)
for i = 1:2
    errorbar(theta, response_mean_ndf0(i, :), response_ste_ndf0(i, :))
    hold on
end
%% fit tuning curve (moving bar)

load('DS161017.mat', 'response_max_norm')
drug = 1;
response_max_norm = response_max_norm{drug};
response_max_norm{2} = [response_max_norm{2}; response_max_norm{3}; response_max_norm{4}];
response_max_norm(3:4) = [];
response_max_norm_ll{2} = response_max_norm;

% load('DS170629.mat', 'response_max_norm')
% drug = 1;
% response_max_norm_ll{1} = response_max_norm{drug};

load('DS170629.mat', 'response_max_norm')
response_max_norm_1 = response_max_norm;
load('DS180413.mat', 'response_max_norm')
response_max_norm_2 = response_max_norm;
for dir = 1:2
    response_max_norm_ll{1}{dir} = [response_max_norm_1{drug}{dir}; response_max_norm_2{drug}{dir}];
end

response_max_norm_ll{1}{2}([14 15 18 35 43 46], :, :) = [];
response_max_norm_ll{2}{2}([3 43], :, :) = [];

%%
ctr = 4;

for dir = 1:2
    figure
    for i = 1:2
        for cc = 1:size(response_max_norm_ll{i}{dir}, 1)
            xdata = 0:pi/4:7*pi/4;
            ydata = response_max_norm_ll{i}{dir}(cc, :, ctr);
            
            [f, g] = fit_cos_(xdata, ydata);
            Ymax{dir}{i}(cc) = f.ymax;
            Phi{dir}{i}(cc) = f.phi;
            Alpha{dir}{i}(cc) = f.alpha;
            Width{dir}{i}(cc) = acos(2*0.5^(1/f.alpha) - 1)/pi*360;
            B{dir}{i}(cc) = f.b;

            xfit = linspace(0, 2*pi, 100);
            yfit = f.ymax.*((0.5+0.5*cos(xfit+f.phi)).^f.alpha.*f.b+1-f.b);
            subplot(6, 9, cc)
            plot(xdata, ydata, 'b')
            hold on
            plot(xfit, yfit, 'r')
            ylim([0 1.1])
            width = acos(2 * (0.5.^(1/f.alpha) - 0.5))/pi*360;
            Width{dir}{i}(cc) = width;
            title(['width = ', num2str(width)])


        end
        Ymax_mean(dir, i) = mean(Ymax{dir}{i});
        Phi_mean(dir, i) = mean(Phi{dir}{i});
        Alpha_mean(dir, i) = mean(Alpha{dir}{i});
        Alpha_ste(dir, i) = std(Alpha{dir}{i})/sqrt(length(Alpha{dir}{i}));
        B_mean(dir, i) = mean(B{dir}{i});
        B_ste(dir, i) = std(B{dir}{i})/sqrt(length(B{dir}{i}));
        Width_mean(dir, i) = mean(Width{dir}{i});
        Width_ste(dir, i) = std(Width{dir}{i})/sqrt(length(B{dir}{i}));
    end
end


ct = {'superior', 'others'};
marker = 'xo';
figure
for i = 1:2
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
errorbar(0:1, B_mean(i, :), B_ste(i, :), 'Color', 'k', 'Marker', marker(i), 'MarkerSize', 10);
hold on
end
xlim([-0.2 1.2])
ylim([0 1.1])
legend(ct)
xlabel('light level')
ylabel('beta')
title('moving bar (figS5)')


ct = {'superior', 'others'};
marker = 'xo';
figure
for i = 1:2
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
errorbar(0:1, Alpha_mean(i, :), Alpha_ste(i, :), 'Color', 'k', 'Marker', marker(i), 'MarkerSize', 10);
hold on
end
xlim([-0.2 1.2])
% ylim([0 1.1])
legend(ct)
xlabel('light level')
ylabel('alpha')
title('moving bar (figS5)')

ct = {'superior', 'others'};
marker = 'xo';
figure
for i = 1:2
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
errorbar(0:1, Width_mean(i, :), Width_ste(i, :), 'Color', 'k', 'Marker', marker(i), 'MarkerSize', 10);
hold on
end
xlim([-0.2 1.2])
% ylim([0 1.1])
legend(ct)
xlabel('light level')
ylabel('width')
title('moving bar (figS5)')
ylim([0 250])
%%
load('DS150603.mat', 'bgfr_ct')
bgfr_ct_wt = bgfr_ct;
for i = 1:5
    bgfr_ct_wt{i}{2} = cell2mat(bgfr_ct_wt{i}(2:4));
    bgfr_ct_wt{i}(3:4) = [];
    for j = 1:2
        bgfr_ct_wt_mean(i, j) = mean(bgfr_ct_wt{i}{j});
        bgfr_ct_wt_ste(i, j) = std(bgfr_ct_wt{i}{j})/sqrt(length(bgfr_ct_wt{i}{j}));
    end
end
load('DS161208.mat', 'bgfr_ct')
bgfr_ct_1208 = bgfr_ct;
load('DS161212.mat', 'bgfr_ct')
bgfr_ct_1212 = bgfr_ct;
for i = 1:5
    bgfr_ct_1208{i}{2} = cell2mat(bgfr_ct_1208{i}(2:4));
    bgfr_ct_1208{i}(3:4) = [];
    bgfr_ct_1212{i}{2} = cell2mat(bgfr_ct_1212{i}(2:4));
    bgfr_ct_1212{i}(3:4) = [];
    for j = 1:2
        bgfr_ct_ko{i}{j} = [bgfr_ct_1208{i}{j} bgfr_ct_1212{i}{j}];
        bgfr_ct_ko_mean(i, j) = mean(bgfr_ct_ko{i}{j});
        bgfr_ct_ko_ste(i, j) = std(bgfr_ct_ko{i}{j})/sqrt(length(bgfr_ct_ko{i}{j}));
    end
end


figure
for i = 1:2
    errorbar(10.^[0:4], bgfr_ct_wt_mean(:, i), bgfr_ct_wt_ste(:, i));
    hold on
    errorbar(10.^[0:4], bgfr_ct_ko_mean(:, i), bgfr_ct_ko_ste(:, i), '--');
end
xlabel('R*/rod/s')
ylabel('background firing rate (HZ)')
legend('WT superior', 'KO superior', 'WT others', 'KO others')
set(gca, 'XScale', 'log')
xlim([10^-0.5 10^4.5])

%%
load('/Users/xyao/matlab/code-private/DS_new/dsmodel_rod.mat')
precision = mean(mean(cell2mat(Precision_all), 2), 3) - 0.5;
theta_temp = cellfun(@cell2mat, theta_e_all, 'UniformOutput', 0);
theta = mean(mean(cell2mat(theta_temp), 2), 3);
figure
subplot(1,2,1)
% yyaxis left
plot(0:4, (90-theta)/max(90-theta), 'ro-');
hold on
xx = linspace(-1, 5, 100);
yy = (theta(1)-theta(2))/max(90-theta)*xx + (90 - theta(1))/max(90 - theta(1));
plot(xx, yy, 'r--')
% ylim([0 90])

% yyaxis right
plot(0:4, precision/max(precision), 'bo-')
xx = linspace(-1, 5, 100);
yy = (precision(2)-precision(1))/max(precision)*xx + precision(1)/max(precision);
plot(xx, yy, 'b--')
% ylim([0.5 1])
xlim([-0.5 4.5])
ylim([0.7 1.1])

load('/Users/xyao/matlab/code-private/DS_new/dsmodel_cone.mat')
precision = mean(mean(cell2mat(Precision_all), 2), 3) - 0.5;
theta_temp = cellfun(@cell2mat, theta_e_all, 'UniformOutput', 0);
theta = mean(mean(cell2mat(theta_temp), 2), 3);
subplot(1,2,2)
% yyaxis left
plot(0:4, (90-theta)/max(90-theta), 'ro-');
hold on
xx = linspace(-1, 5, 100);
yy = (theta(1)-theta(2))/max(90-theta)*xx + (90 - theta(1))/max(90 - theta(1));
plot(xx, yy, 'r--')
% ylim([0 90])

% yyaxis right
plot(0:4, precision/max(precision), 'bo-')
xx = linspace(-1, 5, 100);
yy = (precision(2)-precision(1))/max(precision)*xx + precision(1)/max(precision);
plot(xx, yy, 'b--')
% ylim([0.5 1])
xlim([-0.5 4.5])
ylim([0.82 1.02])

