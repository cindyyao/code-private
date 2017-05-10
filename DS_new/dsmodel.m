%% model parameters (vector of 4)

clear errorSup errorNarrow errorBroad
% NDF 4
% CRF 
c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% sigma1 = [0.1 0.1 0.1 0.1]*1; % std of incident noise, should be same for all cells
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
Alpha2 = [0.8 0.8 0.8 0.8; 0.8 0.8 0.8 2; 0.8 0.8 2 2; 0.8 2 2 2; 2 2 2 2]; % tuning width 
r_base = [1 1 1 1]*0.2; % background firing
ymax = [1 1 1 1]*5;


% ctr = 5:5:100;
ctr = [5 10 20 40 80 150 300]; % contrast
for i = 1:5
    alpha2 = Alpha2(i, :);
    for ii = 1:length(ctr)
        % run_MLdecoder
        [ERROR(i, ii), THETA_E{i, ii}] = run_poiss_MLdecoder(c50,alpha1,phi,alpha2,r_base,ymax,ctr(ii));
        ii
    end
end

figure
for i = 1:5
    plot(ctr, ERROR(i, :))
    hold on
end
set(gca, 'xscale', 'log')
xlim([4 400])
xlabel('contrast')
ylabel('error')
title('low light level')

%%
% ctr = 5:5:100;
ctr = [5 10 20 40 80 150 300]; % contrast
for ii = 1:length(ctr)
% NDF 4
% CRF 
c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% sigma1 = [0.1 0.1 0.1 0.1]*1; % std of incident noise, should be same for all cells
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
alpha2 = [0.8 2 2 2]; % tuning width 
r_base = [0.2 0.2 0.2 0.2]; % background firing
% independent noise
% sigma2 = [1 1 1 1]*1; % should be same for all cells
ymax = [1 1 1 1]*5;


% run_MLdecoder
errorSup(ii) = run_poiss_MLdecoder(c50,alpha1,phi,alpha2,r_base,ymax,ctr(ii));



%
% NDF 4
% CRF 

% ymax = [1 1 1 1]*5;
% c50 = [1.8 1.8 1.8 1.8];
% alpha1 = [7 7 7 7];
% sigma1 = [0.1 0.1 0.1 0.1]*1; % std of incident noise, should be same for all cells
% phi = [0 1/2 1 3/2]*pi; % preferred direction
alpha2 = [1 1 1 1]*2; % tuning width 
% r_base = [1 1 1 1]*0.2; % background firing
% independent noise
% sigma2 = [1 1 1 1]*1; % should be same for all cells

% run MLdecoder
% mean(theta_e_mean)

errorNarrow(ii) = run_poiss_MLdecoder(c50,alpha1,phi,alpha2,r_base,ymax,ctr(ii));

% ymax = [1 1 1 1]*8;
% c50 = [1 1 1 1]*1.5;
% alpha1 = [7 7 7 7];
% sigma1 = [0.1 0.1 0.1 0.1]*1; % std of incident noise, should be same for all cells
% phi = [0 1/2 1 3/2]*pi; % preferred direction
alpha2 = [1 1 1 1]*0.8; % tuning width 
% r_base = [1 1 1 1]*3; % background firing
% independent noise
% sigma2 = [1 1 1 1]*1; % should be same for all cells

% run MLdecoder
% mean(theta_e_mean)

errorBroad(ii) = run_poiss_MLdecoder(c50,alpha1,phi,alpha2,r_base,ymax,ctr(ii));

end

figure
plot(ctr, errorSup, 'b')
hold on
plot(ctr, errorNarrow, 'r')
plot(ctr, errorBroad, 'k')
legend('one broad', 'all narrow', 'all broad');
set(gca, 'xscale', 'log')
xlim([4 400])
xlabel('contrast')
ylabel('error')
title('high light level')

%%
% NDF 0
% CRF 
ymax = [35 25 25 25];
c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
sigma1 = [0.1 0.1 0.1 0.1]; % std of incident noise, should be same for all cells
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
alpha2 = [1 1 1 1]*2; % tuning width 
sigma3 = [1 1 1 1]*1;
r_base = [1 0 0 0]; % background firing
% independent noise
sigma2 = [1 1 1 1]; % should be same for all cells


%% plot model CRF and tuning curve
% CRF
x = linspace(0.7,3, 100);
for cc = 1:4
    y(cc, :) = ymax(cc).*(x.^alpha1(cc)./(x.^alpha1(cc) + c50(cc).^alpha1(cc)));
end
figure
plot(x', y')

% tuning curves
clear x y
x = linspace(0, 2*pi, 100);
for cc = 1:4
    y(cc, :) = ymax(cc).*(0.5 + 0.5 * sin(x + phi(cc))).^alpha2(cc);
%     y(cc, :) = ymax(cc).*exp(-(mod(x-phi(cc)+pi,2*pi)-pi).^2./(2*sigma3(cc).^2));
end
figure
plot(x', y')

% all
ctr = log10(80);
for cc = 1:4
    y(cc, :) = ymax(cc).*(ctr.^alpha1(cc)./(ctr.^alpha1(cc) + c50(cc).^alpha1(cc))).*(0.5 + 0.5 * sin(x + phi(cc))).^alpha2(cc)+r_base(cc);
end
figure
plot(x', y')

%%
repeat_n = 1000;
x = log10(80);
THETA = [0:1/18:35/18]*pi; % direction

Output = zeros(length(THETA), repeat_n, 4);
for dir = 1:length(THETA)
    for repeat = 1:repeat_n
        theta = THETA(dir);
        input = ymax.*(x.^alpha1./(x.^alpha1 + c50.^alpha1) + sigma1.*random('norm', 0, 1, [1, 4]));
        gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2;
        output = max(max(input.*gain+r_base, 0) + sigma2.*random('norm', 0, 1, [1, 4]), 0);
        Output(dir, repeat, :) = output;
    end
end


Output_mean = squeeze(mean(Output, 2));
Output_ste = squeeze(std(Output, [], 2)/sqrt(repeat_n));
Output_std = squeeze(std(Output, [], 2));
Output_var = squeeze(var(Output, [], 2));


figure
for cc = 1:4
    errorbar(THETA, Output_mean(:, cc), Output_std(:, cc))
    hold on
end

figure
scatter(Output_mean(:), Output_var(:))


%% 170314
clear error theta_e
ctr = [5 10 20 40 80 150 300]; % contrast
% y = linspace(0.2, 0.99, 20);
% Alpha2 = 1./(-log2(y));

theta = unique(sort([linspace(10,50,4)/180*pi linspace(50,80,10)/180*pi linspace(80,170,8)/180*pi]));
Alpha2 = log(0.5)./log(0.5+0.5*cos(theta));


for i = 1:length(ctr)
    for j = 1:length(Alpha2)
        % NDF 4
        % CRF 
        c50 = [1.8 1.8 1.8 1.8];
        alpha1 = [7 7 7 7];
        % sigma1 = [0.1 0.1 0.1 0.1]*1; % std of incident noise, should be same for all cells
        % tuning curve
        phi = [0 1/2 1 3/2]*pi; % preferred direction
        alpha2 = [1 1 1 1]*Alpha2(j); % tuning width 
        r_base = [1 1 1 1]; % background firing
        % independent noise
        % sigma2 = [1 1 1 1]*1; % should be same for all cells
        ymax = [1 1 1 1]*25; % ymax = 5 for rod light level


        % run_MLdecoder
        [error(i,j), theta_e_temp{i,j}] = run_poiss_MLdecoder(c50,alpha1,phi,alpha2,r_base,ymax,ctr(i));
        [i,j]
    end
end 

ERROR = error;
figure
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

y = linspace(0.01, 0.99, 20);
alpha2 = 1./(-log2(y));
theta = linspace(10,170,20)/180*pi;
alpha2 = log(0.5)./log(0.5+0.5*cos(theta));
theta = 0:0.1:2*pi;
phi = pi;
clear tuning
for i = 1:length(alpha2)
    for j = 1:length(theta)
        tuning(j, i) = (0.5 + 0.5 * cos(theta(j) + phi)).^alpha2(i);
    end
end
figure
plot(repmat(theta', 1, length(alpha2))*180/pi, tuning)

% plot(1./(2.^(1./alpha)))
% 
% 
for i = 1:size(theta_e, 1)
    for j = 1:size(theta_e, 2)
        theta_e_small{i, j} = theta_e{i, j}(1:200, :);
        error_small(i, j) = mean(theta_e_small{i, j}(:));
    end
end


%% compare 3 model
% NDF 4
% CRF 
c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% sigma1 = [0.1 0.1 0.1 0.1]*1; % std of incident noise, should be same for all cells
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
alpha2 = [0.8 1.6 1.6 1.6]; % tuning width 
r_base = [5 0.5 0.5 0.5]; % background firing
% independent noise
sigma2 = [1 1 1 1]*1; % should be same for all cells
ymax = [20 10 10 10];

%% 

%%ctr = [5 10 20 40 80 150 300]; % contrast
% score = (90 - Error)/180 + 0.5;
score = Error;

detectImprove = (Precision_ratio(2, :) - Precision_ratio(3, :))./Precision_ratio(3, :);
angleImprove = (score(3, :) - score(2, :))./score(3, :);
figure
plot(ctr, detectImprove)
hold on
plot(ctr, angleImprove)
plot(ctr, detectImprove + angleImprove)
plot(ctr, zeros(length(ctr)), 'k--')
set(gca, 'xscale', 'log')
legend('detection', 'angle prediction', 'sum')
xlabel('contrast')
ylabel('% improvement')
xlim([4,400])


clear ERROR
for i = 1:length(ctr)
    for j = 1:length(Alpha2)
        ERROR(i, j) = mean(mean(theta_e{i, j}(200:900, :)));
    end
end

figure
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

load('dsDetection.mat')
contrast = [5 10 20 40 80 150 300];
figure
subplot(1, 2, 1)
for i = 1:3
    plot(contrast, Precision_null(i, :))
    hold on
end
set(gca, 'xscale', 'log')
xlabel('contrast')
ylabel('% correct')
legend('all broad', 'one broad', 'all narrow')
title('null likelihood')
xlim([4 400])

subplot(1, 2, 2)
for i = 1:3
    plot(contrast, Precision_ratio(i, :))
    hold on
end
set(gca, 'xscale', 'log')
xlabel('contrast')
ylabel('% correct')
legend('all broad', 'one broad', 'all narrow')
title('likelihood ratio')
xlim([4 400])

% likelihood ratio performance
Error_score = (90 - Error) / 90;
detection_score = (Precision_ratio - 0.5) * 2;


detectImprove = (detection_score(2, :) - detection_score(3, :))./detection_score(3, :);
angleImprove = (Error_score(2, :) - Error_score(3, :))./Error_score(3, :);
figure
plot(contrast, detectImprove)
hold on
plot(contrast, angleImprove)
plot(contrast, detectImprove + angleImprove)
plot(contrast, zeros(length(contrast)), 'k--')
set(gca, 'xscale', 'log')
legend('detection', 'angle prediction', 'sum')
xlabel('contrast')
ylabel('% improvement')
xlim([4,400])

%% percentage improvement
error_total_area = sum(Error(1, :) - Error(5, :));
detection_total_area = sum(Precision_ratio(1, :) - Precision_ratio(5, :));

for i = 1:4
    error_improv_rate(5-i) = sum(Error(i, :) - Error(i + 1, :))/error_total_area;
    detection_improv_rate(5-i) = sum(Precision_ratio(i, :) - Precision_ratio(i + 1, :))/detection_total_area;
end

