%% model parameters (vector of 4)
% NDF 4 control
% CRF 
% ymax = [1 1 1.5 1]*20;
% c50 = [1.8 1.8 1.8 1.8];
% alpha1 = [7 7 7 7];
% sigma1 = [0.1 0.1 0.1 0.1]*3; % std of incident noise, should be same for all cells
% % tuning curve
% phi = [0 1/2 1 3/2]*pi; % preferred direction
% alpha2 = [3 3 1.5 3]; % tuning width 
% r_base = [0 0 10 0]*1; % background firing
% % independent noise
% sigma2 = [1 1 1 1]*1; % should be same for all cells

c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
sigma1 = [0.1 0.1 0.1 0.1]*3; % std of incident noise, should be same for all cells
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
alpha2 = [0.5 1 1 1]*2; % tuning width 
r_base = [10 1 1 1]*1; % background firing
% independent noise
sigma2 = [1 1 1 1]*2; % should be same for all cells
ymax = [2 1 1 1]*10;


error_control = run_MLdecoder(c50, alpha1, sigma1, phi, alpha2, r_base, sigma2, ymax);
%% NDF 4 same

% ymax = [1 1 1 1]*20;
% c50 = [1.8 1.8 1.8 1.8];
% alpha1 = [7 7 7 7];
% sigma1 = [0.1 0.1 0.1 0.1]*3; % std of incident noise, should be same for all cells
% % tuning curve
% phi = [0 1/2 1 3/2]*pi; % preferred direction
% alpha2 = [3 3 3 3]; % tuning width 
% r_base = [0 0 0 0]*1; % background firing
% % independent noise
% sigma2 = [1 1 1 1]*1; % should be same for all cells


% alpha2 = [1 1 1 1]*3; % tuning width 
% r_base = [1 1 1 1]*1; % background firing
% ymax = [1 1 1 1]*10;

alpha2 = [1 1 1 1]*2; % tuning width 
r_base = [1 1 1 1]*1; % background firing
ymax = [1 1 1 1]*10;


error_replaced = run_MLdecoder(c50, alpha1, sigma1, phi, alpha2, r_base, sigma2, ymax);

%%
figure
errorbar(1, mean(error_control(:)), std(error_control(:))/sqrt(length(error_control(:))), 'bo')
hold on
errorbar(2, mean(error_replaced(:)), std(error_replaced(:))/sqrt(length(error_control(:))), 'ro')
xtick = {'', 'control', '"KO"', ''};
set(gca,'XTicklabel',xtick)
ylabel('mean error(degree)')
%% model parameters (vector of 4)
% NDF 4 control

ERROR = [];

for cc = linspace(log10(1), log10(1000), 20)

    c = 10^cc;

    c50 = [1.6 1.8 1.8 1.8];
    alpha1 = [7 7 7 7];
    % tuning curve
    phi = [0 1/2 1 3/2]*pi; % preferred direction
    alpha2 = [0.8 1 1 1]*3; % tuning width 
    r_base = [10 1 1 1]*1; % background firing
    % independent noise
    ymax = [2 1 1 1]*10;


    error_control = run_poiss_MLdecoder(c50, alpha1, phi, alpha2, r_base, ymax, c);

    %
    c50 = [1.8 1.8 1.8 1.8];
    alpha2 = [1 1 1 1]*3; % tuning width 
    r_base = [1 1 1 1]*1; % background firing
    ymax = [1 1 1 1]*10;


    error_replaced = run_poiss_MLdecoder(c50, alpha1, phi, alpha2, r_base, ymax, c);

    %
    c50 = [1.6 1.6 1.6 1.6];
    alpha2 = [0.8 0.8 0.8 0.8]*3; % tuning width 
    r_base = [10 10 10 10]*1; % background firing
    ymax = [2 2 2 2]*10;


    error_replaced2 = run_poiss_MLdecoder(c50, alpha1, phi, alpha2, r_base, ymax, c);


    ERROR = [ERROR; error_control error_replaced error_replaced2];
    close all
end

ctr =  linspace(log10(1), log10(1000), 20);
ctr_l = 10.^ctr;
figure
plot(repmat(ctr_l', 1, 3), ERROR)
set(gca, 'xscale', 'log')
legend('control', '"ko"', 'all superior')
xlabel('contrast %')
ylabel('mean error (deg)')
title('Poisson model')

