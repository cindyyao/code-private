repeat = 10;
rmax = 4.9;

%% compare different number of broad tuning

contrast = [40]; % contrast
a = 1.8; b = 0.8;
Alpha2 = [a a a a; a a b a; a b b a; b b b a; b b b b]; % tuning width 


c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
r_base = [1 1 1 1]*0.02; % background firing
repeat_n = 1000;
a = 0.95; b = 0.6;
Beta = [a a a a; a a b a; a b b a; b b b a; b b b b];
a = 1; b = 1.4;
Ymax = [a a a a; a a b a; a b b a; b b b a; b b b b]*rmax;
[Precision_rod, theta_e_rod] = deal(cell(1, repeat));
for rp = 1:repeat
    clear Precision_ratio
    for k = 1:size(Beta, 1)
        ymax = Ymax(k, :); % ymax = 5 for rod light level
        alpha2 = Alpha2(k, :);
        for ctr = 1:length(contrast)
            [~, Precision(k, ctr, :)] = run_poiss_Detection_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,Beta(k,:));
            [error(k,ctr), theta_e{k,ctr}] = run_poiss_MLdecoder_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,Beta(k,:));
    %         [width ctr]
            disp(['repeat ' num2str(rp) ' detection task: beta ' num2str(k) ' contrast ' num2str(ctr) ' finished.'])
        end
    end
    Precision_rod{rp} = Precision;
    ERROR_rod{repeat} = error;
    theta_e_rod{rp} = theta_e;

    save('dsmodel_revision2.mat', 'Precision_rod', 'ERROR_rod', 'theta_e_rod')
end

%%
repeat = 10;
rmax = 25;

%% compare different number of broad tuning

contrast = [40]; % contrast
a = 1.8; b = 0.8;
Alpha2 = [a a a a; a a b a; a b b a; b b b a; b b b b]; % tuning width 


c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
r_base = [1 1 1 1]*0.02; % background firing
repeat_n = 1000;
a = 0.95; b = 0.6;
Beta = [a a a a; a a b a; a b b a; b b b a; b b b b];
a = 1; b = 1.4;
Ymax = [a a a a; a a b a; a b b a; b b b a; b b b b]*rmax;
[Precision_cone, theta_e_cone] = deal(cell(1, repeat));
for rp = 1:repeat
    clear Precision_ratio
    for k = 1:size(Beta, 1)
        ymax = Ymax(k, :); % ymax = 5 for rod light level
        alpha2 = Alpha2(k, :);
        for ctr = 1:length(contrast)
            [~, Precision(k, ctr, :)] = run_poiss_Detection_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,Beta(k,:));
            [error(k,ctr), theta_e{k,ctr}] = run_poiss_MLdecoder_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,Beta(k,:));
    %         [width ctr]
            disp(['repeat ' num2str(rp) ' detection task: beta ' num2str(k) ' contrast ' num2str(ctr) ' finished.'])
        end
    end
    Precision_cone{rp} = Precision;
    ERROR_cone{repeat} = error;
    theta_e_cone{rp} = theta_e;

    save('dsmodel_revision2.mat', 'Precision_cone', 'ERROR_cone', 'theta_e_cone', '-append')
end

%%

temp(1, 1, :) = ERROR_rod;
error_rod_all = mean(cell2mat(temp), 3);

temp(1, 1, :) = Precision_rod;
Precision_rod_all = mean(cell2mat(temp), 3);

error_all = error_rod_all;
precision_all = Precision_rod_all - 0.5;

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

clear temp
temp(1, 1, :) = ERROR_cone;
error_cone_all = mean(cell2mat(temp), 3);

temp(1, 1, :) = Precision_cone;
Precision_cone_all = mean(cell2mat(temp), 3);

error_all = error_cone_all;
precision_all = Precision_cone_all - 0.5;

c1 = 1;
subplot(2, 1, 2)
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
