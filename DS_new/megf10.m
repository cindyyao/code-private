load dsStat.mat
tuning_curves{1} = dsStatsWT.allDS_tuningCurves/8;
tuning_curves{2} = dsStatsKO.allDS_tuningCurves/8;
xdata = dsStatsWT.dirs/180*pi;
for i = 1:2
    figure
    for cc = 1:size(tuning_curves{i}, 1)
        ydata = tuning_curves{i}(cc, :);
        [f, g] = fit_cos(xdata, ydata);
        Ymax{i}(cc) = f.ymax;
        Phi{i}(cc) = f.phi;
        Alpha{i}(cc) = f.alpha;
        Width{i}(cc) = acos(2*0.5^(1/f.alpha) - 1)/pi*360;
        B{i}(cc) = f.b;

        xfit = linspace(0, 2*pi, 100);
        yfit = f.ymax.*(0.5+0.5*cos(xfit+f.phi)).^f.alpha+f.b;
        subplot(8, 9, cc)
        plot(xdata, ydata, 'b')
        hold on
        plot(xfit, yfit, 'r')
%                 ylim([0 1.1])
        width = acos(2 * (0.5.^(1/f.alpha) - 0.5))/pi*360;
        title(['width = ', num2str(width)])
    end
end

%%
alpha2 = [1 1 1 1]*mean(Alpha{1});
phi = [0 1/2 1 3/2]*pi; % preferred direction
r_base = [1 1 1 1]*mean(B{1}); % background firing
repeat_n = 200;
% ymax = [1 1 1 1]*mean(Ymax{1});
ymax = [1 1 1 1]*1000;

[error_wt_mean, error_wt_std] = run_poiss_MLdecoder_megf10(phi,alpha2,r_base,ymax,repeat_n);
% save('dsmodel_megf10.mat', 'error_wt')

alpha2 = [1 1 1 1]*mean(Alpha{2});
phi = [0 1/2 1 3/2]*pi; % preferred direction
r_base = [1 1 1 1]*mean(B{2}); % background firing
repeat_n = 2000;
ymax = [1 1 1 1]*mean(Ymax{2});
[error_ko_mean, error_ko_std] = run_poiss_MLdecoder_megf10(phi,alpha2,r_base,ymax,repeat_n);
% save('dsmodel_megf10.mat', 'error_ko', '-append')

xx = linspace(0, 360, 100);
figure
errorbar(xx, error_wt_mean, error_wt_std, 'r')
figure
errorbar(xx, error_ko_mean, error_ko_std, 'r')
% legend('WT', 'KO')
xlabel('direction (degree)')
ylabel('mean error')
ylim([7 13])
xlim([-10 370])