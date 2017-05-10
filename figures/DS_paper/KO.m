%% RF size
load('DS161122.mat')
rf_wt = rf_area_clean{5}{1};
load('DS161208.mat')
rf_ko = rf_area_clean{1};

rf = [rf_wt' rf_ko'];

ONOFF = {'ON', 'OFF'};

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = ONOFF;
model_series = cellfun(@mean, rf);
model_error = cellfun(@std, rf)./sqrt(cellfun(@length, rf));
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('RF size (mm^2)')
legend('WT', 'KO');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

%% dim flash response threshold
load('DS160324.mat')
thresh_wt = xthreshold_0304and0324{1};

load('DS161208.mat')
thresh_ko = xthreshold{1};

thresh = [{thresh_wt}; {thresh_ko}];
model_series = 10.^cellfun(@mean, thresh);
model_uerror = 10.^(cellfun(@mean, thresh) + cellfun(@std, thresh)./sqrt(cellfun(@length, thresh))) - 10.^cellfun(@mean, thresh);
model_lerror = 10.^cellfun(@mean, thresh) - 10.^(cellfun(@mean, thresh) - cellfun(@std, thresh)./sqrt(cellfun(@length, thresh)));

h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other
xtick = {'WT', 'KO'};

set(gca,'yscale','log','XTicklabel',xtick)
ylabel('RF size (mm^2)')
% legend('WT', 'KO');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_lerror(:,i), model_uerror(:,i), 'k', 'linestyle', 'none');
end

%% direction tuning

load('DS150603.mat')
rho_mean_wt = rho_dg_mean;
rho_ste_wt = rho_dg_ste;
dsi_wt = dsi_dg_wt;
load('DS161208.mat')
rho_mean_ko = rho_dg_mean;
rho_ste_ko = rho_dg_ste;
dsi_ko = dsi_dg;

xx = [0:pi/4:7*pi/4] - pi;
xx = xx/pi*180+22.5;
ct = 1;
ll = 1;
figure
errorbar(xx, rho_mean_wt{ct}(ll, :), rho_ste_wt{ct}(ll, :), 'k');
hold on
errorbar(xx, rho_mean_ko{ct}(ll, :), rho_ste_ko{ct}(ll, :), 'k--');
xlabel('direction (degree)')
ylabel('normalized average response')
legend('WT', 'KO')


for d = 1:5
    dsi_superior_wt_mean(d) = mean(dsi_wt{d}{1});
    dsi_superior_wt_ste(d) = std(dsi_wt{d}{1})/sqrt(length(dsi_wt{d}{1}));
    
    dsi_superior_ko_mean(d) = mean(dsi_ko{d}{1});
    dsi_superior_ko_ste(d) = std(dsi_ko{d}{1})/sqrt(length(dsi_ko{d}{1}));
    
    dsi_temp = cell2mat(dsi_wt{d}(2:4)');
    dsi_others_wt_mean(d) = mean(dsi_temp);
    dsi_others_wt_ste(d) = std(dsi_temp)/sqrt(length(dsi_temp));

    dsi_temp = cell2mat(dsi_ko{d}(2:4)');
    dsi_others_ko_mean(d) = mean(dsi_temp);
    dsi_others_ko_ste(d) = std(dsi_temp)/sqrt(length(dsi_temp));
end

figure
x = 0:4;
errorbar(x, dsi_superior_wt_mean, dsi_superior_wt_ste, 'k')
hold on
errorbar(x, dsi_superior_ko_mean, dsi_superior_ko_ste, 'k--')
errorbar(x, dsi_others_wt_mean, dsi_others_wt_ste, 'color', [1 1 1]*0.8)
errorbar(x, dsi_others_ko_mean, dsi_others_ko_ste, 'color', [1 1 1]*0.8, 'LineStyle', '--')
xlabel('log(light intensity)')
ylabel('DSI')
legend('WT superior', 'KO superior', 'WT others', 'KO others')
xlim([-1 5])
