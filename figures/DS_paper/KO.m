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

% load('DS150603.mat')
% rho_mean_wt = rho_dg_mean;
% rho_ste_wt = rho_dg_ste;
% dsi_wt = dsi_dg_wt;
% rho_dg_wt = rho_dg;
% load('DS161208.mat')
% dsi_ko_1208 = dsi_dg;
% rho_dg_ko_1208 = rho_dg;
% gfp_rho_1208 = gfp_rho;
% gfp_n_rho_1208 = gfp_n_rho;
% gfp_dsi_1208 = gfp_dsi{1};
% gfp_n_dsi_1208 = gfp_n_dsi{1};
% dsiMin_1208 = disMin_oo{1};
load('DS150603.mat')
rho_mean_wt_0603 = rho_dg_mean_bgnd;
rho_ste_wt_0603 = rho_dg_ste_bgnd;
dsi_wt_0603 = dsi_dg_wt_bgnd;
rho_dg_wt_0603 = rho_dg_bgnd;

load('DS160130.mat')
rho_mean_wt_0130 = rho_dg_mean_bgnd;
rho_ste_wt_0130 = rho_dg_ste_bgnd;
dsi_wt_0130 = dsi_dg_wt_bgnd;
rho_dg_wt_0130 = rho_dg_bgnd;

load('DS161208.mat')
dsi_ko_1208 = dsi_dg_bgnd;
rho_dg_ko_1208 = rho_dg_bgnd;
gfp_rho_1208 = gfp_rho_bgnd;
gfp_n_rho_1208 = gfp_n_rho_bgnd;
gfp_dsi_1208 = gfp_dsi_bgnd{1};
gfp_n_dsi_1208 = gfp_n_dsi_bgnd{1};
dsiMin_1208 = disMin_oo{1};

% dsi_nocorr_1208 = dsi_dg{1}{1}(idx_ko);
% dsi_corr_1208 = dsi_dg{1}{1}(~idx_ko);
id_dir{1}(9) = [];
[~, i] = intersect(id_dir{1}, unique(corr_cells_test));
idx_ko_gfp(i) = 0;
dsi_nocorr_1208 = dsi_dg{1}{1}(idx_ko_gfp);
dsi_corr_1208 = dsi_dg{1}{1}(~idx_ko_gfp);


% load('DS161212.mat')
% dsi_ko_1212 = dsi_dg;
% rho_dg_ko_1212 = rho_dg;
% gfp_rho_1212 = gfp_rho;
% gfp_n_rho_1212 = gfp_n_rho;
% gfp_dsi_1212 = gfp_dsi{1};
% gfp_n_dsi_1212 = gfp_n_dsi{1};
% dsiMin_1212 = disMin_oo{1};

load('DS161212.mat')
dsi_ko_1212 = dsi_dg_bgnd;
rho_dg_ko_1212 = rho_dg_bgnd;
gfp_rho_1212 = gfp_rho_bgnd;
gfp_n_rho_1212 = gfp_n_rho_bgnd;
gfp_dsi_1212 = gfp_dsi_bgnd{1};
gfp_n_dsi_1212 = gfp_n_dsi_bgnd{1};
dsiMin_1212 = disMin_oo{1};


% dsi_nocorr_1212 = dsi_dg{1}{1}(idx_ko);
% dsi_corr_1212 = dsi_dg{1}{1}(~idx_ko);
[~, i] = intersect(id_dir{1}, unique(corr_cells_test));
idx_ko_gfp(i) = 0;

dsi_nocorr_1212 = dsi_dg{1}{1}(idx_ko_gfp);
dsi_corr_1212 = dsi_dg{1}{1}(~idx_ko_gfp);

gfp_rho = [gfp_rho_1208; gfp_rho_1212];
gfp_n_rho = [gfp_n_rho_1208; gfp_n_rho_1212];
gfp_dsi = [gfp_dsi_1208; gfp_dsi_1212];
gfp_n_dsi = [gfp_n_dsi_1208; gfp_n_dsi_1212];
dsiMin = [dsiMin_1208 dsiMin_1212];
dsi_nocorr = [dsi_nocorr_1208; dsi_nocorr_1212];
dsi_corr = [dsi_corr_1208; dsi_corr_1212];
% dsi_ko = [dsi_ko_1208{1}{1}; dsi_ko_1212{1}{1}];
% min distance
XX = 0:10:150;
figure
a = hist(dsiMin, XX);
bar(XX, a, 1)
xlabel('um')
ylabel('number of cells')
xlim([0 150])
% ylim([0 8])

gfp = sum(dsiMin < 50);
gfp_n = sum(dsiMin >= 50);
figure
bar([1 2], [gfp gfp_n], 1)
xlim([0 3])

for ll = 1:5
    for ct = 1:4
        rho_dg_ko{ll}{ct} = [rho_dg_ko_1208{ll}{ct}; rho_dg_ko_1212{ll}{ct}];
        rho_dg_wt{ll}{ct} = [rho_dg_wt_0603{ll}{ct}; rho_dg_wt_0130{ll}{ct}];
        rho_mean_ko{ll}(ct, :) = mean(rho_dg_ko{ll}{ct});
        rho_mean_wt{ll}(ct, :) = mean(rho_dg_wt{ll}{ct});
        rho_ste_ko{ll}(ct, :) = std(rho_dg_ko{ll}{ct})/sqrt(size(rho_dg_ko{ll}{ct}, 1));
        rho_ste_wt{ll}(ct, :) = std(rho_dg_wt{ll}{ct})/sqrt(size(rho_dg_wt{ll}{ct}, 1));
        dsi_ko{ll}{ct} = [dsi_ko_1208{ll}{ct}; dsi_ko_1212{ll}{ct}];
        dsi_wt{ll}{ct} = [dsi_wt_0603{ll}{ct}; dsi_wt_0130{ll}{ct}];
        dsi_mean_ko(ll, ct) = mean(dsi_ko{ll}{ct});
        dsi_mean_wt(ll, ct) = mean(dsi_wt{ll}{ct});
        dsi_ste_ko(ll, ct) = std(dsi_ko{ll}{ct})/sqrt(length(dsi_ko{ll}{ct}));
        dsi_ste_wt(ll, ct) = std(dsi_wt{ll}{ct})/sqrt(length(dsi_wt{ll}{ct}));
    end
end

dsi_others_wt = cell2mat(dsi_wt{1}(2:4)');
dsi_others_ko = cell2mat(dsi_ko{1}(2:4)');


% clean the dsi by excluding outliers
thresh_low = mean(dsi_others_wt) - std(dsi_others_wt) * 2;
thresh_high = mean(dsi_others_wt) + std(dsi_others_wt) * 2;
dsi_others_wt_clean = dsi_others_wt((dsi_others_wt < thresh_high) & (dsi_others_wt > thresh_low));

thresh_low = mean(dsi_others_ko) - std(dsi_others_ko) * 2;
thresh_high = mean(dsi_others_ko) + std(dsi_others_ko) * 2;
dsi_others_ko_clean = dsi_others_ko((dsi_others_ko < thresh_high) & (dsi_others_ko > thresh_low));

thresh_low = mean(gfp_dsi) - std(gfp_dsi) * 2;
thresh_high = mean(gfp_dsi) + std(gfp_dsi) * 2;
gfp_dsi_clean = gfp_dsi((gfp_dsi < thresh_high) & (gfp_dsi > thresh_low));

thresh_low = mean(gfp_n_dsi) - std(gfp_n_dsi) * 2;
thresh_high = mean(gfp_n_dsi) + std(gfp_n_dsi) * 2;
gfp_n_dsi_clean = gfp_n_dsi((gfp_n_dsi < thresh_high) & (gfp_n_dsi > thresh_low));

thresh_low = mean(dsi_nocorr) - std(dsi_nocorr) * 2;
thresh_high = mean(dsi_nocorr) + std(dsi_nocorr) * 2;
dsi_nocorr_clean = dsi_nocorr((dsi_nocorr < thresh_high) & (dsi_nocorr > thresh_low));

thresh_low = mean(dsi_corr) - std(dsi_corr) * 2;
thresh_high = mean(dsi_corr) + std(dsi_corr) * 2;
dsi_corr_clean = dsi_corr((dsi_corr < thresh_high) & (dsi_corr > thresh_low));

xx = [0:pi/4:7*pi/4] - pi;
xx = xx/pi*180+22.5;
ct = 1;
ll = 1;
figure
subplot(1,2,1)
errorbar(xx, rho_mean_wt{ll}(ct, :), rho_ste_wt{ll}(ct, :), 'k');
hold on
errorbar(xx, rho_mean_ko{ll}(ct, :), rho_ste_ko{ll}(ct, :), 'k--');
xlabel('direction (degree)')
ylabel('normalized average response')
legend('WT', 'KO')
ylim([0 1])


ll = 1;
subplot(1,2,2)
errorbar(xx, rho_dg_wt_others_mean{ll}, rho_dg_wt_others_ste{ll}, 'k');
hold on
errorbar(xx, rho_dg_ko_others_mean{ll}, rho_dg_ko_others_ste{ll}, 'k--');
xlabel('direction (degree)')
ylabel('normalized average response')
legend('WT', 'KO')
ylim([0 1])


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
ylim([0 1])



for ll = 1:5
    rho_dg_wt_others{ll} = cell2mat(rho_dg_wt{ll}(2:4)');
    rho_dg_wt_others_mean{ll} = mean(rho_dg_wt_others{ll});
    rho_dg_wt_others_ste{ll} = std(rho_dg_wt_others{ll})/sqrt(size(rho_dg_wt_others{ll}, 1));
end

for ll = 1:5
    rho_dg_ko_others{ll} = cell2mat(rho_dg_ko{ll}(2:4)');
    rho_dg_ko_others_mean{ll} = mean(rho_dg_ko_others{ll});
    rho_dg_ko_others_ste{ll} = std(rho_dg_ko_others{ll})/sqrt(size(rho_dg_ko_others{ll}, 1));
end



figure
errorbar(xx, mean(gfp_rho), std(gfp_rho)/sqrt(size(gfp_rho, 1)), 'b')
hold on
errorbar(xx, mean(gfp_n_rho), std(gfp_n_rho)/sqrt(size(gfp_n_rho, 1)), 'r')
ylim([0 1])

[~, p(1)] = ttest2(gfp_dsi_clean, gfp_n_dsi_clean);

[~, p(2)] = ttest2(gfp_dsi_clean, dsi_others_ko);

[~, p(3)] = ttest2(gfp_n_dsi_clean, dsi_others_ko);

[~, p(4)] = ttest2(dsi_nocorr_clean, dsi_others_ko);

[~, p(5)] = ttest2(dsi_nocorr_clean, dsi_corr_clean);

% bar plot
dsi_mean(1) = mean(gfp_dsi_clean);
dsi_mean(2) = mean(gfp_n_dsi_clean);
dsi_mean(3) = mean(dsi_others_ko);

dsi_ste(1) = std(gfp_dsi_clean)/sqrt(length(gfp_dsi_clean));
dsi_ste(2) = std(gfp_n_dsi_clean)/sqrt(length(gfp_n_dsi_clean));
dsi_ste(3) = std(dsi_others_ko)/sqrt(length(dsi_others_ko));


Ticks = {'GFP+', 'GFP-', 'KO others'};

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);

xtick = Ticks;
model_series = dsi_mean';
model_error = dsi_ste';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('DSI')
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end
ylim([0 1])
%%
for ll = 1:5
    gfp_dsi_mean(ll) = mean(gfp_dsi{ll});
    gfp_dsi_ste(ll) = std(gfp_dsi{ll})/sqrt(size(gfp_dsi{ll}, 1));
    gfp_n_dsi_mean(ll) = mean(gfp_n_dsi{ll});
    gfp_n_dsi_ste(ll) = std(gfp_n_dsi{ll})/sqrt(size(gfp_n_dsi{ll}, 1));
end
figure
errorbar([0:4], gfp_dsi_mean, gfp_dsi_ste, 'b')
hold on
errorbar([0:4], gfp_n_dsi_mean, gfp_n_dsi_ste, 'r')
ylim([0 1])

%%

dsi_ko = [dsi_ko_1208{1}{1}; dsi_ko_1212{1}{1}];

figure
subplot(3,1,1)
bar(0.025:0.05:1, histcounts(dsi_dg_wt{1}{1}, 0:0.05:1), 1)
ylabel('# of cells')
legend('WT superior')
% title('WT superior')
subplot(3,1,2)
bar(0.025:0.05:1, histcounts(dsi_ko, 0:0.05:1), 1)
ylabel('# of cells')
% title('KO superior')
subplot(3,1,3)
bar(0.025:0.05:1, histcounts(cell2mat(dsi_dg_wt{1}(2:4)'), 0:0.05:1), 1)
ylabel('# of cells')
xlabel('DSI')
% title('WT others')


figure
subplot(3,1,1)
bar(0.05:0.1:1, histcounts(dsi_dg_wt{1}{1}, 0:0.1:1), 1)
ylabel('# of cells')
title('WT superior')
subplot(3,1,2)
bar(0.05:0.1:1, histcounts(dsi_ko, 0:0.1:1), 1)
ylabel('# of cells')
title('KO superior')
subplot(3,1,3)
bar(0.05:0.1:1, histcounts(cell2mat(dsi_dg_wt{1}(2:4)'), 0:0.1:1), 1)
ylabel('# of cells')
xlabel('DSI')
title('WT others')


figure
subplot(3,1,1)
bar(0.0125:0.025:1, histcounts(dsi_dg_wt{1}{1}, 0:0.025:1), 1)
ylabel('# of cells')
title('WT superior')
subplot(3,1,2)
bar(0.0125:0.025:1, histcounts(dsi_ko, 0:0.025:1), 1)
ylabel('# of cells')
title('KO superior')
subplot(3,1,3)
bar(0.0125:0.025:1, histcounts(cell2mat(dsi_dg_wt{1}(2:4)'), 0:0.025:1), 1)
ylabel('# of cells')
xlabel('DSI')
title('WT others')