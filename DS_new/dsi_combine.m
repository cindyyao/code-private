load('DS150603.mat', 'dsi_dg')
dsi_150603 = dsi_dg;
load('DS160130.mat', 'dsi_dg')
dsi_160130 = dsi_dg;
ct = {'superior', 'anterior', 'inferior', 'posterior'};

for ll = 1:5
    for ct = 1:4
        dsi_all{ll}{ct} = [dsi_150603{ll}{ct}; dsi_160130{ll}{ct}];
        dsi_mean{ll}(ct) = mean(dsi_all{ll}{ct});
        dsi_ste{ll}(ct) = std(dsi_all{ll}{ct})/sqrt(length(dsi_all{ll}{ct}));
    end
end
dsi_mean = cell2mat(dsi_mean');
dsi_ste = cell2mat(dsi_ste');

% DSI
% FigHandle = figure;
% set(FigHandle, 'Position', [100, 100, 1000, 500]);
% 
xtick = ct;
model_series = dsi_mean';
model_error = dsi_ste';
% h = bar(model_series);
% set(h,'BarWidth',1); % The bars will now touch each other
% 
% set(gca,'XTicklabel',xtick)
% ylabel('DSI')
% legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
% hold on;
% 
% numgroups = size(model_series, 1); 
% numbars = size(model_series, 2); 
% 
% groupwidth = min(0.8, numbars/(numbars+1.5));
% 
% for i = 1:numbars
% % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
% x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
% errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
% end

% dsi curve
marker = 'xosd';
figure
subplot(2,1,1)
for i = 1:4
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
errorbar(0:4, model_series(i,:), model_error(i,:), 'Color', 'k', 'Marker', marker(i), 'MarkerSize', 10);
hold on
end
xlim([-0.5 4.5])
ylim([0 1.1])
legend('superior', 'anterior', 'inferior', 'posterior')


dsi_all = dsi_150603;
for ct = 1:4
    dsi_ll{ct} = [];
    group{ct} = [];
    for ll = 1:5
        dsi_ll{ct} = [dsi_ll{ct} dsi_all{ll}{ct}'];
        group{ct} = [group{ct} ones(1,length(dsi_all{ll}{ct}))*ll];
    end
    p(ct) = anova1(dsi_ll{ct}, group{ct});
end
%% width combine
load('DS150603.mat', 'Width')
Width_150603 = Width;
load('DS160130.mat', 'Width')
Width_160130 = Width;
for ll = 1:5
    for ct = 1:4
        Width_all{ct}{ll} = [Width_150603{ct}{ll} Width_160130{ct}{ll}];
        WidthMean(ct, ll) = mean(Width_all{ct}{ll});
        WidthSte(ct, ll) = std(Width_all{ct}{ll})/sqrt(length(Width_all{ct}{ll}));
    end
end

marker = 'xosd';
subplot(2,1,2)
for dir = 1:4
    errorbar(0:4, WidthMean(dir, :), WidthSte(dir, :), 'Color', 'k', 'Marker', marker(dir), 'MarkerSize', 10)
    hold on
end
legend('superior', 'anterior', 'inferior', 'posterior')
ylim([0 250])
xlabel('light level')
ylabel('tuning width (degree)')
xlim([-0.5 4.5])

%% 
for ct = 1:4
    width_ll{ct} = [];
    group{ct} = [];
    for ll = 1:5
        width_ll{ct} = [width_ll{ct} Width_all{ct}{ll}];
        group{ct} = [group{ct} ones(1,length(Width_all{ct}{ll}))*ll];
    end
    p(ct) = anova1(width_ll{ct}, group{ct});
end


