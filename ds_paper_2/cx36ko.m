%%%%%% RF size
%% Cx36 ko
load('DS161212.mat', 'rf_all', 'id_dir', 'fs_idx_new', 'idx_dir')
for ll = 1:3
    rf_onoff{ll} = cell(length(rf_all{ll}), 1);
    for cc = 1:length(rf_all{ll})
        if ~isempty(rf_all{ll}{cc})
            rf_onoff{ll}{cc} = sum(rf_all{ll}{cc}, 3);
        end
    end
end

%% fit RF with Gaussian                
PixelArea = (30*4)^2/10^6;
% fit and compute rf area
for ll = 1:3
    for dir = 1:4
        clear rf_area_temp
        rf_area_temp = [];
        for cc = 1:length(id_dir{dir})
            if ~fs_idx_new(idx_dir{dir}(cc), ll)
                data = rf_onoff{ll}{idx_dir{dir}(cc)};
                params = fit_2d_gaussian(data);
                rf_area_temp = [rf_area_temp params.xd * params.yd * pi * PixelArea];
            end
        end
        rf_area_ko{ll}{dir} = rf_area_temp;
    end
end

%% WT
load('DS161122.mat', 'rf_all', 'id_dir', 'fs_idx_new', 'idx_dir')
for ll = 1:5
    rf_onoff{ll} = cell(length(rf_all{ll}), 1);
    for cc = 1:length(rf_all{ll})
        if ~isempty(rf_all{ll}{cc})
            rf_onoff{ll}{cc} = sum(rf_all{ll}{cc}, 3);
        end
    end
end

%% fit RF with Gaussian                
PixelArea = (30*4)^2/10^6;
% fit and compute rf area
for ll = 1:5
    for dir = 1:4
        clear rf_area_temp
        rf_area_temp = [];
        for cc = 1:length(id_dir{dir})
            if ~fs_idx_new(idx_dir{dir}(cc), ll)
                data = rf_onoff{ll}{idx_dir{dir}(cc)};
                params = fit_2d_gaussian(data);
                rf_area_temp = [rf_area_temp params.xd * params.yd * pi * PixelArea];
            end
        end
        rf_area_wt{ll}{dir} = rf_area_temp;
    end
end

rf_plot{1,1} = rf_area_wt{1}{1};
rf_plot{1,2} = rf_area_ko{1}{1};
rf_plot{2,1} = cell2mat(rf_area_wt{1}(2:4));
rf_plot{2,2} = cell2mat(rf_area_ko{1}(2:4));

% exclude outliers
stdn = 2;
for ds = 1:2
    for dir = 1:2
        notdone = 1;
        rf_area_temp = rf_plot{ds, dir};
        while notdone
            a = length(rf_area_temp);
            rf_area_temp(rf_area_temp > std(rf_area_temp)*stdn + mean(rf_area_temp)) = [];
            b = length(rf_area_temp);
            if a == b
                notdone = 0;
                rf_plot{ds, dir} = rf_area_temp;
            end
        end
    end
end


%%
ct = {'superior', 'others'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = ct;
model_series = cellfun(@mean, rf_plot);
model_error = cellfun(@std, rf_plot)./sqrt(cellfun(@length, rf_plot));
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('RF area (mm^2)')
legend('WT','FACx');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

color = 'brgk';
marker = 'od';
c = 2;
figure
for i = 1:2
    for g = 1:2
        plot((i*c+g)*ones(1, length(rf_plot{i, g})), rf_plot{i, g}, ['k' marker(g)], 'markersize', 7)
        hold on
    end
    errorbar([i*c+1 i*c+2], model_series(i, :), model_error(i, :), 'ks-', 'markersize', 15);
end
xlim([2.5 6.5])
legend('WT', 'FACx')
ylabel('RF area (mm^2)')
%% absolute sensitivity

load('DS160304.mat', 'xthreshold_0304')
xthreshold_wt = xthreshold_0304;
load('DS161208.mat', 'xthreshold')
xthreshold_ko = xthreshold;

clear xthreshold
xthreshold{1,1} = xthreshold_wt{1};
xthreshold{1,2} = xthreshold_ko{1};
xthreshold{2,1} = cell2mat(xthreshold_wt(2:3));
xthreshold{2,2} = cell2mat(xthreshold_ko(2:4));

% exclude outliers
for ct = 1:2
    notdone = 1;
    xthreshold_temp = xthreshold{ct, 1};
    while notdone
        a = length(xthreshold_temp);
        xthreshold_temp(xthreshold_temp > std(xthreshold_temp)*2+mean(xthreshold_temp)) = [];
        b = length(xthreshold_temp);
        if a == b
            notdone = 0;
            xthreshold{ct, 1} = xthreshold_temp;
        end
    end
end

model_series = cellfun(@mean, xthreshold);
model_error = cellfun(@std, xthreshold)./sqrt(cellfun(@length, xthreshold));
color = 'brgk';
marker = 'od';

c = 2;
figure
for i = 1:2
    for g = 1:2
        plot((i*c+g)*ones(1, length(xthreshold{i, g})), xthreshold{i, g}, ['k' marker(g)], 'markersize', 7)
        hold on
    end
    errorbar([i*c+1 i*c+2], model_series(i, :), model_error(i, :), 'ks-', 'markersize', 15);
end
xlim([2.5 6.5])
legend('WT', 'FACx')
ylabel('log(R*/rod)')

%% Fig (RF examples)
% WT superior: 2016-11-22-0 cell_id:7011 (ON+OFF)
% WT anterior: 2016-11-22-0 cell_id:2179  (ON+OFF)
% KO superior: 2016-12-12-0 cell_id:362 (ON+OFF)
% KO other: 2016-12-12-0 cell_id:5449 (ON+OFF)