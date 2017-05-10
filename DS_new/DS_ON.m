clear all
color = 'brgck';
%% direction tuning
% 2015-06-03-0
load('DS150603.mat')
id = {[1607 3902 4998 6737], [6347], []};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);

dirn = 3;
D = 5;
T = 2;
p_direction = DG_cut{D}.angle{T}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

%subtypes
clear rho_dg dsi_dg rho_dg_mean rho_dg_ste dsi_dg_mean dsi_dg_ste
for d = 1:5
    subplot(2, 3, d)
    for i = 1:dirn
        rho_dg{d}{i} = [];
        dsi_dg{d}{i} = [];
        if ~isempty(idx{i})
            for cc = 1:length(idx{i})
                if ~dg_idx(idx{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx{i}(cc), :))>0
                [xsort, seq] = sort(xx(idx{i}(cc), :));
                y_temp = DG_cut{d}.rho{T}(idx{i}(cc), :);
                plot(xsort, y_temp(seq), color(i))
                ylim([0 1])
                hold on
                rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
                dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx{i}(cc))];
                end
            end
        end
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    xlim([-pi pi])
end

clearvars -EXCEPT rho_dg dsi_dg dirn D T color
% 2015-06-18-0
load('DS150618.mat')
id = {[6766], [], [1952]};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);

p_direction = DG_cut{D}.angle{T}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

%subtypes
for d = 1:5
    subplot(2, 3, d)
    for i = 1:dirn
        if ~isempty(idx{i})
            for cc = 1:length(idx{i})
                if ~dg_idx(idx{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx{i}(cc), :))>0
                [xsort, seq] = sort(xx(idx{i}(cc), :));
                y_temp = DG_cut{d}.rho{T}(idx{i}(cc), :);
                plot(xsort, y_temp(seq), color(i))
                ylim([0 1])
                hold on
                rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
                dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx{i}(cc))];
                end
            end
        end
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    xlim([-pi pi])
end

clearvars -EXCEPT rho_dg dsi_dg dirn D T color
% 2015-07-03-0
load('DS150703-1.mat')
id = {[261 4576 4923 5087 5239 6017], [2596 3391 4352 7397], [287 1653 6123]};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);

p_direction = DG_cut{D}.angle{T}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

%subtypes
for d = 1:5
    subplot(2, 3, d)
    for i = 1:dirn
        if ~isempty(idx{i})
            for cc = 1:length(idx{i})
                if ~dg_idx(idx{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx{i}(cc), :))>0
                [xsort, seq] = sort(xx(idx{i}(cc), :));
                y_temp = DG_cut{d}.rho{T}(idx{i}(cc), :);
                plot(xsort, y_temp(seq), color(i))
                ylim([0 1])
                hold on
                rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
                dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx{i}(cc))];
                end
            end
        end
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    xlim([-pi pi])
end

clearvars -EXCEPT rho_dg dsi_dg dirn D T color

% 2016-01-30-0
load('DS160130.mat')
id = {[2582 3123 6917 7306], [4502], [392 4863]};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);

p_direction = DG_cut{D}.angle{T}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

%subtypes
for d = 1:5
    subplot(2, 3, d)
    for i = 1:dirn
        if ~isempty(idx{i})
            for cc = 1:length(idx{i})
                if ~dg_idx(idx{i}(cc), d) && sum(DG_cut{d}.rho{T}(idx{i}(cc), :))>0
                [xsort, seq] = sort(xx(idx{i}(cc), :));
                y_temp = DG_cut{d}.rho{T}(idx{i}(cc), :);
                plot(xsort, y_temp(seq), color(i))
                ylim([0 1])
                hold on
                rho_dg{d}{i} = [rho_dg{d}{i}; y_temp(seq)];
                dsi_dg{d}{i} = [dsi_dg{d}{i}; DG_cut{d}.dsindex{T}(idx{i}(cc))];
                end
            end
        end
    end
    xlabel('direction (rad)')
    ylabel('normalized response')
    xlim([-pi pi])
end

clearvars -EXCEPT rho_dg dsi_dg dirn D T color xsort

% average
ct = {'superior', 'anterior', 'inferior'};
ll = {'NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0'};
for d = 1:5
    for i = 1:dirn
        rho_dg_mean{d}(i, :) = mean(rho_dg{d}{i});
        rho_dg_ste{d}(i, :) = std(rho_dg{d}{i})/sqrt(size(rho_dg{d}{i}, 1));
        dsi_dg_mean{d}(i) = mean(dsi_dg{d}{i});
        dsi_dg_ste{d}(i) = std(dsi_dg{d}{i})/sqrt(length(dsi_dg{d}{i}));
    end
end
dsi_dg_mean = cell2mat(dsi_dg_mean');
dsi_dg_ste = cell2mat(dsi_dg_ste');

% plot average (cell type)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for d = 1:5
    subplot(2, 3, d)
    for i = 1:dirn
        errorbar(xsort, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(i));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
    title(ll{d})
end
legend(ct)

for i = 1:5
    for dir = 1:3
        dsi_mean(dir, i) = mean(dsi_dg{i}{dir});
        dsi_ste(dir, i) = std(dsi_dg{i}{dir})/sqrt(length(dsi_dg{i}{dir}));
    end
end

% plot average (light level)
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:dirn
    subplot(2, 2, i)
    for d = 1:5
        errorbar(xsort/pi*180, rho_dg_mean{d}(i, :), rho_dg_ste{d}(i, :), color(d));
        hold on
    end
    xlabel('direction (rad)')
    ylabel('normalized average response')
    title(ct{i})
    xlim([-200 200])
end
legend(ll)

marker = 'xo*d';
subplot(2, 2, 4)
for dir = 1:3
    errorbar([0:4], dsi_mean(dir, :), dsi_ste(dir, :), ['k-' marker(dir)], 'MarkerSize', 10);
    hold on
end
ylim([0 1.1])
xlim([-0.5 4.5])
legend(ct)

% average across direction
for i = 1:5
    rho_dg_all_mean(i, :) = mean(cell2mat(rho_dg{i}'));
    rho_dg_all_ste(i, :) = std(cell2mat(rho_dg{i}'))/sqrt(size(cell2mat(rho_dg{i}'), 1));
end

load('DS150603.mat')

D = 5;
T = 1;
p_direction = DG_cut{D}.angle{T}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

%subtypes
sub_i = 2;
for d = 1:5
    rho_dg_oo{d} = [];
    dsi_dg_oo{d} = [];
    for cc = 1:length(idx_sub{sub_i})
        if ~dg_idx(idx_sub{sub_i}(cc), d) && sum(DG_cut{d}.rho{T}(idx_sub{sub_i}(cc), :))>0
        [xsort, seq] = sort(xx(idx_sub{sub_i}(cc), :));
        y_temp = DG_cut{d}.rho{T}(idx_sub{sub_i}(cc), :);
        rho_dg_oo{d} = [rho_dg_oo{d}; y_temp(seq)];
        dsi_dg_oo{d} = [dsi_dg_oo{d}; DG_cut{d}.dsindex{T}(idx_sub{sub_i}(cc))];
        end
    end
end



figure
subplot(1,2,1)
errorbar(xsort/pi*180, rho_dg_all_mean(5, :), rho_dg_all_ste(5, :), 'b')
hold on
errorbar(xsort/pi*180, mean(rho_dg_oo{5}), std(rho_dg_oo{5})/sqrt(size(rho_dg_oo{5}, 1)), 'r');
legend('ON DSGC', 'ON-OFF DSGC')
subplot(1,2,2)
errorbar(1,mean(dsi_dg_oo{5}), std(dsi_dg_oo{5})/sqrt(length(dsi_dg_oo{5})), 'ro');
hold on
errorbar(2, mean(cell2mat(dsi_dg{5}')), std(cell2mat(dsi_dg{5}'))/sqrt(length(cell2mat(dsi_dg{5}'))), 'bo')
ylim([0.7 1])
legend('ON-OFF DSGC', 'ON DSGC')

% DSI
% FigHandle = figure;
% set(FigHandle, 'Position', [100, 100, 1000, 500]);
% 
% xtick = ct;
% model_series = [dsi_dg_mean(1,1) dsi_dg_mean(2,1) dsi_dg_mean(3,1) dsi_dg_mean(4,1) dsi_dg_mean(5,1); dsi_dg_mean(1,2) dsi_dg_mean(2,2) dsi_dg_mean(3,2) dsi_dg_mean(4,2) dsi_dg_mean(5,2);dsi_dg_mean(1,3) dsi_dg_mean(2,3) dsi_dg_mean(3,3) dsi_dg_mean(4,3) dsi_dg_mean(5,3); dsi_dg_mean(1,4) dsi_dg_mean(2,4) dsi_dg_mean(3,4) dsi_dg_mean(4,4) dsi_dg_mean(5,4)];   
% model_error = [dsi_dg_ste(1,1) dsi_dg_ste(2,1) dsi_dg_ste(3,1) dsi_dg_ste(4,1) dsi_dg_ste(5,1); dsi_dg_ste(1,2) dsi_dg_ste(2,2) dsi_dg_ste(3,2) dsi_dg_ste(4,2) dsi_dg_ste(5,2);dsi_dg_ste(1,3) dsi_dg_ste(2,3) dsi_dg_ste(3,3) dsi_dg_ste(4,3) dsi_dg_ste(5,3); dsi_dg_ste(1,4) dsi_dg_ste(2,4) dsi_dg_ste(3,4) dsi_dg_ste(4,4) dsi_dg_ste(5,4)];
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
% 
% 
% 


%% dim flashes
Pc_on = [];
% 2016-02-19-0
load('DS160219.mat')
idx_on = [5 17 21 28];
Pc_on = [Pc_on; Pc(idx_on, :)];
% 2016-03-24-0
load('DS160324_SpikeDistance.mat')
idx_on = [6 16 30 44];
Pc_on = [Pc_on; Pc(idx_on, :)];

figure
plot(log10(Irel), Pc_on)

% compare across cell type
for ct = 1:4
    Pc_dir{ct} = Pc(idx_dir_flash{ct}, :);
    Pc_dir_mean(ct, :) = mean(Pc_dir{ct}, 1);
    Pc_dir_ste(ct, :) = std(Pc_dir{ct}, [], 1)/sqrt(length(idx_dir_flash{ct}));
end

Pc_dir_mean(5,:) = mean(Pc_on);
Pc_dir_ste(5,:) = std(Pc_on, [], 1)/sqrt(size(Pc_on, 1));

figure
for ct = 1:5
    errorbar(log10(Irel), Pc_dir_mean(ct, :), Pc_dir_ste(ct, :), 'color', color(ct));
    hold on
end
legend('superior', 'anterior', 'inferior', 'posterior', 'ON')
xlabel('log(R*/rod)')
ylabel('probability')

%% response polarity
% 2015-06-03-0
load('DS150603.mat')
id = {[1607 3902 4998 6737], [6347], []};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);

for i = 1:5
    ratio_on{i} = cell(3,1);
    for dir = 1:3
        if ~isempty(id{dir})
            ratio_on{i}{dir} = ratio{i}(idx{dir}, :);
        end
    end
end

% 2015-06-18-0
load('DS150618.mat')
id = {[6766], [], [1952 3182]};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);

for i = 1:5
    for dir = 1:3
        if ~isempty(id{dir})
            ratio_on{i}{dir} = [ratio_on{i}{dir}; ratio{i}(idx{dir}, :)];
        end
    end
end

% 2015-07-03-0
load('DS150703-1.mat')
id = {[261 4576 4923 5087 5239 6017], [2596 3391 4352 7397], [287 1653 6123]};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);
for i = 1:5
    for dir = 1:3
        if ~isempty(id{dir})
            ratio_on{i}{dir} = [ratio_on{i}{dir}; ratio{i}(idx{dir}, :)];
        end
    end
end




% 2016-01-30-0
load('DS160130.mat')
id = {[2582 3123 6917 7306], [4502], [392 4863]};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);
for i = 1:5
    for dir = 1:3
        if ~isempty(id{dir})
            ratio_on{i}{dir} = [ratio_on{i}{dir}; ratio{i}(idx{dir}, :)];
        end
    end
end

for i = 1:5
    for dir = 1:3
        ratio_on{i}{dir} = exciseRows_empty(ratio_on{i}{dir});
        ratio_on_mean(i, dir, :) = mean(ratio_on{i}{dir});
        ratio_on_ste(i, dir, :) = std(ratio_on{i}{dir})/sqrt(size(ratio_on{i}{dir}, 1));
    end
    ratio_on_all{i} = cell2mat(ratio_on{i});
    ratio_on_all_mean(i, :) = mean(ratio_on_all{i});
    ratio_on_all_ste(i, :) = std(ratio_on_all{i})/sqrt(size(ratio_on_all{i}, 1));
end

% plot average f2/f1
for T = 1:2
ct = {'superior', 'anterior', 'inferior'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = ct;
model_series = ratio_on_mean(:,:,T)';
model_error = ratio_on_ste(:,:,T)';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('F2/F1')
legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end
if T == 1
title('high speed')
else
title('low speed')
end
end

% plot average f2/f1
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = {'high speed', 'low speed'};
model_series = ratio_on_all_mean';
model_error = ratio_on_all_ste';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('F2/F1')
legend('NDF4','NDF3', 'NDF2', 'NDF1', 'NDF0', 'location', 'northeast');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

title('ON DSGC')

% compare
% 2015-06-03-0
load('DS150603.mat')
figure
for ct = 1:4
    errorbar(10.^(0:4), ratio_dir_mean(:, ct, 2), ratio_dir_ste(:, ct, 2), '--')
    hold on 
end
errorbar(10.^(0:4), ratio_oo_mean(:, 2), ratio_oo_ste(:, 2))
errorbar(10.^(0:4), ratio_on_all_mean(:, 2), ratio_on_all_ste(:, 2))
set(gca, 'XScale', 'Log')
legend('S', 'A', 'I', 'P', 'on-off','on')
ylim([0 3])
xlim(10.^[-0.5 4.5])
ylabel('F2/F1')
xlabel('R*/rod/s')

%% full field pulses

% 2015-06-03-0
load('DS150603.mat')
id = {[1607 3902 4998 6737], [6347], []};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);
for i = 1:5
    for dir = 1:3
        ratio_ffp_on{i}{dir} = ratio_ffp{i}(idx{dir});
    end
end

% 2015-06-18-0
load('DS150618.mat')
id = {[6766], [], [1952]};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);
for i = 1:5
    for dir = 1:3
        ratio_ffp_on{i}{dir} = [ratio_ffp_on{i}{dir} ratio_ffp{i}(idx{dir})];
    end
end

% 2015-07-03-0
load('DS150703-1.mat')
id = {[261 4576 4923 5239 6017], [2596 3391 4352 7397], [287 1653 6123]};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);
for i = 1:5
    for dir = 1:3
        ratio_ffp_on{i}{dir} = [ratio_ffp_on{i}{dir} ratio_ffp{i}(idx{dir})];
    end
end

% 2016-01-30-0
load('DS160130.mat')
id = {[2582 3123 6917 7306], [4502], [392 4863]};
[~, idx] = cellfun(@(id) intersect(ds_id, id), id, 'UniformOutput', false);
for i = 1:5
    for dir = 1:3
        ratio_ffp_on{i}{dir} = [ratio_ffp_on{i}{dir} ratio_ffp{i}(idx{dir})];
    end
end

for i = 1:5
    for dir = 1:3
        ratio_temp = ratio_ffp_on{i}{dir};
        ratio_temp(isnan(ratio_temp)) = [];
        ratio_ffp_avg(i,dir) = mean(ratio_temp);
        ratio_ffp_ste(i,dir) = std(ratio_temp)/sqrt(length(ratio_temp));
    end
    ratio_ffp_on_all{i} = cell2mat(ratio_ffp_on{i});
    ratio_ffp_on_all{i}(isnan(ratio_ffp_on_all{i})) = [];
    ratio_ffp_all_avg(i) = mean(ratio_ffp_on_all{i});
    ratio_ffp_all_ste(i) = std(ratio_ffp_on_all{i})/sqrt(length(ratio_ffp_on_all{i}));

end

figure
for dir = 1:3
    errorbar([0:4], ratio_ffp_avg(:, dir), ratio_ffp_ste(:,dir))
    hold on
end
title('ffp off on ratio')
xlabel('light level')
ylabel('off on ratio')
xlim([-0.5 4.5])

figure
errorbar([0:4], ratio_ffp_all_avg, ratio_ffp_all_ste)
title('ffp off on ratio')
xlabel('light level')
ylabel('off on ratio')
xlim([-0.5 4.5])

% compare with ON-OFF DSGCs
clear ratio_ffp_avg ratio_ffp_ste
load('DS150603.mat')
for d = 1:5
    ratio_temp = ratio_ffp{d}(idx_sub{2});
    ratio_temp(isnan(ratio_temp)) = [];
    ratio_ffp_avg(d) = mean(ratio_temp);
    ratio_ffp_ste(d) = std(ratio_temp)/sqrt(length(ratio_temp));
end
figure
errorbar([0:4], ratio_ffp_avg, ratio_ffp_ste)
hold on
errorbar([0:4], ratio_ffp_all_avg, ratio_ffp_all_ste)
legend('on-off', 'on')
title('ffp off on ratio')
xlabel('light level')
ylabel('off on ratio')
xlim([-0.5 4.5])
