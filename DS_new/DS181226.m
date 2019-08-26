cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

datadg = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-12-26-0/data005-sorted/data005-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2018-12-26-0/stimuli/s05.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

params_idx = [4 5]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

load('DS181226.mat')
datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-12-26-0/data001-data003-map/data001-data003-map', opt);
time_points = [2400];
datamfs(1:2) = split_datarun(datarun, time_points);
datamfs{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2018-12-26-0/stimuli/s01.mat';
datamfs{1} = load_stim_mfs(datamfs{1});

datamfs{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2018-12-26-0/stimuli/s03.mat';
datamfs{2} = load_stim_mfs(datamfs{2});

%% dg
n = 1;
i = 1;
[raster_dg, dg, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
dg{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);


delta_p = 1; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, dg{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(dg{i});
rep = datadg.stimulus.repetitions;
%%
L = 1;
mag_pca = MAG_all_norm_dg{L};
mag_pca = mag_pca';
[id_sub, idx_sub] = deal(cell(2, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(mag_pca);
pc1 = 1; pc2 = 3;
plot(scores(:, pc1), scores(:, pc2), 'o')
hold on
for i = 1:2
    [x, y] = ginput;
    plot(x, y)
    IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
    [~, idx_sub{i}] = find(IN' == 1);
    id_sub{i} = ds_id(idx_sub{i});
end
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
title('NDF 0')

figure
plot(scores(idx_sub{1}, pc1), scores(idx_sub{1}, pc2), 'ro', scores(idx_sub{2}, pc1), scores(idx_sub{2}, pc2), 'bo')
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
title('NDF 0')

v = 4*datadg.stimulus.params.SPATIAL_PERIOD./datadg.stimulus.params.TEMPORAL_PERIOD;
figure
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{1})), 'r')
hold on
semilogx(v, exciseColumn(MAG_all_norm_dg{L}(:, idx_sub{2})), 'b')
xlabel('micron/second')
ylabel('Response')
xlim([v(end) v(1)])

t = 5;
figure
compass(dg{1}.U{t}(idx_sub{1}), dg{1}.V{t}(idx_sub{1}), 'r')
hold on
compass(dg{1}.U{t}(idx_sub{2}), dg{1}.V{t}(idx_sub{2}), 'b')

%%
d = 1;
t = 5;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(dg{d}.U{t}, dg{d}.V{t});
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(dg{d}.U{t}, dg{d}.V{t}, x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = I;
    id_dir{i} = ds_id(idx_dir{i});
end


%% plot rfs
cell_type = {'superior', 'anterior', 'inferior', 'posterior'};
field_width = 17; field_height = 17;
field_width_sta = 40; field_height_sta = 40;
subregion = 0;
stop = 0.5; %second
for i = 1:2
%     fs_raster{i} = get_fs_raster(datafs{i}, ds_id, 'stop', 0.5);
    fs_raster{i} = get_fs_raster(datamfs{i}, ds_id);
    for cc = 1:length(ds_id)
        if fs_idx(cc)
            fs_raster{i}{cc} = [];
        end
    end
    fs_spike{i} = fs_get_spike(fs_raster{i});
    [rf_all{i}, rf_std{i}] = get_fs_rf(fs_spike{i}, field_width, field_height,subregion);
end

%%
for cc = 1:length(ds_id)
    figure(1)
    set(gcf, 'Position', [1 1 1000 600])
    id = ds_id(cc);
    for i = 1:2
        if ~isempty(rf_all{i}{cc})
%             rf = padarray(rf_all{i}{cc},[7,7]);
            rf = rf_all{i}{cc};
    
            subplot(2,3,3*(i-1)+1)
            imagesc(sum(rf,3))
            colormap gray
            axis image
            axis off

            subplot(2,3,3*(i-1)+2)
            imagesc(rf(:,:,1))
            colormap gray
            axis image
            axis off

            subplot(2,3,3*i)
            imagesc(rf(:,:,2))
            colormap gray
            axis image
            axis off

        end
    end
    print_close(1,[14 10],num2str(id))
%     pause
%     close(1)
end
%% 
for ll = 1:2
    trigger = datamfs{ll}.triggers(1:2:end);
    list = datamfs{ll}.stimulus.trial_list;
    repeat = datamfs{ll}.stimulus.repetitions;
    for i = 1:max(list)
        index(i, :) = find(list == i);
    end
    raster_onoff{ll} = cell(length(ds_id), 1);
    for i = 1:length(ds_id)
        if ~fs_idx(i)
            idx = get_cell_indices(datamfs{ll}, ds_id(i));
            raster = get_raster(datamfs{ll}.spikes{idx}, trigger, 'plot', false);
            raster_onoff{ll}{i} = raster(index);
        end
    end
end

%
bin_size = 0.02;
xx = [bin_size/2:bin_size:2-bin_size/2];
XX = [bin_size/2:bin_size:1-bin_size/2];
for ll = 1:2
    histg_onoff{ll} = cell(length(ds_id), 1);
    for cc = 1:length(raster_onoff{ll})
        if ~isempty(raster_onoff{ll}{cc})
            for p = 1:289%size(raster_onoff{ll}{1}, 1)
                raster_all_onoff{ll}{cc}{p} = sort(cell2mat(raster_onoff{ll}{cc}(p,:)'));
                histg_onoff{ll}{cc}{p} = hist(raster_all_onoff{ll}{cc}{p}, xx);
                for onoff = 1:2
                    raster_all{ll}{cc}{onoff}{p} = sort(cell2mat(fs_raster{ll}{cc}(p,:,onoff)'));
                    histg_all{ll}{cc}{onoff}{p} = hist(raster_all{ll}{cc}{onoff}{p}, XX);
                end
            end
        end
    end
end

h = subplot(field_height, field_width, 1); p1 = get(h, 'pos');
h = subplot(field_height, field_width, 2); p2 = get(h, 'pos');
width = (p2(1) - p1(1))*0.8;

h = subplot(field_height, field_width, 1); p1 = get(h, 'pos');
h = subplot(field_height, field_width, field_width+1); p2 = get(h, 'pos');
height = p1(2) - p2(2);
for cc = 2:2%length(ds_id)
    for ll = 3:3
        H = figure(1);
        set(H, 'Position', [1 1 1080 1080])
        if ~fs_idx(cc,ll)
            for p = 1:289
                h = subplot(field_width, field_height, p);
                l = get(h, 'pos');
                l(3) = width; l(4) = height;
                set(h, 'pos', l);
%                 plot(xx, histg_onoff{ll}{cc}{p})
                n = 1/bin_size;
                plot(xx(1:n), histg_onoff{ll}{cc}{p}(1:n), 'b')
                hold on
                plot(xx(n+1:end), histg_onoff{ll}{cc}{p}(n+1:end), 'r')
                ylim([0 max(max(cell2mat(histg_onoff{ll}{cc}')))])
                axis off
            end
        end
%     pause 
%     close all
        name = [num2str(ds_id(cc)) '_' LL{ll}];
%         print_close(1, [14 14], name);
    end
end


%% fit RF with Gaussian                
PixelArea = (30*4)^2/10^6;
% fit and compute rf area
for ds = 1:2
    for dir = 1:4
        clear rf_area_temp
        for onoff = 1:2
            rf_area_temp = [];
            for cc = 1:length(id_dir{dir})
                if ~fs_idx(idx_dir{dir}(cc))
                    data = rf_all{ds}{idx_dir{dir}(cc)}(:, :, onoff);
    %                     data = rf_wt{dir}{cc}{onoff};
                    if sum(sum(data > mean(data(:))+4*std(data(:))))>0
    %                         figure(100)
    %                         imagesc(data)
    %                         colormap gray
    %                         pause

                        params = fit_2d_gaussian(data);
    %                     Gaussian_params{ll}{dir}{cc}{onoff} = params;
                        rf_area_temp = [rf_area_temp params.xd * params.yd * pi * PixelArea];
                    end
                end
            end
            rf_area{ds}{dir}{onoff} = rf_area_temp(rf_area_temp < 0.3);
%             rf_area{ds}{dir}{onoff} = rf_area_temp;
        end
    end
end

% exclude outliers
stdn = 1.7;
for ds = 1:2
    for dir = 1:4
        for onoff = 1:2
            notdone = 1;
            rf_area_temp = rf_area{ds}{dir}{onoff};
            while notdone
                a = length(rf_area_temp);
                rf_area_temp(rf_area_temp > std(rf_area_temp)*stdn + mean(rf_area_temp)) = [];
                b = length(rf_area_temp);
                if a == b
                    notdone = 0;
                    rf_area_clean{ds}{dir}{onoff} = rf_area_temp;
                end
            end
            rf_area_clean_mean{ds}{onoff}(dir) = mean(rf_area_clean{ds}{dir}{onoff});
            rf_area_clean_ste{ds}{onoff}(dir) = std(rf_area_clean{ds}{dir}{onoff})/sqrt(length(rf_area_clean{ds}{dir}{onoff}));
        end
    end
end

% plot 
color = 'brgkc';
figure
for onoff = 1:2
    subplot(1,2,onoff)
    for dir = 1:4
        ll = 1;
        n = length(rf_area_clean{ll}{dir}{onoff});
        h(ll) = plot((dir-1)*5+ll*ones(n,1), rf_area_clean{ll}{dir}{onoff}, [color(ll) 'o']);
        hold on
        errorbar((dir-1)*5+ll+0.5, rf_area_clean_mean{ll}{onoff}(dir), rf_area_clean_ste{ll}{onoff}(dir), [color(ll) 'd'])
        ll = 2;
        n = length(rf_area_clean{ll}{dir}{onoff});
        h(ll) = plot((dir-1)*5+ll*ones(n,1), rf_area_clean{ll}{dir}{onoff}, [color(ll) 'o']);
        errorbar((dir-1)*5+ll+0.5, rf_area_clean_mean{ll}{onoff}(dir), rf_area_clean_ste{ll}{onoff}(dir), [color(ll) 'd'])
    end
%     set(gca, 'yscale', 'log')
    legend([h(1), h(2)], 'control', 'SR')
    if onoff == 1
        title('ON')
    else
        title('OFF')
    end
    ylabel('RF area (mm^2)')
    set(gca, 'xtick', [])
%         ylim([0 0.3])

end


        
%% combine datasets
load('DS181226.mat')
rf_area_clean_all{1} = rf_area_clean;
load('DS190107.mat')
rf_area_clean_all{2} = rf_area_clean;
clear rf_area_clean
for ds = 1:2
    for ct = 1:4
        for onoff = 1:2
            rf_area_clean{ds}{ct}{onoff} = [rf_area_clean_all{1}{ds}{ct}{onoff} rf_area_clean_all{2}{ds}{ct}{onoff}];
            rf_area_clean_mean{ds}{onoff}(ct) = mean(rf_area_clean{ds}{ct}{onoff});
            rf_area_clean_ste{ds}{onoff}(ct) = std(rf_area_clean{ds}{ct}{onoff})/sqrt(length(rf_area_clean{ds}{ct}{onoff}));
        end
    end
end


% plot 
color = 'brgkc';
figure
for onoff = 1:2
    subplot(1,2,onoff)
    for dir = 1:4
        for ll = 1:2
            n = length(rf_area_clean{ll}{dir}{onoff});
            h(ll) = plot((dir-1)*5+ll*ones(n,1), rf_area_clean{ll}{dir}{onoff}, [color(ll) 'o']);
            hold on
            errorbar((dir-1)*5+ll+0.5, rf_area_clean_mean{ll}{onoff}(dir), rf_area_clean_ste{ll}{onoff}(dir), [color(ll) 'd'])
        end
    end
%     set(gca, 'yscale', 'log')
    legend([h(1), h(2)], 'control', 'SR')
    if onoff == 1
        title('ON')
    else
        title('OFF')
    end
    ylabel('RF area (mm^2)')
    set(gca, 'xtick', [])
         ylim([0 0.2])

end

%%
for ll = 1:2
    for onoff = 1:2
        rf_area_clean_combine{ll}{onoff}{1} = rf_area_clean{ll}{1}{onoff};
        temp = cat(1, rf_area_clean{ll}{2:4});
        rf_area_clean_combine{ll}{onoff}{2} = cell2mat(temp(:, onoff)');
    end
end

stdn = 1.7;
for ds = 1:2
    for dir = 1:2
        for onoff = 1:2
            notdone = 1;
            rf_area_temp = rf_area_clean_combine{ds}{onoff}{dir};
            while notdone
                a = length(rf_area_temp);
                rf_area_temp(rf_area_temp > std(rf_area_temp)*stdn + mean(rf_area_temp)) = [];
%                 rf_area_temp(rf_area_temp < -std(rf_area_temp)*stdn + mean(rf_area_temp)) = [];
                b = length(rf_area_temp);
                if a == b
                    notdone = 0;
                    rf_area_clean_combine{ds}{onoff}{dir} = rf_area_temp;
                end
            end
            rf_area_clean_combine_mean{ds}{onoff}(dir) = mean(rf_area_clean_combine{ds}{onoff}{dir});
            rf_area_clean_combine_ste{ds}{onoff}(dir) = std(rf_area_clean_combine{ds}{onoff}{dir})/sqrt(length(rf_area_clean{ds}{onoff}{dir}));
        end
    end
end

figure
TITLE = {'ON', 'OFF'};
for onoff = 1:2
    subplot(1,2,onoff)
    marker = 'od';
    for ds = 1:2
        for ct = 1:2
            plot((ds-1)*2+ct*ones(1, length(rf_area_clean_combine{ds}{onoff}{ct})), rf_area_clean_combine{ds}{onoff}{ct}, ['k' marker(ct)], 'markersize', 7)
            hold on
        end
        errorbar((ds-1)*2+[1 2], cellfun(@mean, rf_area_clean_combine{ds}{onoff}), cellfun(@std, rf_area_clean_combine{ds}{onoff})./sqrt(cellfun(@length, rf_area_clean_combine{ds}{onoff})), 'ks-', 'markersize', 15);
    end
    xlim([0.5 4.5])
    ylabel('log(R*/rod)') 
    legend('superior', 'other')
    title(TITLE{onoff})
    ylim([0 0.15])
end

%%
temp_mean = cat(1, rf_area_clean_combine_mean{:});
temp_ste = cat(1, rf_area_clean_combine_ste{:});
for onoff = 1:2
    rf_mean_onoff{onoff} = cell2mat(temp_mean(:, onoff));
    rf_ste_onoff{onoff} = cell2mat(temp_ste(:, onoff));
end

ct = {'control', 'SR'};
figure
for onoff = 1:2
    subplot(1,2,onoff)
    set(gcf, 'DefaultLineLineWidth', 1.5)
    xtick = ct;
    model_series = rf_mean_onoff{onoff};
    model_error = rf_ste_onoff{onoff};
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('RF area (mm^2)')
    legend('superior','other');
    hold on;

    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 

    groupwidth = min(0.8, numbars/(numbars+1.5));

    for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end
    ylim([0 0.1])
    xlim([.5 2.5])
end

%% center pixel gain
repeat = 5;
LL = {'control', 'SR'};
clear hist_oo_ct

c = [bin_size/2:bin_size:2-bin_size/2];
figure
for ll = 1:2
    for ct = 1:4
        CC = 1;
        for cc = 1:length(id_dir{ct})
            if ~isempty(histg_onoff{ll}{idx_dir{ct}(cc)})
                [~, argmax] = max(cellfun(@max, histg_onoff{ll}{idx_dir{ct}(cc)}));
                hist_oo_ct{ll, ct}(CC, :) = histg_onoff{ll}{idx_dir{ct}(cc)}{argmax}/repeat/bin_size;
                CC = CC + 1;
            end
        end
    end
    subplot(2,1,ll)
    for ct = 1:4
        plot(xx, mean(hist_oo_ct{ll, ct}))
        hold on
    end
     title(LL{ll})
    ylabel('firing rate (Hz)')
    if ll == 1
        legend('superior', 'anterior', 'inferior', 'posterior')
    end
    if ll == 2
        xlabel('time (s)')
    end
end

figure
color = 'kr';
for ll = 1:2
    hist_oo_ct_combine{ll, 1} = hist_oo_ct{ll, 1};
    hist_oo_ct_combine{ll, 2} = cell2mat(hist_oo_ct(ll, 2:4)');
    subplot(2,1,ll)
    for ct = 1:2
        hist_temp = hist_oo_ct_combine{ll, ct};
        plot(xx, mean(hist_temp), color(ct))
        hold on
        shadedplot(xx, mean(hist_temp)+std(hist_temp)/sqrt(size(hist_temp, 1)), mean(hist_temp)-std(hist_temp)/sqrt(size(hist_temp, 1)), color(ct), color(ct));
        hold on
    end
    title(LL{ll})
    ylabel('firing rate (Hz)')
    if ll == 1
        legend('superior', 'others')
    end
    if ll == 2
        xlabel('time (s)')
    end
end


% mean spike number
center_gain = cellfun(@sum, cellfun(@mean, hist_oo_ct_combine(1, :), 'UniformOutput', 0));

%% non-parametric RF (total spike number)
ds = 2;
idx{1} = idx_dir{1};
idx{2} = cell2mat(idx_dir(2:4));
for ds = 1:2
    for ct = 1:2
        CC = 1;
        for cc = 1:length(idx{ct})
            if ~isempty(rf_all{ds}{cc})
                for onoff = 1:2
                    temp = rf_all{ds}{cc}(:, :, onoff);
                    temp = temp(:);
                    temp(temp < mean(temp) + 1.5*std(temp)) = 0;
                    
                    rf_spikes{ds}{ct}(CC, onoff) = sum(temp);
                end
                CC = CC + 1;
            end
        end
    end
end








