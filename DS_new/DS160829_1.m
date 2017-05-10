cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);
load('DS160829.mat')

datadg = load_data('/Volumes/lab/analysis/2016-08-29-0/data003-sorted/data003-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-08-29-0/stimuli/s03.mat';
datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);

% [NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
% ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
% params_idx = [2 5]; % which parameters to use for classification

% [ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

datafs{1} = load_data('/Volumes/lab/analysis/2016-08-29-0/data000-map/data000-map', opt);
datafs{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-08-29-0/stimuli/s00.mat';
datafs{1} = load_stim_matlab(datafs{1});

datafs{2} = load_data('/Volumes/lab/analysis/2016-08-29-0/data001-map/data001-map', opt);
datafs{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-08-29-0/stimuli/s01.mat';
datafs{2} = load_stim_matlab(datafs{2});

datafs{3} = load_data('/Volumes/lab/analysis/2016-08-29-0/data004-map/data004-map', opt);
datafs{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-08-29-0/stimuli/s04.mat';
datafs{3} = load_stim_matlab(datafs{3});

ds_id_fs = intersect(datafs{3}.cell_ids, ds_id);

%% dg
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);


delta_p = 4; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(DG{i});
rep = datadg.stimulus.repetitions;

%% plot individual cells
[DG_cut, raster_dg_cut, raster_p_sum_cut] = cut_dg(DG, raster_dg, raster_p_sum, 2, [5]);
for cc = 64:64 %length(ds_id)
    plot_ds_raster(DG, raster_dg, cc, ds_id(cc), '', 1, 1, 0)
end
%% classification based on speed tunning
L = 1;
mag_pca = MAG_all_norm_dg{L};
mag_pca = mag_pca';
[id_sub, idx_sub] = deal(cell(2, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(mag_pca);
pc1 = 1; pc2 = 2;
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

t = 4;
figure
compass(DG{1}.U{t}(idx_sub{1}), DG{1}.V{t}(idx_sub{1}), 'r')
hold on
compass(DG{1}.U{t}(idx_sub{2}), DG{1}.V{t}(idx_sub{2}), 'b')

%%
d = 1;
t = 5;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(DG{d}.U{t}(idx_sub{2}), DG{d}.V{t}(idx_sub{2}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}(idx_sub{2}), DG{d}.V{t}(idx_sub{2}), x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = idx_sub{2}(I);
    id_dir{i} = ds_id(idx_dir{i});
end

d = 1;
t = 3;
h = figure;
dirn = 3;
set(h, 'Position', [1 1 1080 500])
compass(DG{d}.U{t}(idx_sub{1}), DG{d}.V{t}(idx_sub{1}));
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}(idx_sub{1}), DG{d}.V{t}(idx_sub{1}), x, y);
    [~, I] = find(IN == 1);
    idx_dir_on{i} = idx_sub{1}(I);
    id_dir_on{i} = ds_id(idx_dir_on{i});
end

%% plot rfs
cell_type = {'superior', 'anterior', 'inferior', 'posterior'};
field_width = 20; field_height = 20;
field_width_sta = 30; field_height_sta = 30;
subregion = 0;
stop = 0.5; %second
for i = 1:3
%     fs_raster{i} = get_fs_raster(datafs{i}, ds_id, 'stop', 0.5);
    fs_raster{i} = get_fs_raster(datafs{i}, ds_id);
    for cc = 1:length(ds_id)
        if fs_idx(cc,i)
            fs_raster{i}{cc} = [];
        end
    end
    fs_spike{i} = fs_get_spike(fs_raster{i});
    [rf_all{i}, rf_std{i}] = get_fs_rf(fs_spike{i}, field_width, field_height,subregion);
end

for cc = 1:10%length(ds_id)
    set(gcf, 'Position', [1 1 1000 1000])
    id = ds_id(cc);
    for i = 1:3
        if ~isempty(rf_all{i}{cc})
            rf = padarray(rf_all{i}{cc},[7,7]);
    
            subplot(3,3,3*(i-1)+1)
            imagesc(sum(rf,3))
            colormap gray
            axis image

            subplot(3,3,3*(i-1)+2)
            imagesc(rf(:,:,1))
            colormap gray
            title('on')
            axis image

            subplot(3,3,3*(i-1)+3)
            imagesc(rf(:,:,2))
            colormap gray
            title('off')
            axis image
        end
    end
    print_close(1,[12 12],num2str(id))
end

%%
% fs_spike_temp = cellfun(@(fs_spike) sum(fs_spike,2), fs_spike, 'UniformOutput', false);
% fs_spike_temp = cellfun(@(fs_spike_temp) sum(fs_spike_temp,3), fs_spike_temp, 'UniformOutput', false);
% [~, I] = cellfun(@(fs_spike_temp) sort(fs_spike_temp, 'descend'), fs_spike_temp, 'UniformOutput', false);trigger = datafs.triggers(2:2:end);
for ll = 1:3
    trigger = datafs{ll}.triggers(2:2:end);
    list = datafs{ll}.stimulus.trial_list;
    repeat = datafs{ll}.stimulus.repetitions;
    for i = 1:max(list)
        index(i, :) = find(list == i);
    end
    raster_onoff{ll} = cell(length(ds_id_fs), 1);
    for i = 1:length(ds_id_fs)
        if ~fs_idx(i,ll)
            idx = get_cell_indices(datafs{ll}, ds_id_fs(i));
            raster = get_raster(datafs{ll}.spikes{idx}, trigger, 'plot', false);
            raster_onoff{ll}{i} = raster(index);
        end
    end
end

%%
LL = {'NDF4', 'NDF2', 'NDF0'};
bin_size = 0.01;
xx = [bin_size/2:bin_size:2-bin_size/2];
XX = [bin_size/2:bin_size:1-bin_size/2];
for ll = 1:3
    for cc = 1:length(raster_onoff{ll})
        if ~isempty(raster_onoff{ll}{cc})
            for p = 1:400%size(raster_onoff{ll}{1}, 1)
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

h = subplot(20, 20, 1); p1 = get(h, 'pos');
h = subplot(20, 20, 2); p2 = get(h, 'pos');
width = (p2(1) - p1(1))*0.8;

h = subplot(20, 20, 1); p1 = get(h, 'pos');
h = subplot(20, 20, 21); p2 = get(h, 'pos');
height = p1(2) - p2(2);
for ll = 1:3
    for cc = 1:length(ds_id)
        H = figure(1);
        set(H, 'Position', [1 1 1080 1080])
        if ~fs_idx(cc,ll)
            for p = 1:400
                h = subplot(20, 20, p);
                l = get(h, 'pos');
                l(3) = width; l(4) = height;
                set(h, 'pos', l);
                plot(xx, histg_onoff{ll}{cc}{p})
                ylim([0 max(max(cell2mat(histg_onoff{ll}{cc}')))])
                axis off
            end
        end
%     pause 
%     close all
        name = [num2str(ds_id(cc)) '_' LL{ll}];
        print_close(1, [14 14], name);
    end
end

        
%% exclude noisy pixels
threshold = 0.3;
for ll = 1:3
    for cc = 1:length(ds_id)
        if ~fs_idx(cc,ll)
            for onoff = 1:2
                rf_temp = rf_all{ll}{cc}(:,:,onoff);
                max_temp = max(rf_temp(:));
                rf_temp(rf_temp < max_temp*threshold) = 0;
                rf_clean{ll}{cc}(:,:,onoff) = rf_temp;
            end
        end
    end
end

for cc = 1:length(ds_id)
    set(gcf, 'Position', [1 1 1000 1000])
    id = ds_id(cc);
    for i = 1:3
        if ~isempty(rf_clean{i}{cc})
            rf = padarray(rf_clean{i}{cc},[7,7]);
    
            subplot(3,3,3*(i-1)+1)
            imagesc(sum(rf,3))
            colormap gray
            axis image

            subplot(3,3,3*(i-1)+2)
            imagesc(rf(:,:,1))
            colormap gray
            title('on')
            axis image

            subplot(3,3,3*(i-1)+3)
            imagesc(rf(:,:,2))
            colormap gray
            title('off')
            axis image
        end
    end
    print_close(1,[12 12],num2str(id))
end

%% superior RF size estimation
LL = {'NDF2', 'NDF0'};
bin_n = 30;
threshold = 0.3;
PixelArea = (20*4)^2/10^6; %mm^2

% pulse filter
for j = 11:11; %second
    clear rf_wt rf_wt_area rf_wt_area_mean rf_wt_area_ste rf_area rf_area_clean
    tau = 0.01*j;
    tt = -3*tau:bin_size:3*tau;
    filter = exp(-tt.^2/(2*tau^2));
    npixel = 3;
    for dir = 1:4
        
    for ll = 2:3
        for cc = 1:length(id_dir{dir})
            if ~fs_idx(idx_dir{dir}(cc), ll)
                for onoff = 1:2
                    % use the npixel brightest pixels to calculate a psth template
                    [~, idx] = sort(cellfun(@max,histg_all{ll}{idx_dir{dir}(cc)}{onoff}), 'descend');
                    idx = idx(1:npixel);
                    temp = histg_all{ll}{idx_dir{dir}(cc)}{onoff}(idx);
                    for i = 1:npixel
                        temp{i} = temp{i}/norm(temp{i});
                    end
                    temp_mean = mean(cell2mat(temp'));
                    temp_mean_norm = temp_mean/norm(temp_mean);
%                     temp_mean_norm = conv(temp_mean_norm, filter, 'same');

    %                 plot(temp_mean_norm)
    %                 pause
                    for p = 1:400
                        rf_wt{ll}{dir}{cc}{onoff}(p) = histg_all{ll}{idx_dir{dir}(cc)}{onoff}{p}*temp_mean_norm';
                    end
                    rf_wt{ll}{dir}{cc}{onoff} = reshape(rf_wt{ll}{dir}{cc}{onoff}, field_width, field_height);
                    rf_wt_area{dir}{ll-1}{onoff}{cc} = sum(sum(rf_wt{ll}{dir}{cc}{onoff} > max(rf_wt{ll}{dir}{cc}{onoff}(:))*threshold))*PixelArea;
                end
            end
        end
        for onoff = 1:2
            rf_wt_area{dir}{ll-1}{onoff} = cell2mat(rf_wt_area{dir}{ll-1}{onoff});
            rf_wt_area_mean(dir, ll-1, onoff) = mean(rf_wt_area{dir}{ll-1}{onoff});
            rf_wt_area_ste(dir, ll-1, onoff) = std(rf_wt_area{dir}{ll-1}{onoff})/sqrt(length(rf_wt_area{dir}{ll-1}{onoff}));
        end
    end
    end
%     subplot(3,3,j)
%     imagesc(rf_wt{3}{1}{2}')
%     colormap gray


% fit RF with Gaussian

for dir = 1:4;
    for ll = 2:3
        clear rf_area_temp
        for cc = 1:length(id_dir{dir})
            if ~isempty(rf_wt{ll}{dir}{cc})
                for onoff = 1:2
                    data = rf_wt{ll}{dir}{cc}{onoff};
                    params = fit_2d_gaussian(data);
    %                 pause
                    Gaussian_params{ll}{dir}{cc}{onoff} = params;
                    rf_area_temp{cc}(onoff) = params.xd * params.yd * pi;
                end
            end
        end
        rf_area{ll}{dir} = cell2mat(rf_area_temp');
    end
end

% exclude outliers
for ll = 2:3
    for dir = 1:4
        for onoff = 1:2
            notdone = 1;
            rf_area_temp = rf_area{ll}{dir}(:,onoff);
            while notdone
                a = length(rf_area_temp);
                rf_area_temp(rf_area_temp > std(rf_area_temp)*2 + mean(rf_area_temp)) = [];
                b = length(rf_area_temp);
                if a == b
                    notdone = 0;
                    rf_area_clean{ll}{dir}{onoff} = rf_area_temp;
                end
            end
        end
    end
end

rf_area_clean_all{j} = rf_area_clean;
end

% plot 
for onoff = 1:2
    figure
    for dir = 1:4
        for ll = 2:3
            n = size(rf_area_clean{ll}{dir}, 1);
            plot((dir-1)*5+ll*ones(n,1), rf_area_clean{ll}{dir}{onoff}, 'o')
            hold on
        end
    end
end

for i = 1:length(rf_area_clean_all)
    for ll = 2:3
        for dir = 1:4
            for onoff = 1:2
                rf_area_mean{ll}{onoff}(dir,i) = mean(rf_area_clean_all{i}{ll}{dir}{onoff})*PixelArea;
                rf_area_ste{ll}{onoff}(dir,i) = std(rf_area_clean_all{i}{ll}{dir}{onoff})/sqrt(length(rf_area_clean_all{i}{ll}{dir}{onoff}))*PixelArea;
            end
        end
    end
end

color = 'brgk';
figure
for ll = 2:3
    for onoff = 1:2
        subplot(2,2,(ll-2)*2+onoff)
        dir = 1;
        errorbar(rf_area_mean{ll}{onoff}(dir,1:10), rf_area_ste{ll}{onoff}(dir,1:10), 'color', color(dir))
        hold on
        for dir = 2:4
            errorbar(rf_area_mean{ll}{onoff}(dir,1:10), rf_area_ste{ll}{onoff}(dir,1:10), 'color', color(dir))
            hold on
        end
        xlabel('filter size')
        ylabel('RF area (mm^2)')
        legend(cell_type)
        xlim([0 11])
    end
end

%% added @ 2016-11-30
% fit RF with Gaussian
                
PixelArea = (15*4)^2/10^6;
% fit and compute rf area
for ll = 1:3
    for dir = 1:4
        clear rf_area_temp
        for onoff = 1:2
            rf_area_temp = [];
            for cc = 1:length(id_dir{dir})
                if ~isempty(rf_all{ll}{idx_dir{dir}(cc)})
                    data = rf_all{ll}{idx_dir{dir}(cc)}(:, :, onoff);
                    if sum(sum(data > mean(data(:))+3*std(data(:))))>0
%                         figure(100)
%                         imagesc(data)
%                         colormap gray

                        params = fit_2d_gaussian(data);
    %                     Gaussian_params{ll}{dir}{cc}{onoff} = params;
                        rf_area_temp = [rf_area_temp params.xd * params.yd * pi * PixelArea];
                        
%                         params.xd * params.yd * pi * PixelArea
%                         id_dir{dir}(cc)
%                         pause

                    end
                end
            end
            rf_area{ll}{dir}{onoff} = rf_area_temp;
        end
    end
end


% exclude outliers
stdn = 2;
for ll = 1:3
    for dir = 1:4
        for onoff = 1:2
            notdone = 1;
            rf_area_temp = rf_area{ll}{dir}{onoff};
            while notdone
                a = length(rf_area_temp);
                rf_area_temp(rf_area_temp > mean(rf_area_temp) + std(rf_area_temp)*stdn) = [];
                rf_area_temp(rf_area_temp < mean(rf_area_temp) - std(rf_area_temp)*stdn) = [];
                b = length(rf_area_temp);
                if a == b
                    notdone = 0;
                    rf_area_clean{ll}{dir}{onoff} = rf_area_temp;
                end
            end
            rf_area_clean_mean{onoff}(ll, dir) = mean(rf_area_clean{ll}{dir}{onoff});
            rf_area_clean_ste{onoff}(ll, dir) = std(rf_area_clean{ll}{dir}{onoff})/sqrt(length(rf_area_clean{ll}{dir}{onoff}));
        end
    end
end




% plot 
color = 'brgkc';
figure
for onoff = 1:2
    subplot(1,2,onoff)
    for dir = 1:4
        for ll = 2:3
            n = size(rf_area_clean{ll}{dir}{onoff}, 1);
%             n = length(rf_area{ll}{dir}{onoff});
            h{ll} = plot((dir-1)*5+ll*ones(n,1), rf_area_clean{ll}{dir}{onoff}, [color(ll) 'o']);
%             h{ll} = plot((dir-1)*5+ll*ones(n,1), rf_area{ll}{dir}{onoff}, [color(ll) 'o']);
            hold on
        end
    end
%     set(gca, 'yscale', 'log')
    legend([h{2}(1), h{3}(1)], 'NDF 2', 'NDF 0')
    if onoff == 1
        title('ON')
    else
        title('OFF')
    end
    ylabel('RF area (mm^2)')
    set(gca, 'xtick', [])
%     ylim([0 0.4])
end

figure
for onoff = 1:2
    subplot(1,2,onoff)
    for dir = 1:4
        errorbar([2 4], rf_area_clean_mean{onoff}(2:3, dir), rf_area_clean_ste{onoff}(2:3, dir))
        hold on
    end
%     ylim([0 0.05])
    xlabel('log(background intensity)')
    ylabel('RF area (mm^2)')
    legend(cell_type)
end

%%
dir = 4;
onoff = 2;
y = [];group = [];
for ll = 2:3
    y = [y rf_area_clean{ll}{dir}{onoff}];
    group = [group ones(1, length(rf_area_clean{ll}{dir}{onoff}))*ll];
end

anova1(y, group)

%
onoff = 2;
y = [];group = [];
for ll = 3:4
    y = [y rf_area_clean_all{ll}{onoff}];
    group = [group ones(1, length(rf_area_clean_all{ll}{onoff}))*ll];
end

anova1(y, group)
%%
stdn = 2;
for ll = 2:3
    for onoff = 1:2
        rf_area_temp = [];
        for dir = 2:4
            rf_area_temp = [rf_area_temp rf_area{ll}{dir}{onoff}];
        end
        notdone = 1;
        while notdone
            a = length(rf_area_temp);
            rf_area_temp(rf_area_temp > mean(rf_area_temp) + std(rf_area_temp)*stdn) = [];
            rf_area_temp(rf_area_temp < mean(rf_area_temp) - std(rf_area_temp)*stdn) = [];
            b = length(rf_area_temp);
            if a == b
                notdone = 0;
                rf_area_clean_all{ll}{onoff} = rf_area_temp;
            end
        end
        rf_area_clean_all_mean{onoff}(ll) = mean(rf_area_clean_all{ll}{onoff});
        rf_area_clean_all_ste{onoff}(ll) = std(rf_area_clean_all{ll}{onoff})/sqrt(length(rf_area_clean_all{ll}{onoff}));
    end
end


figure
for onoff = 1:2
    errorbar([2 4], rf_area_clean_all_mean{onoff}(2:3), rf_area_clean_all_ste{onoff}(2:3))
    hold on
%     ylim([0 0.05])
    xlabel('log(background intensity)')
    ylabel('RF area (mm^2)')
    legend('ON', 'OFF')
end

onoff = 1;
y = [];group = [];
for ll = 2:3
    y = [y rf_area_clean_all{ll}{onoff}];
    group = [group ones(1, length(rf_area_clean_all{ll}{onoff}))*ll];
end

anova1(y, group)
