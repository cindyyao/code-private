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

datafs{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-08-29-0/data000-map/data000-map', opt);
datafs{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-08-29-0/stimuli/s00.mat';
datafs{1} = load_stim_matlab(datafs{1});

datafs{2} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-08-29-0/data001-map/data001-map', opt);
datafs{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-08-29-0/stimuli/s01.mat';
datafs{2} = load_stim_matlab(datafs{2});

datafs{3} = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-08-29-0/data004-map/data004-map', opt);
datafs{3}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-08-29-0/stimuli/s04.mat';
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

for cc = 14:14%length(ds_id)
    figure
    set(gcf, 'Position', [1 1 1000 1000])
    id = ds_id(cc);
    for i = 1:3
        if ~isempty(rf_all{i}{cc})
%             rf = padarray(rf_all{i}{cc},[7,7]);
            rf = rf_all{i}{cc};
    
            subplot(3,3,3*(i-1)+1)
            imagesc(sum(rf,3))
%             colormap gray
            axis image
            axis off

            subplot(3,3,3*(i-1)+2)
            imagesc(rf(:,:,1))
%             colormap gray
%             title('on')
            axis image
            axis off

            subplot(3,3,3*(i-1)+3)
            imagesc(rf(:,:,2))
%             colormap gray
%             title('off')
            axis image
            axis off

        end
    end
%     print_close(1,[12 12],num2str(id))
end

max_spike = cell(4,1);
figure
for ct = 1:4
    subplot(2,2,ct)
    for cc = 1:length(id_dir_fs{ct})
        max_spike{ct}(cc) = max(rf_all{idx_dir_fs{ct}(cc)}(:));
    end
    hist(max_spike{ct})
    xlabel('spike #')
    ylabel('cell #')
    title(cell_type{ct})
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

%% get 12s spontaneous activity
% data000 use last 12s (note: color [0.5 0.5 0.5], not [0.1 0.1 0.1])
for cc = 1:length(ds_id)
    if ~fs_idx(cc,1)
        idx = get_cell_indices(datafs{1},ds_id(cc));
        spikes_temp = datafs{1}.spikes{idx};
        bgnd_firing{1}{cc} = spikes_temp(spikes_temp>datafs{1}.duration-12) - (datafs{1}.duration-12);
    end
end

% data001, data004: use first 12s (note: color [0.1 0.1 0.1])
for ll = 2:3
    for cc = 1:length(ds_id)
        if ~fs_idx(cc,ll)
            idx = get_cell_indices(datafs{ll},ds_id(cc));
            spikes_temp = datafs{ll}.spikes{idx};
            bgnd_firing{ll}{cc} = spikes_temp(spikes_temp < 52 & spikes_temp > 40)-40;
        end
    end
end

for dir = 1:4
    for cc = 1:length(id_dir{dir})
        figure(1)
        for ll = 1:3
            if ~fs_idx(idx_dir{dir}(cc), ll)
                subplot(3,1,ll)
                stem(bgnd_firing{ll}{idx_dir{dir}(cc)}, ones(length(bgnd_firing{ll}{idx_dir{dir}(cc)}),1))
                xlim([0 12])
            end
        end
        dir
        pause
        close(1)
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

%% fit RF with Gaussian
PixelArea = (20*4)^2/10^6;
% fit and compute rf area
for ll = 2:3
    for dir = 1:3
        clear rf_area_temp
        for cc = 1:length(id_dir{dir+1})
            if ~isempty(rf_all{ll}{idx_dir{dir+1}(cc)})
                for onoff = 1:2
                    data = rf_all{ll}{idx_dir{dir+1}(cc)}(:, :, onoff);
                    params = fit_2d_gaussian(data);
                    Gaussian_params{ll}{dir}{cc}{onoff} = params;
                    rf_area_temp{cc}(onoff) = params.xd * params.yd * pi * PixelArea;
                end
            end
        end
        rf_area{ll}{dir+1} = cell2mat(rf_area_temp');
    end
end

for ll = 2:3
    for dir = 1:4
        clear rf_area_temp
        for cc = 1:length(id_dir{dir})
            if ~isempty(rf_all{ll}{idx_dir{dir}(cc)})
                for onoff = 1:2
                    data = rf_all{ll}{idx_dir{dir}(cc)}(:, :, onoff);
                    params = fit_2d_gaussian(data);
                    Gaussian_params{ll}{dir}{cc}{onoff} = params;
                    rf_area_temp{cc}(onoff) = params.xd * params.yd * pi * PixelArea;
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


% plot 
color = 'brgkc';
figure
for onoff = 1:2
    subplot(1,2,onoff)
    for dir = 1:4
        for ll = 2:3
            n = size(rf_area_clean{ll}{dir}, 1);
            h{ll-1} = plot((dir-1)*5+ll*ones(n,1), rf_area_clean{ll}{dir}{onoff}, [color(ll-1) 'o']);
            hold on
        end
    end
    set(gca, 'yscale', 'log')
    legend([h{1}(1), h{2}(1)], 'NDF 2', 'NDF 0')
    if onoff == 1
        title('ON')
    else
        title('OFF')
    end
    ylabel('RF area (mm^2)')
    set(gca, 'xtick', [])
end

%% superior RF size estimation
LL = {'NDF2', 'NDF0'};
dir = 1;
bin_n = 30;
threshold = 0.3;
PixelArea = (20*4)^2/10^6; %mm^2

% pulse filter
clear rf_wt rf_wt_area rf_wt_area_mean rf_wt_area_ste
for j = 5:5; %second
    tau = 0.01*j;
    tt = -3*tau:bin_size:3*tau;
    filter = exp(-tt.^2/(2*tau^2));
    npixel = 5;
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
%                     figure
%                     subplot(1,2,1)
%                     plot([bin_size:bin_size:1], temp_mean_norm)
%                     xlabel('second')
%                     ylabel('response')
                    temp_mean_norm = conv(temp_mean_norm, filter, 'same');
%                     subplot(1,2,2)
%                     plot([bin_size:bin_size:1], temp_mean_norm)
%                     xlabel('second')
%                     ylabel('response')

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
end

% step filter
for ll = 2:3
    for cc = 1:length(id_dir{dir})
        if ~fs_idx(idx_dir{dir}(cc), ll)
            for onoff = 1:2
                for p = 1:400
                    rf_step{ll}{cc}{onoff}(p) = sum(histg_all{ll}{idx_dir{dir}(cc)}{onoff}{p}(0.1/bin_size+1:0.4/bin_size));
                end
                rf_step{ll}{cc}{onoff} = reshape(rf_step{ll}{cc}{onoff}, field_width, field_height);
            end
        end
    end
end


for cc = 1:length(id_dir{dir})
    h = figure(1);
    set(h, 'position', [1 1 1980 1080])
    
    for ll = 1:2
        if ~isempty(rf_all{ll+1}{idx_dir{dir}(cc)})
            % original
            rf_on = rf_all{ll+1}{idx_dir{dir}(cc)}(:,:,1);
            subplot(4,8,(ll-1)*16+1)
            imagesc(rf_on)
            colormap gray
            if ll == 1
                title('spike #')
            end
            ylabel([LL{ll} ' ON'])
            subplot(4,8,(ll-1)*16+2)
            hist(rf_on(:), bin_n)
            xlim([0 max(rf_on(:))])
            rf_off = rf_all{ll+1}{idx_dir{dir}(cc)}(:,:,2);
            subplot(4,8,(ll-1)*16+9)
            imagesc(rf_off)
            colormap gray
            ylabel([LL{ll} ' OFF'])
            subplot(4,8,(ll-1)*16+10)
            hist(rf_off(:), bin_n)
            xlim([0 max(rf_off(:))])
            
            % pulse filter
            rf_on = rf_wt{ll+1}{dir}{cc}{1}';
            subplot(4,8,(ll-1)*16+4)
            imagesc(rf_on)
            colormap gray
            if ll == 1
                title('pulse filter')
            end
            subplot(4,8,(ll-1)*16+5)
            hist(rf_on(:), bin_n)
            xlim([0 max(rf_on(:))])
            rf_off = rf_wt{ll+1}{dir}{cc}{2}';
            subplot(4,8,(ll-1)*16+12)
            imagesc(rf_off)
            colormap gray
            subplot(4,8,(ll-1)*16+13)
            hist(rf_off(:), bin_n)
            xlim([0 max(rf_off(:))])
            
            % step filter
            rf_on = rf_step{ll+1}{cc}{1}';
            subplot(4,8,(ll-1)*16+7)
            imagesc(rf_on)
            colormap gray
            if ll == 1
                title('step filter')
            end
            subplot(4,8,(ll-1)*16+8)
            hist(rf_on(:), bin_n)
            xlim([0 max(rf_on(:))])
            rf_off = rf_step{ll+1}{cc}{2}';
            subplot(4,8,(ll-1)*16+15)
            imagesc(rf_off)
            colormap gray
            subplot(4,8,(ll-1)*16+16)
            hist(rf_off(:), bin_n)
            xlim([0 max(rf_off(:))])

            
            
        end
    end
    print_close(1,[24,12], num2str(id_dir{dir}(cc)))
end

% figure
% plot(meshgrid(1:2,1:length(id_dir{1}))', rf_wt_area{1}', 'ko-')
% hold on
% plot(meshgrid(3:4,1:length(id_dir{1}))', rf_wt_area{2}', 'ko-')
% xlim([0 5])
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

for dir = 2:2
    for cc = 1:length(id_dir{dir})
        figure(1)
        for ll = 1:2
            if ~isempty(rf_wt{ll+1}{dir}{cc})
                for onoff = 1:2
                    subplot(2,2,(ll-1)*2+onoff)
                    imagesc(rf_wt{ll+1}{dir}{cc}{onoff})
                    colormap gray
                end
            end
        end
        pause
        close(1)
    end
end
            