cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

datadg = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-01-11-0/data007/data007', opt);
datadg.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2018-01-11-0/stimuli/s07.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

params_idx = [3 2]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);


datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2018-01-11-0/datamfs-map/datamfs-map', opt);
time_points = [2566 5059];
datamfs(1:3) = split_datarun(datarun, time_points);
datamfs{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2018-01-11-0/stimuli/s01.mat';
datamfs{1} = load_stim_mfs(datamfs{1});

datamfs{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2018-01-11-0/stimuli/s03.mat';
datamfs{2} = load_stim_mfs(datamfs{2});

datamfs{3}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2018-01-11-0/stimuli/s05.mat';
datamfs{3} = load_stim_mfs(datamfs{3});
load('DS180111.mat')
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

t = 2;
figure
compass(dg{1}.U{t}(idx_sub{1}), dg{1}.V{t}(idx_sub{1}), 'r')
hold on
compass(dg{1}.U{t}(idx_sub{2}), dg{1}.V{t}(idx_sub{2}), 'b')

%%
d = 1;
t = 1;
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
for i = 1:3
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
    set(gcf, 'Position', [1 1 1000 1000])
    id = ds_id(cc);
    for i = 1:3
        if ~isempty(rf_all{i}{cc})
%             rf = padarray(rf_all{i}{cc},[7,7]);
            rf = rf_all{i}{cc};
    
            subplot(3,3,3*(i-1)+1)
            imagesc(sum(rf,3))
            colormap gray
            axis image
            axis off

            subplot(3,3,3*(i-1)+2)
            imagesc(rf(:,:,1))
            colormap gray
            axis image
            axis off

            subplot(3,3,3*i)
            imagesc(rf(:,:,2))
            colormap gray
            axis image
            axis off

        end
    end
    print_close(1,[24 12],num2str(id))
%     pause
%     close(1)
end

%% fit RF with Gaussian                
PixelArea = (30*4)^2/10^6;
% fit and compute rf area
for ds = 1:3
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
stdn = 2;
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
%                     temp(temp < mean(temp) + 1.5*std(temp)) = 0;
                    
                    rf_spikes{ds}{ct}(CC, onoff) = sum(temp);
                end
                CC = CC + 1;
            end
        end
    end
end
