%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load drifting grating data
datadg = load_data('/Volumes/lab/analysis/2016-10-19-0/data008/data008', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-10-19-0/stimuli/s08.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

% identify DS cells
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [1 2]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);


datafs{1} = load_data('/Volumes/lab/analysis/2016-10-19-0/data000-map/data000-map', opt);
datafs{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-10-19-0/stimuli/s00.mat';
datafs{1} = load_mfs_161019(datafs{1});
datafs{2} = load_data('/Volumes/lab/analysis/2016-10-19-0/data001-map/data001-map', opt);
datafs{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-10-19-0/stimuli/s01.mat';
datafs{2} = load_mfs_161019(datafs{2});
datafs{3} = load_data('/Volumes/lab/analysis/2016-10-19-0/data003-map/data003-map', opt);
datafs{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-10-19-0/stimuli/s03.mat';
datafs{3} = load_mfs_161019(datafs{3});
datafs{4} = load_data('/Volumes/lab/analysis/2016-10-19-0/data007-map/data007-map', opt);
datafs{4}.names.stimulus_path = '/Volumes/lab/analysis/2016-10-19-0/stimuli/s07.mat';
datafs{4} = load_stim_mfs(datafs{4});



%%
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);

delta_p = 1; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(DG{i});
rep = datadg.stimulus.repetitions;


d = 1;
t = 2;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(DG{d}.U{t}, DG{d}.V{t});
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}, DG{d}.V{t}, x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = I;
    id_dir{i} = ds_id(idx_dir{i});
end
%%
cell_type = {'superior', 'anterior', 'inferior', 'posterior'};
field_width = 15; field_height = 15;
% field_width_sta = 30; field_height_sta = 30;
subregion = 0;
stop = 0.5; %second
for i = 1:4
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

for cc = 1:length(ds_id)
    set(gcf, 'Position', [1 1 1000 500])
    id = ds_id(cc);
    for i = 1:4
        if ~isempty(rf_all{i}{cc})
            rf = rf_all{i}{cc};
            
            subplot(2,4,i)
            imagesc(rf(:,:,1))
            colormap gray
            title('on')
            axis image

            subplot(2,4,4+i)
            imagesc(rf(:,:,2))
            colormap gray
            title('off')
            axis image
        end
    end
    pause
    print_close(1,[12 8],num2str(id))
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
for ll = 1:4
    trigger = datafs{ll}.triggers(1:2:end);
    list = datafs{ll}.stimulus.trial_list;
    repeat = datafs{ll}.stimulus.repetitions;
    clear index
    for i = 1:max(list)
        index(i, :) = find(list == i);
    end
    raster_onoff{ll} = cell(length(ds_id), 1);
    for i = 1:length(ds_id)
        if ~fs_idx(i,ll)
            idx = get_cell_indices(datafs{ll}, ds_id(i));
            raster = get_raster(datafs{ll}.spikes{idx}, trigger, 'plot', false);
            raster_onoff{ll}{i} = raster(index);
        end
    end
end

%%
LL = {'NDF4', 'NDF3', 'NDF2', 'NDF0'};
bin_size = 0.05;
xx = [bin_size/2:bin_size:2-bin_size/2];
XX = [bin_size/2:bin_size:1-bin_size/2];
for ll = 1:4
    for cc = 1:length(raster_onoff{ll})
        if ~isempty(raster_onoff{ll}{cc})
            for p = 1:225%size(raster_onoff{ll}{1}, 1)
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

h = subplot(15, 15, 1); p1 = get(h, 'pos');
h = subplot(15, 15, 2); p2 = get(h, 'pos');
width = (p2(1) - p1(1))*0.8;

h = subplot(15, 15, 1); p1 = get(h, 'pos');
h = subplot(15, 15, 16); p2 = get(h, 'pos');
height = p1(2) - p2(2);
for cc = 1:length(ds_id)
    for ll = 1:4
        H = figure(1);
        set(H, 'Position', [1 1 1080 1080])
        if ~fs_idx(cc,ll)
            for p = 1:225
                h = subplot(15, 15, p);
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

%% fit RF with Gaussian                
PixelArea = (40*4)^2/10^6;
% fit and compute rf area
for ll = 1:4
    for dir = 1:4
        clear rf_area_temp
        for onoff = 1:2
            rf_area_temp = [];
            for cc = 1:length(id_dir{dir})
                if ~fs_idx(idx_dir{dir}(cc), ll)
                    data = rf_all{ll}{idx_dir{dir}(cc)}(:, :, onoff);
%                     data = rf_wt{ll}{dir}{cc}{onoff};
                    if sum(sum(data > mean(data(:))+3*std(data(:))))>0
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
            rf_area{ll}{dir}{onoff} = rf_area_temp;
        end
    end
end


% exclude outliers
stdn = 2;
for ll = 1:4
    for dir = 1:4
        for onoff = 1:2
            notdone = 1;
            rf_area_temp = rf_area{ll}{dir}{onoff};
            while notdone
                a = length(rf_area_temp);
                rf_area_temp(rf_area_temp > std(rf_area_temp)*stdn + mean(rf_area_temp)) = [];
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
        for ll = 1:4
            n = size(rf_area_clean{ll}{dir}{onoff}, 1);
            h{ll} = plot((dir-1)*5+ll*ones(n,1), rf_area_clean{ll}{dir}{onoff}, [color(ll) 'o']);
%             n = length(rf_area{ll}{dir}{onoff});
%             h{ll} = plot((dir-1)*5+ll*ones(n,1), rf_area{ll}{dir}{onoff}, [color(ll) 'o']);
            hold on
        end
    end
%     set(gca, 'yscale', 'log')
    legend([h{1}(1), h{2}(1), h{3}(1), h{4}(1)], 'NDF 4', 'NDF 3', 'NDF 2', 'NDF 0')
    if onoff == 1
        title('ON')
    else
        title('OFF')
    end
    ylabel('RF area (mm^2)')
    set(gca, 'xtick', [])
%         ylim([0 0.5])

end

figure
for onoff = 1:2
    subplot(1,2,onoff)
    for dir = 1:4
        errorbar([1 2 3 5], rf_area_clean_mean{onoff}(:, dir), rf_area_clean_ste{onoff}(:, dir), 'color', color(dir))
        hold on
    end
    ylim([0 0.3])
    xlabel('log(background intensity)')
    ylabel('RF area (mm^2)')
    legend(cell_type)

end

%% 

stdn = 2;
for ll = 1:4
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

% for ll = 1:5
%     for onoff = 1:2
%         temp = [];
%         for dir = 1:3
%             temp = [temp rf_area_clean{ll}{dir}{onoff}];
%         end
%         rf_area_clean_all_mean{onoff}(ll) = mean(temp);
%         rf_area_clean_all_ste{onoff}(ll) = std(temp)/sqrt(length(temp));
%         rf_area_clean_all{onoff}{ll} = temp;
%     end
% end

figure
for onoff = 1:2
    errorbar([1 2 3 5], rf_area_clean_all_mean{onoff}, rf_area_clean_all_ste{onoff})
    hold on
%     ylim([0 0.05])
    xlabel('log(background intensity)')
    ylabel('RF area (mm^2)')
    legend('ON', 'OFF')
end
