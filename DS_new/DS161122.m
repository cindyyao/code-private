%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);
path = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2016-11-22-0/';
% load drifting grating data
datadg = load_data(strcat(path, 'data007-sorted/data007-sorted'), opt);
datadg.names.stimulus_path = strcat(path, 'stimuli/s07.txt');
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);
datadg = load_ei(datadg, ds_id, 'array_type', 519);
% identify DS cells
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [2 3]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);

datafs{1} = load_data(strcat(path, 'data000-map/data000-map'), opt);
datafs{1}.names.stimulus_path = strcat(path, 'stimuli/s00.mat');
datafs{1} = load_stim_mfs(datafs{1});
datafs{2} = load_data(strcat(path, 'data001-map/data001-map'), opt);
datafs{2}.names.stimulus_path = strcat(path, 'stimuli/s01.mat');
datafs{2} = load_stim_mfs(datafs{2});
datafs{3} = load_data(strcat(path, 'data002-map/data002-map'), opt);
datafs{3}.names.stimulus_path = strcat(path, 'stimuli/s02.mat');
datafs{3} = load_stim_mfs(datafs{3});
datafs{4} = load_data(strcat(path, 'data003-map/data003-map'), opt);
datafs{4}.names.stimulus_path = strcat(path, 'stimuli/s03.mat');
datafs{4} = load_stim_mfs(datafs{4});
datafs{5} = load_data(strcat(path, 'data006-map/data006-map'), opt);
datafs{5}.names.stimulus_path = strcat(path, 'stimuli/s06.mat');
datafs{5} = load_stim_mfs(datafs{5});

load('DS161122.mat')
%% dg
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);


delta_p = 3; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(DG{i});
rep = datadg.stimulus.repetitions;

%%
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
t = 3;
figure
compass(DG{1}.U{t}(idx_sub{1}), DG{1}.V{t}(idx_sub{1}), 'r')
hold on
compass(DG{1}.U{t}(idx_sub{2}), DG{1}.V{t}(idx_sub{2}), 'b')

%% plot rfs
cell_type = {'superior', 'anterior', 'inferior', 'posterior'};
field_width = 17; field_height = 17;
field_width_sta = 40; field_height_sta = 40;
subregion = 0;
stop = 0.5; %second
for i = 1:5
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
%%
for cc = 1:length(ds_id)
    figure(1)
    set(gcf, 'Position', [1 1 1000 500])
    id = ds_id(cc);
    for i = 1:5
        if ~isempty(rf_all{i}{cc})
%             rf = padarray(rf_all{i}{cc},[7,7]);
            rf = rf_all{i}{cc};
    
            subplot(3,5,i)
            imagesc(sum(rf,3))
            colormap gray
            axis image
            axis off

            subplot(3,5,5+i)
            imagesc(rf(:,:,1))
            colormap gray
            
            data = rf(:, :, 1);
            if sum(sum(data > mean(data(:))+5*std(data(:)))) == 0
                title('excluded')
            end
            
            axis image
            axis off

            subplot(3,5,10+i)
            imagesc(rf(:,:,2))
            colormap gray
            
            data = rf(:, :, 2);
            if sum(sum(data > mean(data(:))+5*std(data(:)))) == 0
                title('excluded')
            end
            
            axis image
            axis off

        end
    end
%     print_close(1,[24 12],num2str(id))
    pause
    close(1)
end

%
% fs_spike_temp = cellfun(@(fs_spike) sum(fs_spike,2), fs_spike, 'UniformOutput', false);
% fs_spike_temp = cellfun(@(fs_spike_temp) sum(fs_spike_temp,3), fs_spike_temp, 'UniformOutput', false);
% [~, I] = cellfun(@(fs_spike_temp) sort(fs_spike_temp, 'descend'), fs_spike_temp, 'UniformOutput', false);trigger = datafs.triggers(2:2:end);
for ll = 1:5
    trigger = datafs{ll}.triggers(1:2:end);
    list = datafs{ll}.stimulus.trial_list;
    repeat = datafs{ll}.stimulus.repetitions;
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

%
LL = {'NDF4', 'NDF3', 'NDF2', 'NDF1', 'NDF0'};
bin_size = 0.02;
xx = [bin_size/2:bin_size:2-bin_size/2];
XX = [bin_size/2:bin_size:1-bin_size/2];
for ll = 1:5
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

%% gaussian filter
PixelArea = (30*4)^2/10^6;
threshold = 0.3;

clear rf_wt rf_wt_area rf_wt_area_mean rf_wt_area_ste rf_area rf_area_clean
tau = 0.05;
tt = -3*tau:bin_size:3*tau;
filter = exp(-tt.^2/(2*tau^2));

filter = filter/norm(filter);
npixel = 5;
for dir = 1:4

for ll = 1:5
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
                
%                 figure
%                 subplot(1, 2, 1); plot(temp_mean_norm)
                
                temp_mean_norm = conv(temp_mean_norm, filter, 'same');

%                 subplot(1, 2, 2); plot(temp_mean_norm)
%                 pause
                
                for p = 1:289
                    rf_wt{ll}{dir}{cc}{onoff}(p) = histg_all{ll}{idx_dir{dir}(cc)}{onoff}{p}*temp_mean_norm';
                end
                rf_wt{ll}{dir}{cc}{onoff} = reshape(rf_wt{ll}{dir}{cc}{onoff}, field_width, field_height);
                rf_wt_area{dir}{ll}{onoff}{cc} = sum(sum(rf_wt{ll}{dir}{cc}{onoff} > max(rf_wt{ll}{dir}{cc}{onoff}(:))*threshold))*PixelArea;
            end
        end
    end
    for onoff = 1:2
        rf_wt_area{dir}{ll}{onoff} = cell2mat(rf_wt_area{dir}{ll}{onoff});
        rf_wt_area_mean(dir, ll, onoff) = mean(rf_wt_area{dir}{ll}{onoff});
        rf_wt_area_ste(dir, ll, onoff) = std(rf_wt_area{dir}{ll}{onoff})/sqrt(length(rf_wt_area{dir}{ll}{onoff}));
    end
end
end


%% fit RF with Gaussian                
PixelArea = (30*4)^2/10^6;
% fit and compute rf area
for ll = 1:5
    for dir = 1:4
        clear rf_area_temp
        for onoff = 1:3
            rf_area_temp = [];
%             rf_center_temp = [];
            for cc = 1:length(id_dir{dir})
                if ~fs_idx(idx_dir{dir}(cc), ll)
                    if onoff < 3
                        data = rf_all{ll}{idx_dir{dir}(cc)}(:, :, onoff);
                    else
                        data = sum(rf_all{ll}{idx_dir{dir}(cc)}, 3);
                    end
                    if sum(sum(data > mean(data(:))+5*std(data(:))))>0
%                         figure(100)
%                         imagesc(data)
%                         colormap gray
%                         pause
                        
                        params = fit_2d_gaussian(data);
    %                     Gaussian_params{ll}{dir}{cc}{onoff} = params;
                        rf_area_temp = [rf_area_temp params.xd * params.yd * pi * PixelArea];
                        params_all{ll}{dir}{onoff}{cc} = params;
                    end
                end
            end
            
            rf_area{ll}{dir}{onoff} = rf_area_temp;
        end
    end
end


% exclude outliers
stdn = 3;
for ll = 1:5
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
        for ll = 1:5
            n = size(rf_area_clean{ll}{dir}{onoff}, 1);
            h{ll} = plot((dir-1)*5+ll*ones(n,1), rf_area_clean{ll}{dir}{onoff}, [color(ll) 'o']);
%             n = length(rf_area{ll}{dir}{onoff});
%             h{ll} = plot((dir-1)*5+ll*ones(n,1), rf_area{ll}{dir}{onoff}, [color(ll) 'o']);
            hold on
        end
    end
%     set(gca, 'yscale', 'log')
    legend([h{1}(1), h{2}(1), h{3}(1), h{4}(1), h{5}(1)], 'NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0')
    if onoff == 1
        title('ON')
    else
        title('OFF')
    end
    ylabel('RF area (mm^2)')
    set(gca, 'xtick', [])
%         ylim([0 0.3])

end

figure
for onoff = 1:2
    subplot(1,2,onoff)
    for dir = 1:4
        errorbar(0:4, rf_area_clean_mean{onoff}(:, dir), rf_area_clean_ste{onoff}(:, dir), 'color', color(dir))
        hold on
    end
    ylim([0 0.14])
    xlim([-1 5])
    xlabel('log(background intensity)')
    ylabel('RF area (mm^2)')
    legend(cell_type)

end

ONOFF = {'ON', 'OFF'};
figure
for onoff = 1:2
    dir = 1;
    errorbar(0:4, rf_area_clean_mean{onoff}(:, dir), rf_area_clean_ste{onoff}(:, dir), 'color', color(onoff))
    hold on
    ylim([0 0.14])
    xlim([-1 5])
    xlabel('log(background intensity)')
    ylabel('RF area (mm^2)')
    legend(ONOFF)

end

%% mosaic
ll = 5;
PIXSIZE = 30;

color = 'brk';
figure
for ct = 1:4
    subplot(2, 2, ct)
    for cc = 1:length(params_all{ll}{ct}{onoff})
        for onoff = 3:3
            if ~isempty(params_all{ll}{ct}{onoff}{cc})
                params = params_all{ll}{ct}{onoff}{cc};
                drawEllipse(params.x0*PIXSIZE, params.y0*PIXSIZE, params.xd*PIXSIZE, params.yd*PIXSIZE, params.angle, 'color', color(onoff))
                hold on
            end
            if ct == 2 && cc == 23
                drawEllipse(params.x0*PIXSIZE, params.y0*PIXSIZE, params.xd*PIXSIZE, params.yd*PIXSIZE, params.angle, 'color', 'r')
            end
        end
    end
    xlim([100, 500])
    ylim([100, 500])
    daspect([1 1 1])
%     legend('ON', 'OFF', 'ALL')
end
%%
dir = 4;
onoff = 1;
y = [];group = [];
for ll = 3:5
    y = [y rf_area_clean{ll}{dir}{onoff}];
    group = [group ones(1, length(rf_area_clean{ll}{dir}{onoff}))*ll];
end

anova1(y, group)
%%
stdn = 2;
for ll = 1:5
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
    errorbar(0:4, rf_area_clean_all_mean{onoff}, rf_area_clean_all_ste{onoff})
    hold on
%     ylim([0 0.05])
    xlim([-1 5])
    xlabel('log(background intensity)')
    ylabel('RF area (mm^2)')
    legend('ON', 'OFF')
end

%%
onoff = 1;
y = [];group = [];
for ll = 3:5
    y = [y rf_area_clean_all{ll}{onoff}];
    group = [group ones(1, length(rf_area_clean_all{ll}{onoff}))*ll];
end

anova1(y, group)

%% neighboring pairs
pos = datadg.ei.position;
mode = 'neg';
neighbors = [];
ct = 1;
for cc1 = 1:length(id_dir{ct})
    for cc2 = cc1+1:length(id_dir{ct})
        id1 = id_dir{ct}(cc1);
        idx1 = get_cell_indices(datadg, id1);
        ei1 = datadg.ei.eis{idx1};
        com1 = ei_com_xy(ei1, pos, 30*3, mode);
        id2 = id_dir{ct}(cc2);
        idx2 = get_cell_indices(datadg, id2);
        ei2 = datadg.ei.eis{idx2};
        com2 = ei_com_xy(ei2, pos, 30*3, mode);
        if pdist([com1;com2]) < 150
            neighbors = [neighbors; id1 id2];
        end
    end
end

ct = 4;
celltype = {'superior', 'anterior', 'inferior', 'posterior'};
coms = [];
for cc = 1:length(id_dir{ct})
    id = id_dir{ct}(cc);
    idx = get_cell_indices(datadg, id);
    ei = datadg.ei.eis{idx};
    com = ei_com_xy(ei, pos, 30*3, mode);
    coms = [coms; com];
end

corner_i = [4 126 195 264 386 455 4];
corner_position = datadg.ei.position(corner_i, :);
figure
for cc = 1:length(id_dir{ct})
    plot(coms(cc, 1), coms(cc, 2),'ko')
    hold on
%     text(coms(cc, 1)+5, coms(cc, 2)+5, num2str(id_dir{1}(cc)), 'FontSize', 10)
    
end

%%
duration = 3000;
bin_size = 0.0005;
max_lag = 20;
ct = 1;
N = 10000;
ll = {'NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0'};
xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
peak = zeros(size(neighbors, 1), 5);
for cp = 1:size(neighbors, 1)
    if mod(cp, 5) == 1
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 2000 2000])
        fig_i = 1;
    end
    id1 = neighbors(cp, 1);
    id2 = neighbors(cp, 2);
    for i = 1:5
        if ~fs_idx(find(ds_id == id1), i)
            idx1 = get_cell_indices(datafs{i}, id1);
            idx2 = get_cell_indices(datafs{i}, id2);
            spikes1 = datafs{i}.spikes{idx1};
            spikes1_TF= ceil(spikes1/bin_size);
            spikes1 = zeros(duration/bin_size, 1);
            spikes1(spikes1_TF) = 1;

            spikes2 = datafs{i}.spikes{idx2};
            spikes2_TF= ceil(spikes2/bin_size);
            spikes2 = zeros(duration/bin_size, 1);
            spikes2(spikes2_TF) = 1;

            A = xcorr(spikes1, spikes2, max_lag, 'coeff');
            peak(cp, i) = max(A);
% %     A = xcorr(spikes1, spikes2, max_lag);
%     [h, filteredA] = find_smallest_h(A);
%     [bootstat,bootsam] = bootstrp(N,@find_smallest_h_hist,rude(round(filteredA), 1:max_lag*2+1), max_lag);
%     p = sum(bootstat > h)/N;
            subplot(5, 5, (fig_i-1)*5+i)
%     if p < 0.01
%         bar(xx, A, 'r')
%     else
            bar(xx, A, 'b')
%     end
            if fig_i == 1
                title(ll{i})
            end
            xlim([-0.01 0.01])
        end
    end
    fig_i = fig_i + 1;
end

figure
for i = 1:size(peak, 1)
    plot([1 10 100 1000 10000], peak(i, :))
    hold on
end
set(gca, 'xscale', 'log')
xlabel('R*/rod/s')
ylabel('correlation coefficient')

%% compare responses to center flash between on-off and on DSGC
ll = 1;
CC = 1;

for ct = 1:3
    for cc = 1:length(id_dir_on{ct})
        [~, argmax] = max(cellfun(@max, histg_onoff{ll}{idx_dir_on{ct}(cc)}));
        hist_on(CC, :) = histg_onoff{ll}{idx_dir_on{ct}(cc)}{argmax};
        CC = CC+1;
    end
end

CC = 1;

for ct = 1:4
    for cc = 1:length(id_dir{ct})
        if ~isempty(histg_onoff{ll}{idx_dir{ct}(cc)})
            [~, argmax] = max(cellfun(@max, histg_onoff{ll}{idx_dir{ct}(cc)}));
            hist_oo(CC, :) = histg_onoff{ll}{idx_dir{ct}(cc)}{argmax};
            CC = CC+1;
        end
    end
end

repeat = 5;
hist_on = hist_on/repeat/bin_size;
hist_oo = hist_oo/repeat/bin_size;

xx = [bin_size/2:bin_size:2-bin_size/2];
figure
shadedplot(xx, mean(hist_oo)+std(hist_oo)/sqrt(size(hist_oo, 1)), mean(hist_oo)-std(hist_oo)/sqrt(size(hist_oo, 1)), 'k', 'k');
hold on
plot(xx, mean(hist_oo), 'k')
shadedplot(xx, mean(hist_on)+std(hist_on)/sqrt(size(hist_on, 1)), mean(hist_on)-std(hist_on)/sqrt(size(hist_on, 1)), 'r', 'r');
hold on
plot(xx, mean(hist_on), 'r')

xlabel('time (second)')
ylabel('firing rate (Hz)')
legend('ON-OFF', 'ON')


xx = [bin_size/2:bin_size:2-bin_size/2];
figure
shadedplot(xx, mean(hist_oo)+std(hist_oo), mean(hist_oo)-std(hist_oo), 'k', 'k');
hold on
plot(xx, mean(hist_oo), 'k')
shadedplot(xx, mean(hist_on)+std(hist_on), mean(hist_on)-std(hist_on), 'r', 'r');
hold on
plot(xx, mean(hist_on), 'r')

xlabel('time (second)')
ylabel('firing rate (Hz)')
legend('ON-OFF', 'ON')


%% compare responses to center flash between superior on-off and others
repeat = 5;
LL = {'NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0'};
clear hist_oo_ct

c = [bin_size/2:bin_size:2-bin_size/2];
figure
for ll = 1:5
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
    subplot(5,1,ll)
    for ct = 1:4
        plot(xx, mean(hist_oo_ct{ll, ct}))
        hold on
    end
    title(LL{ll})
    ylabel('firing rate (Hz)')
    if ll == 1
        legend('superior', 'anterior', 'inferior', 'posterior')
    end
    if ll == 5
        xlabel('time (s)')
    end
end

figure
for ll = 1:5
    hist_oo_ct_combine{ll, 1} = hist_oo_ct{ll, 1};
    hist_oo_ct_combine{ll, 2} = cell2mat(hist_oo_ct(ll, 2:4)');
    subplot(5,1,ll)
    for ct = 1:2
        plot(xx, mean(hist_oo_ct_combine{ll, ct}))
        hold on
    end
    title(LL{ll})
    ylabel('firing rate (Hz)')
    if ll == 1
        legend('superior', 'others')
    end
    if ll == 5
        xlabel('time (s)')
    end
end


% mean spike number
center_gain = cellfun(@sum, cellfun(@mean, hist_oo_ct(1, :), 'UniformOutput', 0));
%% fit RF with Gaussian                
PixelArea = (30*4)^2/10^6;
% fit and compute rf area
for ll = 1:5
    for dir = 1:4
        rf_area_temp = [];
        for cc = 1:length(id_dir{dir})
            if ~fs_idx(idx_dir{dir}(cc), ll)
                data = sum(rf_all{ll}{idx_dir{dir}(cc)}, 3);
                if sum(sum(data > mean(data(:))+5*std(data(:))))>0
                    params = fit_2d_gaussian(data);
                    rf_area_temp = [rf_area_temp params.xd * params.yd * pi * PixelArea];
                end
            end
        end
        rf_area_all{ll}{dir} = rf_area_temp;
    end
end


% exclude outliers
stdn = 2;
for ll = 1:5
    for dir = 1:4
        notdone = 1;
        rf_area_temp = rf_area_all{ll}{dir};
        while notdone
            a = length(rf_area_temp);
            rf_area_temp(rf_area_temp > std(rf_area_temp)*stdn + mean(rf_area_temp)) = [];
            b = length(rf_area_temp);
            if a == b
                notdone = 0;
                rf_area_all_clean{ll}{dir} = rf_area_temp;
            end
        end
        rf_area_all_clean_mean(ll, dir) = mean(rf_area_all_clean{ll}{dir});
        rf_area_all_clean_ste(ll, dir) = std(rf_area_all_clean{ll}{dir})/sqrt(length(rf_area_all_clean{ll}{dir}));

    end
end

area_type = rf_area_all_clean_mean(1, :);
load('DS160324.mat', 'xthreshold')