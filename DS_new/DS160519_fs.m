%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

datadg = load_data('/Volumes/lab/analysis/2016-05-19-0/data007-sorted/data007-sorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-05-19-0/stimuli/s07.mat';
datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);

datafs{1} = load_data('/Volumes/lab/analysis/2016-05-19-0/data000-001-map/data000-001-map', opt);
load('/Volumes/lab/analysis/2016-05-19-0/stimuli/s00.mat');
stim{1} = stim_out;
load('/Volumes/lab/analysis/2016-05-19-0/stimuli/s01.mat');
stim{2} = stim_out;
datafs{1}.stimulus = [];
datafs{1}.stimulus.repetitions = stim{1}.repeats;
datafs{1}.stimulus.trial_list = [stim{1}.trial_list(1:845) stim{2}.trial_list];
triggers_temp = [datafs{1}.triggers(2:2:1690); datafs{1}.triggers(2135:2:end)];
datafs{1}.stimulus.triggers = triggers_temp';
triggers_temp = [datafs{1}.triggers(1:1691); datafs{1}.triggers(2135:end)];
datafs{1}.triggers = triggers_temp';


datafs{2} = load_data('/Volumes/lab/analysis/2016-05-19-0/data002-004-map/data002-004-map', opt);
load('/Volumes/lab/analysis/2016-05-19-0/stimuli/s02.mat');
stim{1} = stim_out;
load('/Volumes/lab/analysis/2016-05-19-0/stimuli/s03.mat');
stim{2} = stim_out;
load('/Volumes/lab/analysis/2016-05-19-0/stimuli/s04.mat');
stim{3} = stim_out;

datafs{2}.stimulus = [];
datafs{2}.stimulus.repetitions = stim{1}.repeats;
datafs{2}.stimulus.trial_list = [stim{1}.trial_list(1:845) stim{2}.trial_list(1:676) stim{3}.trial_list];
triggers_temp = [datafs{2}.triggers(2:2:1690); datafs{2}.triggers(1747:2:3097); datafs{2}.triggers(3285:2:end)];
datafs{2}.stimulus.triggers = triggers_temp';
triggers_temp = [datafs{2}.triggers(1:1691); datafs{2}.triggers(1747:3098); datafs{2}.triggers(3285:end)];
datafs{2}.triggers = triggers_temp';

datafs{3} = load_data('/Volumes/lab/analysis/2016-05-19-0/data008-010-map/data008-010-map', opt);
load('/Volumes/lab/analysis/2016-05-19-0/stimuli/s08.mat');
stim{1} = stim_out;
load('/Volumes/lab/analysis/2016-05-19-0/stimuli/s09.mat');
stim{2} = stim_out;
load('/Volumes/lab/analysis/2016-05-19-0/stimuli/s10.mat');
stim{3} = stim_out;
datafs{3}.stimulus = [];
datafs{3}.stimulus.repetitions = stim{1}.repeats;
datafs{3}.stimulus.trial_list = [stim{1}.trial_list(1:676) stim{2}.trial_list(1:676) stim{3}.trial_list];
triggers_temp = [datafs{3}.triggers(2:2:1352); datafs{3}.triggers(1560:2:2910); datafs{3}.triggers(2944:2:end)];
datafs{3}.stimulus.triggers = triggers_temp';
triggers_temp = [datafs{3}.triggers(1:1353); datafs{3}.triggers(1560:2911); datafs{3}.triggers(2944:end)];
datafs{3}.triggers = triggers_temp';

%%

[NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [1 2]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

%% 
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg,ds_id,0,1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));
raster_dg{i} = get_ds_raster(datadg, ds_id);

delta_p = 2; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

[raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
MAG_all_norm_dg{i} = normalize_MAG(DG{i});
rep = datadg.stimulus.repetitions;

%% classify ds direction
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

%% plot individual cell summary
for cc = 1:length(ds_id)
    id = ds_id(cc);
    plot_ds_raster(DG, raster_dg, cc, id, '', 1, 1, 1)
end
%% map camera and display coordinates
cd /Volumes/Suk/Analysis/xyao/2016-05-19-0/images/
im_s = imread('WN.jpg'); % load stimulus picture taken by camera
im_array = imread('array.jpg'); % load array image taken by camera
cd /Users/xyao/matlab/code-private/DS_new/
stixel_size = 30; % frame shown in WN.jpg
movie_path = '/Volumes/lab/acquisition/movie-xml/BW-30-6-0.48-11111-20x20-60.35.xml';

mov = get_movie(movie_path, 0, 1);
mov_frame = matrix_scaled_up(squeeze(mov(:,:,1)), stixel_size);
clear movingPoints fixedPoints
cpselect(im_s, mov_frame) % select 4 control points

%% register two images
x_start = 105;
stixel_size_ = 15;
tform = fitgeotrans(movingPoints, fixedPoints, 'projective');
registered = imwarp(im_s, tform,'OutputView',imref2d(size(mov_frame)));
figure 
imshow(registered);
figure
imshowpair(mov_frame,registered,'blend');

% transform array image into display coordinates
registered_array = imwarp(im_array, tform, 'OutputView', imref2d(size(mov_frame)));
figure
imshow(registered_array);

% get array location in display coordinates

%              EI & Display                              
%
%               195 (1)                         
%                 / \                              
%               /     \                            
%   264 (6)    |       |    126 (2)               
%   386 (5)    |       |    4   (3)       
%               \     /                           
%                 \ /                              
%                455 (4)                          
array_location_display = ginput;

% get array location in ei coordinates
load('electrode_position.mat')
elec_corner = [195 126 4 455];
array_location_ei = position(elec_corner,:);
Tform = maketform('projective', array_location_ei, array_location_display);
test = tformfwd(Tform, array_location_ei)-array_location_display % should be equal or close to zeros

% distance = 1; % range used to calculate EI com
% center_ei = get_ei_com(datafs, ds_id_fs, distance);
% center_ei_display = (tformfwd(Tform, center_ei) - x_start)/stixel_size_;

%% plot rfs
load('DS160519.mat')
cell_type = {'superior', 'anterior', 'inferior', 'posterior'};
field_width = 13; field_height = 13;
field_width_sta = 40; field_height_sta = 40;
subregion = 1;

for i = 1:3
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
    figure
    set(gcf, 'Position', [1 1 1000 1000])
    id = ds_id(cc);
    for i = 1:3
        if ~isempty(rf_all{i}{cc})
            rf = padarray(rf_all{i}{cc},[7,7]);
    %     ei = center_ei_display(cc,:);
    %     ei = ei + 7;
    
            subplot(3,3,3*(i-1)+1)
            imagesc(sum(rf,3))
            colormap gray
        %     hold on 
        %     plot(ei(1),ei(2),'ro')
            axis image

            subplot(3,3,3*(i-1)+2)
            imagesc(rf(:,:,1))
            colormap gray
        %     hold on
        %     plot(ei(1),ei(2),'ro')
            title('on')
            axis image

            subplot(3,3,3*(i-1)+3)
            imagesc(rf(:,:,2))
            colormap gray
        %     hold on
        %     plot(ei(1),ei(2),'ro')
            title('off')
            axis image
        end
    end
    print_close(1,[12 12],num2str(id))
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

for ct = 1:4
    for ll = 1:3
        for cc = 1:length(id_dir_rf{ct})
            idx = find(ds_id == id_dir_rf{ct}(cc));
            for onoff = 1:2
                a = rf_all{ll}{idx}(:,:,onoff);
                sig_stixel = a>mean(a(:))+2*std(a(:));
                sig_stixel = reshape(sig_stixel(1:13, 1:13)', 169, 1);
                sig_stixel_all{ct}{cc}{ll}{onoff} = sig_stixel;
                sig_spike = fs_spike{ll}{idx}(sig_stixel, :, onoff);
                rf_size{ct}{ll}(cc, onoff) = sum(sum(a>mean(a(:))+2*std(a(:))))/4;
                sig_spike_all{ct}{cc}{ll}{onoff} = sig_spike;
            end
        end
        rf_size_mean(ct, ll, :) = mean(rf_size{ct}{ll});
        rf_size_ste(ct, ll, :) = std(rf_size{ct}{ll})/sqrt(length(id_dir_rf{ct}));
    end
end

% plot rf size
ct = {'superior', 'anterior', 'inferior', 'posterior'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = ct;
model_series = squeeze(rf_size_mean(:,:,1));
model_error = squeeze(rf_size_ste(:,:,1));
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('rf size (pixel #)')
legend('NDF 4','NDF 2', 'NDF 0');
title('ON')
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

% off
ct = {'superior', 'anterior', 'inferior', 'posterior'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = ct;
model_series = squeeze(rf_size_mean(:,:,2));
model_error = squeeze(rf_size_ste(:,:,2));
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('rf size (pixel #)')
legend('NDF 4','NDF 2', 'NDF 0');
title('OFF')
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

%% 
id = 3273;
idx = find(ds_id == id);
figure
for ll = 1:3
    for onoff = 1:2
        subplot(3, 2, (ll-1)*2+onoff)
        imagesc(rf_all{ll}{idx}(3:15, 7:19, onoff))
        colormap gray
        axis image
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
    end
end

%%

for ct = 1:4
    for ll = 1:3
        for cc = 1:length(id_dir_rf{ct})
            for onoff = 1:2
                temp = sig_spike_all{ct}{cc}{ll}{onoff};
                snr{ct}(cc, ll) = mean(mean(temp, 2)./std(temp, [], 2));
            end
        end
    end
    snr{ct} = nan2empty(snr{ct});
    snr_mean(ct, :) = mean(snr{ct});
    snr_ste(ct, :) = std(snr{ct})/sqrt(size(snr{ct}, 1));
end

ct = {'superior', 'anterior', 'inferior', 'posterior'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = ct;
model_series = snr_mean;
model_error = snr_ste;
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('SNR')
legend('NDF 4','NDF 2', 'NDF 0');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

%%
for ct = 1:4
    for onoff = 1:2
        for cc = 1:length(id_dir_rf{ct})
            sig_stixel_all_int{ct}{cc}{onoff} = sig_stixel_all{ct}{cc}{1}{onoff} & sig_stixel_all{ct}{cc}{2}{onoff} & sig_stixel_all{ct}{cc}{3}{onoff};
            if sum(sig_stixel_all_int{ct}{cc}{onoff}) ~= 0
                idx = find(ds_id == id_dir_rf{ct}(cc));
                for ll = 1:3
                    spike_temp = fs_spike{ll}{idx}(sig_stixel_all_int{ct}{cc}{onoff}, :, onoff);
                    spike_count{ct}{onoff}(cc, ll) = mean(spike_temp(:));
                end
            end
        end
        spike_count{ct}{onoff} = zero2empty(spike_count{ct}{onoff});
        spike_count_mean{onoff}(ct, :) = mean(spike_count{ct}{onoff});
        spike_count_ste{onoff}(ct, :) = std(spike_count{ct}{onoff})/sqrt(size(spike_count{ct}{onoff}, 1));
    end
end

ct = {'superior', 'anterior', 'inferior', 'posterior'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = ct;
model_series = spike_count_mean{1};
model_error = spike_count_ste{1};
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('RF gain (spike count)')
legend('NDF 4','NDF 2', 'NDF 0');
title('ON')
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

% off
ct = {'superior', 'anterior', 'inferior', 'posterior'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = ct;
model_series = spike_count_mean{2};
model_error = spike_count_ste{2};
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('RF gain (spike count)')
legend('NDF 4','NDF 2', 'NDF 0');
title('OFF')
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

%% added on 2016-11-30

% get center location
for cc = 1:length(ds_id)
    figure(1)
    id = ds_id(cc);
    i = 3;
    while isempty(rf_all{i}{cc})
        i = i - 1;
        if i == 0
            break
        end
    end
    if i>0
        rf = rf_all{i}{cc};

        imagesc(sum(rf,3))
        colormap gray
        axis image
        axis off
        center(cc, :) = ginput;
        close(1)
    end
end
center = floor(center);
center = max(center, ones(length(ds_id),2)*7);
center = min(center, ones(length(ds_id),2)*20);

% fit RF with Gaussian

for ll = 1:3
    rf_all_center{ll} = cell(1, length(ds_id));
    for cc = 1:length(ds_id)
        if ~isempty(rf_all{ll}{cc})
            rf_all_center{ll}{cc} = rf_all{ll}{cc}(center(cc, 2)-6:center(cc, 2)+6, center(cc, 1)-6:center(cc, 1)+6, :);
        end
    end
end
                
PixelArea = (15*4)^2/10^6;
% fit and compute rf area
for ll = 1:3
    for dir = 1:4
        clear rf_area_temp
        for onoff = 1:2
            rf_area_temp = [];
            for cc = 1:length(id_dir{dir})
                if ~isempty(rf_all_center{ll}{idx_dir{dir}(cc)})
                    data = rf_all_center{ll}{idx_dir{dir}(cc)}(:, :, onoff);
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
        for ll = 1:3
            n = size(rf_area_clean{ll}{dir}{onoff}, 1);
%             n = length(rf_area{ll}{dir}{onoff});
            h{ll} = plot((dir-1)*5+ll*ones(n,1), rf_area_clean{ll}{dir}{onoff}, [color(ll) 'o']);
%             h{ll} = plot((dir-1)*5+ll*ones(n,1), rf_area{ll}{dir}{onoff}, [color(ll) 'o']);
            hold on
        end
    end
%     set(gca, 'yscale', 'log')
    legend([h{1}(1), h{2}(1), h{3}(1)], 'NDF 4', 'NDF 2', 'NDF 0')
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
        errorbar([1 3 5], rf_area_clean_mean{onoff}(:, dir), rf_area_clean_ste{onoff}(:, dir))
        hold on
    end
%     ylim([0 0.05])
    xlabel('log(background intensity)')
    ylabel('RF area (mm^2)')
    legend(cell_type)
end
%% 

stdn = 2;
for ll = 1:3
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
    errorbar([1 3 5], rf_area_clean_all_mean{onoff}, rf_area_clean_all_ste{onoff})
    hold on
%     ylim([0 0.05])
    xlabel('log(background intensity)')
    ylabel('RF area (mm^2)')
    legend('ON', 'OFF')
end

%%
dir = 1;
onoff = 2;
y = [];group = [];
for ll = 1:5
    y = [y rf_area_clean{ll}{dir}{onoff}];
    group = [group ones(1, length(rf_area_clean{ll}{dir}{onoff}))*ll];
end

anova1(y, group)

%
onoff = 1;
y = [];group = [];
for ll = 2:3
    y = [y rf_area_clean_all{ll}{onoff}];
    group = [group ones(1, length(rf_area_clean_all{ll}{onoff}))*ll];
end

anova1(y, group)
