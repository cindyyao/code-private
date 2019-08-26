%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

datadg = load_data('/Volumes/lab/analysis/2017-03-31-0/data000/data000', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2017-03-31-0/stimuli/s00.txt';
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);
% load('DS170315.mat');

% identify DS cells
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);
params_idx = [2 3]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg,ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg);


datarun = load_data('/Volumes/lab/analysis/2017-03-31-0/data001-009-map/data001-009-map', opt);
time_points = [422 1372 1792 1892 2292 3242 5142 5542];

datamb_temp(1:9) = split_datarun(datarun, time_points);
datamb{1} = datamb_temp{2};
datamb{1}.names.stimulus_path = '/Volumes/lab/analysis/2017-03-31-0/stimuli/s02.txt';
datamb{1} = load_stim(datamb{1}, 'user_defined_trigger_set', [1:2:560]);
datamb{1}.stimulus.triggers = datamb{1}.stimulus.triggers';

datamb{2} = datamb_temp{6};
datamb{2}.names.stimulus_path = '/Volumes/lab/analysis/2017-03-31-0/stimuli/s06.txt';
datamb{2} = load_stim(datamb{2}, 'user_defined_trigger_set', [1:2:560]);
datamb{2}.stimulus.triggers = datamb{2}.stimulus.triggers';

datamb{3} = datamb_temp{9};
datamb{3}.names.stimulus_path = '/Volumes/lab/analysis/2017-03-31-0/stimuli/s09.txt';
datamb{3} = load_stim(datamb{3}, 'user_defined_trigger_set', [1:2:560]);
datamb{3}.stimulus.triggers = datamb{3}.stimulus.triggers';
load('DS170331.mat')

datawash{1} = datamb_temp{3};
datawash{1}.names.stimulus_path = '/Volumes/lab/analysis/2017-03-31-0/stimuli/s03.txt';
datawash{1} = load_stim(datawash{1}, 'user_defined_trigger_set', [1:2:length(datawash{1}.triggers)]);
datawash{1}.stimulus.triggers = datawash{1}.stimulus.triggers';


datawash{2} = datamb_temp{7};
datawash{2}.names.stimulus_path = '/Volumes/lab/analysis/2017-03-31-0/stimuli/s07.txt';
datawash{2} = load_stim(datawash{2}, 'user_defined_trigger_set', [1:2:length(datawash{2}.triggers)]);
datawash{2}.stimulus.triggers = datawash{2}.stimulus.triggers';

%% dg
% ds_id = datadg.cell_ids;
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

%% mb
n = 3;
for i = 1:n
    duration(i) = datamb{i}.triggers(end);
end
bin_size = 0.025; %sec
[raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx, raster_mb_all] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},ds_id,duration(i),bin_size);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb{i} = get_mb_raster(datamb{i}, ds_id, duration(i));
    raster_mb_all{i} = combine_repeats(raster_mb{i});
    for j = 1:length(raster_mb{i})
        if(mb_idx(j))
            raster_mb{i}{j} = [];
            raster_mb_all{i}{j} = [];
        end
    end
    trial_dur{i} = get_mb_trial_dur(datamb{i}, 400, 400, 0.5);
end

ctr_p = 7; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_mb = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}, raster_p_sum_all{i}] = get_pdirection_raster(raster_mb{i}, MB{1}.angle{ctr_p});
    [raster_n_sum{i}, n_idx{i}, raster_n_sum_all{i}] = get_ndirection_raster(raster_mb{i}, MB{1}.angle{ctr_p});
    MAG_all_norm_mb{i} = normalize_MAG(MB{i});
    rep = datamb{i}.stimulus.repetitions;
end

%%
d = 1;
t = 1;
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



%% mb wash
bin_size = 0.025; %sec
[raster_mb_wash, MB_wash, raster_mb_all_wash] = deal(cell(n, 1));
for i = 1:2    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datawash{i},ds_id,duration(i),bin_size);
    MB_wash{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datawash{i}));
    raster_mb_wash{i} = get_mb_raster(datawash{i}, ds_id, duration(i));
    raster_mb_all_wash{i} = combine_repeats(raster_mb_wash{i});
    for j = 1:length(raster_mb_wash{i})
        spikes_wash{i}{j} = [];
        if(mb_idx(j))
            raster_mb_wash{i}{j} = [];
            raster_mb_all_wash{i}{j} = [];
        else
            spikes_wash{i}{j} = squeeze(cellfun(@length, raster_mb_wash{i}{j}))';
        end
    end
end

for dir = 1:2
    CC = 1;
    for cc = 1:length(id_dir{dir})
        if ~mb_idx(idx_dir{dir}(cc))
            wash_response{dir}{1}(:, CC) = spikes_wash{1}{idx_dir{dir}(cc)}(1:20, 2);
            wash_mean{dir}(1,CC) = mean(wash_response{dir}{1}(:, CC));
            wash_response{dir}{2}(:, CC) = spikes_wash{1}{idx_dir{dir}(cc)}(end-19:end, 2);
            wash_mean{dir}(2,CC) = mean(wash_response{dir}{2}(:, CC));
            wash_response{dir}{3}(:, CC) = spikes_wash{2}{idx_dir{dir}(cc)}(1:20, 2);
            wash_mean{dir}(3,CC) = mean(wash_response{dir}{3}(:, CC));
            wash_response{dir}{4}(:, CC) = spikes_wash{2}{idx_dir{dir}(cc)}(end-19:end, 2);
            wash_mean{dir}(4,CC) = mean(wash_response{dir}{4}(:, CC));
            CC = CC + 1;
        end
    end
    for i = 1:4
        wash_response{dir}{i} = reshape(wash_response{dir}{i}, prod(size(wash_response{dir}{i})), 1);
    end
    wash_response{dir} = cell2mat(wash_response{dir});
end

ds_cell_type = {'superior', 'anterior'};
figure
for dir = 1:2
    subplot(1,2,dir)
    plot(wash_mean{dir})
    hold on
    for i = 1:size(wash_response{dir}, 1)
        scatter([1 2 3 4], wash_response{dir}(i, :))
    end
    
    title(ds_cell_type{dir})
    ylabel('spike#')
end

conditions = {'wash in start', 'wash in end', 'wash out start', 'wash out end'};
for dir = 1:2
    figure
    for i = 1:4
        subplot(4, 1, i)
        if dir == 1
            hist(wash_response{dir}(:, i), 0.5:1:30)
        else
            hist(wash_response{dir}(:, i), 0.5:1:10)
        end
        if i == 1
            title(ds_cell_type{dir})
        end
        ylabel(conditions{i})
        xlabel('# of spikes')
    end
end
    
%%
filter1 = ones(3, 1)/3;
filter2 = ones(3, 1)/3;

for dir = 1:2
    for cc = 1:length(id_dir{dir})
        figure(1)
        subplot(2, 1, 1)
        temp = conv2(spikes_wash{1}{idx_dir{dir}(cc)}, filter1, 'valid');
        plot(temp)
        hold on
        plot([24 24], [0 max(temp(:))], 'k')
        xlim([0 300])
        subplot(2, 1, 2)
        temp = conv2(spikes_wash{2}{idx_dir{dir}(cc)}, filter2, 'valid');
        plot(temp)
        hold on
        plot([19 19], [0 max(temp(:))], 'k')
        print_close(figure(1), [8 8], num2str(id_dir{dir}(cc)))
    end
end
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
for dir = 1:4
    id_dir_mb{dir} = id_dir{dir}(~mb_idx(idx_dir{dir}));
    idx_dir_mb{dir} = idx_dir{dir}(~mb_idx(idx_dir{dir}));
end

%% plot cell summary
for dir = 1:1
    for cc = 2:length(id_dir{dir})
        plot_mb_raster_ctr(MB, raster_mb, trial_dur, idx_dir{dir}(cc), id_dir{dir}(cc), '', 3, 7, 1)
    end
end
%% get spontaneous activity
for drug = 1:3
    for dir = 1:4
        for cc = 1:length(id_dir_mb{dir})
            idx = get_cell_indices(datamb{drug},id_dir_mb{dir}(cc));
            spikes_temp = datamb{drug}.spikes{idx};
            bgnd_firing{dir}(drug, cc) = length(spikes_temp(spikes_temp < 950 & spikes_temp > 900))/50;
        end
    end
end
%% max window
% on-off DSGC
clear Max_i
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'wash'};
step_size = 60;
for drug = 1:3
    for dir = [1 2 4]
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                for ctr = 7:-1:1
                    a = raster_p_sum{drug}{idx_dir{dir}(cc)}{ctr};
                    hist_temp = hist(a, xx);
                    if drug == 1 && ctr == 7
                        [max_p, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
                        Max_i{dir}(CC) = max_i;
                    else
                        max_p = conv(hist_temp, ones(1,step_size), 'valid');
                        max_p = max_p(Max_i{dir}(CC));
                    end
                    response_pmax{drug}{dir}(CC, ctr) = max_p/datamb{drug}.stimulus.repetitions - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
                end
                CC = CC + 1;
            end
        end
        response_pmax_norm{drug}{dir} = response_pmax{drug}{dir}./repmat(max(response_pmax{1}{dir}, [], 2), 1, size(response_pmax{drug}{dir},2));
%         response_pmax_norm{drug}{dir} = response_pmax{drug}{dir}./repmat(max(response_pmax{drug}{dir}, [], 2), 1, size(response_pmax{drug}{dir},2));
    end
    response_pmax_all{drug} = cell2mat(response_pmax{drug}(1)');
    response_pmax_norm_all{drug} = cell2mat(response_pmax_norm{drug}(1)');
%     response_pmax_all{drug} = cell2mat(response_pmax{drug}');
%     response_pmax_norm_all{drug} = cell2mat(response_pmax_norm{drug}');

end


ctr_x = [10 20 40 80 150 300 400];
color = 'brgkc';
figure
set(gcf, 'Position', [1 1 900 800])

for dir = [1 2 4]
    subplot(2,2,dir)
    for drug = 1:3
        errorbar(ctr_x, mean(exciseRows_empty(response_pmax_norm{drug}{dir}), 1), std(exciseRows_empty(response_pmax_norm{drug}{dir}), [], 1)/sqrt(size(exciseRows_empty(response_pmax_norm{drug}{dir}), 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
%     xlim([3 400])
end

