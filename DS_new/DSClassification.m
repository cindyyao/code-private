function [ds_id, id_dir, id_dir_on, idx_dir, idx_dir_on] = DSClassification(data_path, stim_path, varargin)

p = inputParser;

% specify list of optional parameters
p.addParameter('params_idx', [2 3]); % which parameters to use for classification
p.addParameter('delta_p', 3); %  which params to use to calculate prefer direction indices 
p.addParameter('classify_onoff', true);
p.addParameter('PCs', [1 2]);
p.addParameter('classify_onoff_direction', true);
p.addParameter('classify_on_direction', true);
p.addParameter('manual', false);


% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);
datadg = load_data(data_path, opt);
datadg.names.stimulus_path = stim_path;
if strcmp(stim_path(end-3:end), '.mat')
    datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);
else
    datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);
end

%% identify DSGCs
[NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg, datadg.cell_ids, 0, 1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb, datadg);
params_idx = params.params_idx; % which parameters to use for classification

[ds_id, nonds_id, id_init] = classify_ds(datadg, ds_struct, params_idx, 'manual', params.manual);


if params.classify_onoff
    %% classify ON-OFF vs ON DSGCs based on speed tunning
    [NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg,ds_id,0,1);
    DG = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));

    delta_p = params.delta_p; % choose which params to use to calculate prefer direction indices 
    MAG_all_norm_dg = normalize_MAG(DG);

    mag_pca = MAG_all_norm_dg;
    mag_pca = mag_pca';
    [id_sub, idx_sub] = deal(cell(2, 1));

    FigHandle = figure;
    set(FigHandle, 'Position', [1 1 380 400])

    [~,scores,~,~] = princomp(mag_pca);
    pc1 = params.PCs(1); pc2 = params.PCs(2);
    plot(scores(:, pc1), scores(:, pc2), 'o')
    hold on
    for i = 1:2
        [x, y] = ginput;
        plot(x, y)
        IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
        [~, idx_sub{i}] = find(IN' == 1);
        id_sub{i} = ds_id(idx_sub{i});
    end
    xlabel(['PC ' num2str(pc1)])
    ylabel(['PC ' num2str(pc2)])

    figure
    plot(scores(idx_sub{1}, pc1), scores(idx_sub{1}, pc2), 'ro', scores(idx_sub{2}, pc1), scores(idx_sub{2}, pc2), 'bo')
    xlabel(['PC ' num2str(pc1)])
    ylabel(['PC ' num2str(pc2)])

    v = 4*datadg.stimulus.params.SPATIAL_PERIOD./datadg.stimulus.params.TEMPORAL_PERIOD;
    figure
    semilogx(v, exciseColumn(MAG_all_norm_dg(:, idx_sub{1})), 'r')
    hold on
    semilogx(v, exciseColumn(MAG_all_norm_dg(:, idx_sub{2})), 'b')
    xlabel('micron/second')
    ylabel('Response')
    xlim([v(end) v(1)])

    t = delta_p;
    figure
    compass(DG.U{t}(idx_sub{1}), DG.V{t}(idx_sub{1}), 'r')
    hold on
    compass(DG.U{t}(idx_sub{2}), DG.V{t}(idx_sub{2}), 'b')
else
    idx_sub{1} = [];
    idx_sub{2} = 1:length(ds_id);
end

%% Classify DSGCs by direction

idx_dir = [];  id_dir = [];
if params.classify_onoff_direction
    % ON-OFF
    t = delta_p;
    h = figure;
    dirn = 4;
    set(h, 'Position', [1 1 1080 500])
    compass(DG.U{t}(idx_sub{2}), DG.V{t}(idx_sub{2}));
    color = 'brgkc';

    for i = 1:dirn
        hold on
        [x, y] = ginput;
        plot(x, y, color(i));

        IN = inpolygon(DG.U{t}(idx_sub{2}), DG.V{t}(idx_sub{2}), x, y);
        [~, I] = find(IN == 1);
        idx_dir{i} = idx_sub{2}(I);
        id_dir{i} = ds_id(idx_dir{i});
    end
end

% ON
idx_dir_on = [];  id_dir_on = [];
if params.classify_on_direction
    t = delta_p;
    h = figure;
    dirn = 3;
    set(h, 'Position', [1 1 1080 500])
    compass(DG.U{t}(idx_sub{1}), DG.V{t}(idx_sub{1}));
    color = 'brgkc';

    for i = 1:dirn
        hold on
        [x, y] = ginput;
        plot(x, y, color(i));

        IN = inpolygon(DG.U{t}(idx_sub{1}), DG.V{t}(idx_sub{1}), x, y);
        [~, I] = find(IN == 1);
        idx_dir_on{i} = idx_sub{1}(I);
        id_dir_on{i} = ds_id(idx_dir_on{i});
    end
end
