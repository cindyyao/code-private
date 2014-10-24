%% 2013-02-21-0

opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun4 = load_data('/Analysis/xyao/2013-02-21-0/data004/data004', opt);
datarun2 = load_data('/Analysis/xyao/2013-02-21-0/data002/data002', opt);
datarun0 = load_data('/Analysis/xyao/2013-02-21-0/data000/data000', opt);

cell_type = {'ON brisk transient', 'ON transient', 'OFF brisk transient', 'OFF sustained','OFF transient'};

%% mapped from ndf 4
[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map3(datarun4, ...
    datarun2, datarun0, cell_type);

% to ndf 0
cell_type_n = 4;
original = get_cell_ids(datarun4, cell_type{cell_type_n});
cell_spec = cell_id{cell_type_n}(:, 1);
overlap = intersect(cell_spec, original);

plot_rf_summaries(datarun4, cell_spec, 'fit_color', 'b', 'foa', 0);
plot_rf_summaries(datarun4, original, 'fit_color', 'r', 'clear', false)
plot_rf_summaries(datarun4, overlap, 'fit_color', 'c', 'clear', false)

% to ndf 2
cell_type_n = 4;
original = get_cell_ids(datarun2, cell_type{cell_type_n});
cell_spec = cell_id{cell_type_n}(:, 2);
overlap = intersect(cell_spec, original);

plot_rf_summaries(datarun2, cell_spec, 'fit_color', 'b', 'foa', 0);
plot_rf_summaries(datarun2, original, 'fit_color', 'r', 'clear', false)
plot_rf_summaries(datarun2, overlap, 'fit_color', 'c', 'clear', false)

%% mapped from ndf 0
[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map3(datarun0, ...
    datarun2, datarun4, cell_type);

% to ndf 4
cell_type_n = 4;
original = get_cell_ids(datarun0, cell_type{cell_type_n});
cell_spec = cell_id{cell_type_n}(:, 1);
overlap = intersect(cell_spec, original);

plot_rf_summaries(datarun0, cell_spec, 'fit_color', 'b', 'foa', 0);
plot_rf_summaries(datarun0, original, 'fit_color', 'r', 'clear', false)
plot_rf_summaries(datarun0, overlap, 'fit_color', 'c', 'clear', false)

% to ndf 2
cell_type_n = 4;
original = get_cell_ids(datarun2, cell_type{cell_type_n});
cell_spec = cell_id{cell_type_n}(:, 2);
overlap = intersect(cell_spec, original);

plot_rf_summaries(datarun2, cell_spec, 'fit_color', 'b', 'foa', 0);
plot_rf_summaries(datarun2, original, 'fit_color', 'r', 'clear', false)
plot_rf_summaries(datarun2, overlap, 'fit_color', 'c', 'clear', false)


%% mapped from ndf 2

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map3(datarun0, ...
    datarun4, datarun2, cell_type);

% to ndf 4
cell_type_n = 4;
original = get_cell_ids(datarun0, cell_type{cell_type_n});
cell_spec = cell_id{cell_type_n}(:, 1);
overlap = intersect(cell_spec, original);

plot_rf_summaries(datarun0, cell_spec, 'fit_color', 'b', 'foa', 0);
plot_rf_summaries(datarun0, original, 'fit_color', 'r', 'clear', false)
plot_rf_summaries(datarun0, overlap, 'fit_color', 'c', 'clear', false)


% to ndf 0
cell_type_n = 4;
original = get_cell_ids(datarun4, cell_type{cell_type_n});
cell_spec = cell_id{cell_type_n}(:, 2);
overlap = intersect(cell_spec, original);

plot_rf_summaries(datarun4, cell_spec, 'fit_color', 'b', 'foa', 0);
plot_rf_summaries(datarun4, original, 'fit_color', 'r', 'clear', false)
plot_rf_summaries(datarun4, overlap, 'fit_color', 'c', 'clear', false)

%%

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map3(datarun4, ...
    datarun2, datarun0, cell_type);

cell_type_n = 4;
cell_spec = cell_id{cell_type_n}(:, 2);

marks_params.thresh = 4.5;
datarun2 = get_sta_summaries(datarun2, cell_spec, 'marks_params', marks_params);
plot_time_courses(datarun2, cell_spec, 1, 1)

%% mapped from ndf 4 to ndf 0

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun4, datarun0, cell_type);

cell_type_n = 4;
original = get_cell_ids(datarun4, cell_type{cell_type_n});
cell_spec = cell_id{cell_type_n}(:, 1);
overlap = intersect(cell_spec, original);

plot_rf_summaries(datarun4, cell_spec, 'fit_color', 'b', 'foa', 0);
plot_rf_summaries(datarun4, original, 'fit_color', 'r', 'clear', false)
plot_rf_summaries(datarun4, overlap, 'fit_color', 'c', 'clear', false)

%% mapped from ndf 0 to ndf 4

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun0, datarun4, cell_type);

for i = 1:4
cell_type_n = i;
original = get_cell_ids(datarun0, cell_type{cell_type_n});
cell_spec = cell_id{cell_type_n}(:, 1);
overlap = intersect(cell_spec, original);

plot_rf_summaries(datarun0, cell_spec, 'fit_color', 'b', 'foa', 0);
plot_rf_summaries(datarun0, original, 'fit_color', 'r', 'clear', false)
plot_rf_summaries(datarun0, overlap, 'fit_color', 'c', 'clear', false)
end

%% mapped from ndf 0 to ndf 2

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun2, datarun4, cell_type);

for i = 1:5
cell_type_n = i;
original = get_cell_ids(datarun2, cell_type{cell_type_n});
cell_spec = cell_id{cell_type_n}(:, 1);
overlap = intersect(cell_spec, original);

plot_rf_summaries(datarun2, cell_spec, 'fit_color', 'b', 'foa', 0);
plot_rf_summaries(datarun2, original, 'fit_color', 'r', 'clear', false)
plot_rf_summaries(datarun2, overlap, 'fit_color', 'c', 'clear', false)
end

%% mapped from ndf 2 to ndf 0

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun4, datarun2, cell_type);

for i = 1:4
cell_type_n = i;
original = get_cell_ids(datarun4, cell_type{cell_type_n});
cell_spec = cell_id{cell_type_n}(:, 1);
overlap = intersect(cell_spec, original);

plot_rf_summaries(datarun4, cell_spec, 'fit_color', 'b', 'foa', 0);
plot_rf_summaries(datarun4, original, 'fit_color', 'r', 'clear', false)
plot_rf_summaries(datarun4, overlap, 'fit_color', 'c', 'clear', false)
end

%% mapped from ndf 4 to ndf 2

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun2, datarun0, cell_type);

for i = 1:4
cell_type_n = i;
original = get_cell_ids(datarun2, cell_type{cell_type_n});
cell_spec = cell_id{cell_type_n}(:, 1);
overlap = intersect(cell_spec, original);

plot_rf_summaries(datarun2, cell_spec, 'fit_color', 'b', 'foa', 0);
plot_rf_summaries(datarun2, original, 'fit_color', 'r', 'clear', false)
plot_rf_summaries(datarun2, overlap, 'fit_color', 'c', 'clear', false)
end

%% mapped from ndf 2 to ndf 4

[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun0, datarun2, cell_type);

for i = 1:4
cell_type_n = i;
original = get_cell_ids(datarun0, cell_type{cell_type_n});
cell_spec = cell_id{cell_type_n}(:, 1);
overlap = intersect(cell_spec, original);

plot_rf_summaries(datarun0, cell_spec, 'fit_color', 'b', 'foa', 0);
plot_rf_summaries(datarun0, original, 'fit_color', 'r', 'clear', false)
plot_rf_summaries(datarun0, overlap, 'fit_color', 'c', 'clear', false)
end
