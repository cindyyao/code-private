wn_map_path = '/Analysis/xyao/2014-07-07-0/data007-map/data007-map';
wn_path = '/Analysis/xyao/2014-07-07-0/data007/data007';
cd /Users/xyao/matlab/code-private/
%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);
datawn_map = load_data(wn_map_path, opt);
datawn_map = load_sta(datawn_map);
datawn = load_data(wn_path, opt);
datawn = load_sta(datawn);
datawn = get_rf_coms(datawn, 'all');


%%
% get array location in display coordinates

%                 EI                               DISPLAY
%
%               195 (1)                           126(2) 4(3)  
%                 / \                               ______
%               /     \                            /      \
%   264 (6)    |       |    126 (2)               /        \
%   386 (5)    |       |    4   (3)       195(1)  \        / 455(4)
%               \     /                            \      /
%                 \ /                               ------
%                455 (4)                         264(6) 386(5)  
array_location_display = [209 416; 261 320; 374 317; 422 425];

% get array location in ei coordinates
elec_corner = [195 126 4 455];
array_location_ei = datawn.ei.position(elec_corner,:);
Tform = maketform('projective', array_location_ei, array_location_display);
test = tformfwd(Tform, array_location_ei)-array_location_display % should be equal or close to zeros

%% get cell location from ei
on1 = get_cell_ids(datawn, 'ON type 1');
on2 = get_cell_ids(datawn, 'ON type 2');

off1 = get_cell_ids(datawn, 'OFF type 1');
off2 = get_cell_ids(datawn, 'OFF type 2');
off3 = get_cell_ids(datawn, 'OFF type 3');
off4 = get_cell_ids(datawn, 'OFF type 4');
off5 = get_cell_ids(datawn, 'OFF type 5');

robust_id = [on1 on2 off1 off2 off3 off4 off5];
robust_idx = get_cell_indices(datawn, robust_id);

distance = 1; % range used to calculate EI com
stixel_size_ = datawn.stimulus.stixel_width;
center_ei = get_ei_com(datawn, robust_id, distance);
center_ei_display = tformfwd(Tform, center_ei);
center_sta = cell2mat(datawn.stas.rf_coms(robust_idx)) * stixel_size_;
center_offset = center_sta - center_ei_display;
[center_offset_pol(:,1), center_offset_pol(:,2)] = cart2pol(center_offset(:,1),center_offset(:,2));

% plot ei-sta offset
figure
quiver(center_ei_display(:,1),center_ei_display(:,2),center_offset(:,1),center_offset(:,2),0)
id_text = cellfun(@num2str, num2cell(robust_id), 'Uni', false);
text(center_ei_display(:,1), center_ei_display(:,2), id_text);

% exclude outliers
figure
hist(center_offset_pol(:,2), 30)
xlabel('offset')
[x,~] = ginput;
outlier_i = find(center_offset_pol(:,2) > x);
center_offset_pol(outlier_i,:) = [];
center_offset(outlier_i,:) = [];
robust_id(outlier_i) = [];
robust_idx(outlier_i) = [];
center_ei_display(outlier_i,:) = [];
center_sta(outlier_i,:) = [];

% replot offset map
figure
quiver(center_ei_display(:,1),center_ei_display(:,2),center_offset(:,1),center_offset(:,2),0)
id_text = cellfun(@num2str, num2cell(robust_id), 'Uni', false);
text(center_ei_display(:,1), center_ei_display(:,2), id_text);

% fit plane
dataX(:,1:2) = center_ei_display;
dataX(:,3) = center_offset(:,1);
dataY(:,1:2) = center_ei_display;
dataY(:,3) = center_offset(:,2);

[normal_x, meanX, pcaExplainedX, sseX] = fit_plane(dataX);
[normal_y, meanY, pcaExplainedY, sseY] = fit_plane(dataY);

%% get DS cell id
dg_path = '/Analysis/xyao/2014-07-07-0/data005-map/data005-map';
dg_stimulus_path = '/Analysis/xyao/2014-07-07-0/stimuli/s05';
datadg = load_data(dg_path, opt);
datadg.names.stimulus_path = dg_stimulus_path;
datadg = load_stim(datadg, 'user_defined_trigger_interval', 10);

[NumSpikesCell, StimComb] = get_spikescellstim(datadg,datadg.cell_ids,0);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);
params_idx = input('indices of parameters for classification: (eg. [1 2])\n'); 
[ds_id, nonds_id] = classify_ds(datadg, ds_struct, params_idx);
ds_idx = get_cell_indices(datadg, ds_id);


%% get DS cell ei location (display coordinate)
center_ei_ds = get_ei_com(datadg, ds_id, 1);
center_ei_ds = tformfwd(Tform, center_ei_ds);

offset_x_ds = meanX*normal_x*ones(length(ds_id),1) - center_ei_ds*normal_x(1:2);
offset_y_ds = meanY*normal_y*ones(length(ds_id),1) - center_ei_ds*normal_y(1:2);
offset_ds = [offset_x_ds offset_y_ds];
% center_sta_ds_est = (center_ei_ds + offset_ds)/stixel_size_;
center_sta_ds_est = center_ei_ds + offset_ds;

%% manually correct center location if needed
center_sta_ds_est = correct_ei_center(datawn_map, ds_id, center_sta_ds_est, stixel_size_);

%% generate circular masks
radius = 30;
display_width = 800;
display_height = 600;

n = size(center_sta_ds_est,1);
radius_all = ones(n,1)*radius;
masks = make_circular_masks(center_sta_ds_est, radius_all, display_width, display_height);
% maploc = '/Volumes/lab/analysis/2015-07-03-0/stimuli/targeted_flash/';
maploc = '/Analysis/xyao/test/stimuli/targeted_flash/';
for i = 1:length(ds_id)
    name = ['map-' num2str(i-1, '%04.0f') '.txt'];
    dlmwrite([maploc name], masks{i},  'delimiter', '\t', 'newline', 'pc')
end
