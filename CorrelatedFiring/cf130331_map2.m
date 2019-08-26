opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
datarun{5} = load_data('/Analysis/xyao/2013-03-31-0/data009-map2/data009-map2', opt);
datarun{6} = load_data('/Analysis/xyao/2013-03-31-0/data010-map2/data010-map2', opt);
datarun{6} = load_sta(datarun{6});
marks_params.thresh = 4.0;
datarun{6} = get_sta_summaries(datarun{6}, 'all', 'marks_params', marks_params);
datarun{6} = get_sta_fits_from_vision(datarun{6});

% map ei from WN to gray screen
cell_type = {'ON brisk transient', 'ON transient', 'OFF transient1', 'OFF transient2'};
n = length(cell_type);
cell_id_m2 = cell(n, 1);
cell_idx_m2 = cell(n, 1);
for i = 1:n
    id_temp = get_cell_ids(datarun{6}, cell_type{i});
    for j = 1:length(id_temp)
        if isempty(find(datarun{5}.cell_ids == id_temp(j)))
            id_temp(j) = 0;
        end
    end
    cell_id_m2{i} = id_temp(id_temp>0);
    cell_idx_m2{i}(:, 1) = get_cell_indices(datarun{5}, cell_id_m2{i})';
    cell_idx_m2{i}(:, 2) = get_cell_indices(datarun{6}, cell_id_m2{i})';
end

cell_id_matched_m2 = [cell_id_m2{1} cell_id_m2{2} cell_id_m2{3} cell_id_m2{4}];
cell_idx_matched_m2(:, 1) = get_cell_indices(datarun{5}, cell_id_matched_m2)';
cell_idx_matched_m2(:, 2) = get_cell_indices(datarun{6}, cell_id_matched_m2)';

%plot results of mapping
figure; 
for i = 1:n
    subplot(2, 4, i)
    plot_rf_summaries(datarun{6}, cell_type{i});
    hold on
    plot_rf_summaries(datarun{6}, cell_id_m2{i}, 'fit_color', 'r');
    title(cell_type{i})
end

% set parameters
extension_factor = 1.5;

duration = 1200; % unit: s
bin_size = 0.0005; % unit: s
edges = 0:bin_size:duration;

maxlag = round(0.1/bin_size); 

    
%% get cell pairs for each cell type 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cps_m = cell(n, 1);
for ct = 1:n
    % get center locations of all cells in this mosaic
    center_pts = rf_centers(datarun{6}, cell_id_m2{ct}, 'com');
    % compute distance to all neighbors
    distances = ipdm(center_pts);
    % calculate mean rf radius for this cell type
    idx = get_cell_indices(datarun{6}, cell_type{ct});
    radii = zeros(length(idx), 1);
    for cc = 1:length(idx)
        radii(cc) = max(datarun{6}.stas.fits{idx(cc)}.sd);
    end
    radius = mean(radii);
    % get neighbor cell pairs indices
    cps_temp = [];
    for cc1 = 1:length(cell_id_m2{ct})-1
        for cc2 = cc1+1:length(cell_id_m2{ct})
            if distances(cc1, cc2) < extension_factor*radius*2
                cps_temp = [cps_temp; cc1 cc2 distances(cc1, cc2)];
            end
        end
    end
    cps_m{ct} = cps_temp;
end



%% get spike histgram count

for cc = 1:length(cell_id_matched_m2)
    idx = get_cell_indices(datarun{5}, cell_id_matched_m2(cc));
    spikes = datarun{5}.spikes{idx};
    spikeshc = histc(spikes, edges);
%     spikeshc = sparse(spikeshc);
    datarun{5}.spikeshc{idx} = spikeshc;
end

        
%% plot correlation function
% get cell type index
ct = 4;

figure
for cp = 1:10
    cpi = 21:30;
% get cell indices
idx1 = cell_idx_m2{ct}(cps_m{ct}(cpi(cp), 1), 1);
idx2 = cell_idx_m2{ct}(cps_m{ct}(cpi(cp), 2), 1);
id1 = cell_id_m2{ct}(cps_m{ct}(cpi(cp), 1));
id2 = cell_id_m2{ct}(cps_m{ct}(cpi(cp), 2));

cf = xcorr(datarun{5}.spikeshc{idx1}, datarun{5}.spikeshc{idx2}, maxlag);

subplot(5, 8, 4*cp-3:4*cp-1)
bar([-maxlag:maxlag]*bin_size*1000, cf, 'hist')
if cp == 9
    xlabel('time lag (ms)')
end
xlim([-maxlag*bin_size*1000 maxlag*bin_size*1000])
ylim([0 max(cf)*1.2])
ylabel('correlation')
if cp == 1
    title(cell_type{ct});
end
subplot(5, 8, 4*cp)
plot_rf_summaries(datarun{6}, cell_type{ct});
hold on
plot_rf_summaries(datarun{6}, [id1 id2], 'fit_color', 'r');
end

