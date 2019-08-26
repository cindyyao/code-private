opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
datarun{1} = load_data('/Analysis/xyao/2013-03-31-0/data009/data009', opt);
datarun{2} = load_data('/Analysis/xyao/2013-03-31-0/data010/data010', opt);
datarun{2} = load_sta(datarun{2});
marks_params.thresh = 4.0;
datarun{2} = get_sta_summaries(datarun{2}, 'all', 'marks_params', marks_params);
datarun{2} = get_sta_fits_from_vision(datarun{2});

% map ei from WN to gray screen
cell_type = {'ON brisk transient', 'ON transient', 'OFF transient1','OFF transient2', ...
    'OFF sustained', 'OFF slow'};
n = length(cell_type);
[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun{1}, ...
    datarun{2}, cell_type);


%plot results of mapping
figure; 
for i = 1:n
    subplot(2, 4, i)
    plot_rf_summaries(datarun{2}, cell_type{i});
    hold on
    plot_rf_summaries(datarun{2}, cell_id{i}(:, 2), 'fit_color', 'r');
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

cps = cell(n, 1);
for ct = 1:n
    % get center locations of all cells in this mosaic
    center_pts = rf_centers(datarun{2}, cell_id{ct}(:, 2), 'com');
    % compute distance to all neighbors
    distances = ipdm(center_pts);
    % calculate mean rf radius for this cell type
    idx = get_cell_indices(datarun{2}, cell_type{ct});
    radii = zeros(length(idx), 1);
    for cc = 1:length(idx)
        radii(cc) = max(datarun{2}.stas.fits{idx(cc)}.sd);
    end
    radius = mean(radii);
    % get neighbor cell pairs indices
    cps_temp = [];
    for cc1 = 1:size(cell_id{ct}, 1)-1
        for cc2 = cc1+1:size(cell_id{ct}, 1)
            if distances(cc1, cc2) < extension_factor*radius*2
                cps_temp = [cps_temp; cc1 cc2 distances(cc1, cc2)];
            end
        end
    end
    cps{ct} = cps_temp;
end



%% get spike histgram count

for cc = 1:size(cell_idx_matched, 1)
    idx = cell_idx_matched(cc, 1);
    spikes = datarun{1}.spikes{idx};
    spikeshc = histc(spikes, edges);
%     spikeshc = sparse(spikeshc);
    datarun{1}.spikeshc{idx} = spikeshc;
end

        
%% plot correlation function
% get cell type index
ct = 1;

figure
for cp = 1:1
    cpi = 2:2;
% get cell indices
idx1 = cell_idx{ct}(cps{ct}(cpi(cp), 1), 1);
idx2 = cell_idx{ct}(cps{ct}(cpi(cp), 2), 1);
id1 = cell_id{ct}(cps{ct}(cpi(cp), 1), 2);
id2 = cell_id{ct}(cps{ct}(cpi(cp), 2), 2);

cf = xcorr(datarun{1}.spikeshc{idx1}, datarun{1}.spikeshc{idx2}, maxlag);

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
plot_rf_summaries(datarun{2}, cell_type{ct});
hold on
plot_rf_summaries(datarun{2}, [id1 id2], 'fit_color', 'r');
end


%% get cell pairs for each cell type pair
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

cps = cell(n, n);
for ct1 = 1:n-1
    for ct2 = ct1+1:n
    % get center locations of all cells in these two mosaics
    center_pts1 = rf_centers(datarun{2}, cell_id{ct1}(:, 2), 'com');
    center_pts2 = rf_centers(datarun{2}, cell_id{ct2}(:, 2), 'com');
    
    % compute distance to all neighbors
    distances = ipdm(center_pts1, center_pts2);
    % calculate mean rf radius for these two cell types
    idx = get_cell_indices(datarun{2}, cell_type{ct1});
    radii = zeros(length(idx), 1);
    for cc = 1:length(idx)
        radii(cc) = max(datarun{2}.stas.fits{idx(cc)}.sd);
    end
    radius1 = mean(radii);
    idx = get_cell_indices(datarun{2}, cell_type{ct2});
    radii = zeros(length(idx), 1);
    for cc = 1:length(idx)
        radii(cc) = max(datarun{2}.stas.fits{idx(cc)}.sd);
    end
    radius2 = mean(radii);
    % get neighbor cell pairs indices
    cps_temp = [];
    for cc1 = 1:size(cell_id{ct1}, 1)
        for cc2 = 1:size(cell_id{ct2}, 1)
            if distances(cc1, cc2) < extension_factor*(radius1 + radius2);
                cps_temp = [cps_temp; cc1 cc2 distances(cc1, cc2)];
            end
        end
    end
    cps{ct1, ct2} = cps_temp;
    end
end

        
    
%% plot correlation function
% get cell type index
ct1 = 5;
ct2 = 6;

figure
for cp = 1:10
    cpi = 1:10;
% get cell indices
idx1 = cell_idx{ct1}(cps{ct1, ct2}(cpi(cp), 1), 1);
idx2 = cell_idx{ct2}(cps{ct1, ct2}(cpi(cp), 2), 1);
id1 = cell_id{ct1}(cps{ct1, ct2}(cpi(cp), 1), 2);
id2 = cell_id{ct2}(cps{ct1, ct2}(cpi(cp), 2), 2);

cf = xcorr(datarun{1}.spikeshc{idx1}, datarun{1}.spikeshc{idx2}, maxlag);

subplot(5, 8, 4*cp-3:4*cp-1)
bar([-maxlag:maxlag]*bin_size*1000, cf, 'hist')
if cp == 9
    xlabel('time lag (ms)')
end
xlim([-maxlag*bin_size*1000 maxlag*bin_size*1000])
ylim([0 max(cf)*1.2])
ylabel('correlation')
if cp == 1
    
    title([cell_type{ct1} ' and ' cell_type{ct2}]);
end
subplot(5, 8, 4*cp)
% plot_rf_summaries(datarun{2}, cell_type{ct1});
% hold on
% plot_rf_summaries(datarun{2}, cell_type{ct2}, 'fit_color', 'b');
plot_rf_summaries(datarun{2}, [id1 id2]);
end

