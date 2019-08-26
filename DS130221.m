%%
% cd /Users/xyao/Desktop/dscellcalculation 
%within the ds run, how many and which cells are actually ds?
datarun=load_data('/Volumes/lab/analysis/2013-02-21-0/data008-sorted/data008-sorted');
datarun=load_neurons(datarun);
datarun = load_ei(datarun, id_dir{1}, 'array_type', 519);
% datarun=load_ei(datarun, id_sds);
% datarun_path='/Volumes/backup/data/2016-06-16-0/';
datarun.names.stimulus_path = '/Volumes/lab/analysis/2013-02-21-0/stimuli/s08';
datarun = load_stim(datarun, 'user_defined_trigger_interval', 10);
 
[NumSpikesCell, ~, StimComb]=get_spikescellstim(datarun,datarun.cell_ids,0, 1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb); 

params_idx = [1 2]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datarun, ds_struct, params_idx);


[NumSpikesCell, ~, StimComb]=get_spikescellstim(datarun,ds_id,0, 1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb); 

%% dg
% ds_id = datadg.cell_ids;
n = 1;
i = 1;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
[NumSpikesCell, ~,StimComb] = get_spikescellstim(datarun,ds_id,0,1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
raster_dg{i} = get_ds_raster(datarun, ds_id);

%% classify into 4 directions
d = 1;
t = 1;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(ds_struct.U{t}, ds_struct.V{t});
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(ds_struct.U{t}, ds_struct.V{t}, x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = I;
    id_dir{i} = ds_id(idx_dir{i});
end

%% plot cell summary
for dir = 1:4
    for cc = 1:length(id_dir{dir})
        plot_ds_raster(DG, raster_dg, idx_dir{dir}(cc), id_dir{dir}(cc), '', 1, 1, 1)
    end
end


%% neighboring pairs
pos = datarun.ei.position;
mode = 'neg';
neighbors = [];
for cc1 = 1:length(id_sds)
    for cc2 = cc1+1:length(id_sds)
        id1 = id_sds(cc1);
        idx1 = get_cell_indices(datarun, id1);
        ei1 = datarun.ei.eis{idx1};
        com1 = ei_com_xy(ei1, pos, 30*3, mode);
        id2 = id_sds(cc2);
        idx2 = get_cell_indices(datarun, id2);
        ei2 = datarun.ei.eis{idx2};
        com2 = ei_com_xy(ei2, pos, 30*3, mode);
        if pdist([com1;com2]) < 150
            neighbors = [neighbors; id1 id2];
        end
    end
end

coms = [];
for cc = 1:length(id_sds)
    id = id_sds(cc);
    idx = get_cell_indices(datarun, id);
    ei = datarun.ei.eis{idx};
    com = ei_com_xy(ei, pos, 30*3, mode);
    coms = [coms; com];
end

corner_i = [4 126 195 264 386 455 4];
corner_position = datarun.ei.position(corner_i, :);
figure
for cc = 1:length(id_sds)
    plot(coms(cc, 1), coms(cc, 2),'ko')
    hold on
%     text(coms(cc, 1)+5, coms(cc, 2)+5, num2str(id_dir{1}(cc)), 'FontSize', 10)
    
end

%%
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);
datafs{1} = load_data('/Volumes/lab/analysis/2013-02-21-0/data000-map/data000-map', opt);
datafs{2} = load_data('/Volumes/lab/analysis/2013-02-21-0/data002-map/data002-map', opt);
datafs{3} = load_data('/Volumes/lab/analysis/2013-02-21-0/data004-map/data004-map', opt);


duration = 3600;
bin_size = 0.0005;
max_lag = 20;
ct = 1;
N = 10000;
ll = {'NDF 4', 'NDF 3', 'NDF 2', 'NDF 1', 'NDF 0'};
xx = -max_lag*bin_size:bin_size:max_lag*bin_size;
peak = zeros(size(neighbors, 1), 3);
for cp = 1:size(neighbors, 1)
    if mod(cp, 5) == 1
        FigHandle = figure;
        set(FigHandle, 'Position', [1 1 2000 2000])
        fig_i = 1;
    end
    id1 = neighbors(cp, 1);
    id2 = neighbors(cp, 2);
    for i = 1:3
%         if ~fs_idx(find(ds_id == id1), i)
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
            subplot(5, 3, (fig_i-1)*3+i)
%     if p < 0.01
%         bar(xx, A, 'r')
%     else
            bar(xx, A, 'b')
%     end
            if fig_i == 1
                title(ll{i})
            end
            xlim([-0.01 0.01])
%         end
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
