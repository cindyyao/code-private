cd /Users/xyao/matlab/code-private/
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);


datadg{1} = load_data('/Volumes/lab/analysis/2016-09-15-0/data007/data007', opt);
datadg{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-09-15-0/stimuli/s07';
datadg{1} = load_stim(datadg{1}, 'user_defined_trigger_interval', 10);

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg{1},datadg{1}.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);
params_idx = [1 2]; % which parameters to use for classification

[ds_id, nonds_id] = classify_ds(datadg{1}, ds_struct, params_idx);
ds_idx = get_cell_indices(datadg{1}, ds_id);
[NumSpikesCell, ~, StimComb] = get_spikescellstim(datadg{1},ds_id,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);

datadg{2} = load_data('/Volumes/lab/analysis/2016-09-15-0/data009-map/data009-map', opt);
datadg{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-09-15-0/stimuli/s09';
datadg{2} = load_stim(datadg{2}, 'user_defined_trigger_interval', 10);

%% dg
n = 2;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:2
    [NumSpikesCell,~, StimComb] = get_spikescellstim(datadg{i},ds_id,0,1);
    DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster_dg{i} = get_ds_raster(datadg{i}, ds_id);
end


delta_p = 2; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_dg = cell(n, 1);
max_r = cell(n, 1);
norm_max_r = cell(n, 1);

for i = 1:2
    [raster_p_sum{i}, p_idx{i}] = get_pdirection_raster(raster_dg{i}, DG{i}.angle{delta_p});
    MAG_all_norm_dg{i} = normalize_MAG(DG{i});
    rep = datadg{i}.stimulus.repetitions;
end

%% plot cell summary
% map_idx = map_ei(datadg{1}, datadg{2});
tp = 2;
for cc = 2:length(ds_id)
    if sum(DG{2}.RHO{tp}(cc, :))>0
        fh = figure(1);
        set(fh, 'Position', [1 1 1080 500])
        [idx, xx, yy] = subplot_idx_12(1, 2);
        tt = DG{1}.theta{1}(1, :);
        for j = 1:2 
            cell_idx = cc;
            h = subplot(xx, yy, idx(j, 13:16)) ;
            u_temp = DG{j}.U{tp}(cell_idx);
            v_temp = DG{j}.V{tp}(cell_idx);
            alim = max(sqrt(u_temp^2+v_temp^2), 3);
            P = polar(0, alim);
            set(P, 'Visible', 'off')
            hold on
            compass(DG{j}.U{tp}(cell_idx), DG{j}.V{tp}(cell_idx), 'r');
            polar(tt, DG{j}.rho{tp}(cell_idx, :), 'b');
            polar_theta_off(h)
            for i = 1:12
                subplot(xx, yy, idx(j, i)); plot_raster(squeeze(raster_dg{j}{cell_idx}(1, tp, i, :)), 0, 8)

                if mod(idx(j, i), yy) == 1
                    ylabel('trial number')
                end
                if idx(j, i) > yy*(xx-1)
                    xlabel('time (s)')
                end
            end
        end
        name = [num2str(ds_id(cc))];
        screen_size = [13 8];
        set(figure(1), 'paperpositionmode', 'auto');
        set(gcf, 'PaperUnits', 'inch');
        set(figure(1), 'PaperSize', screen_size);
        print(figure(1), '-dpdf', name)
        close

    end
end
%%
t = 2;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(DG{1}.U{t}, DG{1}.V{t});
color = 'brgkc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{1}.U{t}, DG{1}.V{t}, x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = I;
    id_dir{i} = ds_id(idx_dir{i});
end

%% AP5 suppression
color = 'brgkcmy';
dirn = 4;
D = 1;
T = 2;


p_direction = DG{D}.angle{T}';
xx = 0:pi/6:11*pi/6;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 12);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

for drug = 1:2
    for i = 1:4
        RHO_dg{drug}{i} = [];
        rho_dg{drug}{i} = [];
        for cc = 1:length(idx_dir{i})
            if sum(DG{2}.RHO{T}(idx_dir{i}(cc), :))>0
                [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
                Y_temp = DG{drug}.RHO{T}(idx_dir{i}(cc), :);
                y_temp = DG{drug}.rho{T}(idx_dir{i}(cc), :);
                RHO_dg{drug}{i} = [RHO_dg{drug}{i}; Y_temp(seq)];
                rho_dg{drug}{i} = [rho_dg{drug}{i}; y_temp(seq)];
            end
        end
        RHO_mean{drug}(i, :) = mean(RHO_dg{drug}{i},1);
        RHO_ste{drug}(i, :) = std(RHO_dg{drug}{i},[],1)/sqrt(size(RHO_dg{drug}{i}, 1));
        rho_mean{drug}(i, :) = mean(rho_dg{drug}{i},1);
        rho_ste{drug}(i, :) = std(rho_dg{drug}{i},[],1)/sqrt(size(rho_dg{drug}{i}, 1));
    end
end

Drug = {'control', 'AP5'};
h = figure;
set(h, 'Position', [1 1 1520,1080])
for i = 1:4
    subplot(2, 2, i)
    for drug = 1:2
        errorbar(xsort/pi*180, RHO_mean{drug}(i, :), RHO_ste{drug}(i, :), color(drug));
        hold on
    end
    xlabel('degrees')
    ylabel('spike count')
    legend(Drug)
end

Drug = {'control', 'AP5'};
h = figure;
set(h, 'Position', [1 1 1520,1080])
for i = 1:4
    subplot(2, 2, i)
    for drug = 1:2
        errorbar(xsort/pi*180, rho_mean{drug}(i, :), rho_ste{drug}(i, :), color(drug));
        hold on
    end
    xlabel('degrees')
    ylabel('normalized response')
    legend(Drug)
end
