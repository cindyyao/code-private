
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);


datadg{1} = load_data('/Volumes/lab/analysis/2016-02-18-0/data004-sorted/data004-sorted', opt);
datadg{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-02-18-0/stimuli/s04.mat';
datadg{1} = load_stim_matlab(datadg{1}, 'user_defined_trigger_interval', 10);

datamb{1} = load_data('/Volumes/lab/analysis/2016-02-18-0/data000-map/data000-map', opt);
datamb{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-02-18-0/stimuli/s00.mat';
datamb{1} = load_stim_matlab(datamb{1});
datamb{2} = load_data('/Volumes/lab/analysis/2016-02-18-0/data002-map/data002-map', opt);
datamb{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-02-18-0/stimuli/s02.mat';
datamb{2} = load_stim_matlab(datamb{2});
datamb{3} = load_data('/Volumes/lab/analysis/2016-02-18-0/data003-map/data003-map', opt);
datamb{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-02-18-0/stimuli/s03.mat';
datamb{3} = load_stim_matlab(datamb{3});
datamb{4} = load_data('/Volumes/lab/analysis/2016-02-18-0/data005-map/data005-map', opt);
datamb{4}.names.stimulus_path = '/Volumes/lab/analysis/2016-02-18-0/stimuli/s05.mat';
datamb{4} = load_stim_matlab(datamb{4});
datamb{5} = load_data('/Volumes/lab/analysis/2016-02-18-0/data006-map/data006-map', opt);
datamb{5}.names.stimulus_path = '/Volumes/lab/analysis/2016-02-18-0/stimuli/s06.mat';
datamb{5} = load_stim_matlab(datamb{5});
datamb{6} = load_data('/Volumes/lab/analysis/2016-02-18-0/data007-map/data007-map', opt);
datamb{6}.names.stimulus_path = '/Volumes/lab/analysis/2016-02-18-0/stimuli/s07.mat';
datamb{6} = load_stim_matlab(datamb{6});

datawn = load_data('/Volumes/lab/analysis/2016-02-18-0/data008-map/data008-map', opt);
%%
n = 5;
i = 1;
[NumSpikesCell,~, StimComb] = get_spikescellstim(datadg{i},datadg{i}.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb,datadg{i});

params_idx = [1 2]; % which parameters to use for classification
[ds_id, nonds_id] = classify_ds(datadg{i}, ds_struct, params_idx);

[NumSpikesCell, StimComb] = get_spikescellstim(datadg{i},ds_id,0);
DG{1} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg{i}));


% i = 5;
% [NumSpikesCell,~,StimComb] = get_spikescellstim_mb(datamb{i},datamb{i}.cell_ids,704,1);
% ds_struct = mbcellanalysis(NumSpikesCell, StimComb,datamb{i});
% 
% params_idx = [1 1]; % which parameters to use for classification
% [ds_id, nonds_id] = classify_ds_(datamb{i}, ds_struct, params_idx, [1 1]);

n = 6;
duration = 704;
[raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, ~, StimComb] = get_spikescellstim_mb(datamb{i},ds_id,duration,1);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb{i} = get_mb_raster(datamb{i}, ds_id, duration);
    for j = 1:length(raster_mb{i})
        if(mb_idx(j, i))
            raster_mb{i}{j} = [];
        end
    end
    trial_dur{i} = get_mb_trial_dur(datamb{i});
end

%% plot cell summary
for cc = 1:length(ds_id)
    plot_mb_raster_ctr(MB(1:3), raster_mb(1:3), trial_dur(1:3), cc, ds_id(cc), 'NDF3', 3, 3, 1)
    plot_mb_raster_ctr(MB(4:6), raster_mb(4:6), trial_dur(4:6), cc, ds_id(cc), 'NDF1', 3, 3, 1)
end

%% classify DSGC into subtypes (directions)
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

%% Direction tuning (moving bar)
% all ds cells
MB_NDF = reshape(MB, 3, 2);
raster_mb_NDF = reshape(raster_mb, 3, 2);

dirn = 4;
D = 4;
T = 1;
BW = 1;
CL = 1;

p_direction = MB{D}.angle{T,BW,CL}';
xx = 0:pi/4:7*pi/4;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 8);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;


%subtypes
clear rho_mb_mean rho_mb_ste dsi_mb_mean dsi_mb_ste rho_mb dsi_mb

for d = 1:2
    figure
    for drug = 1:3
        for cl = 1:3
            subplot(3, 3, (drug-1)*3+cl)
            for i = 1:dirn
                rho_mb{drug, d}{i}{cl} = [];
                dsi_mb{drug, d}{i}{cl} = [];
                for cc = 1:length(idx_dir{i})
                    if ~mb_idx(idx_dir{i}(cc), (d-1)*3+drug) && sum(MB_NDF{drug, d}.rho{T, BW,cl}(idx_dir{i}(cc), :))>0
                    [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
                    y_temp = MB_NDF{drug, d}.rho{T,BW,cl}(idx_dir{i}(cc), :);
                    plot(xsort, y_temp(seq), color(i))
                    ylim([0 1])
        %             pause
                    hold on
                    rho_mb{drug, d}{i}{cl} = [rho_mb{drug, d}{i}{cl}; y_temp(seq)];
                    dsi_mb{drug, d}{i}{cl} = [dsi_mb{drug, d}{i}{cl}; MB_NDF{drug,d}.dsindex{T,BW,cl}(idx_dir{i}(cc))];
                    end
                end
                if ~isempty(rho_mb{drug, d}{i}{cl})
                    rho_mb_mean{cl}{drug, d}(i, :) = mean(rho_mb{drug, d}{i}{cl},1);
                    rho_mb_ste{cl}{drug, d}(i, :) = std(rho_mb{drug, d}{i}{cl},[],1)/sqrt(size(rho_mb{drug, d}{i}{cl}, 1));
                    dsi_mb_mean{cl}{drug, d}(i) = mean(dsi_mb{drug, d}{i}{cl});
                    dsi_mb_ste{cl}{drug, d}(i) = std(dsi_mb{drug, d}{i}{cl})/sqrt(length(dsi_mb{drug, d}{i}{cl}));
                end
            end
            xlabel('direction (rad)')
            ylabel('normalized response')
%             title(ll{d})
            xlim([-pi pi])
            dsi_mb_mean_all{cl} = cell2mat(dsi_mb_mean{cl}(:,d)');
            dsi_mb_ste_all{cl} = cell2mat(dsi_mb_ste{cl}(:,d)');            

        end
    end
end

% plot average (cell type)
ct = {'superior', 'anterior', 'inferior', 'posterior'};
for d = 1:2
    h = figure;
    set(h, 'Position', [1 1 1520,1080])
    for drug = 1:3
        for cl = 1:3
            subplot(3, 3, (drug-1)*3+cl)
            for i = 1:dirn
                if i<=size(rho_mb_mean{cl}{drug, d},1)
                    errorbar(xsort/pi*180, rho_mb_mean{cl}{drug, d}(i, :), rho_mb_ste{cl}{drug, d}(i, :), color(i));
                    hold on
                end
            end
            xlabel('degrees')
            ylabel('normalized mean spike number')
%                 title(ll{d});
            
        end
    end
    legend(ct)
end

% plot average (cell type)
Drug = {'control', 'drug', 'wash'};
Ctr = {'high contrast', 'medium contrast', 'low contrast'};
for d = 1:2
    h = figure;
    set(h, 'Position', [1 1 1520,1080])
    for i = 1:dirn
        for cl = 1:3
            subplot(3, 4, (cl-1)*4+i)
            for drug = 1:3
                if i<=size(rho_mb_mean{cl}{drug, d},1)
                    errorbar(xsort/pi*180, rho_mb_mean{cl}{drug, d}(i, :), rho_mb_ste{cl}{drug, d}(i, :), color(drug));
                    hold on
                end
            end
            xlabel('degrees')
            
            if cl == 1
                title(ct{i})
            end
            if i == 1
                ylabel(Ctr{cl})
            end
        end
        if i == 1
            legend(Drug)
        end
    end
end

%% non-ds cells
nonds_id = get_cell_ids(datawn, 'OFF sustained');
n = 6;
duration = 704;

for i = 1:n    
    [NumSpikesCell, ~, StimComb] = get_spikescellstim_mb(datamb{i},nonds_id,duration,1);
    MB_n{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb_n{i} = get_mb_raster(datamb{i}, nonds_id, duration);
    trial_dur{i} = get_mb_trial_dur(datamb{i});
end

for cc = 11:11 %length(ds_id)
    plot_mb_raster_ctr(MB_n(1:3), raster_mb_n(1:3), trial_dur(1:3), cc, nonds_id(cc), 'NDF3', 3, 3, 1)
    plot_mb_raster_ctr(MB_n(4:6), raster_mb_n(4:6), trial_dur(4:6), cc, nonds_id(cc), 'NDF1', 3, 3, 1)
end

