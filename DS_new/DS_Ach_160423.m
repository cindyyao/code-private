cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);


datadg = load_data('/Volumes/lab/analysis/2016-04-23-0/data002-DSsorted/data002-DSsorted', opt);
datadg.names.stimulus_path = '/Volumes/lab/analysis/2016-04-23-0/stimuli/s02.mat';
datadg = load_stim_matlab(datadg, 'user_defined_trigger_interval', 10);

datarun = load_data('/Volumes/lab/analysis/2016-04-23-0/data003-005-map/data003-005-map', opt);
time_points = [2130 4270];
datamb(1:3) = split_datarun(datarun, time_points);
datamb{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-04-23-0/stimuli/s03.mat';
datamb{1} = load_stim_matlab(datamb{1});
datamb{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-04-23-0/stimuli/s04.mat';
datamb{2} = load_stim_matlab(datamb{2});
datamb{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-04-23-0/stimuli/s05.mat';
datamb{3} = load_stim_matlab(datamb{3});

datawn = load_data('/Volumes/lab/analysis/2016-04-23-0/data001-map/data001-map', opt);

%% 
load('DS160423.mat')
i = 1;
[NumSpikesCell, ~,StimComb] = get_spikescellstim(datadg,ds_id,0,1);
DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datadg));

n = 3;
duration = [1781.4 2131.4 1781.4]; %sec
bin_size = 0.025; %sec
[raster_mb, MB, trial_dur, raster_p_sum_mb, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},ds_id,duration(i),bin_size);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb{i} = get_mb_raster(datamb{i}, ds_id, duration(i));
    for j = 1:length(raster_mb{i})
        if(mb_idx(j))
            raster_mb{i}{j} = [];
        end
    end
    trial_dur{i} = get_mb_trial_dur(datamb{i});
end

ctr_p = 1; % choose which params to use to calculate prefer direction indices 
MAG_all_norm_mb = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}, raster_p_sum_all{i}] = get_pdirection_raster(raster_mb{i}, MB{1}.angle{ctr_p});
    MAG_all_norm_mb{i} = normalize_MAG(MB{i});
    rep = datamb{i}.stimulus.repetitions;
end

%% plot cell summary
for cc = 1:length(ds_id)
    plot_mb_raster_ctr_12(MB, raster_mb, trial_dur, cc, ds_id(cc), 'NDF0', 3, 6, 0)
end
cc = 4;
plot_mb_raster_one_12(MB(1), raster_mb(1), trial_dur(1), cc, ds_id(cc), 'NDF0', 1, 1, 0)
plot_mb_raster_ctr_12(MB(1), raster_mb(1), trial_dur(1), cc, ds_id(cc), 'NDF0', 1, 1, 0)
%% Direction tuning (moving bar)
% all ds cells

color = 'brgkcm';
dirn = 4;
D = 1;
T = 1;
BW = 1;
CL = 1;

p_direction = MB{D}.angle{T,BW,CL}';
xx = 0:pi/6:11*pi/6;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 12);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;


%subtypes
clear rho_mb_mean rho_mb_ste dsi_mb_mean dsi_mb_ste rho_mb dsi_mb

% figure
for drug = 1:3
    for cl = 1:6
%         subplot(3, 6, (drug-1)*6+cl)
        for i = 1:dirn
            rho_mb{drug}{i}{cl} = [];
            RHO_mb{drug}{i}{cl} = [];
            dsi_mb{drug}{i}{cl} = [];
            for cc = 1:length(idx_dir{i})
                if ~mb_idx(idx_dir{i}(cc))% && sum(MB_NDF{drug, d}.RHO{T, BW,cl}(idx_dir{i}(cc), :))>0
                [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
                Y_temp = MB{drug}.RHO{T,BW,cl}(idx_dir{i}(cc), :);
                y_temp = MB{drug}.rho{T,BW,cl}(idx_dir{i}(cc), :);
%                 y_temp = MB{drug}.RHO{T,BW,cl}(idx_dir{i}(cc), :)/bin_size;
%                 plot(xsort, y_temp(seq), color(i))
%                     ylim([0 1])
    %             pause
%                 hold on
                    if sum(y_temp ~= 0)
                        rho_mb{drug}{i}{cl} = [rho_mb{drug}{i}{cl}; y_temp(seq)];
                        RHO_mb{drug}{i}{cl} = [RHO_mb{drug}{i}{cl}; Y_temp(seq)];
                        dsi_mb{drug}{i}{cl} = [dsi_mb{drug}{i}{cl}; MB{drug}.dsindex{T,BW,cl}(idx_dir{i}(cc))];
                    end
                end
            end
            if ~isempty(rho_mb{drug}{i}{cl})
                rho_mb_mean{cl}{drug}(i, :) = mean(rho_mb{drug}{i}{cl},1);
                rho_mb_ste{cl}{drug}(i, :) = std(rho_mb{drug}{i}{cl},[],1)/sqrt(size(rho_mb{drug}{i}{cl}, 1));
                RHO_mb_mean{cl}{drug}(i, :) = mean(RHO_mb{drug}{i}{cl},1);
                RHO_mb_ste{cl}{drug}(i, :) = std(RHO_mb{drug}{i}{cl},[],1)/sqrt(size(RHO_mb{drug}{i}{cl}, 1));
                dsi_mb_mean{cl}{drug}(i) = mean(dsi_mb{drug}{i}{cl});
                dsi_mb_ste{cl}{drug}(i) = std(dsi_mb{drug}{i}{cl})/sqrt(length(dsi_mb{drug}{i}{cl}));
            end
        end
%         xlabel('direction (rad)')
%         ylabel('spike number')
%             title(ll{d})
%         xlim([-pi pi])
        dsi_mb_mean_all{cl} = cell2mat(dsi_mb_mean{cl}');
        dsi_mb_ste_all{cl} = cell2mat(dsi_mb_ste{cl}');            

    end
end

% plot average (cell type)
ct = {'superior', 'anterior', 'inferior', 'posterior'};
h = figure;
set(h, 'Position', [1 1 1520,1080])
for drug = 1:3
    for cl = 1:6
        subplot(3, 6, (drug-1)*6+cl)
        for i = 1:dirn
            if i<=size(rho_mb_mean{cl}{drug},1)
                errorbar(xsort/pi*180, rho_mb_mean{cl}{drug}(i, :), rho_mb_ste{cl}{drug}(i, :), color(i));
                hold on
            end
        end
        xlabel('degrees')
        ylabel('spike number')
%                 title(ll{d});

    end
end
legend(ct)

% plot average (cell type)
Drug = {'control', 'drug', 'wash'};
% Ctr = {'high contrast', 'medium contrast', 'low contrast'};
h = figure;
set(h, 'Position', [1 1 1520,1080])
for i = 1:dirn
    for cl = 1:6
        subplot(4, 6, (i-1)*6+cl)
        for drug = 1:2
            if i<=size(rho_mb_mean{cl}{drug},1)
                errorbar(xsort/pi*180, rho_mb_mean{cl}{drug}(i, :), rho_mb_ste{cl}{drug}(i, :), color(drug));
                hold on
            end
        end
        xlabel('degrees')

        if cl == 1
            title(ct{i})
        end
%         if i == 1
%             ylabel(Ctr{cl})
%         end
    end
    if i == 1
        legend(Drug)
    end
end

% control
ctr = {'80%', '40%', '20%', '10%', '5%', '2.5%'};
ct = {'superior', 'anterior', 'inferior', 'posterior'};
h = figure;
set(h, 'Position', [1 1 1520,1080])
drug = 1;
for i = 1:dirn
    subplot(3, 2, i)
    for cl = 1:6
        if i<=size(rho_mb_mean{cl}{drug},1)
            errorbar(xsort/pi*180, rho_mb_mean{cl}{drug}(i, :)/max(rho_mb_mean{cl}{drug}(i, :)), rho_mb_ste{cl}{drug}(i, :), color(cl));
            hold on
        end
    end
    xlabel('degrees')
    ylabel('spike number')
    title(ct{i})
%                 title(ll{d});

end
legend(ctr)
subplot(3,2,5)

% unnormalized
ctr = {'80%', '40%', '20%', '10%', '5%', '2.5%'};
ct = {'superior', 'anterior', 'inferior', 'posterior'};
h = figure;
set(h, 'Position', [1 1 1520,1080])
drug = 1;
for i = 1:dirn
    subplot(3, 2, i)
    for cl = 1:6
        if i<=size(RHO_mb_mean{cl}{drug},1)
            errorbar(xsort/pi*180, RHO_mb_mean{cl}{drug}(i, :), RHO_mb_ste{cl}{drug}(i, :), color(cl));
            hold on
        end
    end
    xlabel('degrees')
    ylabel('spike number')
    title(ct{i})
%                 title(ll{d});

end
legend(ctr)
subplot(3,2,5)

%% tuning width
drug = 1;
for ctr = 1:6
    C = 1;
    for ct = 1:4
        figure
        CC = 1;
        for cc = 1:size(rho_mb{1}{ct}{ctr}, 1)
            if ~(sum(rho_mb{1}{ct}{ctr}(cc, :)) == 0)
                [f{ct}{ctr}{CC}, G{ct}{ctr}{CC}] = fit_gaussian(xsort, rho_mb{1}{ct}{ctr}(cc, :));
                xfit = linspace(xsort(1), xsort(end), 100);
                yfit = f{ct}{ctr}{CC}.a*exp(-((xfit-f{ct}{ctr}{CC}.b)/f{ct}{ctr}{CC}.c).^2)+f{ct}{ctr}{CC}.d;
                subplot(5,7,cc)
                plot(xsort, rho_mb{1}{ct}{ctr}(cc, :), 'b')
                hold on
                plot(xfit, yfit, 'r')
                tuning_width{ct}{ctr}(CC) = f{ct}{ctr}{CC}.c*0.83*180*2/pi;
                tuning_width_all{ctr}(C) = f{ct}{ctr}{CC}.c*0.83*180*2/pi;
                dc{ct}{ctr}(CC) = f{ct}{ctr}{CC}.d;
                CC = CC + 1;
                C = C + 1;
            end
        end
        tuning_width_mean{ct}(ctr) = mean(tuning_width{ct}{ctr});
        tuning_width_ste{ct}(ctr) = std(tuning_width{ct}{ctr})/sqrt(length(tuning_width{ct}{ctr}));
        dc_mean{ct}(ctr) = mean(dc{ct}{ctr});
        dc_ste{ct}(ctr) = std(dc{ct}{ctr})/sqrt(length(dc{ct}{ctr}));
    end
    tuning_width_mean_all(ctr) = mean(tuning_width_all{ctr});
    tuning_width_ste_all(ctr) = std(tuning_width_all{ctr})/sqrt(length(tuning_width_all{ctr}));

end

aa = 5;
ctr_x = [80 40 20 10 5 2.5];
figure
for ct = 1:4
    errorbar(ctr_x(1:aa), tuning_width_mean{ct}(1:aa), tuning_width_ste{ct}(1:aa), 'color', color(ct));
    hold on
end
set(gca, 'Xscale', 'log')
ylim([0 200])
xlim([3 100])
legend(dscell_type)
xlabel('contrast %')
ylabel('degree')

figure
errorbar(ctr_x, tuning_width_mean_all, tuning_width_ste_all)
set(gca, 'Xscale', 'log')
xlabel('contrast %')
ylabel('degree')
title('tuning width')
ylim([0 200])



%% tuning width
% circular standard deviation


%% contrast response function

% DS cell
color = 'brgkc';
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'Hex', 'wash'};
for drug = 1:3
    for cc = 1:length(raster_p_sum{1})
        if ~isempty(raster_p_sum{drug}{cc})
            pd_spikes{drug}(cc,:) = cellfun('length', raster_p_sum{drug}{cc})/datamb{drug}.stimulus.repetitions;
        end
    end
    
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                spike_temp = cellfun('length', raster_p_sum{drug}{idx_dir{dir}(cc)})/datamb{drug}.stimulus.repetitions;
                if sum(spike_temp) ~= 0
                    pd_dir_spikes{drug}{dir}(CC,:) = cellfun('length', raster_p_sum{drug}{idx_dir{dir}(cc)})/datamb{drug}.stimulus.repetitions;
                    pd_dir_spikes_nor{drug}{dir}(CC,:) = pd_dir_spikes{drug}{dir}(CC,:)/max(pd_dir_spikes{1}{dir}(CC,:));
                    CC = CC + 1;
                end
            end
        end
        pd_dir_spikes_nor{drug}{dir} = nan2empty(pd_dir_spikes_nor{drug}{dir});
    end
end

drug = 1;
figure
for dir = 1:4
    errorbar([80 40 20 10 5 2.5], mean(pd_dir_spikes_nor{drug}{dir}), std(pd_dir_spikes_nor{drug}{dir})/sqrt(size(pd_dir_spikes_nor{drug}{dir}, 1)), 'color', color(dir));
    hold on
end
legend(dscell_type)
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('normalized response')
title(condition{drug})

% non-DSRGC
clear nonds_id
nonds_type = {'ON transient', 'ON sustained', 'OFF transient', 'OFF transient slow', 'OFF transient large'};
for i = 1:length(nonds_type)
    nds_temp = get_cell_ids(datawn, nonds_type{i});
    nonds_id{i} = intersect(nds_temp, datarun.cell_ids);
end
nonds_id_all = sort(cell2mat(nonds_id));
for i = 1:length(nonds_type)
    [~,nonds_idx{i}] = intersect(nonds_id_all, nonds_id{i});
end

for i = 1:n
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{i},nonds_id_all,duration(i),bin_size);
    MB_n{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb, datamb{i}));
    raster_mb_nds{i} = get_mb_raster(datamb{i}, nonds_id_all, duration(i));
end

ctr_p = 1; % choose which params to use to calculate prefer direction indices 

for i = 1:n
    [raster_p_sum_nds{i}, p_idx{i}, raster_p_sum_nds_all{i}] = get_pdirection_raster(raster_mb_nds{i}, MB_n{i}.angle{ctr_p});
end

for drug = 1:3
%     for cc = 1:length(raster_p_sum_nds{1})
%         if ~isempty(raster_p_sum_nds{drug}{cc})
%             pd_spikes_nds{drug}(cc,:) = cellfun('length', raster_p_sum_nds{drug}{cc})/datamb{drug}.stimulus.repetitions;
%         end
%     end
    
    for ct = 1:length(nonds_type)
        CC = 1; 
        for cc = 1:length(nonds_id{ct})
            if ~isempty(raster_p_sum_nds{drug}{nonds_idx{ct}(cc)})
                spike_temp = cellfun('length', raster_p_sum_nds{drug}{nonds_idx{ct}(cc)})/datamb{drug}.stimulus.repetitions;
                if sum(spike_temp) ~= 0
                    pd_ct_spikes{drug}{ct}(CC,:) = cellfun('length', raster_p_sum_nds{drug}{nonds_idx{ct}(cc)})/datamb{drug}.stimulus.repetitions;
                    pd_ct_spikes_nor{drug}{ct}(CC,:) = pd_ct_spikes{drug}{ct}(CC,:)/max(pd_ct_spikes{1}{ct}(CC,:));
                    CC = CC + 1;
                end
            end         
        end
        pd_ct_spikes_nor{drug}{ct} = nan2empty(pd_ct_spikes_nor{drug}{ct});
    end
end

drug = 1;

figure
for ct = 1:5
    errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes_nor{drug}{ct}), std(pd_ct_spikes_nor{drug}{ct})/sqrt(size(pd_ct_spikes_nor{drug}{ct}, 1)), 'color', color(ct));
    hold on
end
legend(nonds_type)
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('normalized response')
title(condition{drug})

figure
for ct = 1:5
    h1 = errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes_nor{drug}{ct}), std(pd_ct_spikes_nor{drug}{ct})/sqrt(size(pd_ct_spikes_nor{drug}{ct}, 1)), 'b');
    hold on
end
for ct = 1:4
    h2 = errorbar([80 40 20 10 5 2.5], mean(pd_dir_spikes_nor{drug}{ct}), std(pd_dir_spikes_nor{drug}{ct})/sqrt(size(pd_dir_spikes_nor{drug}{ct}, 1)), 'r');
    hold on
end
legend([h1 h2], {'nDS cell', 'DS cell'})
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('normalized response')
title(condition{drug})

%%
figure
errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes_nor{1}{1}), std(pd_ct_spikes_nor{1}{1})/sqrt(size(pd_ct_spikes_nor{1}{1}, 1)), 'color', color(1));
hold on
errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes_nor{2}{1}), std(pd_ct_spikes_nor{2}{1})/sqrt(size(pd_ct_spikes_nor{2}{1}, 1)), 'color', color(1), 'LineStyle','--');
errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes_nor{1}{4}), std(pd_ct_spikes_nor{1}{4})/sqrt(size(pd_ct_spikes_nor{1}{4}, 1)), 'color', color(2));
errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes_nor{2}{4}), std(pd_ct_spikes_nor{2}{4})/sqrt(size(pd_ct_spikes_nor{2}{4}, 1)), 'color', color(2), 'LineStyle','--');
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('normalized response')
legend('ON transient control', 'ON transient HEX', 'OFF transient slow control', 'OFF transient slow HEX')

figure
errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes_nor{1}{2}), std(pd_ct_spikes_nor{1}{2})/sqrt(size(pd_ct_spikes_nor{1}{2}, 1)), 'color', color(1));
hold on
errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes_nor{2}{2}), std(pd_ct_spikes_nor{2}{2})/sqrt(size(pd_ct_spikes_nor{2}{2}, 1)), 'color', color(1), 'LineStyle','--');
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('normalized response')
legend('ON sustained control', 'ON sustained HEX')

%%
for i = 1:3
    pd_temp{i}{1} = pd_dir_spikes_nor{i}{1};
    pd_temp{i}{2} = cell2mat(pd_dir_spikes_nor{i}(2:4)');
end
figure
errorbar([80 40 20 10 5 2.5], mean(pd_temp{1}{1}), std(pd_temp{1}{1})/sqrt(size(pd_temp{1}{1}, 1)), 'color', color(1));
hold on
errorbar([80 40 20 10 5 2.5], mean(pd_temp{2}{1}), std(pd_temp{2}{1})/sqrt(size(pd_temp{2}{1}, 1)), 'color', color(1), 'LineStyle','--');
errorbar([80 40 20 10 5 2.5], mean(pd_temp{1}{2}), std(pd_temp{1}{2})/sqrt(size(pd_temp{1}{2}, 1)), 'color', color(2));
errorbar([80 40 20 10 5 2.5], mean(pd_temp{2}{2}), std(pd_temp{2}{2})/sqrt(size(pd_temp{2}{2}, 1)), 'color', color(2), 'LineStyle','--');
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('normalized response')
legend('superior control', 'superior HEX', 'other DS control', 'other DS HEX');



%% plot cell summary
cc = find(nonds_id_all == 797);
% for cc = 1:length(ds_id)
    plot_mb_raster_ctr_12(MB_n, raster_mb_nds, trial_dur, cc, nonds_id_all, 'NDF0', 3, 6, 0)
% end

%%
bin_size = 0.1;
nonds_idx_temp{2} = nonds_idx{2}([1 2 3 4 5 6 8 9 11 13]);
% nonds_idx_temp{2} = nonds_idx{2};

ct = 2;
for cc = 1:length(nonds_idx_temp{ct})
    for i = 1:3
        h = hist(raster_p_sum_nds{i}{nonds_idx_temp{ct}(cc)}{1}, 0:bin_size:5);
        t = find(h > 2*mean(h),1);
        base = mean(h(1:t-2));
        for ctr = 1:6
            h = hist(raster_p_sum_nds{i}{nonds_idx_temp{ct}(cc)}{ctr}, 0:bin_size:5);
            r_n = sum(h - base);
            pd_r_nds{i}{ct}(cc,ctr) = r_n;
        end
        pd_r_nds_nor{i}{ct}(cc,:) = pd_r_nds{i}{ct}(cc,:)/max(pd_r_nds{i}{ct}(cc,:));
    end
end

figure
for cc = 1:length(nonds_idx_temp{2})
    subplot(4,4,cc)
    hist(raster_p_sum_nds{1}{nonds_idx_temp{2}(cc)}{1}, 0:bin_size:5);
end

xx = [80 40 20 10 5 2.5];
figure
errorbar(xx, mean(pd_r_nds_nor{1}{2}), std(pd_r_nds_nor{1}{2})/sqrt(size(pd_r_nds_nor{1}{2},1)), 'r')
hold on
drug = 1;
for ct = [1 4]
    errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes_nor{drug}{ct}), std(pd_ct_spikes_nor{drug}{ct})/sqrt(size(pd_ct_spikes_nor{drug}{ct}, 1)), 'color', color(ct));
    hold on
end

errorbar(xx, mean(pd_r_nds_nor{2}{2}), std(pd_r_nds_nor{2}{2})/sqrt(size(pd_r_nds_nor{2}{2},1)), 'r--')
drug = 2;
for ct = [1 4]
    errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes_nor{drug}{ct}), std(pd_ct_spikes_nor{drug}{ct})/sqrt(size(pd_ct_spikes_nor{drug}{ct}, 1)), 'color', color(ct), 'LineStyle', '--');
    hold on
end

set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('normalized response')
legend(nonds_type([1 4 2]))

%%
figure
drug = 1;
for dir = 1:4
    errorbar([80 40 20 10 5 2.5], mean(pd_dir_spikes_nor{drug}{dir}), std(pd_dir_spikes_nor{drug}{dir})/sqrt(size(pd_dir_spikes_nor{drug}{dir}, 1)), 'color', color(dir));
    hold on
end
drug = 2;
for dir = 1:4
    errorbar([80 40 20 10 5 2.5], mean(pd_dir_spikes_nor{drug}{dir}), std(pd_dir_spikes_nor{drug}{dir})/sqrt(size(pd_dir_spikes_nor{drug}{dir}, 1)), 'color', color(dir), 'LineStyle', '--');
    hold on
end

legend(dscell_type)
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('normalized response')
title(condition{drug})

figure
drug = 1;
for ct = [1 2 4]
    errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes_nor{drug}{ct}), std(pd_ct_spikes_nor{drug}{ct})/sqrt(size(pd_ct_spikes_nor{drug}{ct}, 1)), 'color', color(ct));
    hold on
end

drug = 2;
for ct = [1 2 4]
    errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes_nor{drug}{ct}), std(pd_ct_spikes_nor{drug}{ct})/sqrt(size(pd_ct_spikes_nor{drug}{ct}, 1)), 'color', color(ct), 'LineStyle', '--');
    hold on
end
legend(nonds_type([1 2 4]))
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('normalized response')
title(condition{drug})
%%
cc = 13;
CC = find(nonds_id_all == nonds_id{1}(cc));
plot_mb_raster_ctr_12(MB_n, raster_mb_nds, trial_dur, CC, nonds_id{1}(cc), 'NDF0', 1, 3, 0)

id = 35;
CC = find(ds_id == id);
plot_mb_raster_ctr_12(MB, raster_mb, trial_dur, CC, id, 'NDF0', 1, 3, 0)

%%
ct = 4;
figure
errorbar([80 40 20 10 5 2.5], mean(pd_dir_spikes{1}{ct}), std(pd_dir_spikes{1}{ct})/sqrt(size(pd_dir_spikes{1}{ct}, 1)), 'k');
hold on
errorbar([80 40 20 10 5 2.5], mean(pd_dir_spikes{2}{ct}), std(pd_dir_spikes{2}{ct})/sqrt(size(pd_dir_spikes{2}{ct}, 1)), 'r');
errorbar([80 40 20 10 5 2.5], mean(pd_dir_spikes{3}{ct}), std(pd_dir_spikes{3}{ct})/sqrt(size(pd_dir_spikes{3}{ct}, 1)), 'b');
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('normalized response')
legend('control', 'HEX', 'wash');


ct = 4;
figure
errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes{1}{ct}), std(pd_ct_spikes{1}{ct})/sqrt(size(pd_ct_spikes{1}{ct}, 1)), 'k');
hold on
errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes{2}{ct}), std(pd_ct_spikes{2}{ct})/sqrt(size(pd_ct_spikes{2}{ct}, 1)), 'r');
errorbar([80 40 20 10 5 2.5], mean(pd_ct_spikes{3}{ct}), std(pd_ct_spikes{3}{ct})/sqrt(size(pd_ct_spikes{3}{ct}, 1)), 'b');
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('normalized response')
legend('control', 'HEX', 'wash');

%% fit ds
for ct = 1:4
    for drug = 1:3
        figure
        set(gcf, 'Position', [1 1 900 800])
        for cc = 1:size(pd_dir_spikes{drug}{ct}, 1)

            ydata = pd_dir_spikes{drug}{ct}(cc, :);
            xdata = log10([80 40 20 10 5 2.5]);
            [f, G] = fit_mm(xdata, ydata, 'upper', [100, 100, log10(300)]);
            fit_all{drug}{ct}{cc} = f;
            G_all{drug}{ct}{cc} = G;

            x = linspace(min(xdata), max(xdata), 100);
            y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a);

            subplot(5,5,cc)
            plot(xdata, ydata)
            hold on
            plot(x, y)

            sigma{drug}{ct}(cc) = 10^f.sigma;
            ymax{drug}{ct}(cc) = f.ymax;
        end
    end
end

% exclude sigma outlier (2 std)
for drug = 1:3
    for ct = 1:4
        notdone = 1;
        sigma_temp = sigma{drug}{ct};
        while notdone
            a = length(sigma_temp);
            sigma_temp(sigma_temp > std(sigma_temp)*2+mean(sigma_temp)) = [];
            b = length(sigma_temp);
            if a == b
                notdone = 0;
                sigma{drug}{ct} = sigma_temp;
            end
        end
        sigma_mean(drug, ct) = mean(sigma{drug}{ct});
        sigma_ste(drug,ct) = std(sigma{drug}{ct})/sqrt(length(sigma{drug}{ct}));
    end
end

% exclude ymax outlier (2 std)
for drug = 1:3
    for ct = 1:4
        notdone = 1;
        ymax_temp = ymax{drug}{ct};
        while notdone
            a = length(ymax_temp);
            ymax_temp(ymax_temp > std(ymax_temp)*2+mean(ymax_temp)) = [];
            b = length(ymax_temp);
            if a == b
                notdone = 0;
                ymax{drug}{ct} = ymax_temp;
            end
        end
        ymax_mean(drug, ct) = mean(ymax{drug}{ct});
        ymax_ste(drug,ct) = std(ymax{drug}{ct})/sqrt(length(ymax{drug}{ct}));
    end
end

figure
for drug = 1:3
    for ct = 1:4
      h(drug) = scatter(ct*ones(1,length(sigma{drug}{ct}))+0.2*(drug-2), sigma{drug}{ct}, color(drug));
      hold on
    end
end
set(gca, 'XTickLabel', [{''} dscell_type])
ylabel('half-maximum contrast')
legend([h(1), h(2), h(3)],'control', 'Hex', 'wash')

figure
for drug = 1:3
    for ct = 1:4
      h(drug) = scatter(ct*ones(1,length(ymax{drug}{ct}))+0.2*(drug-2), ymax{drug}{ct}, color(drug));
      hold on
    end
end
set(gca, 'XTickLabel', [{''} dscell_type])
ylabel('ymax')
legend([h(1), h(2), h(3)],'control', 'Hex', 'wash')

% plot sigma mean
ct = {'superior', 'anterior', 'inferior', 'posterior'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = ct;
model_series = sigma_mean';
model_error = sigma_ste';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('half-maximum contrast')
legend('control','HEX', 'wash');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

% plot ymax mean
ct = {'superior', 'anterior', 'inferior', 'posterior'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = ct;
model_series = ymax_mean';
model_error = ymax_ste';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('ymax')
legend('control','HEX', 'wash');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

%% fit nds
for ct = [1 2 4]
    for drug = 1:3
        figure
        for cc = 1:size(pd_ct_spikes{drug}{ct}, 1)
            ydata = pd_ct_spikes{drug}{ct}(cc, :);
            dc = min(ydata);
            ydata_fit = ydata - dc;
            xdata = log10([80 40 20 10 5 2.5]);
            [f, G] = fit_mm(xdata, ydata_fit, 'upper', [100, 100, log10(100)]);
            fit_all{drug}{ct}{cc} = f;
            G_all{drug}{ct}{cc} = G;

            x = linspace(min(xdata), max(xdata), 100);
            y = f.ymax*x.^f.a./(x.^f.a + f.sigma^f.a) + dc;

            subplot(5,5,cc)
            plot(xdata, ydata)
            hold on
            plot(x, y)

            sigma_nds{drug}{ct}(cc) = 10^f.sigma;
            ymax_nds{drug}{ct}(cc) = f.ymax;
        end
    end
end

% exclude outlier (2 std)
for drug = 1:3
    for ct = 1:4
        notdone = 1;
        sigma_temp = sigma_nds{drug}{ct};
        while notdone
            a = length(sigma_temp);
            sigma_temp(sigma_temp > std(sigma_temp)*2+mean(sigma_temp)) = [];
            b = length(sigma_temp);
            if a == b
                notdone = 0;
                sigma_nds{drug}{ct} = sigma_temp;
            end
        end
        sigma_mean_nds(drug, ct) = mean(sigma_nds{drug}{ct});
        sigma_ste_nds(drug,ct) = std(sigma_nds{drug}{ct})/sqrt(length(sigma_nds{drug}{ct}));
    end
end

figure
for drug = 1:3
    for ct = 1:4
      scatter(ct*ones(1,length(sigma_nds{drug}{ct}))+0.1*drug, sigma_nds{drug}{ct}, color(drug))
      hold on
    end
end

ct = {'superior', 'anterior', 'inferior', 'posterior'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = ct;
model_series = sigma_mean_nds';
model_error = sigma_ste_nds';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('F2/F1')
legend('control','HEX', 'wash');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end


%
sigma_mean_all = [sigma_mean sigma_mean_nds(:,1)];
sigma_ste_all = [sigma_ste sigma_ste_nds(:,1)];
ct = {'superior', 'anterior', 'inferior', 'posterior', 'ON transient'};
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
xtick = ct;
model_series = sigma_mean_all';
model_error = sigma_ste_all';
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('half-maximum contrast')
legend('control','Hex', 'wash');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

for ct = 1:4
    for drug = 1:3
        for drug2 = 1:3
            [~,P{ct}(drug,drug2)] = ttest2(sigma{drug}{ct}, sigma{drug2}{ct});
        end
    end
end

%% spike timing precision
clear spike_1st_std_all spike_1st_std_all_mean spike_1st_std_ste spike_1st_std
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
drug = 1;
for ctr = 1:6
    for dir = 1:4
        spike_1st_std{ctr}{dir} = [];
        for cc = 1:length(id_dir{dir});
            idx = idx_dir{dir}(cc);
            
            if ~isempty(raster_p_sum_all{drug}{idx})
            
                spike = squeeze(raster_p_sum_all{drug}{idx}(1,1,ctr,:));
                spike_1st = [];
                for trial = 1:length(spike)
                    if ~isempty(spike{trial})
                        spike_1st = [spike_1st spike{trial}(1)];
                    end
                end
                if length(spike_1st)>=3
                    spike_1st_all{dir}{cc}{ctr} = spike_1st;
                    spike_1st_std{ctr}{dir} = [spike_1st_std{ctr}{dir} robust_std(spike_1st)];
                end
            end
        end
    end

    spike_1st_std_mean(:, ctr) = cellfun(@mean, spike_1st_std{ctr})';
    spike_1st_std_ste(:, ctr) = cellfun(@std,spike_1st_std{ctr})'/sqrt(length(spike_1st_std{ctr}{1}));

%     spike_1st_std{dir}(isnan(spike_1st_std{dir})) = 0;
%     spike_1st_std{dir}(sum(spike_1st_std{dir}, 2) == 0, :) = [];
    spike_1st_std{ctr}{1} = [];
end

spike_1st_std_all = cellfun(@cell2mat, spike_1st_std, 'UniformOutput', false);
spike_1st_std_all_mean = cellfun(@mean, spike_1st_std_all);
spike_1st_std_all_ste = cellfun(@std,spike_1st_std_all)./sqrt(cellfun(@length, spike_1st_std_all));
figure
errorbar(ctr_x(1:end-1), spike_1st_std_all_mean(1:end-1), spike_1st_std_all_ste(1:end-1)) 
set(gca, 'Xscale', 'log')
% xlim([3 400])

marker = 'sodx';
figure
for dir = 1:4
    errorbar(ctr_x, spike_1st_std_mean(dir, :), spike_1st_std_ste(dir, :), 'Marker', marker(dir), 'color', 'k', 'MarkerSize', 8)
    hold on
end

set(gca, 'Xscale', 'log')
legend(dscell_type)
xlabel('contrast %')
ylabel('std(second)')
title('spike timing precision')
% ylim([0 0.5])

%% spike timing precision (nds cells)
clear spike_1st_std_all spike_1st_std_all_mean spike_1st_std_ste spike_1st_std
drug = 1;
for ctr = 1:6
    for dir = 1:1
        spike_1st_std{ctr}{dir} = [];
        for cc = 1:length(nonds_id{1});
            idx = nonds_idx{1}(cc);
            
            if ~isempty(raster_p_sum_nds_all{drug}{idx})
            
                spike = squeeze(raster_p_sum_nds_all{drug}{idx}(1,1,ctr,:));
                spike_1st = [];
                for trial = 1:length(spike)
                    if ~isempty(spike{trial})
                        spike_1st = [spike_1st spike{trial}(1)];
                    end
                end
                if length(spike_1st)>=3
                    spike_1st_all{dir}{cc}{ctr} = spike_1st;
                    spike_1st_std{ctr}{dir} = [spike_1st_std{ctr}{dir} robust_std(spike_1st)];
                end
            end
        end
    end

    spike_1st_std_mean(:, ctr) = cellfun(@mean, spike_1st_std{ctr})';
    spike_1st_std_ste(:, ctr) = cellfun(@std,spike_1st_std{ctr})'/sqrt(length(spike_1st_std{ctr}{1}));

%     spike_1st_std{dir}(isnan(spike_1st_std{dir})) = 0;
%     spike_1st_std{dir}(sum(spike_1st_std{dir}, 2) == 0, :) = [];
%     spike_1st_std{ctr}{1} = [];
end

spike_1st_std_all = cellfun(@cell2mat, spike_1st_std, 'UniformOutput', false);
spike_1st_std_all_mean = cellfun(@mean, spike_1st_std_all);
spike_1st_std_all_ste = cellfun(@std,spike_1st_std_all)./sqrt(cellfun(@length, spike_1st_std_all));
figure
errorbar(ctr_x(1:end-1), spike_1st_std_all_mean(1:end-1), spike_1st_std_all_ste(1:end-1)) 
set(gca, 'Xscale', 'log')
% xlim([3 400])

marker = 'sodx';
figure
for dir = 1:4
    errorbar(ctr_x, spike_1st_std_mean(dir, :), spike_1st_std_ste(dir, :), 'Marker', marker(dir), 'color', 'k', 'MarkerSize', 8)
    hold on
end

set(gca, 'Xscale', 'log')
legend(dscell_type)
xlabel('contrast %')
ylabel('std(second)')
title('spike timing precision')
% ylim([0 0.5])

