% generate figures for 2015 Neurobiology Retreat poster
% 9/25/2015
cd /Users/xyao/matlab/code-private/figures/NeuroRetreat2015/
%% load data
% 2014-07-07-0
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);

% load data
datarun1{1} = load_data('/Analysis/xyao/2014-07-07-0/data001-map/data001-map', opt);
datarun1{1}.names.stimulus_path = '/Analysis/xyao/2014-07-07-0/stimuli/s01';
datarun1{1} = load_stim(datarun1{1}, 'user_defined_trigger_interval', 10);

datarun1{2} = load_data('/Analysis/xyao/2014-07-07-0/data002-map/data002-map', opt);
datarun1{2}.names.stimulus_path = '/Analysis/xyao/2014-07-07-0/stimuli/s02';
datarun1{2} = load_stim(datarun1{2}, 'user_defined_trigger_interval', 10);

datarun1{3} = load_data('/Analysis/xyao/2014-07-07-0/data004-map/data004-map', opt);
datarun1{3}.names.stimulus_path = '/Analysis/xyao/2014-07-07-0/stimuli/s04';
datarun1{3} = load_stim(datarun1{3}, 'user_defined_trigger_interval', 10);

datarun1{4} = load_data('/Analysis/xyao/2014-07-07-0/data005-map/data005-map', opt);
datarun1{4}.names.stimulus_path = '/Analysis/xyao/2014-07-07-0/stimuli/s05';
datarun1{4} = load_stim(datarun1{4}, 'user_defined_trigger_interval', 10);

datarun2{1} = load_data('/Analysis/xyao/2014-03-20-0/data000/data000', opt);
datarun2{1}.names.stimulus_path = '/Analysis/xyao/2014-03-20-0/stimuli/s00';
datarun2{1} = load_stim(datarun2{1}, 'user_defined_trigger_interval', 10);

datarun2{2} = load_data('/Analysis/xyao/2014-03-20-0/data001-map/data001-map', opt);
datarun2{2}.names.stimulus_path = '/Analysis/xyao/2014-03-20-0/stimuli/s01';
datarun2{2} = load_stim(datarun2{2}, 'user_defined_trigger_interval', 10);

datarun2{3} = load_data('/Analysis/xyao/2014-03-20-0/data002-map/data002-map', opt);
datarun2{3}.names.stimulus_path = '/Analysis/xyao/2014-03-20-0/stimuli/s02';
datarun2{3} = load_stim(datarun2{3}, 'user_defined_trigger_interval', 10);

datarun2{4} = load_data('/Analysis/xyao/2014-03-20-0/data003-map/data003-map', opt);
datarun2{4}.names.stimulus_path = '/Analysis/xyao/2014-03-20-0/stimuli/s03';
datarun2{4} = load_stim(datarun2{4}, 'user_defined_trigger_interval', 10);

datarun2{5} = load_data('/Analysis/xyao/2014-03-20-0/data004-map/data004-map', opt);
datarun2{5}.names.stimulus_path = '/Analysis/xyao/2014-03-20-0/stimuli/s04';
datarun2{5} = load_stim(datarun2{5}, 'user_defined_trigger_interval', 10);

datarun2{6} = load_data('/Analysis/xyao/2014-03-20-0/data005-map/data005-map', opt);
datarun2{6}.names.stimulus_path = '/Analysis/xyao/2014-03-20-0/stimuli/s05';
datarun2{6} = load_stim(datarun2{6}, 'user_defined_trigger_interval', 10);

datadg{1} = load_data('/Analysis/xyao/2015-07-03-0/data000-map/data000-map', opt);
datadg{1}.names.stimulus_path = '/Analysis/xyao/2015-07-03-0/stimuli/s00.mat';
datadg{1} = load_stim_matlab(datadg{1}, 'user_defined_trigger_interval', 10);
datadg{2} = load_data('/Analysis/xyao/2015-07-03-0/data003-map/data003-map', opt);
datadg{2}.names.stimulus_path = '/Analysis/xyao/2015-07-03-0/stimuli/s03.mat';
datadg{2} = load_stim_matlab(datadg{2}, 'user_defined_trigger_interval', 10);
datadg{3} = load_data('/Analysis/xyao/2015-07-03-0/data006-map/data006-map', opt);
datadg{3}.names.stimulus_path = '/Analysis/xyao/2015-07-03-0/stimuli/s06.mat';
datadg{3} = load_stim_matlab(datadg{3}, 'user_defined_trigger_interval', 10);
datadg{4} = load_data('/Analysis/xyao/2015-07-03-0/data009-map/data009-map', opt);
datadg{4}.names.stimulus_path = '/Analysis/xyao/2015-07-03-0/stimuli/s09.mat';
datadg{4} = load_stim_matlab(datadg{4}, 'user_defined_trigger_interval', 10);
datadg{5} = load_data('/Analysis/xyao/2015-07-03-0/data012-map/data012-map', opt);
datadg{5}.names.stimulus_path = '/Analysis/xyao/2015-07-03-0/stimuli/s12.mat';
datadg{5} = load_stim_matlab(datadg{5}, 'user_defined_trigger_interval', 10);

%% classification
load('fig1.mat')
Title = {'NDF3 SP60', 'NDF3 SP240', 'NDF0 SP60', 'NDF0 SP240'};
% DS cell example 
% 2014-07-07-0 id:3437
plot_ds_raster_(DS1, raster, 15, 3437, 7, Title, 2, 2, 0)

% non-DS cell example
% 2014-07-07-0 id:19
plot_ds_raster_(nDS1, raster_temp, 1, 19, 7, Title, 2, 2, 0)

% tuning curve comparison
dstuning = DS1{4}.rho{1}(15, :);
ndstuning = nDS1{4}.rho{1}(1, :);
figure
x = 0:45:315;
plot(x, dstuning)
hold on
plot(x, ndstuning, 'r')
xlim([0 315])

% vector sum plot
% blue: 2014-03-20-0
% red: 2015-07-03-0
[NumSpikesCell, StimComb] = get_spikescellstim(datarun2{6},datarun2{6}.cell_ids,0);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);
figure(100)
subplot(1, 2, 1)
plot(ds_struct.mag{2, 1}, ds_struct.mag{3, 1}, 'bo')
xlim([0 2.5]); ylim([0 2.5]);

[NumSpikesCell, StimComb] = get_spikescellstim(datadg{5},datadg{5}.cell_ids,0);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);
figure(100)
subplot(1, 2, 2)
plot(ds_struct.mag{4, 1}, ds_struct.mag{5, 1}, 'ro')
xlim([0 2.5]); ylim([0 2.5]);

% polar plot
figure
subplot(1, 2, 1)
compass(DS2{6}.U{3}, DS2{6}.V{3})
subplot(1, 2, 2)
compass(DS3{5}.U{4}, DS3{5}.V{4})

% PCA
load('fig2.mat')
% 2014-03-20-0
figure
subplot(1, 2, 1)
plot(scores2(:, 1), scores2(:, 3), 'o')
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
ylim([-0.4 0.3])
title('2014-03-20-0')
% 2015-07-03-0
subplot(1, 2, 2)
plot(scores3(:, 1), scores3(:, 3), 'o')
xlabel('1st Principal Component')
ylabel('3rd Principal Component')
title('2015-07-03-0')

% speed tuning curves
% 2014-03-20-0
figure
subplot(2, 2, 1)
for cc = 1:length(idx2{1})
    if ~isempty(raster2{6}{idx2{1}(cc)})
        m = MAG_all_norm2{6}(:, idx2{1}(cc));
        semilogx(v2, m, 'b')
        hold on
    end
end
xlabel('speed')
ylabel('Response')
xlim([v2(end) v3(1)])
title('ON')
subplot(2, 2, 2)
for cc = 1:length(idx2{2})
    if ~isempty(raster2{6}{idx2{2}(cc)})
        m = MAG_all_norm2{6}(:, idx2{2}(cc));
        semilogx(v2, m, 'b')
        hold on
    end
end
xlabel('speed')
ylabel('Response')
xlim([v2(end) v3(1)])
title('ON-OFF')

% 2015-07-03-0
% figure
subplot(2, 2, 3)
semilogx(v3, exciseColumn(MAG_all_norm3{2}(:, idx3{1})), 'r')
xlabel('micron/second')
ylabel('Response')
title('ON')
xlim([v2(end) v3(1)])
subplot(2, 2, 4)
semilogx(v3, exciseColumn(MAG_all_norm3{2}(:, idx3{2})), 'r')
xlabel('micron/second')
ylabel('Response')
title('ON-OFF')
xlim([v2(end) v3(1)])

% polar plot
% 2014-03-20-0
figure
subplot(2, 2, 1)
compass(DS2{6}.U{3}(idx2{1}), DS2{6}.V{3}(idx2{1}))
subplot(2, 2, 2)
compass(DS2{6}.U{3}(idx2{2}), DS2{6}.V{3}(idx2{2}))
% 2015-07-03-0
% figure
subplot(2, 2, 3)
compass(DS3{5}.U{4}(idx3{1}), DS3{5}.V{4}(idx3{1}))
subplot(2, 2, 4)
compass(DS3{5}.U{4}(idx3{2}), DS3{5}.V{4}(idx3{2}))

%% frequency doubling
%2014-03-20-0 2015-07-03-0
%on vs onoff
figure
subplot(1, 2, 1)
errorbar(v2(1:7), ratio2_mean{1,1}, ratio2_ste{1,1}, 'b')
hold on
errorbar(v2(1:7), ratio2_mean{1,2}, ratio2_ste{1,2}, 'r')
errorbar(v3(1:7), ratio3_mean{1,1}, ratio3_ste{1,1}, 'color', [0 0 0.5])
errorbar(v3(1:7), ratio3_mean{1,2}, ratio3_ste{1,2}, 'color', [0.5 0 0])

xlim([v2(7) v3(1)])
set(gca, 'XScale', 'log')

subplot(1, 2, 2)
errorbar(v2(1:7), ratio2_mean{2,1}, ratio2_ste{2,1}, 'b')
hold on
errorbar(v2(1:7), ratio2_mean{2,2}, ratio2_ste{2,2}, 'r')
errorbar(v3(1:7), ratio3_mean{2,1}, ratio3_ste{2,1}, 'color', [0 0 0.5])
errorbar(v3(1:7), ratio3_mean{2,2}, ratio3_ste{2,2}, 'color', [0.5 0 0])

xlim([v2(7) v3(1)])
set(gca, 'XScale', 'log')

%individual direction
%2014-03-20-0
figure
for ct = 1:4
    subplot(2, 2, ct)
    errorbar(v2(1:7), ratio2_dir_mean{1,1}{ct}, ratio2_dir_ste{1,1}{ct}, 'b')
    hold on
    errorbar(v2(1:7), ratio2_dir_mean{1,2}{ct}, ratio2_dir_ste{1,2}{ct}, 'r')
    set(gca, 'XScale', 'log')
    xlim([v2(7) v2(1)])
    xlabel('micron/sec')
    ylabel('F2/F1')
end
    
figure
for ct = 1:2
    subplot(1, 2, ct)
    errorbar(v2(1:7), ratio2_dir_mean{2,1}{ct}, ratio2_dir_ste{2,1}{ct}, 'b')
    hold on
    errorbar(v2(1:7), ratio2_dir_mean{2,2}{ct}, ratio2_dir_ste{2,2}{ct}, 'r')
    set(gca, 'XScale', 'log')
    xlim([v2(7) v2(1)])
    xlabel('micron/sec')
    ylabel('F2/F1')
end

%2014-07-03-0
celltype_oo = {'posterior', 'inferior', 'anterior', 'superior'};
celltype_on = {'inferior', 'anterior', 'superior'};
figure
for ct = 1:4
    subplot(2, 2, ct)
    errorbar(v3(1:7), ratio3_dir_mean{1,1}{ct}, ratio3_dir_ste{1,1}{ct}, 'b')
    hold on
    errorbar(v3(1:7), ratio3_dir_mean{1,2}{ct}, ratio3_dir_ste{1,2}{ct}, 'r')
    set(gca, 'XScale', 'log')
    xlim([v3(7) v3(1)])
    title(celltype_oo{ct})
    xlabel('micron/sec')
    ylabel('F2/F1')
end
    
figure
for ct = 1:3
    subplot(2, 2, ct)
    errorbar(v3(1:7), ratio3_dir_mean{2,1}{ct}, ratio3_dir_ste{2,1}{ct}, 'b')
    hold on
    errorbar(v3(1:7), ratio3_dir_mean{2,2}{ct}, ratio3_dir_ste{2,2}{ct}, 'r')
    set(gca, 'XScale', 'log')
    xlim([v3(7) v3(1)])
    title(celltype_on{ct})
    xlabel('micron/sec')
    ylabel('F2/F1')

end

%% speed tuning
