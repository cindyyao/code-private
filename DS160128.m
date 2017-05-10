cd /Users/xyao/matlab/code-private/

opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
datarun{1} = load_data('/Volumes/lab/analysis/2016-01-28-0/data000-map/data000-map', opt);
datarun{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-01-28-0/stimuli/s00';
datarun{1} = load_stim(datarun{1}, 'user_defined_trigger_interval', 10);

datarun{2} = load_data('/Volumes/lab/analysis/2016-01-28-0/data001/data001', opt);
datarun{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-01-28-0/stimuli/s01';
datarun{2} = load_stim(datarun{2}, 'user_defined_trigger_interval', 10);

datarun{3} = load_data('/Volumes/lab/analysis/2016-01-28-0/data004/data004', opt);
datarun{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-01-28-0/stimuli/s04';
datarun{3} = load_stim(datarun{3}, 'user_defined_trigger_interval', 10);

%%

[NumSpikesCell, ~,StimComb] = get_spikescellstim(datarun{3},datarun{3}.cell_ids,0,1);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);

figure
hist(ds_struct.mag{1}, 50)
[x,~] = ginput;
ds_idx = find(ds_struct.mag{1} > x);
ds_id = datarun{3}.cell_ids(ds_idx);

n = 3;
[raster_dg, DG, ~, raster_p_sum, p_idx] = deal(cell(n, 1));

for i = 1:3
    [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim(datarun{i},ds_id,0, 0.05);
    DG{i} = sort_direction(dscellanalysis(NumSpikesCell, StimComb));
    raster_dg{i} = get_ds_raster(datarun{i}, ds_id);
end

i = 3;
for cc = 1:length(raster_dg{i})
    raster_dg_12{i}{cc} = raster_dg{i}{cc}(:, :, 1:3:end, :);
end

%% classify DSGC into subtypes (directions)
d = 3;
t = 1;
h = figure;
dirn = 4;
set(h, 'Position', [1 1 1080 500])
compass(DG{d}.U{t}, DG{d}.V{t});
color = 'bkrgc';

for i = 1:dirn
    hold on
    [x, y] = ginput;
    plot(x, y, color(i));

    IN = inpolygon(DG{d}.U{t}, DG{d}.V{t}, x, y);
    [~, I] = find(IN == 1);
    idx_dir{i} = I;
    id_dir{i} = ds_id(idx_dir{i});
end

%% plot ds raster
% id_quad = [1801 4743 482 361];
id_quad = [3754 4126 4666 7396];
for i = 1:1
cc = find(ds_id == id_quad(i));
% for cc = 1:length(ds_id)
    plot_ds_raster_12(DG(3), raster_dg_12(3), cc, ds_id(cc), '', 0)
% end
end
%% direction decoding
ll = 3;
% id_quad = [1801 4743 482 361];
[~,~,idx] = intersect(id_quad, ds_id, 'stable');
raster_quad = raster_dg{ll}(idx);
for i = 1:length(id_quad)
    SpikeN_quad(:,:,i) = squeeze(cellfun('length',raster_quad{i}));
end

for t = 1:datarun{ll}.stimulus.repetitions
    test_trial = squeeze(SpikeN_quad(:,t,:));
    training_trial = SpikeN_quad;
    training_trial(:,t,:) = [];
    training_avg = squeeze(mean(training_trial,2));
    for theta = 1:length(datarun{ll}.stimulus.params.DIRECTION)
        test_trial_r = repmat(test_trial(theta, :), length(datarun{ll}.stimulus.params.DIRECTION),1);
        p = exp(-(test_trial_r - training_avg).^2./(2*training_avg))./sqrt(2*pi*training_avg);
%         p = (training_avg.*exp(-training_avg./test_trial_r)).^test_trial_r./factorial(test_trial_r);
        p = prod(p, 2);
        [~,i] = max(p);
        error = abs(i-theta)*10;
        if error > 180
            error = 360-error;
        end
        theta_e(t,theta) = error;
    end
end
theta_e_o = theta_e;
% rsme = sqrt(sum(theta_e(2:end,:).^2)/11);
figure
errorbar(0:10:350, mean(theta_e(2:end,:)), std(theta_e(2:end,:))/sqrt(11))
% errorbar(0:10:350, mean(theta_e), std(theta_e)/sqrt(12))
% ylim([-5 30])
xlim([-10 360])
xlabel('degree')
ylabel('error(deg)')

% tuning curve
figure
tuning_avg = squeeze(mean(SpikeN_quad, 2));
plot(0:10:350, tuning_avg)
xlabel('degree')
ylabel('spike #')
xlim([-10 360])
mean(mean(theta_e(2:end,:)))
%%
ll = 1;
% id_quad = [1801 4743 482 361];
[~,~,idx] = intersect(id_quad, ds_id, 'stable');
raster_quad = raster_dg{ll}(idx);
for i = 1:length(id_quad)
    SpikeN_quad(:,:,i) = squeeze(cellfun('length',raster_quad{i}));
end


% a = SpikeN_quad(:,:,4);
% Spike_s_new = circshift(a,18);
% a = SpikeN_quad(:,:,2);
% Spike_s_new = circshift(a,-9);
% a = SpikeN_quad(:,:,3);
% Spike_s_new = circshift(a,9);

Spike_s = SpikeN_quad(:,:,1);
base = max(min(Spike_s));
% base = repmat(min(Spike_s), 36,1);
% base = 13;
Spike_s_new = Spike_s - base;
ratio = max(Spike_s)./max(Spike_s_new);
Spike_s_new = round(Spike_s_new.*repmat(ratio, 36,1));
Spike_s_new(Spike_s_new < 0) = 0;

SpikeN_quad_r = SpikeN_quad;
SpikeN_quad_r(:,:,1) = Spike_s_new;

figure
tuning_avg = squeeze(mean(SpikeN_quad_r, 2));
plot(0:10:350, tuning_avg)
xlabel('degree')
ylabel('spike #')



for t = 1:datarun{ll}.stimulus.repetitions
    test_trial = squeeze(SpikeN_quad_r(:,t,:));
    training_trial = SpikeN_quad_r;
    training_trial(:,t,:) = [];
    training_avg = squeeze(mean(training_trial,2));
    for theta = 1:length(datarun{ll}.stimulus.params.DIRECTION)
        test_trial_r = repmat(test_trial(theta, :), length(datarun{ll}.stimulus.params.DIRECTION),1);
        p = exp(-(test_trial_r - training_avg).^2./(2*training_avg))./sqrt(2*pi*training_avg);
%         p = (training_avg.*exp(-training_avg./test_trial_r)).^test_trial_r./factorial(test_trial_r);
        p = prod(p, 2);
        [~,i] = max(p);
        error = abs(i-theta)*10;
        if error > 180
            error = 360-error;
        end
        theta_e(t,theta) = error;
    end
end

% rsme = sqrt(sum(theta_e(2:end,:).^2)/11);
figure
errorbar(0:10:350, mean(theta_e(2:end,:)), std(theta_e(2:end,:))/sqrt(11))
hold on
errorbar(0:10:350, mean(theta_e_o(2:end,:)), std(theta_e_o(2:end,:))/sqrt(11))
xlabel('degree')
ylabel('error(deg)')
xlim([-10 360])
% errorbar(0:10:350, mean(theta_e), std(theta_e)/sqrt(12))
% ylim([-5 30])

mean(mean(theta_e(2:end,:)))

figure
plot(0:10:350, mean(theta_e_o(2:end,:)) - mean(theta_e(2:end,:)))
hold on
plot([0 350],[0 0])
xlabel('degree')
ylabel('difference of error (deg)')
xlim([-10 360])