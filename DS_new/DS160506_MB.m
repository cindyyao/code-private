cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);


datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-05-06-0/data000-001-map/data000-001-map', opt);
datamb(1:2) = split_datarun(datarun, 1890);
datamb{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-05-06-0/stimuli/s00.mat';
datamb{1} = load_stim_matlab(datamb{1});
datamb{2}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-05-06-0/stimuli/s01.mat';
datamb{2} = load_stim_matlab(datamb{2});

datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-05-06-0/data002-003-map/data002-003-map', opt);
datamb(3:4) = split_datarun(datarun, 1880);
datamb{3}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-05-06-0/stimuli/s02.mat';
datamb{3} = load_stim_matlab(datamb{3});
datamb{4}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-05-06-0/stimuli/s03.mat';
datamb{4} = load_stim_matlab(datamb{4});

datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-05-06-0/data005-006-map/data005-006-map', opt);
datamb(5:6) = split_datarun(datarun, 1880);
datamb{5}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-05-06-0/stimuli/s05.mat';
datamb{5} = load_stim_matlab(datamb{5});
datamb{6}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2016-05-06-0/stimuli/s06.mat';
datamb{6} = load_stim_matlab(datamb{6});
% 
% datawn = load_data('/Volumes/lab/analysis/2016-05-06-0/data004-map/data004-map', opt);
% datawn = load_sta(datawn);
% datawn = get_rfs(datawn, 'all');
% 
% load('DS160506.mat')

% color = 'brgkcmy';
% plot_rf_summaries(datawn, 'DS superior', 'fit_color', color(1), 'clear', false)
% plot_rf_summaries(datawn, 'DS anterior', 'fit_color', color(2), 'clear', false)
% plot_rf_summaries(datawn, 'DS inferior', 'fit_color', color(3), 'clear', false)
% plot_rf_summaries(datawn, 'DS posterior', 'fit_color', color(4), 'clear', false)
% 
% figure
% 
% plot_rf_summaries(datawn, [3020 4219 4836 4771], 'fit_color', color(4), 'clear', false)
datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2016-05-06-0/data007-sorted/data007-sorted', opt);
datarun = load_ei(datarun, id_dir{1}, 'array_id', 1551);
%% load data
cd /Users/xyao/matlab/code-private/DS_new/
opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);


datarun = load_data('/Volumes/janacek/Analysis/2016-05-06-0/data000-001-map/data000-001-map', opt);
datamb(1:2) = split_datarun(datarun, 1890);
datamb{1}.names.stimulus_path = '/Volumes/janacek/Analysis/2016-05-06-0/stimuli/s00.mat';
datamb{1} = load_stim_matlab(datamb{1});
datamb{2}.names.stimulus_path = '/Volumes/janacek/Analysis/2016-05-06-0/stimuli/s01.mat';
datamb{2} = load_stim_matlab(datamb{2});

datarun = load_data('/Volumes/janacek/Analysis/2016-05-06-0/data002-003-map/data002-003-map', opt);
datamb(3:4) = split_datarun(datarun, 1880);
datamb{3}.names.stimulus_path = '/Volumes/janacek/Analysis/2016-05-06-0/stimuli/s02.mat';
datamb{3} = load_stim_matlab(datamb{3});
datamb{4}.names.stimulus_path = '/Volumes/janacek/Analysis/2016-05-06-0/stimuli/s03.mat';
datamb{4} = load_stim_matlab(datamb{4});

datarun = load_data('/Volumes/janacek/Analysis/2016-05-06-0/data005-006-map/data005-006-map', opt);
datamb(5:6) = split_datarun(datarun, 1880);
datamb{5}.names.stimulus_path = '/Volumes/janacek/Analysis/2016-05-06-0/stimuli/s05.mat';
datamb{5} = load_stim_matlab(datamb{5});
datamb{6}.names.stimulus_path = '/Volumes/janacek/Analysis/2016-05-06-0/stimuli/s06.mat';
datamb{6} = load_stim_matlab(datamb{6});

datawn = load_data('/Volumes/janacek/Analysis/2016-05-06-0/data004-map/data004-map', opt);
datawn = load_sta(datawn);
datawn = get_rfs(datawn, 'all');

load('DS160506.mat')

%% 
n = 6;
[raster_mb, MB, ~, raster_p_sum, p_idx] = deal(cell(n, 1));
for i = 1:n    
    [NumSpikesCell, ~,StimComb] = get_spikescellstim_mb(datamb{i},ds_id,1860,0.1);
    MB{i} = sort_direction(mbcellanalysis(NumSpikesCell, StimComb,datamb{i}));
    raster_mb{i} = get_mb_raster(datamb{i}, ds_id, 1860);
    raster_mb_all{i} = combine_repeats(raster_mb{i});
%     raster_mb{i} = get_mb_raster(datamb{i}, ds_id, datamb{i}.stimulus.triggers(865));
    for j = 1:length(raster_mb{i})
        if(idx_mb(j, ceil(i/2)))
            raster_mb{i}{j} = [];
            raster_mb_all{i}{j} = [];
        end
    end
end

close all


for i = 1:n;
    for cc = 1:length(raster_mb{i})
        raster_mb_12{i}{cc} = raster_mb{i}{cc}(:, :, 1:3:end, :,:);
        raster_mb_all_12{i}{cc} = squeeze(raster_mb_all{i}{cc}(:, :, 1:3:end));
    end
    trial_dur{i} = get_mb_trial_dur(datamb{i}, 400, 400, 0.5);
end

ctr_p = 1; % choose which params to use to calculate prefer direction indices 
ds = 5;
MAG_all_norm_mb = cell(n, 1);

for i = 1:n
    [raster_p_sum{i}, p_idx{i}, raster_p_sum_all{i}] = get_pdirection_raster(raster_mb{i}, MB{ds}.angle{ctr_p});
    [raster_n_sum{i}, n_idx{i}, raster_n_sum_all{i}] = get_ndirection_raster(raster_mb{i}, MB{ds}.angle{ctr_p});
    MAG_all_norm_mb{i} = normalize_MAG(MB{i});
    rep = datamb{i}.stimulus.repetitions;
end
%% plot ds raster
for cc = 1:1%length(ds_id)
    plot_mb_raster_ctr_12(MB([1 3 5 2 4 6]), raster_mb_12([1 3 5 2 4 6]), trial_dur, cc, ds_id(cc), '', 2,3,1)
end

%% direction decoding
LL = {'NDF 4', 'NDF 3', 'NDF 0'};
ID_Quad = [814 1743 322 3020; ...
    1891 2510 2328 5749; ...
    3153 3091 3362 4024; ...
%     7066 3605 4219 5134; ...
    4655 3813 1054 1457; ...
    7353 5658 4657 6422];

for q = 1:5
    id_quad = ID_Quad(q,:);
    figure
    set(gcf, 'Position', [1 1 1500 1000])
    for ll = 1:3
        light_level = [1 3 5];
        rep = datamb{light_level(ll)}.stimulus.repetitions;
        [~,~,idx] = intersect(id_quad, ds_id, 'stable'); 
        raster_quad = raster_mb{light_level(ll)}(idx);
        for i = 1:length(id_quad)
            SpikeN_quad(:,:,i) = squeeze(cellfun('length',raster_quad{i}));
        end

        for t = 1:rep 
            test_trial = squeeze(SpikeN_quad(:,t,:));
            training_trial = SpikeN_quad;
            training_trial(:,t,:) = [];
            training_avg = squeeze(mean(training_trial,2));
            training_std = max(squeeze(std(training_trial,[],2))-0.1, 0)+0.1;
            for theta = 1:length(datamb{light_level(ll)}.stimulus.params.DIRECTION)
                test_trial_r = repmat(test_trial(theta, :), length(datamb{light_level(ll)}.stimulus.params.DIRECTION),1);
                p = exp(-(test_trial_r - training_avg).^2./(2*training_std.^2))./sqrt(2*pi*training_std.^2);
        %         p = exp(-(test_trial_r - training_avg).^2./(2*training_avg))./sqrt(2*pi*training_avg);
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
        theta_e_o{q,ll} = theta_e;
        subplot(2,3,ll+3)
        % rsme = sqrt(sum(theta_e(2:end,:).^2)/11);
        % errorbar(0:10:350, mean(theta_e(2:end,:)), std(theta_e(2:end,:))/sqrt(rep-1))
        errorbar(0:10:350, mean(theta_e), std(theta_e)/sqrt(rep))
        % ylim([-5 30])
        xlim([-10 360])
        xlabel('degree')
        ylabel('error(deg)')
%         ylim([0 150])

        % tuning curve
        XX = [0:10:350];
        subplot(2,3,ll)
        tuning_avg = squeeze(mean(SpikeN_quad, 2));
        tuning_ste = squeeze(std(SpikeN_quad, [], 2)/sqrt(rep));
        errorbar(repmat(XX, 4,1)', tuning_avg, tuning_ste)
        xlabel('degree')
        ylabel('spike #')
        xlim([-10 360])
        error_mean(q,ll) = mean(mean(theta_e));
        title(LL{ll})
%         ylim([0 30])
    end
end

figure
x = [0 1 4];
for i = 1:5
    plot(x, error_mean(i,:),'color',[0.5 0.5 0.5], 'Marker', 'o')
    hold on
end
errorbar(x, mean(error_mean), std(error_mean)/sqrt(5), 'color', 'k', 'Linewidth', 2, 'Marker', 's'); 
xlim([-0.5 4.5])
xlabel('Log(luminance)(R*/rod/s)')
ylabel('mean error (degree)')
%%
for q = 3:3
    id_quad = ID_Quad(q,:);
    figure
    set(gcf, 'Position', [1 1 1500 1000])
    for ll = 1:3
        [~,~,idx] = intersect(id_quad, ds_id, 'stable');
        raster_quad = raster_mb{light_level(ll)}(idx);
        for i = 1:length(id_quad)
            SpikeN_quad(:,:,i) = squeeze(cellfun('length',raster_quad{i}));
        end


%         a = SpikeN_quad(:,:,4);
%         Spike_s_new = circshift(a,9);
        a = SpikeN_quad(:,:,2);
        Spike_s_new = circshift(a,-9);
%         a = SpikeN_quad(:,:,3);
%         Spike_s_new = circshift(a,18);

%         Spike_s = SpikeN_quad(:,:,1);
%         base = max(min(Spike_s));
% %         base = median(min(Spike_s));
%         Spike_s_new = Spike_s - base;
%         ratio = max(Spike_s)./max(Spike_s_new);
%         Spike_s_new = round(Spike_s_new.*repmat(ratio, 36,1));
%         Spike_s_new(Spike_s_new < 0) = 0;

        SpikeN_quad_r = SpikeN_quad;
        SpikeN_quad_r(:,:,1) = Spike_s_new;

        subplot(2,3,ll)
        tuning_avg = squeeze(mean(SpikeN_quad_r, 2));
        plot(0:10:350, tuning_avg)
        xlabel('degree')
        ylabel('spike #')
        ylim([0 15])



        for t = 1:datamb{light_level(ll)}.stimulus.repetitions
            test_trial = squeeze(SpikeN_quad_r(:,t,:));
            training_trial = SpikeN_quad_r;
            training_trial(:,t,:) = [];
            training_avg = squeeze(mean(training_trial,2));
            training_std = max(squeeze(std(training_trial,[],2))-0.1, 0)+0.1;
            for theta = 1:length(datamb{light_level(ll)}.stimulus.params.DIRECTION)
                test_trial_r = repmat(test_trial(theta, :), length(datamb{light_level(ll)}.stimulus.params.DIRECTION),1);
                p = 1./(sqrt(2*pi)*training_avg).*exp(-(test_trial_r - training_avg).^2./(2*training_std.^2))./sqrt(2*pi*training_std.^2);
        %         p = exp(-(test_trial_r - training_avg).^2./(2*training_avg))./sqrt(2*pi*training_avg);
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
        theta_e_adjusted{q,ll} = theta_e;
        % rsme = sqrt(sum(theta_e(2:end,:).^2)/11);
        subplot(2,3,ll+3)
        errorbar(0:10:350, mean(theta_e_o{q,ll}), std(theta_e_o{q,ll})/sqrt(rep), 'k')
        hold on
        errorbar(0:10:350, mean(theta_e_adjusted{q,ll}), std(theta_e_adjusted{q,ll})/sqrt(rep), 'r')
        xlabel('degree')
        ylabel('error(deg)')
        xlim([-10 360])
        % errorbar(0:10:350, mean(theta_e), std(theta_e)/sqrt(12))
        ylim([0 120])

        error_mean_adjusted(q,ll) = mean(mean(theta_e));

%         figure
%         plot(0:10:350, mean(theta_e_o) - mean(theta_e))
%         hold on
%         plot([0 350],[0 0])
%         xlabel('degree')
%         ylabel('difference of error (deg)')
%         xlim([-10 360])
    end
end

figure
% x = [0 1 4];
% for i = 1:5
%     h1 = plot(x, error_mean(i,:),'ko-.');
%     hold on
%     h2 = plot(x, error_mean_adjusted(i,:),'ro-.');
% end
errorbar(x, mean(error_mean), std(error_mean)/sqrt(5), 'ks-', 'Linewidth', 2); 
hold on
errorbar(x, mean(error_mean_adjusted), std(error_mean_adjusted)/sqrt(5), 'rs-', 'Linewidth', 2); 

xlim([-0.5 4.5])
xlabel('Log(luminance)(R*/rod/s)')
ylabel('mean error (degree)')
legend('control', 'subsituted')

%% direction decoding
for q = 1:5
    id_quad = ID_Quad(q,:);
    figure
    set(gcf, 'Position', [1 1 1500 1000])
    for pp = 1:2:10
    for ll = 1:3
        light_level = [1 3 5];
        rep = datamb{light_level(ll)}.stimulus.repetitions;
        [~,~,idx] = intersect(id_quad, ds_id, 'stable'); 
        raster_quad = raster_mb{light_level(ll)}(idx);
        for i = 1:length(id_quad)
            SpikeN_quad(:,:,i) = round(squeeze(cellfun('length',raster_quad{i}))*pp/10);
        end

        for t = 1:rep 
            test_trial = round(squeeze(SpikeN_quad(:,t,:)));
            training_trial = SpikeN_quad;
            training_trial(:,t,:) = [];
            training_avg = squeeze(mean(training_trial,2));
            training_std = max(squeeze(std(training_trial,[],2))-0.1, 0)+0.1;
            for theta = 1:length(datamb{light_level(ll)}.stimulus.params.DIRECTION)
                test_trial_r = repmat(test_trial(theta, :), length(datamb{light_level(ll)}.stimulus.params.DIRECTION),1);
                p = exp(-(test_trial_r - training_avg).^2./(2*training_std.^2))./sqrt(2*pi*training_std.^2);
        %         p = exp(-(test_trial_r - training_avg).^2./(2*training_avg))./sqrt(2*pi*training_avg);
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
        theta_e_o{q,ll,pp} = theta_e;
        subplot(1,3,ll)
        errorbar(0:10:350, mean(theta_e), std(theta_e)/sqrt(rep))
        hold on
        legend('10%', '30%', '50%', '70%', '90%')
        xlabel('degree')
        ylabel('error(deg)')
        title(LL{ll})

    end
    
    end
end

%% noise
color = 'brgk';
for ll = 1:3;
    light_level = [1 3 5];
    rep = datamb{light_level(ll)}.stimulus.repetitions;
    for ct = 1:4
        [SpikeN_temp{ct}, SpikeN_mean_temp{ct}, SpikeN_std_temp{ct}] = deal(cell(length(id_dir{ct}), 1));
        for cc = 1:length(id_dir{ct})
            raster_temp = raster_mb{light_level(ll)}{idx_dir{ct}(cc)};
            if ~isempty(raster_temp)
                SpikeN_temp{ct}{cc} = squeeze(cellfun('length',raster_temp));
                SpikeN_mean_temp{ct}{cc} = mean(SpikeN_temp{ct}{cc},2);
                SpikeN_std_temp{ct}{cc} = std(SpikeN_temp{ct}{cc},[],2);
            end
        end
        SpikeN_mean_all{ll}{ct} = cell2mat(SpikeN_mean_temp{ct});
        SpikeN_std_all{ll}{ct} = cell2mat(SpikeN_std_temp{ct});
    end
    SpikeN{ll} = SpikeN_temp;
    SpikeN_mean{ll} = SpikeN_mean_temp;
    SpikeN_std{ll} = SpikeN_std_temp;
end

LL = {'NDF 4', 'NDF 3', 'NDF 0'};
for ll = 1:1
    figure
    for ct = 1:4
        for cc = 1:length(id_dir{ct})
            if ~isempty(SpikeN_mean{ll}{ct}{cc})
                h(ll) = plot(SpikeN_mean{ll}{ct}{cc}, SpikeN_std{ll}{ct}{cc}.^2, [color(ct) 'o']);
                hold on
            end
    %         pause
        end
    end
    title(LL{ll})
    xlabel('mean (spike #)')
    ylabel('variance (spike #)')

end
% legend(h, 'NDF 4', 'NDF 3', 'NDF 0')
xlabel('mean (spike #)')
ylabel('variance (spike #)')

% fit with var = a*mean^b
g = fittype('a*x^b');
xx = 0:0.1:20;
for ll = 1:1
    figure
    for ct = 1:4
        x = SpikeN_mean_all{ll}{ct};
        y = SpikeN_std_all{ll}{ct}.^2;
        f = fit(x,y,g, 'upper', [8, 2], 'lower', [0.2 0.2]);
        plot(x,y,[color(ct) 'o'])
        hold on
        plot(xx, f(xx), color(ct));
    end
end

figure
max_r = [20 30 40];
LL = {'NDF 4', 'NDF 3', 'NDF 0'};
for ll = 1:3
    subplot(1,3,ll)
    xx = 0:0.1:max_r(ll);
    x = cell2mat(SpikeN_mean_all{ll}(2:4)');
    y = cell2mat(SpikeN_std_all{ll}(2:4)').^2;
    f = fit(x,y,g, 'upper', [8, 2], 'lower', [0.2 0.2]);
    plot(x,y,'o','color', [1 0.5 0.5])
    hold on
    plot(xx, f(xx), 'r');
    
    x = SpikeN_mean_all{ll}{1};
    y = SpikeN_std_all{ll}{1}.^2;
    f = fit(x,y,g, 'upper', [8, 2], 'lower', [0.2 0.2]);
    plot(x,y,'o','color', [1 1 1]*0.5)
    hold on
    plot(xx, f(xx), 'k');
    xlabel('spike count mean')
    ylabel('spike count variance')
    title(LL{ll})
end


% fit with var = a*mean
g = fittype('a*x');
figure
LL = {'NDF 4', 'NDF 3', 'NDF 0'};
for ll = 1:3
    subplot(1,3,ll)
    xx = 0:0.1:35;
    x = cell2mat(SpikeN_mean_all{ll}(2:4)');
    y = cell2mat(SpikeN_std_all{ll}(2:4)').^2;
    f = fit(x,y,g);
    plot(x,y,'o','color', [1 0.5 0.5])
    hold on
    plot(xx, f(xx), 'r');
    
    x = SpikeN_mean_all{ll}{1};
    y = SpikeN_std_all{ll}{1}.^2;
    f = fit(x,y,g);
    plot(x,y,'o','color', [1 1 1]*0.5)
    hold on
    plot(xx, f(xx), 'k');
    xlabel('spike count mean')
    ylabel('spike count variance')
    title(LL{ll})
    ylim([0 35])
    xlim([0 35])
end

%% correlation coefficient
% id = [3020 4219 4836 4771];
id = [2467 2510];
[~, idx] = intersect(ds_id, id);
ll = 5;
quad = raster_mb{ll}(idx);
for i = 1:2
    quad{i} = cellfun(@length, squeeze(quad{i}));
end

figure
a = 1;
b = 2;
for i = 1:36
    subplot(6,6,i)
    plot(quad{a}(i,:), quad{b}(i,:), 'o')
    spike_a = quad{a}(i,:);
    spike_b = quad{b}(i,:);
%     Covari = cov(spike_a,spike_b);
    temp = corrcoef(spike_a',spike_b');
    coeff(i) = temp(1,2);
    meanspike(i) = sqrt(mean(spike_a)*mean(spike_b));
end
    
figure
plot(coeff, 'b')
hold on
plot(meanspike/max(meanspike), 'r')

%%
figure
for i = 1:4
    idx = [2 5 8 11];
    subplot(4,1,i)
    spike_t = raster_mb_12{3}{1}{1,1,idx(i),1,1};
    X = [spike_t'; spike_t'];
    Y = [zeros(1, length(spike_t)); ones(1, length(spike_t))];
    line(X,Y, 'color', 'k')
end

color = 'brgk';
figure
for i = 1:4
    idx = [3 6 9 12];
    subplot(4,1,i)
    for j = 1:3
        spike_t = raster_mb_12{3}{1}{1,1,idx(i),1,j+3};
        X = [spike_t'; spike_t'];
        Y = [ones(1, length(spike_t))*(j-0.8); ones(1, length(spike_t))*j];
        line(X,Y, 'color', color(i))
    end
    xlim([0.6 1.3])
    ylim([0 3])
end

%% 
for i = 1:4
    x = [0:360]+90*(i-1);
    y = gaussmf(x,[50 180+90*(i-1)]);
    plot(x,y, 'color', color(i))
    hold on
    x = [0:10:360]+90*(i-1);
    y = gaussmf(x,[50 180+90*(i-1)]);
    errorbar(x,y,ones(length(x),1)*0.05, [color(i) 'o'])
end

x = 90:450;
y = gaussmf(x,[50 270]);
plot(x,y)
x = 180:540;
y = gaussmf(x,[50 360]);
plot(x,y)
x = 270:630;
y = gaussmf(x,[50 450]);
plot(x,y)

%% get background firing
T = 1865;
for dir = 1:4
    for cc = 1:length(id_dir{dir})
        for i = 1:6
            if ~(idx_mb(idx_dir{dir}(cc), ceil(i/2))) && ~isempty(find(datamb{i}.cell_ids == id_dir{dir}(cc), 1))
                trial_dur = mean(diff(datamb{i}.stimulus.triggers));
                nullspike{dir}{i}(cc) = min(MB{i}.RHO{1}(idx_dir{dir}(cc), :))/trial_dur;
                idx = get_cell_indices(datamb{i}, id_dir{dir}(cc));
                bkgndspike{dir}{i}(cc) = length(datamb{i}.spikes{idx}(datamb{i}.spikes{idx} > 1865 & datamb{i}.spikes{idx} < 1880))/15;
            end
        end
    end
end
    

%% fit tuning curves
idx_mb_all = idx_mb(:, 1)|idx_mb(:, 2)|idx_mb(:, 3);

for dir = 1:4
    for i = 1:6
        CC = 1;
        for cc = 1:length(id_dir{dir})
            if ~idx_mb_all(idx_dir{dir}(cc)) && ~isempty(find(datamb{i}.cell_ids == id_dir{dir}(cc), 1))
                xdata = MB{i}.theta{1}(idx_dir{dir}(cc), :);
                ydata = MB{i}.RHO{1}(idx_dir{dir}(cc), :);
                [f, g] = fit_cos(xdata, ydata);
                Ymax{dir}{i}(CC) = f.ymax;
                Phi{dir}{i}(CC) = f.phi;
                Alpha{dir}{i}(CC) = f.alpha;
                B{dir}{i}(CC) = f.b;
                CC = CC + 1;
                
                yfit = f.ymax.*(0.5+0.5*cos(xdata+f.phi)).^f.alpha+f.b;
%                 figure(1)
%                 plot(xdata, ydata, 'b')
%                 hold on
%                 plot(xdata, yfit, 'r')
% %                 pause
%                 close(1)

            end
        end
    end
end


fittingDs = struct();
fittingDs.ymax = Ymax;
fittingDs.alpht = Alpha;
fittingDs.bb = B;
fittingDs.phi = Phi;

%% 
for i = 1:4
    for j = 1:6
        temp = MB{j}.RHO{1}(idx_dir{i}, :);
        isMB = sum(temp') ~= 0;
        temp = temp(isMB, :);
        dstuning{j}(i, :) = mean(temp);
    end
end

xdata = MB{1}.theta{1}(1, :);
for i = 1:4
    for j = 1:6
        ydata = dstuning{j}(i, :);
        [f, g] = fit_cos(xdata, ydata);
        Ymax{j}(i) = f.ymax;
        Phi{j}(i) = f.phi;
        Alpha{j}(i) = f.alpha;
        B{j}(i) = f.b;

        yfit = f.ymax.*(0.5+0.5*cos(xdata+f.phi)).^f.alpha+f.b;
        figure(1)
        plot(xdata, ydata, 'b')
        hold on
        plot(xdata, yfit, 'r')
    end
end


%%
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'SR', 'wash'};
step_size = 100;
for ll = 1:6
    CC = 1;
    for cc = 1:length(ds_id)
        if ~isempty(raster_p_sum{ll}{cc})
            a = raster_p_sum{ll}{cc}{1};
            hist_temp = hist(a, xx);
            [max_p, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
            response_p{ll}(CC) = max_p/datamb{ll}.stimulus.repetitions;
            CC = CC + 1;
        end
    end
    response_p_norm{ll} = response_p{ll}./repmat(max(response_p{ll}, [], 2), 1, size(response_p{ll},2));
end
response_p = reshape(response_p, 2,3);
response_p_mean = cellfun(@mean, response_p);
response_p_ste = cellfun(@std, response_p)./sqrt(cellfun(@length, response_p));

color = 'brgk';
figure
for ll = 1:3
    errorbar([2 1], response_p_mean(:,ll), response_p_ste(:,ll))
    hold on
end


FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1000, 500]);
ctr = {'20 %', '80 %'};

xtick = ctr;
model_series = flipud(response_p_mean);   
model_error = flipud(response_p_ste);
h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',xtick)
ylabel('spike count')
legend('NDF4','NDF3', 'NDF0', 'location', 'northeast');
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

%% tuning curves

%subtypes
T = 1;
p_direction = MB{5}.angle{T}';
xx = 0:pi/6:11*pi/6;
xx = repmat(xx, length(ds_id), 1) - repmat(p_direction, 1, 12);
xx(xx>pi) = xx(xx>pi)-2*pi;
xx(xx<-pi) = xx(xx<-pi)+2*pi;

clear rho_mb_mean rho_mb_ste dsi_mb_mean dsi_mb_ste RHO_mb_mean RHO_mb_ste
for d = 1:6
    for i = 1:4
        rho_mb{d}{i} = [];
        RHO_mb{d}{i} = [];
        dsi_mb{d}{i} = [];
        for cc = 1:length(idx_dir{i})
            if ~idx_mb(idx_dir{i}(cc), ceil(d/2)) && sum(MB{d}.rho{T}(idx_dir{i}(cc), :))>0
            [xsort, seq] = sort(xx(idx_dir{i}(cc), :));
            y_temp = MB{d}.rho{T}(idx_dir{i}(cc), 1:3:end);
            Y_TEMP = MB{d}.RHO{T}(idx_dir{i}(cc), 1:3:end);
            rho_mb{d}{i} = [rho_mb{d}{i}; y_temp(seq)];
            RHO_mb{d}{i} = [RHO_mb{d}{i}; Y_TEMP(seq)];
            dsi_mb{d}{i} = [dsi_mb{d}{i}; MB{d}.dsindex{T}(idx_dir{i}(cc))];
            end
        end
        rho_mb_mean{d}(i, :) = mean(rho_mb{d}{i});
        rho_mb_ste{d}(i, :) = std(rho_mb{d}{i})/sqrt(size(rho_mb{d}{i}, 1));
        RHO_mb_mean{d}(i, :) = mean(RHO_mb{d}{i});
        RHO_mb_ste{d}(i, :) = std(RHO_mb{d}{i})/sqrt(size(RHO_mb{d}{i}, 1));
        dsi_mb_mean{d}(i) = mean(dsi_mb{d}{i});
        dsi_mb_ste{d}(i) = std(dsi_mb{d}{i})/sqrt(length(dsi_mb{d}{i}));
    end
    RHO_mb_mean_all{d}(1,:) = RHO_mb_mean{d}(1, :);
    RHO_mb_mean_all{d}(2,:) = mean(cell2mat(RHO_mb{d}(2:4)'));
    RHO_mb_ste_all{d}(1,:) = RHO_mb_ste{d}(1, :);
    RHO_mb_ste_all{d}(2,:) = std(cell2mat(RHO_mb{d}(2:4)'))/sqrt(size(cell2mat(RHO_mb{d}(2:4)'), 1));
end


RHO_mb_mean = flipud(reshape(RHO_mb_mean, 2, 3));
RHO_mb_ste = flipud(reshape(RHO_mb_ste, 2, 3));

RHO_mb_mean_all = flipud(reshape(RHO_mb_mean_all, 2, 3));
RHO_mb_ste_all = flipud(reshape(RHO_mb_ste_all, 2, 3));

theta = linspace(-180, 180, 13);
theta = theta(2:end);
ctr = 1;
color = 'brgk';

figure
for dir = 1:2
    subplot(1,2,dir)
    for ll = 1:3
        errorbar(theta, RHO_mb_mean_all{ctr, ll}(dir, :), RHO_mb_ste_all{ctr, ll}(dir, :), color(ll))
        hold on
    end
    ylim([0 25])
end

%% direction tuning curve (sliding window)

id_dir{2} = cell2mat(id_dir(2:4));
idx_dir{2} = cell2mat(idx_dir(2:4));
id_dir = id_dir(1:2);
idx_dir = idx_dir(1:2);

idx_mb_all = prod(~idx_mb');
for dir = 1:2
    id_dir_mb{dir} = id_dir{dir}(logical(idx_mb_all(idx_dir{dir})));
    idx_dir_mb{dir} = idx_dir{dir}(logical(idx_mb_all(idx_dir{dir})));
end


clear Max_i
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'SR', 'SR+TPMPA', 'wash'};
step_size = 50;
for ll = 1:6
    for dir = 1:2
        for theta = 1:12
            CC = 1;
            for cc = 1:length(idx_dir{dir})
                if ~isempty(raster_mb_all_12{ll}{idx_dir{dir}(cc)})
                    a = raster_mb_all_12{ll}{idx_dir{dir}(cc)}{theta};
                    hist_temp = hist(a, xx);
                    [max_fr, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
                    response_max{ll}{dir}(CC, theta) = max_fr/datamb{ll}.stimulus.repetitions;
                    CC = CC + 1;
                end
            end
            response_max_norm{ll}{dir}(:, theta, :) = squeeze(response_max{ll}{dir}(:, theta, :))./repmat(squeeze(max(response_max{1}{dir}(:, theta, :), [], 1)), 1, size(response_max{ll}{dir},1))';
        end
    end
end



for dir = 1:2
    for cc = 1:length(idx_dir_mb{dir})
        pindex{dir}(cc) = get_pindex(response_max{5}{dir}(cc, :));
        for ll = 1:6
            response_max{ll}{dir}(cc, :) = circshift(response_max{ll}{dir}(cc, :), [0, 6-pindex{dir}(cc)]);
            response_max_norm{ll}{dir}(cc, :) = response_max{ll}{dir}(cc, :)/max(response_max{ll}{dir}(cc, :));
        end
    end
end

for ll = 1:6
    for dir = 1:2
        response_max_mean{ll}{dir} = mean(response_max{ll}{dir});
        response_max_ste{ll}{dir} = std(response_max{ll}{dir})./sqrt(size(response_max{ll}{dir}, 1));
        response_max_norm_mean{ll}{dir} = nanmean(response_max_norm{ll}{dir});
        response_max_norm_ste{ll}{dir} = nanstd(response_max_norm{ll}{dir})./sqrt(size(response_max_norm{ll}{dir}, 1));
    end
end

response_max_mean = flipud(reshape(response_max_mean, 2, 3));
response_max_ste = flipud(reshape(response_max_ste, 2, 3));

theta = linspace(-180, 180, 13);
theta = theta(2:end);
color = 'brgkcmy';
T = 2;

figure
for dir = 1:2
    subplot(1, 2, dir)
    for ll = 1:3
        errorbar(theta, response_max_mean{T, ll}{dir}, response_max_ste{T, ll}{dir}, color(ll))
        hold on
    end
    legend('NDF 4', 'NDF 3', 'NDF 0')
    xlabel('direction (deg)')
    ylabel('spike count')
    ylim([0 20])
end


figure
for dir = 1:2
    for ctr = 1:7
        subplot(2, 7, 7 * (dir - 1) + ctr);
        for drug = 1:2
            errorbar(theta, response_max_norm_mean{drug}{dir}{ctr}, response_max_norm_ste{drug}{dir}{ctr}, 'color', color(drug));
            hold on
        end
        xlim([-3.5 3])
        if dir == 1
            title(contrast{ctr})
            if ctr == 7
                legend('control', 'SR')
            end
        end
        if (ctr == 1)
            ylabel(dscell_type{dir});
        end
        if (dir == 4)
            xlabel('direction')
        end
    end
end

figure
for dir = 1:2
    for ctr = 1:7
        subplot(2, 7, 7 * (dir - 1) + ctr);
        for drug = 1:2
            errorbar(theta, response_max_mean{drug}{dir}{ctr}, response_max_ste{drug}{dir}{ctr}, 'color', color(drug));
            hold on
        end
        xlim([-3.5 3])
        if dir == 1
            title(contrast{ctr})
            if ctr == 7
                legend('control', 'SR')
            end
        end
        if (ctr == 1)
            ylabel(dscell_type{dir});
        end
        if (dir == 4)
            xlabel('direction')
        end
    end
end

figure
for drug = 1:4
    for ctr = 1:7
        subplot(4, 7, 7 * (drug - 1) + ctr);
        for dir = 1:2
            errorbar(theta, response_max_norm_mean{drug}{dir}{ctr}, response_max_norm_ste{drug}{dir}{ctr}, 'color', color(dir));
            hold on
        end
        xlim([-3.5 3])
        if drug == 1
            title(contrast{ctr})
            if ctr == 7
                legend('superior', 'others')
            end
        end
        if (ctr == 1)
            ylabel(condition{drug});
        end
        if (drug == 4)
            xlabel('direction')
        end
        ylim([0 1])
    end
end

for dir = 1:2
figure
for drug = 1:4
    for ctr = 1:7
        subplot(4, 7, 7 * (drug - 1) + ctr);
        
        for cc = 1:size(response_max_norm{1}{1}, 1)
            if max(response_max_norm{drug}{dir}(cc, :, ctr)) == 1
                plot(theta, response_max_norm{drug}{dir}(cc, :, ctr), 'color', color(dir));
                hold on
            end
        end
        xlim([-3.5 3])
        if drug == 1
            title(contrast{ctr})
        end
        if (ctr == 1)
            ylabel(condition{drug});
        end
        if (drug == 4)
            xlabel('direction')
        end
        ylim([0 1])
    end
end
end

%% neighboring pairs
pos = datarun.ei.position;
mode = 'neg';
neighbors = [];
ct = 1;
for cc1 = 1:length(id_dir{ct})
    for cc2 = cc1+1:length(id_dir{ct})
        id1 = id_dir{ct}(cc1);
        idx1 = get_cell_indices(datarun, id1);
        ei1 = datarun.ei.eis{idx1};
        com1 = ei_com_xy(ei1, pos, 30*3, mode);
        id2 = id_dir{ct}(cc2);
        idx2 = get_cell_indices(datarun, id2);
        ei2 = datarun.ei.eis{idx2};
        com2 = ei_com_xy(ei2, pos, 30*3, mode);
        if pdist([com1;com2]) < 150
            neighbors = [neighbors; id1 id2];
        end
    end
end

cn = 0;
ct = 1;
celltype = {'superior', 'anterior', 'inferior', 'posterior'};
coms = [];
for cc = 1:length(id_dir{ct})
    id = id_dir{ct}(cc);
    idx = get_cell_indices(datarun, id);
    ei = datarun.ei.eis{idx};
    com = ei_com_xy(ei, pos, 30*3, mode);
    coms = [coms; com];
end

corner_i = [4 126 195 264 386 455 4];
corner_position = datarun.ei.position(corner_i, :);
figure
for cc = 1:length(id_dir{ct})
    plot(coms(cc, 1), coms(cc, 2),'ko')
    hold on
    text(coms(cc, 1)+5, coms(cc, 2)+5, num2str(id_dir{ct}(cc)), 'FontSize', 10)
    
end

% for cp = 1:size(corr_cells)
%     idx1 = find(id_dir{1} == corr_cells(cp, 1));
%     idx2 = find(id_dir{1} == corr_cells(cp, 2));
%     plot([coms(idx1, 1), coms(idx2, 1)], [coms(idx1, 2), coms(idx2, 2)], 'k');
% end
% 
plot(corner_position(:, 1), corner_position(:, 2), 'color', [.5 .5 .5])
axis off
title(celltype{ct})

% ct = 1;
cp_i = [];
for c1 = 1:length(id_dir{ct})-1
    for c2 = c1+1:length(id_dir{ct})
        if norm([coms(c1, :) - coms(c2, :)]) < 150
            cn = cn + 1;
            cp_i = [cp_i; c1 c2];
        end
    end
end

%% shuffle correction for mb responses
duration = 2.4;
bin_size = 0.00025;
max_lag = 40;
trial_n = 25;
xx = bin_size/2:bin_size:duration;

for ll = 1:length(datamb)
    for cc = 1:length(raster_p_sum_all{ll})
        if ~isempty(raster_p_sum_all{ll}{cc})
            for trial = 1:trial_n
                hist_p{ll}{cc}(trial, :) = hist(raster_p_sum_all{ll}{cc}{trial}, xx);
            end
            hist_p_mean{ll}{cc} = mean(hist_p{ll}{cc});
        end
    end
end

for ll = 1:length(datamb)
    corr_raw{ll} = [];
    for cp = 1:size(cp_i, 1)
        idx1 = idx_dir{1}(cp_i(cp, 1));
        idx2 = idx_dir{1}(cp_i(cp, 2));
        corr_temp = [];
        if ~isempty(hist_p{ll}{idx1}) && ~isempty(hist_p{ll}{idx2})
            for trial = 1:trial_n
                corrected1 = hist_p{ll}{idx1}(trial, :) - hist_p_mean{ll}{idx1};
                corrected2 = hist_p{ll}{idx2}(trial, :) - hist_p_mean{ll}{idx2};
                corr_temp(trial, :) = xcorr(corrected1, corrected2, max_lag, 'coeff');
            end
            corrected_corr{ll}{cp} = corr_temp;
            corrected_corr_mean{ll}{cp} = mean(corr_temp);
        end
    end
end

XX = -max_lag*bin_size:bin_size:max_lag*bin_size;
for i = 1:length(cp_i)
	figure(1)
    for ll = 1:6
        subplot(3,2,ll)
        if ~isempty(corrected_corr_mean{ll}{i})
            bar(XX, corrected_corr_mean{ll}{i})
        end
    end
    pause
end
        
