% cd /Users/xyao/matlab/code-private/DS_new/
% opt = struct('load_params', 1,'load_neurons', 1);%, 'load_ei', 1);
% 
% 
% datarun = load_data('/Volumes/lab/analysis/2016-05-06-0/data000-001-map/data000-001-map', opt);
% datamb(1:2) = split_datarun(datarun, 1890);
% datamb{1}.names.stimulus_path = '/Volumes/lab/analysis/2016-05-06-0/stimuli/s00.mat';
% datamb{1} = load_stim_matlab(datamb{1});
% datamb{2}.names.stimulus_path = '/Volumes/lab/analysis/2016-05-06-0/stimuli/s01.mat';
% datamb{2} = load_stim_matlab(datamb{2});
% 
% datarun = load_data('/Volumes/lab/analysis/2016-05-06-0/data002-003-map/data002-003-map', opt);
% datamb(3:4) = split_datarun(datarun, 1880);
% datamb{3}.names.stimulus_path = '/Volumes/lab/analysis/2016-05-06-0/stimuli/s02.mat';
% datamb{3} = load_stim_matlab(datamb{3});
% datamb{4}.names.stimulus_path = '/Volumes/lab/analysis/2016-05-06-0/stimuli/s03.mat';
% datamb{4} = load_stim_matlab(datamb{4});
% 
% datarun = load_data('/Volumes/lab/analysis/2016-05-06-0/data005-006-map/data005-006-map', opt);
% datamb(5:6) = split_datarun(datarun, 1880);
% datamb{5}.names.stimulus_path = '/Volumes/lab/analysis/2016-05-06-0/stimuli/s05.mat';
% datamb{5} = load_stim_matlab(datamb{5});
% datamb{6}.names.stimulus_path = '/Volumes/lab/analysis/2016-05-06-0/stimuli/s06.mat';
% datamb{6} = load_stim_matlab(datamb{6});
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
%     raster_mb{i} = get_mb_raster(datamb{i}, ds_id, datamb{i}.stimulus.triggers(865));
    for j = 1:length(raster_mb{i})
        if(idx_mb(j, ceil(i/2)))
            raster_mb{i}{j} = [];
        end
    end
end

for i = 1:n;
    for cc = 1:length(raster_mb{i})
        raster_mb_12{i}{cc} = raster_mb{i}{cc}(:, :, 1:3:end, :,:);
    end
    trial_dur{i} = get_mb_trial_dur(datamb{i}, 400, 400, 0.5);
end
close all
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
    end
    SpikeN{ll} = SpikeN_temp;
    SpikeN_mean{ll} = SpikeN_mean_temp;
    SpikeN_std{ll} = SpikeN_std_temp;
end

LL = {'NDF 4', 'NDF 3', 'NDF 0'};
figure
for ll = 1:1
    figure
    for ct = 1:4
        for cc = 6:6%length(id_dir{ct})
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
legend(h, 'NDF 4', 'NDF 3', 'NDF 0')
xlabel('mean (spike #)')
ylabel('variance (spike #)')
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
