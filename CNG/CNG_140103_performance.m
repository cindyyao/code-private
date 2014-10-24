opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
datarun{1} = load_data('/Analysis/xyao/2014-01-14-0/data000-001map/data000-001map', opt);
datarun{2} = load_data('/Analysis/xyao/2014-01-16-0/data000-001map/data000-001map', opt);
datarun{3} = load_data('/Analysis/xyao/2014-03-13-0/data000-001-map/data000-001-map', opt);

% find the start time of stimulus 
for i = 1:3
    trigger_diff{i} = diff(datarun{i}.triggers); 
    start_I{i} = [1; find(trigger_diff{i} >6.1)+1]; % after turning on the output, the duration of each trial is 3sec.
    start_time{i} = datarun{i}.triggers(start_I{i}); % convert the index of triggers to time
end

% delete the wrong triggers
start_time{1}(12) = []; 
start_I{1}(12) = [];
start_time{3}([3 4 5 17]) = [];
start_I{3}([3 4 5 17]) = [];


% calculate stimuli strengths
load('LED_calibration.mat')
load('CNG_1401.mat', 'stimuli_setting')


for i = 1:3
    for j = 1:size(stimuli_setting{i}, 1)
        I1 = find(Intensity == stimuli_setting{i}(j, 1));
        I2 = find(Duration == stimuli_setting{i}(j, 2));
        stimuli{i}(j) = LED_calibration(I1, I2)*NDF_factor(stimuli_setting{i}(j, 3));
    end
    [stimuli{i}, seq{i}] = sort(stimuli{i});
    delta = diff(stimuli{i});
    delta(delta ~= 0) = 1;
    delta = [1 delta];
    k = 0;
    for j = 1:length(stimuli{i})
        if delta(j) == 1
            k = k+1;
            stimuli_nrpt{i}{k} = [];
            stimuli_nrpt{i}{k} = [stimuli_nrpt{i}{k} j];
        else
            stimuli_nrpt{i}{k} = [stimuli_nrpt{i}{k} j];
        end
    end
    start_time{i} = start_time{i}(seq{i});
    start_I{i} = start_I{i}(seq{i});
        
end


% get cell ids
flash_id{1} = datarun{1}.cell_ids;
flash_id{2} = datarun{2}.cell_ids;
flash_id{3} = datarun{3}.cell_ids;

% get raster & psth

begin_time{1} = 0:3:117;
begin_time{2} = 0:3:117;
begin_time{3} = 720:3:957;

bin_size = 0.0125;
tau = 4*bin_size;
tt = -3*tau:bin_size:3*tau;
filter = exp(-tt.^2/(2*tau^2));
circ = (length(filter)-1)/2;


for i = 1:3
    for cn = 1:length(flash_id{i})
        idx = get_cell_indices(datarun{i}, flash_id{i}(cn));
        spike = datarun{i}.spikes{idx};
        for st = 1:length(stimuli{i})
            if i == 3 && st < 17
                trigger = datarun{i}.triggers(start_I{i}(st):start_I{i}(st)+79);
            else
                trigger = datarun{i}.triggers(start_I{i}(st):start_I{i}(st)+39);
            end
            raster{i}{cn}{st} = get_raster(spike, trigger, 'plot', false, 'stop', 3);
            raster_all{i}{cn}{st} = sort(cell2mat(raster{i}{cn}{st})); %spikes from all trials are added together
            for t = 1:length(raster{i}{cn}{st})
                clear psth_temp
                psth_temp = hist(raster{i}{cn}{st}{t}, bin_size/2:bin_size:3-bin_size/2);
                psth_temp = [psth_temp(end-circ+1:end) psth_temp psth_temp(1:circ)];
                psth_flash{i}{cn}{st}{t} = conv(psth_temp, filter, 'valid');
            end
        end
%             spike_number_mean_temp = mean(cellfun(@length, raster{i}{cn}));
        for j = 1:length(stimuli_nrpt{i})
%             spikes_number_mean{i}{cn}(j) = mean(spike_number_mean_temp(stimuli_nrpt{i}{j}));
            raster_all_nrpt{i}{cn}{j} = cell2mat(raster_all{i}{cn}(stimuli_nrpt{i}{j})');
        end
        raster_dark{i}{cn} = get_raster(spike, begin_time{i}, 'plot', false, 'stop', 3);
        for t = 1:length(raster_dark{i}{cn})
            clear psth_temp
            psth_temp = hist(raster_dark{i}{cn}{t}, bin_size/2:bin_size:3-bin_size/2);
            psth_temp = [psth_temp(end-circ+1:end) psth_temp psth_temp(1:circ)];
            psth_dark{i}{cn}{t} = conv(psth_temp, filter, 'valid');
        end

    end
end

%% ideal observer

for i = 1:3
    for cn = 1:length(flash_id{i})
        for st = 1:length(stimuli{i})
            for t = 1:length(raster{i}{cn}{st})
                V = sum(cell2mat(psth_flash{i}{cn}{st}')) - psth_flash{i}{cn}{st}{t};
                template1 = V/norm(V);
                V = sum(cell2mat(psth_dark{i}{cn}(1:length(raster{i}{cn}{st}))')) - psth_dark{i}{cn}{t};
                template0 = V/norm(V);
                template = template1 - template0;
                corr_flash{i}{cn}{st}(t) = psth_flash{i}{cn}{st}{t}*template';
                corr_dark{i}{cn}{st}(t) = psth_dark{i}{cn}{t}*template';
            end
            % two alternative forced choice method 1
%             Pc_temp(st) = (sum(corr_flash{i}{cn}{st}>0) + sum(corr_flash{i}{cn}{st}==0)/2 + ...
%                 sum(corr_dark{i}{cn}{st}<0) + sum(corr_dark{i}{cn}{st}==0)/2)/ ...
%                 (length(corr_flash{i}{cn}{st})+length(corr_dark{i}{cn}{st}));
            % two alternative forced choice method 2
%             Pc_temp(st) = (sum(corr_flash{i}{cn}{st}>0) + sum(corr_flash{i}{cn}{st}==0)/2)/ ...
%                 length(corr_flash{i}{cn}{st});
            % two alternative forced choice (0 spike --> darkness)
              Pc_temp(st) = (sum(corr_flash{i}{cn}{st}>0) + ...
                sum(corr_dark{i}{cn}{st}<0) + sum(corr_dark{i}{cn}{st}==0))/ ...
                (length(corr_flash{i}{cn}{st})+length(corr_dark{i}{cn}{st}));
%             % two interval forced choice
%             Pc_temp(st) = sum((corr_flash{i}{cn}{st}(1:40) - corr_dark{i}{cn}{st}) > 0)/40;

        end
        for st = 1:length(stimuli_nrpt{i})
            Pc{i}{cn}(st) = mean(Pc_temp(stimuli_nrpt{i}{st}));
        end
    end
    i
end
    

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
semilogx(unique(stimuli{1}), mean(cell2mat(Pc{1}')), 'b')
hold on
semilogx(unique(stimuli{2}), mean(cell2mat(Pc{2}')), 'r')
semilogx(unique(stimuli{3})*0.0144, mean(cell2mat(Pc{3}')), 'k') % *0.0144 is because we added 2 ndf on the mirror before this experiment
xlabel('flash strength')
ylabel('performance')
legend('KO', 'rescue', 'WT', 'location', 'northwest')

%% plot raster, psth and performance distribution
bin_n = 10;
n = 60;
XX = 3/n/2:3/n:3-3/n/2;

for cd = 3:3;
    a = ceil(length((stimuli{cd})+1)/2);
    for cn =13:13; %length(flash_id{cd});
        st = 1;
        trial_n = length(raster_dark{cd}{cn});
        b = 1;
        FigHandle = figure;
        set(FigHandle, 'Position', [0, 0, 1920, 1080]);
        while st <= length(stimuli{cd})+1
            if  st == a+1
                b = 4;
            end
            subplot(a, 6, b)
            if st > 1
                trial_n = length(raster{cd}{cn}{st-1});
            end
            for j = 1:trial_n
                if st == 1
                    SpikeTime = raster_dark{cd}{cn}{j};
                else
                    SpikeTime = raster{cd}{cn}{st-1}{j};
                end
                SpikeTime = SpikeTime';
                X = [SpikeTime; SpikeTime];
                Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
                line(X, Y, 'color', 'b');
                axis([0, 3, 0, trial_n]);
                hold on
            end
            if b == 1
                title(num2str(flash_id{cd}(cn)))
            end
            b = b+1;
            subplot(a, 6, b)
            if st > 1
                hist(raster_all{cd}{cn}{st-1}, XX)
            else
                hist(cell2mat(raster_dark{cd}{cn}), XX)
            end
            xlim([0 3])
            
            b = b+1;
            if st > 1
                corrd = corr_dark{cd}{cn}{st-1};
                corrf = corr_flash{cd}{cn}{st-1};
                xx = linspace(min([corrd corrf]), max([corrd corrf]), bin_n);
                clear theta
                theta(:, 1) = hist(corrd, xx);
                theta(:, 2) = hist(corrf, xx);
                subplot(a, 6, b)
                bar(xx, theta)
            else
                subplot(a, 6, b)
                semilogx(unique(stimuli{cd}), Pc{cd}{cn})
                ylim([0.5 1])
            end
            

            b = b+4;
            st = st+1;
        end

%         screen_size = [24 12];
%         set(figure(1), 'paperpositionmode', 'auto');
%         % set(figure(1), 'PaperPosition', [-0.5 -0.25 22 10]);
%         set(gcf, 'PaperUnits', 'inch');
%         set(figure(1), 'PaperSize', screen_size);
%         print(figure(1), '-dpdf', [cell_type{ct} '_' condition{cd} '_' num2str(cn)])
%         close
% 


    end
end

%%
% cc = 1;
% cd = 3;
% figure
% for i = 1:35
%     subplot(5, 7, i)
%     hist(corr_dark{3}{i}{19}, 20)
%     title(num2str(flash_id{cd}(i)))
%     hold on
% end