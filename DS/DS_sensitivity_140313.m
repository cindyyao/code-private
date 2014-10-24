opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
datarun{1} = load_data('/Analysis/xyao/2014-03-13-0/data000-001-map/data000-001-map', opt);
datarun{2} = load_data('/Analysis/xyao/2014-03-13-0/data004-map/data004-map', opt);
datarun{2}.names.stimulus_path = '/Analysis/xyao/2014-03-13-0/stimuli/s04';
datarun{2} = load_stim(datarun{2}, 'user_defined_trigger_interval', 10);
datarun{3} = load_data('/Analysis/xyao/2014-03-13-0/data003/data003', opt);


[NumSpikesCell, StimComb] = get_spikescellstim(datarun{2},datarun{2}.cell_ids,0);
[mag MAG dsindex magmax magave angle rho RHO theta num U V] = dscellanalysis(NumSpikesCell, StimComb);
DS = v2struct(mag, MAG, dsindex, magmax, magave, angle, rho, RHO, theta, num, U, V);

% pull out DS cells

figure
plot(mag{1, 1}, mag{2, 1}, 'o')
title('vector sum plot')
xlabel('TP 32')
ylabel('TP 128')

hold on
[x, y] = ginput;
plot(x, y);

IN = inpolygon(mag{1, 1}, mag{2, 1}, x, y);
[~, I] = find(IN == 1);
id = datarun{2}.cell_ids(I);

id_ds_flash = intersect(id, datarun{1}.cell_ids);

[cell_list_map, failed_cells] = map_ei(datarun{2}, datarun{1}, 'master_cell_type', 'DS');

id_ds_flash = [691 946 2101 3588 5162 574 693];
id_ds_dg = [691 946 2101 3586 5162 436 632];

nds_type = {'non DS on 1', 'non DS on 2','non DS off 1','non DS off 2','non DS off 3','non DS off 4'};
    
[nds_id, nds_idx, nds_id_matched, nds_idx_matched] = cell_map2(datarun{1}, datarun{3}, nds_type);
%%

% sort the sequence of directions
for j = 1:2
    [theta_seq, I_seq] = sort(DS.theta{j}(1, :));
    r = DS.rho{j};
    R = DS.RHO{j};
    rho_seq = r(:, I_seq);
    RHO_seq = R(:, I_seq);
    DS.rho{j} = rho_seq;
    DS.RHO{j} = RHO_seq;
end

%% plot sensitivity curve
n = 120;
XX = 3/n/2:3/n:3-3/n/2;
%     FigHandle = figure;
%     set(FigHandle, 'Position', [0, 0, 1080, 1080]);
% 
for j = 3:3;
    cn = length(id_ds_flash);
    for k = 1:cn
    %         subplot(ceil(sqrt(cn)), ceil(sqrt(cn)), k)
    %         semilogx(unique(stimuli{j}), spikes_number_mean{j}{i}{k})
        clear sens
        for st = 1:length(stimuli_nrpt{j})
            N = hist(raster_all_nrpt{j}{k}{st}, XX)/length(stimuli_nrpt{j}{st})/length(raster{j}{k}{st})*40;
            sens(st) = sqrt(mean((N(1:2/3*length(XX)) - mean(N(2/3*length(XX)+1:end))).^2)) - std(N(2/3*length(XX)+1:end)); 
        end
        sensitivity_all{j}(k, :) = sens;
    %         semilogx(unique(stimuli{j}), sens)
    %         ylim([0 max(sens)])
    %         if k == 1
    %             title([cell_type{i} ' ' condition{j}])
    %         end
    end
end

figure
for cd = 3:3
for j = 1:length(id_ds_flash)
semilogx(unique(stimuli{cd}), sensitivity_all{cd}(j, :), 'b');
%
hold on
pause
end

end
legend(h, condition{1}, condition{2}, 'location', 'northwest')
    
%% performance

load('CNG_140103_pc.mat')

idx_flash = arrayfun(@(x) find(datarun{1}.cell_ids == x), id_ds_flash);
idx_dg = arrayfun(@(x) find(datarun{2}.cell_ids == x), id_ds_dg);

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:length(id_ds_flash)
    semilogx(unique(stimuli{3}), Pc{3}{idx_flash(i)}, 'b')
    hold on
    pause
end
xlabel('flash strength')
ylabel('performance')
% title('flash only')
% legend('KO', 'rescue', 'WT', 'location', 'northwest')

id1_dg = id_ds_dg([1 5]);
id2_dg = id_ds_dg([2 3 4 6 7]);
idx1_flash = idx_flash([1 5]);
idx2_flash = idx_flash([2 3 4 6 7]);


raster_ds = get_ds_raster(datarun{2}, id_ds_dg);

%% compare ds with non ds

figure
set(gcf, 'DefaultLineLineWidth', 1.5)
for i = 1:length(idx1_flash)
    h(1) = semilogx(unique(stimuli{3}), Pc{3}{idx1_flash(i)}, 'b');
    hold on
%     pause
end
for i = 1:length(idx2_flash)
    h(2) = semilogx(unique(stimuli{3}), Pc{3}{idx2_flash(i)}, 'c');
    hold on
%     pause
end

color = 'rgkym';
for i = 1:5
    for j = 1:length(nds_id{i}(:, 1))
    h(i+2) = semilogx(unique(stimuli{3}), Pc{3}{nds_idx{i}(j, 1)}, 'color', color(i));
    hold on
%     pause
    end
end


xlabel('flash strength')
ylabel('performance')
legend(h, 'DS 1', 'DS 2', 'ON 1', 'ON 2', 'OFF 1', 'OFF 2', 'OFF 3', 'location', 'northwest')

%% average pc
figure
set(gcf, 'DefaultLineLineWidth', 1.5)
semilogx(unique(stimuli{3}), mean(cell2mat(Pc{3}(idx1_flash)')), 'b');
hold on
semilogx(unique(stimuli{3}), mean(cell2mat(Pc{3}(idx2_flash)')), 'c');

color = 'rgkym';
for i = 1:5
    semilogx(unique(stimuli{3}), mean(cell2mat(Pc{3}(nds_idx{i}(:, 1))'), 1), 'color', color(i));
end


xlabel('flash strength')
ylabel('performance')
legend('DS 1', 'DS 2', 'ON 1', 'ON 2', 'OFF 1', 'OFF 2', 'OFF 3', 'location', 'northwest')


%% plot ds raster
tt = datarun{2}.stimulus.params.DIRECTION*pi/180;
v = datarun{2}.stimulus.params.SPATIAL_PERIOD./datarun{2}.stimulus.params.TEMPORAL_PERIOD;
x = 3;
y = 6;
a = [8 9 3 2 1 7 13 14 15; 11 12 6 5 4 10 16 17 18];


for cc = 1:length(id_ds_dg)
    FigHandle = figure;
    set(FigHandle, 'Position', get(0, 'ScreenSize'))
    for t = 1:length(v)
        subplot(x, y, a(t, 1)) ; polar(tt, DS.rho{t}(idx_dg(cc), :))
        for i = 2:9
            subplot(x, y, a(t, i)); plot_raster(squeeze(raster_ds{cc}(1, t, i-1, :)), 0, 8)
            if i == 4
                title(id_ds_dg(cc))
            end
        end
%         name = [num2str(id_nDS_all(idx_nds_4(cc))) '_' num2str(t)];
%         screen_size = [24 12];
%         set(figure(1), 'paperpositionmode', 'auto');
%         % set(figure(1), 'PaperPosition', [-0.5 -0.25 22 10]);
%         set(gcf, 'PaperUnits', 'inch');
%         set(figure(1), 'PaperSize', screen_size);
%         print(figure(1), '-dpdf', name)
%         close

    end
    
end
%% flash summary
bin_n = 10;
n = 60;
XX = 3/n/2:3/n:3-3/n/2;

for cd = 3:3;
    a = ceil(length((stimuli{cd})+1)/2);
    for cc =2:length(nds_idx_matched(:, 1));
        cn = nds_idx_matched(cc, 1);
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
                title(['Flash: ' num2str(nds_id_matched(cc, 1)) '  WN: ' num2str(nds_id_matched(cc, 2))])
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

        screen_size = [24 12];
        set(figure(1), 'paperpositionmode', 'auto');
        % set(figure(1), 'PaperPosition', [-0.5 -0.25 22 10]);
        set(gcf, 'PaperUnits', 'inch');
        set(figure(1), 'PaperSize', screen_size);
        print(figure(1), '-dpdf', ['Flash_' num2str(nds_id_matched(cc, 1)) '_WN_' num2str(nds_id_matched(cc, 2))]);
        close



    end
end


