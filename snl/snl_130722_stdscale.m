%% scaled by std/vector 
datarun{1} = datarun2;
datarun{2} = datarun3;
datarun_2 = datarun{2};

display_rate = 60.35;
refresh_rate = 2;
diff_x = -7.2237;
diff_y = 0.3277;

cell_type = {'ON brisk transient', 'ON transient', 'OFF brisk transient', ...
    'OFF transient', 'OFF slow', 'OFF sustained'};
n = length(cell_type);
[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun{1}, ...
    datarun{2}, cell_type);
%%
ct = 5;
cell_num = size(cell_id{ct}, 1);

    for cn = 10:cell_num          
        FigHandle = figure;
        set(FigHandle, 'Position', [100, 100, 1250, 800]);
        
        for i = 1:2
            idx = cell_idx{ct}(cn, i);
            com = datarun{i}.stas.rf_coms{idx};
            gen_signals(:, i) = datarun{i}.stas.snls{idx}.gen_signal;
            spikes(:, i) = datarun{i}.stas.snls{idx}.spikes;
            tc(:, i) = datarun{i}.stas.time_courses{idx};
            r = datarun{i}.stimulus.field_height;
            
            if i == 1
                rc = abs(com(1)-r/2) - abs(com(2)-r/2);
            end
            
            rf = datarun{i}.stas.rfs{idx};
            if rc <= 0
                sf{i} = rf(ceil(com(2)), :);
                diff = diff_x;
                x2 = [1:60]*2/3+diff;
                
            else
                sf{i} = rf(:, ceil(com(1)));
                diff = diff_y;
                x2 = [1:60]*2/3+diff;

            end
        
            [X(:, i), Y(:, i)] = curve_from_binning(gen_signals(:, i), spikes(:, i), 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
          
        end
        Y = Y*display_rate/refresh_rate;
        
        % unscaled

        subplot(3, 4, 1);
        x1 = 1:40;
        plot(x1, sf{1}, 'b');
        hold on
        plot(x2, sf{2}, 'r');
        title([num2str(cell_id{ct}(cn, 1)) ' ' num2str(cell_id{ct}(cn, 2))])
        xlim([0 40])
        
        subplot(3, 4, 5);
        t = -0.033*[15:-1:1];
        plot(t, tc(:, 1), 'b');
        hold on
        plot(t, tc(:, 2), 'r');
        xlabel('time to spike (s)')
        ylabel('STA contrast')
        
        subplot(3, 4, 9);
        plot(X(:, 1), Y(:, 1), 'b');
        hold on
        plot(X(:, 2), Y(:, 2), 'r');
        xlabel('generator signal')
        ylabel('spike rate(HZ)')
        h_legend = legend('NDF 4', 'NDF 0', 'location', 'northwest');
        set(h_legend,'FontSize',10);

        % calculate the scaling factor
        sigma_1 = std(gen_signals(:, 1));
        sigma_2 = std(gen_signals(:, 2));
        f_gs = sigma_1/sigma_2;
        
        % apply this factor to strf and recalculate the sf and tf
        datarun_2.stas.stas{idx} = datarun{2}.stas.stas{idx}*f_gs;
        datarun_2 = get_sta_summaries(datarun_2, datarun_2.cell_ids(idx));
        
        
        gen_signals(:, 2) = gen_signals(:, 2)*f_gs;
        tc(:, 2) = datarun_2.stas.time_courses{idx};
        rf = datarun_2.stas.rfs{idx};
        if rc <= 0
            sf{2} = rf(ceil(com(2)), :);
        else
            sf{2} = rf(:, ceil(com(1)));
        end
        [X(:, 2), Y(:, 2)] = curve_from_binning(gen_signals(:, 2), spikes(:, 2), 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
        Y(:, 2) = Y(:, 2)*display_rate/refresh_rate;
        subplot(3, 4, 2);
        x1 = 1:40;
        plot(x1, sf{1}, 'b');
        hold on
        plot(x2, sf{2}, 'r');
        xlim([0 40])

        subplot(3, 4, 6);
        t = -0.033*[15:-1:1];
        plot(t, tc(:, 1), 'b');
        hold on
        plot(t, tc(:, 2), 'r');
        xlabel('time to spike (s)')
        ylabel('STA contrast')
        
        subplot(3, 4, 10);
        plot(X(:, 1), Y(:, 1), 'b');
        hold on
        plot(X(:, 2), Y(:, 2), 'r');
        xlabel('generator signal')
        ylabel('spike rate(HZ)')
        
        
        % scale spatial filter by std
        
        sf_std{1} = sf{1}/std(sf{1});
        sf_std{2} = sf{2}/std(sf{2});
        
        subplot(3, 4, 3)
        x1 = 1:40;
        plot(x1, sf_std{1}, 'b');
        hold on
        plot(x2, sf_std{2}, 'r');
        xlim([0 40])
        
        
        % scale temporal filter by std
        tc_std(:, 1) = tc(:, 1)/std(tc(:, 1));
        tc_std(:, 2) = tc(:, 2)/std(tc(:, 2));
        
        subplot(3, 4, 7)
        t = -0.033*[15:-1:1];
        plot(t, tc_std(:, 1), 'b');
        hold on
        plot(t, tc_std(:, 2), 'r');
        xlabel('time to spike (s)')
        ylabel('STA contrast')
        
        % nl
        subplot(3, 4, 11);
        plot(X(:, 1), Y(:, 1), 'b');
        hold on
        plot(X(:, 2), Y(:, 2), 'r');
        xlabel('generator signal')
        ylabel('spike rate(HZ)')
        
        
        % scale spatial filter by vector length
        % use linear interpolation to make the dimension of vector consistant
        x_inter = 1:1/3:40;
        sf_v{1} = interp1(x1, sf{1}, x_inter);
        sf_v{1} = sf_v{1}/norm(sf_v{1});
        sf_v{2} = interp1(x2, sf{2}, x_inter+diff);
        sf_v{2} = sf_v{2}/norm(sf_v{2});
       
        subplot(3, 4, 4)
        plot(x_inter, sf_v{1}, 'b');
        hold on
        plot(x_inter+diff, sf_v{2}, 'r');
        xlim([0 40])

        
        tc_v(:, 1) = tc(:, 1)/norm(tc(:, 1));
        tc_v(:, 2) = tc(:, 2)/norm(tc(:, 2));
        subplot(3, 4, 8)
        plot(t, tc_v(:, 1), 'b');
        hold on
        plot(t, tc_v(:, 2), 'r');
        xlabel('time to spike (s)')
        ylabel('STA contrast')
        
        subplot(3, 4, 12)
        plot(X(:, 1), Y(:, 1), 'b');
        hold on
        plot(X(:, 2), Y(:, 2), 'r');
        xlabel('generator signal')
        ylabel('spike rate(HZ)')

        
        
    end


    
