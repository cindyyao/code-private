display_rate = 60.35;
refresh_rate = 2;


ct = 5;
cell_num = size(cell_id{ct}, 1);

    for cn = 11:cell_num
          
        FigHandle = figure;
        set(FigHandle, 'Position', [100, 100, 900, 530]);
        
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
                sf(:, i) = rf(ceil(com(2)), :);
            else
                sf(:, i) = rf(:, ceil(com(1)));
            end
        
            [X(:, i), Y(:, i)] = curve_from_binning(gen_signals(:, i), spikes(:, i), 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
          
        end
        Y = Y*display_rate/refresh_rate;
        
        % unscaled

        subplot(2, 3, 1);
        x1 = 1:40;
        plot(x1, sf(:, 1), 'b');
        hold on
        plot(x1, sf(:, 2), 'r');
        title([num2str(cell_id{ct}(cn, 1)) ' ' num2str(cell_id{ct}(cn, 2))])
        
        subplot(2, 3, 2);
        t = -0.033*[15:-1:1];
        plot(t, tc(:, 1), 'b');
        hold on
        plot(t, tc(:, 2), 'r');
        xlabel('time to spike (s)')
        ylabel('STA contrast')
        
        subplot(2, 3, 3);
        plot(X(:, 1), Y(:, 1), 'b');
        hold on
        plot(X(:, 2), Y(:, 2), 'r');
        xlabel('generator signal')
        ylabel('spike rate(HZ)')
        h_legend = legend('melatonin', 'non-melatonin', 'location', 'northwest');
        set(h_legend,'FontSize',10);

        % scale the linear filter
        
        m_sf = max(sf);
        f_sf = m_sf(1)/m_sf(2);
        sf_scale_l(:, 1) = sf(:, 1);
        sf_scale_l(:, 2) = sf(:, 2)*f_sf;
        
        m_tc = max(tc);
        f_tc = m_tc(1)/m_tc(2);
        tc_scale_l(:, 1) = tc(:, 1);
        tc_scale_l(:, 2) = tc(:, 2)*f_tc;
        
        
        f_l = f_sf*f_tc;
        
        gen_signals_scale_l(:, 1) = gen_signals(:, 1);
        gen_signals_scale_l(:, 2) = gen_signals(:, 2)*f_l;
        
        [X_scale_l(:, 1), Y_scale_l(:, 1)] = curve_from_binning(gen_signals_scale_l(:, 1), spikes(:, 1), 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
        [X_scale_l(:, 2), Y_scale_l(:, 2)] = curve_from_binning(gen_signals_scale_l(:, 2), spikes(:, 2), 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
        Y_scale_l = Y_scale_l*display_rate/refresh_rate;
        
        subplot(2, 3, 4);
        x1 = 1:40;
        plot(x1, sf_scale_l(:, 1), 'b');
        hold on
        plot(x1, sf_scale_l(:, 2), 'r');
        
        subplot(2, 3, 5);
        t = -0.033*[15:-1:1];
        plot(t, tc_scale_l(:, 1), 'b');
        hold on
        plot(t, tc_scale_l(:, 2), 'r');
        xlabel('time to spike (s)')
        ylabel('STA contrast')
        
        subplot(2, 3, 6);
        plot(X_scale_l(:, 1), Y_scale_l(:, 1), 'b');
        hold on
        plot(X_scale_l(:, 2), Y_scale_l(:, 2), 'r');
        xlabel('generator signal')
        ylabel('spike rate(HZ)')
        
%         % scale nonlinearity
%         
%         m_Y = max(Y);
%         [a, I] = min(m_Y);
%         if I == 1
%             [~,II] = min(abs(Y(:, 2) - a));
%             f_nl = X(end, 1)/X(II, 2);
%         else
%             [~,II] = min(abs(Y(:, 1) - a));
%             f_nl = X(II, 1)/X(end, 2);
%         end
%         
%         gen_signals_scale_nl(:, 1) = gen_signals(:, 1);
%         gen_signals_scale_nl(:, 2) = gen_signals(:, 2)*f_nl;
% 
%         
%         
%         
%         
%         
%         f_Y = m_Y(1)/m_Y(2);
%         Y_scale_nl(:, 1) = Y(:, 1);
%         Y_scale_nl(:, 2) = Y(:, 2)*f_Y;
%         
%         subplot(3, 4, 3)
%         x = 1:40;
%         plot(x, sf(:, 1), 'b');
%         hold on
%         plot(x, sf(:, 2)*f_nl/f_tc, 'r');
%         
%         subplot(3, 4, 7);
%         t = -0.033*[15:-1:1];
%         plot(t, tc_scale_l(:, 1), 'b');
%         hold on
%         plot(t, tc_scale_l(:, 2), 'r');
%         xlabel('time to spike (s)')
%         ylabel('STA contrast')
%         
%         subplot(3, 4, 11);
%         plot(X(:, 1), Y(:, 1), 'b');
%         hold on
%         plot(X(:, 2)*f_nl, Y(:, 2), 'r');
%         xlabel('generator signal')
%         ylabel('spike rate(HZ)')
%         
%         subplot(3, 4, 4);
%         x = 1:40;
%         plot(x, sf_scale_l(:, 1), 'b');
%         hold on
%         plot(x, sf_scale_l(:, 2), 'r');
%         
%         subplot(3, 4, 8);
%         t = -0.033*[15:-1:1];
%         plot(t, tc(:, 1), 'b');
%         hold on
%         plot(t, tc(:, 2)*f_nl/f_sf, 'r');
%         xlabel('time to spike (s)')
%         ylabel('STA contrast')
%         
%         subplot(3, 4, 12);
%         plot(X(:, 1), Y(:, 1), 'b');
%         hold on
%         plot(X(:, 2)*f_nl, Y(:, 2), 'r');
%         xlabel('generator signal')
%         ylabel('spike rate(HZ)')
        
    end
    
