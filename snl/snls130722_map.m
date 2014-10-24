opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);

% ndf 4 melatonin

datarun{1} = load_data('/Analysis/xyao/2013-07-22-0/data002-map/data002-map', opt);
datarun{1} = load_java_movie(datarun{1}, '/Volumes/lab/acquisition/movie-xml/BW-15-2-0.48-11111-40x40-60.35.xml');
datarun{1} = get_sta_summaries(datarun{1}, 'all');

datarun{1} = get_snls(datarun{1}, 'all');


% ndf 4 non-melatonin

datarun{2} = load_data('/Analysis/xyao/2013-07-22-0/data005-map/data005-map', opt);
datarun{2} = load_java_movie(datarun{2}, '/Volumes/lab/acquisition/movie-xml/BW-15-2-0.48-11111-40x40-60.35.xml');
datarun{2} = get_sta_summaries(datarun{2}, 'all');

datarun{2} = get_snls(datarun{2}, 'all');



% ndf 0 non-melatonin

datarun{3} = load_data('/Analysis/xyao/2013-07-22-0/data010-map/data010-map', opt);
datarun{3} = load_java_movie(datarun{3}, '/Volumes/lab/acquisition/movie-xml/BW-10-2-0.48-11111-60x60-60.35.xml');
datarun{3} = get_sta_summaries(datarun{3}, 'all');

datarun{3} = get_snls(datarun{3}, 'all');


datarun1 = datarun{1};
datarun2 = datarun{2};
datarun3 = datarun{3};

%% plot snls

cell_type = {'ON brisk transient', 'ON transient', 'OFF brisk transient', ...
    'OFF transient', 'OFF slow'};
n = length(cell_type);
[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun1, ...
    datarun2, cell_type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ndf 4 melatonin
ct = 6;
id = cell_id{ct};
idx = cell_idx{ct};
cell_numb = size(id, 1);

if cell_numb <= ceil(sqrt(cell_numb))*floor(sqrt(cell_numb))
    dimx = floor(sqrt(cell_numb));
    dimy = ceil(sqrt(cell_numb));
else
    dimx = ceil(sqrt(cell_numb));
    dimy = dimx;
end


figure
for i = 1:cell_numb
    gen_signals = datarun1.stas.snls{idx(i, 1)}.gen_signal;
    spikes = datarun1.stas.snls{idx(i, 1)}.spikes;
    fit = datarun1.stas.snls{idx(i, 1)}.fit_params;
    
            
    subplot(dimx, dimy, i)
    plot_snl(gen_signals, spikes, 'fit', fit, 'foa', -1, 'fig_title', num2str(id(i, 1)));  
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ndf 4 non-melatonin


ct = 6;
id = cell_id{ct};
idx = cell_idx{ct};
cell_numb = size(id, 1);

if cell_numb <= ceil(sqrt(cell_numb))*floor(sqrt(cell_numb))
    dimx = floor(sqrt(cell_numb));
    dimy = ceil(sqrt(cell_numb));
else
    dimx = ceil(sqrt(cell_numb));
    dimy = dimx;
end


figure
for i = 1:cell_numb
    gen_signals = datarun2.stas.snls{idx(i, 2)}.gen_signal;
    spikes = datarun2.stas.snls{idx(i, 2)}.spikes;
    fit = datarun2.stas.snls{idx(i, 2)}.fit_params;
    
            
    subplot(dimx, dimy, i)
    plot_snl(gen_signals, spikes, 'fit', fit, 'foa', -1, 'fig_title', num2str(id(i, 2)));  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ndf 0 non-melatonin


n = length(cell_type);
[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun1, ...
    datarun2, cell_type);

ct = 6;
id = cell_id{ct};
idx = cell_idx{ct};
cell_numb = size(id, 1);

if cell_numb <= ceil(sqrt(cell_numb))*floor(sqrt(cell_numb))
    dimx = floor(sqrt(cell_numb));
    dimy = ceil(sqrt(cell_numb));
else
    dimx = ceil(sqrt(cell_numb));
    dimy = dimx;
end


figure
for i = 1:cell_numb
    gen_signals = datarun3.stas.snls{idx(i, 1)}.gen_signal;
    spikes = datarun3.stas.snls{idx(i, 1)}.spikes;
    fit = datarun3.stas.snls{idx(i, 1)}.fit_params;
    
            
    subplot(dimx, dimy, i)
    plot_snl_(gen_signals, spikes, 'fit', fit, 'foa', -1, 'fig_title', num2str(id(i, 1)));  
end



%%
% comparison between melatonin & non-melatonin


ct = 5;
id = cell_id{ct};
idx = cell_idx{ct};
cell_numb = size(id, 1);

if cell_numb <= ceil(sqrt(cell_numb))*floor(sqrt(cell_numb))
    dimx = floor(sqrt(cell_numb));
    dimy = ceil(sqrt(cell_numb));
else
    dimx = ceil(sqrt(cell_numb));
    dimy = dimx;
end


figure
for i = 1:cell_numb
    gen_signals(:, 1) = datarun1.stas.snls{idx(i, 1)}.gen_signal;
    spikes(:, 1) = datarun1.stas.snls{idx(i, 1)}.spikes;

    
    gen_signals(:, 2) = datarun2.stas.snls{idx(i, 2)}.gen_signal;
    spikes(:, 2) = datarun2.stas.snls{idx(i, 2)}.spikes;
    
    [X1, Y1] = curve_from_binning(gen_signals(:, 1),spikes(:, 1),'average_y','mean','average_x','mean','num_bins',20);
    [X2, Y2] = curve_from_binning(gen_signals(:, 2),spikes(:, 2),'average_y','mean','average_x','mean','num_bins',20);
    
%     f = max(Y1)/max(Y2);
%     Y2 = Y2*f;

    subplot(dimx, dimy, i)
    plot(X1,Y1,'b-+')
    hold on
    plot(X2,Y2,'r-+')
    
    if i == 1
        legend('melatonin', 'non-melatonin')
        
    end
    title(sprintf('%d & %d', id(i, 1), id(i, 2)))
            
end


%%
% comparison between ndf 4 & ndf 0

% n = length(cell_type);
% [cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun3, ...
%     datarun2, cell_type);

ct = 5;
id = cell_id{ct};
idx = cell_idx{ct};
cell_numb = size(id, 1);

if cell_numb <= ceil(sqrt(cell_numb))*floor(sqrt(cell_numb))
    dimx = floor(sqrt(cell_numb));
    dimy = ceil(sqrt(cell_numb));
else
    dimx = ceil(sqrt(cell_numb));
    dimy = dimx;
end


figure
for i = 1:cell_numb
    gen_signals(:, 1) = datarun3.stas.snls{idx(i, 1)}.gen_signal;
    spikes(:, 1) = datarun3.stas.snls{idx(i, 1)}.spikes;

    
    gen_signals(:, 2) = datarun2.stas.snls{idx(i, 2)}.gen_signal;
    spikes(:, 2) = datarun2.stas.snls{idx(i, 2)}.spikes;
    
    [X1, Y1] = curve_from_binning(gen_signals(:, 1),spikes(:, 1),'average_y','mean','average_x','mean','num_bins',20);
    [X2, Y2] = curve_from_binning(gen_signals(:, 2),spikes(:, 2),'average_y','mean','average_x','mean','num_bins',20);
    
%     f = max(Y1)/max(Y2);
%     Y2 = Y2*f;

    subplot(dimx, dimy, i)
    plot(X1,Y1,'b-+')
    hold on
    plot(X2,Y2,'r-+')
    
    if i == 1
        legend('ndf 0', 'ndf 4', 'location', 'northwest')
        
    end
    title(sprintf('%d & %d', id(i, 1), id(i, 2)))
            
end

%% summary

datarun{1} = datarun1;
datarun{2} = datarun2;
datarun{3} = datarun3;

figure;

ct = 5;
    for cn = 1:5
        a = [7 8 11 12];
        for i = 1:2
            idx = cell_idx{ct}(a(cn), i);
            com = datarun{i}.stas.rf_coms{idx};
            gen_signals = datarun{i}.stas.snls{idx}.gen_signal;
            spikes = datarun{i}.stas.snls{idx}.spikes;
            tc(:, i) = datarun{i}.stas.time_courses{idx};
            r = datarun{i}.stimulus.field_height;
            rc = abs(com(1)-r/2) - abs(com(2)-r/2);
            rf = datarun{i}.stas.rfs{idx};
            if rc <= 0
                sf(:, i) = rf(ceil(com(2)), :);
            else
                sf(:, i) = rf(:, ceil(com(1)));
            end
        
            [X(:, i), Y(:, i)] = curve_from_binning(gen_signals, spikes, 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
            
            
            
            
            if i == 1
                subplot(3, 5, cn);
                r = 1:40;
                plot(r, sf(:, i));
                hold on
                title(num2str(cell_id{ct}(a(cn), 1)))
                
        
                subplot(3, 5, cn+5); 
                t = -0.033*[15:-1:1];
                plot(t, tc(:, i));
                hold on
                xlabel('time to spike (s)')
                ylabel('STA contrast')
        
                subplot(3, 5, cn+10); plot(X(:, i), Y(:, i));
                hold on
                xlabel('generator signal')
                ylabel('spike rate(HZ)')
                
            else
                subplot(3, 5, cn);
                r = 1:40;
                plot(r, sf(:, i), 'r');
                ylim([min(sf(:)) max(sf(:))]);
                if cn == 1
                    
                   title([cell_type{ct} '   '  num2str(cell_id{ct}(a(cn), 1))]);
                end
        
                subplot(3, 5, cn+5); 
                t = -0.033*[15:-1:1];
                plot(t, tc(:, i), 'r');
                xlim([-0.5 0])
                ylim([min(tc(:)) max(tc(:))]);

        
                subplot(3, 5, cn+10); plot(X(:, i), Y(:, i), 'r');
                ylim([min(Y(:)) max(Y(:))]);
                if cn == 1
                    h_legend = legend('melatonin', 'non-melatonin', 'location', 'northwest');
                    set(h_legend,'FontSize',12);
                end
                

            end
                
        
        end
    end
    
     
    
%%


% ct = 1;
% cell_num = size(cell_id{ct}, 1);
% 
%     for cn = 1:1
%         figure
%                
%         for i = 1:2
%             idx = cell_idx{ct}(cn, i);
%             com = datarun{i}.stas.rf_coms{idx};
%             gen_signals(:, i) = datarun{i}.stas.snls{idx}.gen_signal;
%             spikes(:, i) = datarun{i}.stas.snls{idx}.spikes;
%             tc(:, i) = datarun{i}.stas.time_courses{idx};
%             r = datarun{i}.stimulus.field_height;
%             
%             if i == 1
%                 rc = abs(com(1)-r/2) - abs(com(2)-r/2);
%             end
%             
%             rf = datarun{i}.stas.rfs{idx};
%             if rc <= 0
%                 sf(:, i) = rf(ceil(com(2)), :);
%             else
%                 sf(:, i) = rf(:, ceil(com(1)));
%             end
%         
%             [X(:, i), Y(:, i)] = curve_from_binning(gen_signals(:, i), spikes(:, i), 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
%         end
%         
%         % unscaled
% 
%         subplot(3, 4, 1);
%         x = 1:40;
%         plot(x, sf(:, 1), 'b');
%         hold on
%         plot(x, sf(:, 2), 'r');
%         
%         subplot(3, 4, 5);
%         t = -0.033*[15:-1:1];
%         plot(t, tc(:, 1), 'b');
%         hold on
%         plot(t, tc(:, 2), 'r');
%         xlabel('time to spike (s)')
%         ylabel('STA contrast')
%         
%         subplot(3, 4, 9);
%         plot(X(:, 1), Y(:, 1), 'b');
%         hold on
%         plot(X(:, 2), Y(:, 2), 'r');
%         xlabel('generator signal')
%         ylabel('spike rate(HZ)')
%         
%         % scale the linear filter
%         
%         m_sf = max(sf);
%         f_sf = m_sf(1)/m_sf(2);
%         sf_scale_l(:, 1) = sf(:, 1);
%         sf_scale_l(:, 2) = sf(:, 2)*f_sf;
%         
%         m_tc = max(tc);
%         f_tc = m_tc(1)/m_tc(2);
%         tc_scale_l(:, 1) = tc(:, 1);
%         tc_scale_l(:, 2) = tc(:, 2)*f_tc;
%         
%         
%         f_l = f_sf*f_tc;
%         
%         gen_signals_scale_l(:, 1) = gen_signals(:, 1);
%         gen_signals_scale_l(:, 2) = gen_signals(:, 2)*f_l;
%         
%         [X_scale_l(:, 1), Y_scale_l(:, 1)] = curve_from_binning(gen_signals_scale_l(:, 1), spikes(:, 1), 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
%         [X_scale_l(:, 2), Y_scale_l(:, 2)] = curve_from_binning(gen_signals_scale_l(:, 2), spikes(:, 2), 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
% 
%         subplot(3, 4, 2);
%         x = 1:40;
%         plot(x, sf_scale_l(:, 1), 'b');
%         hold on
%         plot(x, sf_scale_l(:, 2), 'r');
%         
%         subplot(3, 4, 6);
%         t = -0.033*[15:-1:1];
%         plot(t, tc_scale_l(:, 1), 'b');
%         hold on
%         plot(t, tc_scale_l(:, 2), 'r');
%         xlabel('time to spike (s)')
%         ylabel('STA contrast')
%         
%         subplot(3, 4, 10);
%         plot(X_scale_l(:, 1), Y_scale_l(:, 1), 'b');
%         hold on
%         plot(X_scale_l(:, 2), Y_scale_l(:, 2), 'r');
%         xlabel('generator signal')
%         ylabel('spike rate(HZ)')
%         
%         % scale nonlinearity
%         
%         m_Y = max(Y);
%         
%         f_Y = m_Y(1)/m_Y(2);
%         Y_scale_nl(:, 1) = Y(:, 1);
%         Y_scale_nl(:, 2) = Y(:, 2)*f_Y;
%         
%         subplot(3, 4, 3)
%         x = 1:40;
%         plot(x, sf(:, 1), 'b');
%         hold on
%         plot(x, sf(:, 2)/(f_Y*f_tc), 'r');
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
%         plot(X_scale_l(:, 1), Y_scale_nl(:, 1), 'b');
%         hold on
%         plot(X_scale_l(:, 2), Y_scale_nl(:, 2), 'r');
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
%         plot(t, tc(:, 2)/(f_Y*f_sf), 'r');
%         xlabel('time to spike (s)')
%         ylabel('STA contrast')
%         
%         subplot(3, 4, 12);
%         plot(X_scale_l(:, 1), Y_scale_nl(:, 1), 'b');
%         hold on
%         plot(X_scale_l(:, 2), Y_scale_nl(:, 2), 'r');
%         xlabel('generator signal')
%         ylabel('spike rate(HZ)')             
%         
%     end
