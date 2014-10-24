% circadian night
clear datarun

datarun = load_data('/Volumes/lab/Analysis/2013-07-01-0/data008/data008', opt);
datarun = load_java_movie(datarun, '/Volumes/lab/acquisition/movie-xml/BW-10-1-0.48-11111-60x60-60.35.xml');
datarun = get_sta_summaries(datarun, 'all');

datarun = get_snls(datarun, 'all');
save('snls130701.mat', 'datarun')
 

%%

display_rate = 60.35;
refresh_rate = 1;

cell_type = {'ON brisk transient', 'ON transient', 'ON slow transient', 'OFF brisk transient', ...
    'OFF brisk transient2', 'OFF brisk transient large', 'OFF transient', 'OFF slow transient', 'OFF sustained'};
n = length(cell_type);

for i =  1:n
    cell_id{i} = get_cell_ids(datarun, cell_type{i});
    cell_idx{i} = get_cell_indices(datarun, cell_type{i});
end

%%

ct = 9;
cell_num = length(cell_id{ct});

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1300, 900]);


    for cn = 1:5 %cell_num
          
            a = [16:20];
            idx = cell_idx{ct}(a(cn));
            com = datarun.stas.rf_coms{idx};
            gen_signals = datarun.stas.snls{idx}.gen_signal;
            spikes = datarun.stas.snls{idx}.spikes;
            tc = datarun.stas.time_courses{idx};
            r = datarun.stimulus.field_height;
            
            rf = datarun.stas.rfs{idx};
            pix = [1:r] + 0.5;
            [XX YY] = meshgrid(pix);
            xdis = com(1) - XX;
            ydis = com(2) - YY; 

            dis = zeros(r, r);
            for k = 1:r
                for j = 1:r
                    dis(k, j) = norm([xdis(k, j) ydis(k, j)]);
                end
            end
            
            dis_1 = reshape(dis, 1, r^2);
            dis_1 = dis_1*40/r;
            rf_1 = reshape(rf, 1, r^2);
            [Dis1 RF1] = curve_from_binning(dis_1, rf_1, 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:20);
            Dis1 = [-Dis1(end:-1:1); Dis1];
            RF1 = [RF1(end:-1:1); RF1];

            [X, Y] = curve_from_binning(gen_signals, spikes, 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
            Y = Y*display_rate/refresh_rate;

        subplot(4, 5, cn)
        im = norm_image(rf);
        image(im)
        title(num2str(cell_id{ct}(a(cn))))

        
        subplot(4, 5, cn+5);
        plot(Dis1, RF1);
        xlim([-20 20])
        ylim([-max(RF1)*0.2 max(RF1)*1.1])
        xlabel('distance(a.u.)')
        ylabel('STA contrast')
        
        subplot(4, 5, cn+10);
        t = -1/display_rate*refresh_rate*[30:-1:1];
        plot(t, tc);
        xlabel('time to spike (s)')
        ylabel('STA contrast')
        xlim([-0.5 0])
        ylim([min(tc)*1.2 max(tc)*1.2])


        subplot(4, 5, cn+15);
        plot(X, Y);
        xlabel('generator signal')
        ylabel('spike rate(HZ)')
        if cn == 1
            title(cell_type{ct})
        end
        ylim([0 max(Y)*1.2])

%         h_legend = legend('NDF 4', 'NDF 0', 'location', 'northwest');
%         set(h_legend,'FontSize',10);

 
        
        
    end


    
