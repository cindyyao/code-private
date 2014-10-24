%% scaled by std/vector 
datarun{1} = datarun2;
datarun{2} = datarun3;
datarun_2 = datarun{2};

display_rate = 60.35;
refresh_rate = 2;
diff_x = -7.2237;
diff_y = 0.3277;

cell_type = {'ON brisk transient', 'ON transient', 'OFF brisk transient', ...
    'OFF transient', 'OFF slow'};
n = length(cell_type);
[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun{1}, ...
    datarun{2}, cell_type);
%%
ct = 1;
cell_num = size(cell_id{ct}, 1);

    for cn = 1:1 %cell_num          
        FigHandle = figure;
        set(FigHandle, 'Position', [100, 100, 1250, 800]);
        
        for i = 1:2
            idx = cell_idx{ct}(cn, i);
            com = datarun{i}.stas.rf_coms{idx};
            gen_signals(:, i) = datarun{i}.stas.snls{idx}.gen_signal;
            spikes(:, i) = datarun{i}.stas.snls{idx}.spikes;
            tc(:, i) = datarun{i}.stas.time_courses{idx};
            r = datarun{i}.stimulus.field_height;
            
            rf{i} = datarun{i}.stas.rfs{idx};
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
            rf_1 = reshape(rf{i}, 1, r^2);
            [Dis1{i} RF1{i}] = curve_from_binning(dis_1, rf_1, 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:20);
            Dis1{i} = [-Dis1{i}(end:-1:1); Dis1{i}];
            RF1{i} = [RF1{i}(end:-1:1); RF1{i}];
            
            [X(:, i), Y(:, i)] = curve_from_binning(gen_signals(:, i), spikes(:, i), 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
            Y(:, i) = Y(:, i)*display_rate/refresh_rate;

        end
        
        % unscaled
        
        subplot(2, 4, 1)
        im = norm_image(rf{1});
        image(im)
        title([num2str(cell_id{ct}(cn, 1)) '  NDF 4'])

 
        subplot(2, 4, 5)
        im = norm_image(rf{2});
        image(im)
        title([num2str(cell_id{ct}(cn, 2)) '  NDF 0'])

        subplot(2, 4, 2);
        plot(Dis1{1}, RF1{1}, 'b');
        hold on
        plot(Dis1{2}, RF1{2}, 'r');
        xlim([-20 20])
        
        subplot(2, 4, 3);
        t = -0.033*[15:-1:1];
        plot(t, tc(:, 1), 'b');
        hold on
        plot(t, tc(:, 2), 'r');
        xlabel('time to spike (s)')
        ylabel('STA contrast')
        xlim([-0.5 0])

        subplot(2, 4, 4);
        plot(X(:, 1), Y(:, 1), 'b');
        hold on
        plot(X(:, 2), Y(:, 2), 'r');
        xlabel('generator signal')
        ylabel('spike rate(HZ)')
        h_legend = legend('NDF 4', 'NDF 0', 'location', 'northwest');
%         set(h_legend,'FontSize',12);

        % calculate the scaling factor
        sigma_1 = std(gen_signals(:, 1));
        sigma_2 = std(gen_signals(:, 2));
        f_gs = sigma_1/sigma_2;
        
        % apply this factor to strf and recalculate the sf and tf
        datarun_2.stas.stas{idx} = datarun{2}.stas.stas{idx}*f_gs;
        datarun_2 = get_sta_summaries(datarun_2, datarun_2.cell_ids(idx));
        
        
        gen_signals(:, 2) = gen_signals(:, 2)*f_gs;
        tc(:, 2) = datarun_2.stas.time_courses{idx};
        rf{2} = datarun_2.stas.rfs{idx};
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
        rf_1 = reshape(rf{2}, 1, r^2);
        [Dis1{2} RF1{2}] = curve_from_binning(dis_1, rf_1, 'average_y', 'mean','average_x', 'mean', 'bin_edges', 0:20);
        Dis1{2} = [-Dis1{2}(end:-1:1); Dis1{2}];
        RF1{2} = [RF1{2}(end:-1:1); RF1{2}];
        
        RF1_std{1} = RF1{1}/std(RF1{1});
        RF1_std{2} = RF1{2}/std(RF1{2});

        tc_std(:, 1) = tc(:, 1)/std(tc(:, 1));
        tc_std(:, 2) = tc(:, 2)/std(tc(:, 2));

        [X(:, 2), Y(:, 2)] = curve_from_binning(gen_signals(:, 2), spikes(:, 2), 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
        Y(:, 2) = Y(:, 2)*display_rate/refresh_rate;
        
        subplot(2, 4, 6);
        plot(Dis1{1}, RF1_std{1}, 'b');
        hold on
        plot(Dis1{2}, RF1_std{2}, 'r');
        xlim([-20 20])

        subplot(2, 4, 7);
        plot(t, tc_std(:, 1), 'b');
        hold on
        plot(t, tc_std(:, 2), 'r');
        xlabel('time to spike (s)')
        ylabel('STA contrast')
        xlim([-0.5 0])
        ylim([min(tc_std(:))*1.2 max(tc_std(:))*1.2])


        subplot(2, 4, 8);
        plot(X(:, 1), Y(:, 1), 'b');
        hold on
        plot(X(:, 2), Y(:, 2), 'r');
        xlabel('generator signal')
        ylabel('spike rate(HZ)')
          
    end


    
%% NDF 4

n = length(cell_type);

for i =  1:n
    cell_id{i} = get_cell_ids(datarun{2}, cell_type{i});
    cell_idx{i} = get_cell_indices(datarun{2}, cell_type{i});
end

ct = 5;
cell_num = length(cell_id{ct});

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1300, 900]);



    for cn = 1:5 %cell_num
          
            a = [6 7 8 9 11];
            idx = cell_idx{ct}(a(cn));
            com = datarun{2}.stas.rf_coms{idx};
            gen_signals = datarun{2}.stas.snls{idx}.gen_signal;
            spikes = datarun{2}.stas.snls{idx}.spikes;
            tc = datarun{2}.stas.time_courses{idx};
            r = datarun{2}.stimulus.field_height;
            
            rf = datarun{2}.stas.rfs{idx};
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
        title(num2str(cell_id{ct}(a(cn))))
        xlim([-20 20])
        ylim([-max(RF1)*0.2 max(RF1)*1.1])
        xlabel('distance(a.u.)')
        ylabel('STA contrast')
        
        subplot(4, 5, cn+10);
        t = -1/display_rate*refresh_rate*[15:-1:1];
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

