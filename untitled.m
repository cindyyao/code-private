movie_path = '/Volumes/lab/acquisition/movie-xml/BW-40-2.0.48-11111-15x15.xml';
display_frame_rate = 60.35;
mov = get_movie(movie_path, datarun.triggers, 5);

figure
for i = 1:5
    subplot(2, 3, i)
    colormap gray
    imshow(mov(:, :, 1, i))
end


%%

opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);
datarun = load_data('/Volumes/lab/Experiments/Array/Analysis/2013-02-21-0/data004/data004', opt);

idx = 422;
sta = datarun.stas.stas{idx};
im = norm_image(sta);

for i = 1:15
   FigHandle = figure;
   set(FigHandle, 'Position', [1, 1, 310, 300]);
   image(matrix_scaled_up(im(:, :, :, i), 8))
   axis off
end


%%

for ct = 1:1
    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 1150, 700]);
    
    % NDF 4 sta image
    r = datarun{1}.stimulus.field_height;
    sta = sta_c{ct, 1};
    if cell_type{ct}(2) == 'N' 
        [~, b] = max(sta(:));
    else
        [~, b] = min(sta(:));  
    end

    z = ceil(b/r^2);
    sta_1 = sta(:, :, 1, z);
    sta_1 = matrix_scaled_up(sta_1, 8);
    
    subplot(2, 3, 1)
    im = norm_image(sta_1);
    image(im)
    title('Rod')
    axis off
    
    % NDF 0 sta image

    r = datarun{2}.stimulus.field_height;
    sta = sta_c{ct, 2};
    if cell_type{ct}(2) == 'N' 
        [~, b] = max(sta(:));
    else
        [~, b] = min(sta(:));  
    end

    z = ceil(b/r^2);
    sta_1 = sta(:, :, 1, z);
    sta_1 = matrix_scaled_up(sta_1, 8);
    
    subplot(2, 3, 2)
    im = norm_image(sta_1);
    image(im)
    title('Cone')
    axis off
    
    subplot(2, 3, 3)
    plot_rf_summaries(datarun{1}, cell_id{ct}(:, 1), 'fit_color', 'b', 'scale', 4/3)
    hold on
    plot_rf_summaries(datarun{2}, cell_id{ct}(:, 2), 'fit_color', 'r')
    axis off
    title('Receptive Fields Mosaic')
    
%     % bar
%     subplot(2, 3, [7:9])
%     
%     xtick = {'Surround Strengths', 'Center Size', 'Degree of Transience'};
%     model_series = [surround.mean(ct, 1) surround.mean(ct, 2); center.mean(ct, 1)/2 center.mean(ct, 2)/2; dot.mean(ct, 1) dot.mean(ct, 2)];   
%     model_error = [surround.stev(ct, 1) surround.stev(ct, 2); center.stev(ct, 1)/2 center.stev(ct, 2)/2; dot.stev(ct, 1) dot.stev(ct, 2)];
%     h = bar(model_series);
%     set(h,'BarWidth',1); % The bars will now touch each other
% 
%     set(gca,'XTicklabel',xtick)
%     % ylabel('Surround Strengths')
%     legend('NDF 4','NDF 0', 'location', 'northwest');
%     hold on;
%  
%     numgroups = size(model_series, 1); 
%     numbars = size(model_series, 2); 
% 
%     groupwidth = min(0.8, numbars/(numbars+1.5));
%     
%     for i = 1:numbars
%     % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
%     x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
%     errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
%     end
%     
%     title(cell_type{ct})

    % spatial profile
    subplot(2, 3, 4)
    dis = -10:10;
    plot(dis, Spatial_Profile{ct}(:, 1), 'b', dis, Spatial_Profile{ct}(:, 2), 'r');
    xlim([-10 10])
    xlabel('x(pixel)')
    ylabel('D_s(a.u.)')
    title('Spatial Filter')
%    h_legend = legend('NDF 4', 'NDF 0', 'location', 'northwest');
%     set(h_legend,'FontSize',14);


    % temporal profile
    subplot(2, 3, 5)
    t1 = -0.066*[14:-1:0];
    t2 = -0.033*[14:-1:0];
    plot(t1, Temporal_Filter{ct}(1, :), 'b', t2, Temporal_Filter{ct}(2, :), 'r');
    xlim([-0.75 0])
    xlabel('time(s)')
    ylabel('D_t(a.u.)')
    title('Temporal Filter')
    
    % Nonlinearity
    subplot(2, 3, 6)
    errorbar(Nonlinearity_mean(:, 1, ct, 1), Nonlinearity_mean(:, 1, ct, 2), Nonlinearity_stev(:, 1, ct), 'b');
    hold on
    errorbar(Nonlinearity_mean(:, 2, ct, 1), Nonlinearity_mean(:, 2, ct, 2), Nonlinearity_stev(:, 2, ct), 'r');
    title('Average Nonlinearity')
    legend('rod', 'cone', 'location', 'northwest');
    xlabel('Normalized generator signal')
    ylabel('Normalized firing rate')

    xlim([-2.5 2.5])
    
    % individual nonlinearity
%     subplot(2, 4, 7)
%     n = size(cell_id{ct}, 1);
%     for j = 1:n
%         plot(NL{ct}(1, :, j, 1), NL{ct}(2, :, j, 1), 'b', NL{ct}(1, :, j, 2), NL{ct}(2, :, j, 2), 'r');
%         xlim([-2.5 2.5]) 
%         hold on
%     end
%     xlim([-2.5 2.5])    
%     title('Nonlinearity')
end

%%

%% bar graph
% surround strengths
    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 900, 500]);

    xtick = cell_type;
    model_series = [surround.mean(1, 1) surround.mean(1, 2); surround.mean(2, 1) surround.mean(2, 2); surround.mean(3, 1) surround.mean(3, 2); surround.mean(4, 1) surround.mean(4, 2)];   
    model_error = [surround.stev(1, 1) surround.stev(1, 2); surround.stev(2, 1) surround.stev(2, 2);surround.stev(3, 1) surround.stev(3, 2);surround.stev(4, 1) surround.stev(4, 2);];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('Surround Center Volume Ratio')
    legend('Scotopic','Photopic', 'location', 'northwest');
    hold on;
 
    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 

    groupwidth = min(0.8, numbars/(numbars+1.5));
    
    for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end
    
    title('Surround Strengths')

% center size

    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 700,700]);

    xtick = cell_type;
    model_series = [center.mean(1, 1) center.mean(1, 2); center.mean(2, 1) center.mean(2, 2); center.mean(3, 1) center.mean(3, 2); center.mean(4, 1) center.mean(4, 2)];   
    model_error = [center.stev(1, 1) center.stev(1, 2); center.stev(2, 1) center.stev(2, 2); center.stev(3, 1) center.stev(3, 2);center.stev(4, 1) center.stev(4, 2);];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('Center Size(a.u.)')
    legend('Scotopic','Photopic', 'location', 'northwest');
    hold on;
 
    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 

    groupwidth = min(0.8, numbars/(numbars+1.5));
    
    for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end
    
    title('Center Size')

    % DOT
    
    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 900, 500]);

    xtick = cell_type;
    model_series = [dot.mean(1, 1) dot.mean(1, 2); dot.mean(2, 1) dot.mean(2, 2); dot.mean(3, 1) dot.mean(3, 2); dot.mean(4, 1) dot.mean(4, 2)];   
    model_error = [dot.stev(1, 1) dot.stev(1, 2); dot.stev(2, 1) dot.stev(2, 2); dot.stev(3, 1) dot.stev(3, 2);dot.stev(4, 1) dot.stev(4, 2);];
    h = bar(model_series);
    set(h,'BarWidth',1); % The bars will now touch each other

    set(gca,'XTicklabel',xtick)
    ylabel('DOT')
    legend('Scotopic','Photopic', 'location', 'northwest');
    hold on;
 
    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 

    groupwidth = min(0.8, numbars/(numbars+1.5));
    
    for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
    end
    
    ylim([0 1.2])
    title('Degree of Transience')
    legend('Scotopic','Photopic', 'location', 'northwest');

%%

figure
for i = 1:4
    subplot(2, 2, i)
    plot(Nonlinearity_mean(:, 1, i, 1), Nonlinearity_mean(:, 1, i, 2), 'b');
    hold on
    plot(Nonlinearity_mean(:, 2, i, 1), Nonlinearity_mean(:, 2, i, 2), 'r');
    xlim([-2.5 2.5])
    ylim([0 1])
end

%%
figure
x = [2 1 -1 -2 -1 1 2];
y = [0 sqrt(3) sqrt(3) 0 -sqrt(3) -sqrt(3) 0];
plot(x, y)

%%
figure
for i = 1:4
subplot(2, 2, i)
    plot_rf_summaries(datarun{2}, cell_id{i}(:, 2), 'fit_color', 'k')
    axis off
end

