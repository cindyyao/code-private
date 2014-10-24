opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1, 'load_sta', 1);

% ndf 4 melatonin

datarun{1} = load_data('/Analysis/xyao/2013-07-22-0/data002/data002', opt);
datarun{1} = load_java_movie(datarun{1}, '/Volumes/lab/acquisition/movie-xml/BW-15-2-0.48-11111-40x40-60.35.xml');
datarun{1} = get_sta_summaries(datarun{1}, 'all');

datarun{1} = get_snls(datarun{1}, 'all');


% ndf 4 non-melatonin

datarun{2} = load_data('/Analysis/xyao/2013-07-22-0/data005/data005', opt);
datarun{2} = load_java_movie(datarun{2}, '/Volumes/lab/acquisition/movie-xml/BW-15-2-0.48-11111-40x40-60.35.xml');
datarun{2} = get_sta_summaries(datarun{2}, 'all');

datarun{2} = get_snls(datarun{2}, 'all');



% ndf 0 non-melatonin

datarun{3} = load_data('/Analysis/xyao/2013-07-22-0/data010/data010', opt);
datarun{3} = load_java_movie(datarun{3}, '/Volumes/lab/acquisition/movie-xml/BW-10-2-0.48-11111-60x60-60.35.xml');
datarun{3} = get_sta_summaries(datarun{3}, 'all');

datarun{3} = get_snls(datarun{3}, 'all');

%% plot snls

cell_type = {'ON brisk transient', 'ON transient', 'OFF brisk transient', ...
    'OFF transient', 'OFF slow', 'OFF sustained'};
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
[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun3, ...
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comparison between melatonin & non-melatonin


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
        title(cell_type{ct})
    end
    
            
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comparison between ndf 4 & ndf 0

n = length(cell_type);
[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map2(datarun3, ...
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
        legend('ndf 0', 'ndf 4')
        title(cell_type{ct})
    end
    
            
end

