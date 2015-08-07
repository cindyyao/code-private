function[] = polar_plots_all(U,V,num_t, mag)

%Input X and Y vectors of all cells for all stimuli for all speeds
%Plots polar plot for all speeds for all cells

% Sneha Ravi 
% Last revision: 12-18-2012

for f = 1:size(U, 2)
    figure
    axes_handle = [];
    ylim = zeros(size(U,1));
    X = floor(sqrt(size(U,1)));
    for i = 1:size(U,1)
        ax(i) = subplot(X,ceil(size(U,1)/X),i);
        set(ax(i), 'FontName', 'AvantGarde', 'FontSize', 18)
        h1 = compass(U{i,f},V{i,f});
        h = findall(gca, 'type', 'line'); % find all lines in the current figure
        h(h==ax(i)) = []; % remove the line you ploted from the list.
        set(h, 'LineStyle', '--');
        set(h1,'linestyle','-', 'LineWidth', 1); 
        axes_handle = [axes_handle ax(i)];
        ylim(i) = max(mag{i,f});
        titlechar = [num2str(num_t(i))];
        %titlechar = ['Vector Averages Plot for all cells of temporal period: ' num2str(num(i))];
        %title(titlechar, 'Color', 'k', 'FontWeight', 'bold' , 'FontSize', 18 ,'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Right');
        title(titlechar);
    end

    %linkaxes(axes_handle,'xy'); %Link axes to have same x and y limits
    %[C, I] = max(ylim);
    %set(ax(I),'xlimmode','auto');
    %set(ax(I),'ylimmode','auto');
end

end


%Return vector averages for each speed for each cell to next function
