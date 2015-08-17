function center_location = correct_ei_center(datawn, id_wn, center_location, stixel_size)

center = center_location/stixel_size;

for i = 1:length(id_wn)
    if ~isempty(id_wn{i})
%         figure
%         subplot(2,1,1)
%         plot_ei(datawn, id(i), 'alpha', false)
%         title('white noise')
%         subplot(2,1,2)
%         plot_ei(datadg, id(i), 'alpha', false)
%         title('gratings/bars')
%         set(gcf, 'Position', [1000, 1, 500, 1000])
        
        figure
        plot_rf(datawn, id_wn{i}, 'fit', false);
        hold on
        plot(center(i, 1), center(i, 2), 'ro')
        set(gcf, 'Position', [1 1 1000 1000])
        a = ginput;
        if ~isempty(a)
            center(i, :) = a;
        end
        close all
    end
end

center_location = center*stixel_size;

end