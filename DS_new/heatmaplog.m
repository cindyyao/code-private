% heat map
conditions = {'control', 'AP5'};
[xx, yy] = meshgrid(theta, ctr_x);
for drug = 1:2
    figure
    response_max_all_mean{drug} = squeeze(mean(response_max_all{drug}));
%     s = surf(xx, yy, response_max_all_mean{drug}'/max(response_max_all_mean{1}(:)));
    response_log = log10(abs(response_max_all_mean{drug}));
%     if drug == 1
%         response_min = min(response_log(:));
%     end
%     response_log = abs(response_log - response_min);
%     if drug == 1
%         response_max = max(response_log(:));
%     end
%     response_log = response_log/response_max;
    if drug == 1
        r_min = min(response_log(:));
        r_max = max(response_log(:));
    end
    s = surf(xx, yy, response_log');
    s.EdgeColor = 'none';
    set(gca, 'yscale', 'log')
%     xlabel('direction (degree)')
%     ylabel('contrast')
%     zlabel('spike number')
%     title(conditions{drug})
    xlim([min(xx(:)) max(xx(:))])
    ylim([min(yy(:)) max(yy(:))])
    caxis([-3.5 1.5])
    axis off
%     colorbar
end
