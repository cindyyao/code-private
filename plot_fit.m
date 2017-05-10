function [] = plot_fit(datarun, id)  
    idx = get_cell_indices(datarun, id);
    sta = datarun.stas.stas{idx};
    params = get_params(datarun, id);

    sta_fit = sta_fit_function(params);

    % spatial fit
    figure(1)
    subplot(2,1,1)
    hold off
    temp_rf = rf_from_sta(sta);
    imagesc(norm_image(temp_rf))
    hold on
    plot_spatial_sd(params)
    axis image
    drawnow
    title([num2str(id) '  rmse:' num2str(params(21))])

    % temporal fit

    temp_stix = significant_stixels(sta, 'time', 'max', 'select', 'max', 'thresh', 3.5, 'robust_std_method', 1);
    fit_tc = time_course_from_sta(sta_fit, temp_stix);
    norm_factor = max(abs(reshape(fit_tc, 1, [])));
    fit_tc = fit_tc ./ norm_factor;
    subplot(2,1,2)
    plot(fit_tc, '--k')
    hold on

    real_stix = significant_stixels(sta, 'time', 'max', 'select', 'max', 'thresh', 3.5, 'robust_std_method', 1);
    tc = time_course_from_sta(sta, real_stix);
    norm_factor = max(abs(reshape(tc, 1, [])));
    tc = tc ./ norm_factor;
    plot(tc, 'k')
    hold off

    drawnow
