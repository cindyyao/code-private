% script for fitting STAs in matlab

datarun = load_data('/lab/analysis/2013-02-21-0/data004/data004'); % WT mouse photopic
datarun = load_data('/lab/analysis/2013-02-21-0/data000/data000'); % WT mouse scotopic
datarun = load_data('/lab/analysis/2013-02-21-0/data002/data002'); % WT mouse mesopic


datarun = load_data('/lab/analysis/2013-02-14-0/data006/data006'); %GABA vesicular transport KO
datarun = load_data('/lab/analysis/2013-02-14-0/data000/data000'); %GABA vesicular transport KO

datarun = load_neurons(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_params(datarun);

% compute marks, tcs and RFs
marks_params.thresh = 4.0;
datarun = get_sta_summaries(datarun, 'all','marks_params',marks_params);

% plot a random RF
cell_types = {'off transient', 'off brisk transient', 'on transient', 'on brisk transient'};
temp_indices = get_cell_indices(datarun, cell_types);
cell_ids = datarun.cell_ids(temp_indices);

datarun.matlab.sta_fits = cell(length(datarun.cell_ids), 1);

for cc = 1:length(temp_indices)

    temp_sta = datarun.stas.stas{temp_indices(cc)};


    %% hold spatial params fixed
    clear params
    params.fit_center_point_x = false;
    params.fit_center_point_y = false;
    params.fit_center_sd_x = false;
    params.fit_center_sd_y = false;
    params.fit_center_rotation_angle = false;
    params.fit_center_amp_scale = false;
    params.fit_surround = false;
    params.fit_surround_sd_scale = false;
    params.fit_surround_amp_scale = false;
    params.verbose = true;

    params.sig_stixels = datarun.stas.marks{temp_indices(cc)};

    % fit an STA with matlab code.
    temp_fit = fit_sta(temp_sta, params);

    %% use previous fit for temporal params and hold them fixed -- fit center
    clear params
    % set params
    params.initial_scale_one = temp_fit.scale_one;
    params.initial_scale_two = temp_fit.scale_two;
    params.initial_tau_one = temp_fit.tau_one;
    params.initial_tau_two = temp_fit.tau_two;
    params.initial_n_filters = temp_fit.n_filters;

    % fix params
    params.fit_scale_one = false;
    params.fit_scale_two = false;
    params.fit_tau_one = false;
    params.fit_tau_two = false;
    params.fit_n_filters = false;
    params.verbose = false;

    % set spatial surround to 0 and don't fit
    params.initial_surround_sd_scale = 0;
    params.initial_surround_amp_scale = 0;
    params.fit_surround = false;
    params.fit_surround_sd_scale = false;
    params.fit_surround_amp_scale = false;

    params.sig_stixels = datarun.stas.marks{temp_indices(cc)};
    
    temp_fit = fit_sta(temp_sta, params);

    %% Introduce surrround to fit while keeping everything else constant 
    clear params
    % set center parameters and hold fixed
    params.initial_center_point_x = temp_fit.center_point_x;
    params.initial_center_point_y = temp_fit.center_point_y;
    params.initial_center_sd_x = temp_fit.center_sd_x;
    params.initial_center_sd_y = temp_fit.center_sd_y;
    params.initial_center_rotation_angle = temp_fit.center_rotation_angle;
    params.fit_center_point_x = false;
    params.fit_center_point_y = false;
    params.fit_center_sd_x = false;
    params.fit_center_sd_y = false;
    params.fit_center_rotation_angle = false;
    params.fit_center_amp_scale = false;

    %set temporal parameters and hold fixed
    params.initial_scale_one = temp_fit.scale_one;
    params.initial_scale_two = temp_fit.scale_two;
    params.initial_tau_one = temp_fit.tau_one;
    params.initial_tau_two = temp_fit.tau_two;
    params.initial_n_filters = temp_fit.n_filters;
    params.fit_scale_one = false;
    params.fit_scale_two = false;
    params.fit_tau_one = false;
    params.fit_tau_two = false;
    params.fit_n_filters = false;

    % fit surround
    params.initial_surround_sd_scale = 1.5;
    params.initial_surround_amp_scale = 0.1;
    params.fit_surround = true;
    params.fit_surround_sd_scale = true;
    params.fit_surround_amp_scale = true;
    params.verbose = false;

    params.sig_stixels = datarun.stas.marks{temp_indices(cc)};

    temp_fit = fit_sta(temp_sta, params);


    %% Use these above fit results as initial conditions and allow the entire fit to vary
    clear params
    params.initial_center_point_x = temp_fit.center_point_x;
    params.initial_center_point_y = temp_fit.center_point_y;
    params.initial_center_sd_x = temp_fit.center_sd_x;
    params.initial_center_sd_y = temp_fit.center_sd_y;
    params.initial_center_rotation_angle = temp_fit.center_rotation_angle;
    params.fit_center_point_x = true;
    params.fit_center_point_y = true;
    params.fit_center_sd_x = true;
    params.fit_center_sd_y = true;
    params.fit_center_rotation_angle = true;
    params.fit_center_amp_scale = true;

    %set temporal parameters and hold fixed
    params.initial_scale_one = temp_fit.scale_one;
    params.initial_scale_two = temp_fit.scale_two;
    params.initial_tau_one = temp_fit.tau_one;
    params.initial_tau_two = temp_fit.tau_two;
    params.initial_n_filters = temp_fit.n_filters;
    params.fit_scale_one = true;
    params.fit_scale_two = true;
    params.fit_tau_one = true;
    params.fit_tau_two = true;
    params.fit_n_filters = true;

    % fit surround
    params.initial_surround_sd_scale = temp_fit.surround_sd_scale;
    params.initial_surround_amp_scale = temp_fit.surround_amp_scale;
    params.fit_surround = true;
    params.fit_surround_sd_scale = true;
    params.fit_surround_amp_scale = true;
    params.verbose = true;

    params.sig_stixels = datarun.stas.marks{temp_indices(cc)};    
    
    temp_fit = fit_sta(temp_sta, params);

    full_fit = parse_sta_fit(temp_fit);
    figure(1)
    subplot(2,1,1)
    plot_rf(datarun, cell_ids(cc))
    subplot(2,1,2)

    imagesc(norm_image(full_fit(:,:,1,13)))
    axis square
    drawnow

    datarun.matlab.sta_fits{temp_indices(cc)} = temp_fit;

end




