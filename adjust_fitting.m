sig_stixels = significant_stixels(datarun1.stas.stas{1}, 'thresh', 5, 'time', 'max');
tc = time_course_from_sta(datarun1.stas.stas{1}, sig_stixels);
tc = tc/max(abs(max(tc)), abs(min(tc)));

% params
p1 = 2.1060e+05;
p2 = -2.1051e+05;
tau1 = 5.5970;
tau2 = 5.6002;
n1 = 12.8990;



t = 0:14;
fit = p1*(t/tau1).^n1.*exp(-n1*(t/tau1) - 1) + p2*(t/tau2).^n1.*exp(-n1*(t/tau2) - 1);
ff = 1/max(abs(max(fit)), abs(min(fit)));
figure
plot(tc)
hold on
plot(ff*fit(15:-1:1), 'r')



%%

% set thresh for pick out sig stixel.
mark_params.thresh = 4.0;

% fit sta sequentially 
temp_sta = sta_c{4,2};
    
    % hold spatial params fixed
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
        params.verbose = false;
%         params.initial_center_point_x = 14;
%         params.initial_center_point_y = 14;
        params.mark_params = mark_params;

        % fit an STA with matlab code.
        
            temp_fit = fit_sta(temp_sta, params);
      
        
    
    % use previous fit for temporal params and hold them fixed -- fit center
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
            params.verbose = true;

            % set spatial surround to 0 and don't fit
            params.initial_surround_sd_scale = 0;
            params.initial_surround_amp_scale = 0;
            params.fit_surround = false;
            params.fit_surround_sd_scale = false;
            params.fit_surround_amp_scale = false;

            params.mark_params = mark_params;
    
                temp_fit = fit_sta(temp_sta, params);
           
    
    % Introduce surrround to fit while keeping everything else constant 
        clear params
            % set center parameters 
            params.initial_center_point_x = temp_fit.center_point_x;
            params.initial_center_point_y = temp_fit.center_point_y;
            params.initial_center_sd_x = temp_fit.center_sd_x;
            params.initial_center_sd_y = temp_fit.center_sd_y;
            params.initial_center_rotation_angle = temp_fit.center_rotation_angle;
            
  
            %set temporal parameters 
            params.initial_scale_one = temp_fit.scale_one;
            params.initial_scale_two = temp_fit.scale_two;
            params.initial_tau_one = temp_fit.tau_one;
            params.initial_tau_two = temp_fit.tau_two;
            params.initial_n_filters = temp_fit.n_filters;
        
        % hold temporal and center params fixed
        params.fit_center_point_x = false;
        params.fit_center_point_y = false;
        params.fit_center_sd_x = false;
        params.fit_center_sd_y = false;
        params.fit_center_rotation_angle = false;
        params.fit_center_amp_scale = false;
            
        params.fit_scale_one = false;
        params.fit_scale_two = false;
        params.fit_tau_one = false;
        params.fit_tau_two = false;
        params.fit_n_filters = false;
            
        % fit surround
        params.initial_surround_sd_scale = 1.1;
        params.initial_surround_amp_scale = 0.2;
        params.fit_surround = true;
        params.fit_surround_sd_scale = true;
        params.fit_surround_amp_scale = true;
        params.verbose = true;

        params.mark_params = mark_params;

                temp_fit = fit_sta(temp_sta, params);
           
 
  % Use these above fit results as initial conditions, then refit temporal
        clear params
        % set fit params
        params.fit_center_point_x = false;
        params.fit_center_point_y = false;
        params.fit_center_sd_x = false;
        params.fit_center_sd_y = false;
        params.fit_center_rotation_angle = false;
        params.fit_center_amp_scale = false;
        params.fit_surround = true;
        params.fit_surround_sd_scale = false;
        params.fit_surround_amp_scale = false;
        
        % set initial condition
        params.initial_scale_one = temp_fit.scale_one;
        params.initial_scale_two = temp_fit.scale_two;
        params.initial_tau_one = temp_fit.tau_one;
        params.initial_tau_two = temp_fit.tau_two;
        params.initial_n_filters = temp_fit.n_filters;
        
            params.initial_center_point_x = temp_fit.center_point_x;
            params.initial_center_point_y = temp_fit.center_point_y;
            params.initial_center_sd_x = temp_fit.center_sd_x;
            params.initial_center_sd_y = temp_fit.center_sd_y;
            params.initial_center_rotation_angle = temp_fit.center_rotation_angle;
            params.initial_surround_sd_scale = temp_fit.surround_sd_scale;
            params.initial_surround_amp_scale = temp_fit.surround_amp_scale;
        
        params.verbose = true;
        params.mark_params = mark_params;
        
        % fit an STA with matlab code.
                temp_fit = fit_sta(temp_sta, params);
    
  % Use these above fit results as initial conditions and allow the entire fit to vary
    clear params
    params.fit_center_point_x = false;
    params.fit_center_point_y = false;
    params.fit_center_sd_x = false;
    params.fit_center_sd_y = false;
    params.fit_center_rotation_angle = false;
        
    params.fit_scale_one = false;
    params.fit_scale_two = false;
    params.fit_tau_one = false;
    params.fit_tau_two = false;
    params.fit_n_filters = false;
    
        % set center parameters
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
        % set temporal parameters 
        params.initial_scale_one = temp_fit.scale_one;
        params.initial_scale_two = temp_fit.scale_two;
        params.initial_tau_one = temp_fit.tau_one;
        params.initial_tau_two = temp_fit.tau_two;
        params.initial_n_filters = temp_fit.n_filters;
        
%         params.initial_scale_one = 2.1060e+05;
%         params.initial_scale_two = -2.1051e+05;
%         params.initial_tau_one = 5.5970;
%         params.initial_tau_two = 5.6002;
%         params.initial_n_filters = 12.8990;
        
        params.fit_scale_one = true;
        params.fit_scale_two = true;
        params.fit_tau_one = true;
        params.fit_tau_two = true;
        params.fit_n_filters = true;

        % set surround parameters
%        params.initial_surround_sd_scale = temp_fit.surround_sd_scale;
%        params.initial_surround_amp_scale = temp_fit.surround_amp_scale;
       
       params.initial_surround_sd_scale = 1.1;
       params.initial_surround_amp_scale = 0.05;

       params.fit_surround = true;
       params.fit_surround_sd_scale = true;
       params.fit_surround_amp_scale = true;
    
       params.verbose = true;

       params.mark_params = mark_params;
      
                temp_fit = fit_sta(temp_sta, params);

 
  
       fit_info = temp_fit;

