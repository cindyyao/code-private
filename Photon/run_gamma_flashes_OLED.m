function run_gamma_flashes_OLED(params, steps, guns)

if ischar(params)
    
    linear_gamma = linspace(0,1,256);
    mglSetGammaTable(linear_gamma,linear_gamma,linear_gamma);
    
elseif length(params)==9
    
    scale = params(1:3);
    power = params(4:6);
    offset = params(7:9);
    
    set_gamma_from_fit_params(scale, power, offset);
    
else
    
    fprintf('\n\nNeed 9 numbers: scale x 3 , power x 3 , offset x 3\n\n')
    fprintf('GAMMA was NOT SET\n\n')
    return
    
end

coef = linspace(0,1,steps);
seq = get_sequence(steps);

for i = 1:steps   
    mglClearScreen(guns*coef(seq(i)));
    mglFlush
    pause;
end

