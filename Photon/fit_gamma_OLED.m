function [scale, power, offset] = fit_gamma_OLED(data)

xdata = linspace(1e-10,1,size(data,1));
F = @(x, xdata)(min(max(x(1)*xdata.^x(2)+x(3), 0), 1));
cols = 'rgb';

for i=1:size(data,2)
    ydata = data(:,i);
%     y = y-min(y);
    ydata = ydata/max(ydata);
    x0 = [1 1 0];
    [fit_results,~,~,~,output] = lsqcurvefit(F, x0, xdata', ydata);
    scale(i) = fit_results(1);
    power(i) = fit_results(2);
    offset(i) = fit_results(3);    
    
    subplot(2,2,i)
    plot(xdata, ydata, 'color', 'k', 'marker','*');
    hold on
    plot(xdata, min(max(scale(i).*xdata.^power(i)+offset(i), 0), 1),cols(i));
end
