function [F, G] = fit_gaussian(xdata, ydata)

gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
d0 = min(ydata);
c0 = 1.5;
[~, i] = max(ydata);
b0 = xdata(i);
a0 = 1;
startPoints = [a0 b0 c0 d0];

s = fitoptions('METHOD', 'NonlinearLeastSquares', 'Lower', [1-min(ydata) -1 0 min(ydata)], 'Upper', [1 1 3 min(ydata)], 'Startpoint', startPoints);
f = fittype(gaussEqn, 'coefficients', {'a', 'b', 'c', 'd'}, 'independent', 'x', 'options', s);
[F, G] = fit(xdata',ydata',f);

end
