function [F, G] = fit_nr(xdata, ydata, varargin)

p = inputParser;
p.addParamValue('Upper', [0.5, 10000, max(xdata), min(ydata)]);
p.addParamValue('Lower', [min(ydata), 0, 0, 0]);
p.addParamValue('Startpoints', []);

p.parse(varargin{:});
upper = p.Results.Upper;
lower = p.Results.Lower;
if ~isempty(p.Results.Startpoints)
    startPoints = p.Results.Startpoints;
else
    ymax0 = max(ydata);
    a0 = 2;
    [~,i] = min(abs(ymax0/2-ydata));
    sigma0 = xdata(i);
    b0 = min(ydata);
    startPoints = [ymax0 a0 sigma0 b0];
end
mmEqn = 'ymax*x.^a./(x.^a + sigma^a) + b';

s = fitoptions('METHOD', 'NonlinearLeastSquares', 'Lower', lower, 'Upper', upper, 'Startpoint', startPoints);
f = fittype(mmEqn, 'coefficients', {'ymax', 'a', 'sigma', 'b'}, 'independent', 'x', 'options', s);
[F, G] = fit(xdata', ydata', f);
end