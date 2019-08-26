function [F, G] = fit_cos_(xdata, ydata, varargin)

p = inputParser;
% p.addParamValue('Upper', [(max(ydata)-min(ydata))*1.2, 2*pi, 10, max(min(ydata), 0)]);
p.addParamValue('Upper', [max(max(ydata)*1.2, 0), 2*pi, 10, 1]);
p.addParamValue('Lower', [0, 0, 0.3, min((1-min(ydata))*0.9, 1)]);
p.addParamValue('Startpoints', []);

p.parse(varargin{:});
upper = p.Results.Upper;
lower = p.Results.Lower;
if ~isempty(p.Results.Startpoints)
    startPoints = p.Results.Startpoints;
else
    ymax0 = max(ydata)-min(ydata);
    [~, phii] = max(ydata); 
    phi0 = mod(-xdata(phii), 2*pi);
    alpha0 = 3;
    b0 = 1-min(ydata);
    startPoints = [ymax0 phi0 alpha0 b0];
end
mmEqn = 'ymax*((0.5 + 0.5 * cos(x + phi)).^alpha * b + 1 - b)';

s = fitoptions('METHOD', 'NonlinearLeastSquares', 'Lower', lower, 'Upper', upper, 'Startpoint', startPoints);
f = fittype(mmEqn, 'coefficients', {'ymax', 'phi', 'alpha', 'b'}, 'independent', 'x', 'options', s);
[F, G] = fit(xdata', ydata', f);
end