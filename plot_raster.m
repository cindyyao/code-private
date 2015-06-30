function plot_raster(raster, start_time, end_time, varargin)
% plot_raster(raster, start_time, end_time)
% raster: nx1 cells, one cell for a trial
%
% xyao
% 2013-12-16

p = inputParser;
p.addParamValue('color', 'b', @ischar);
p.addParamValue('first_trial', 1, @isnumeric);

p.parse(varargin{:});
params = p.Results;

for j = 1:length(raster)
    SpikeTime = raster{j};
    SpikeTime = SpikeTime';
    X = [SpikeTime; SpikeTime];
    Y = [ones(1, length(SpikeTime))*(j+params.first_trial-1.9); ones(1, length(SpikeTime))*(j+params.first_trial-1)];
    line(X, Y, 'color', params.color);
    axis([start_time, end_time,0,length(raster)+params.first_trial-1])
    hold on
end

