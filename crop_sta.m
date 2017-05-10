function new_sta = crop_sta(sta, varargin)
% Crop sta on its spatial dimension so that the com of sta is in the center
% of the cropped sta.
% There is option to upscale the sta before cropping.
%
% new_sta = crop_sta(sta, varargin)
% 
% INPUTS: 
%   sta             a spatial-temporal-chromatic STA
%
% OPTIONAL INPUTS (HANDLED WITH VARARGIN AND INPUTPARSER):
%
%   scale_factor        10        a factor by which to upscale the sta
%                                 before cropping
%   crop_factor         2         crop sta from mxn to (m/f) x (n/f)
%
% AUTHOR: xyao
% DATE: 2015-10-14

% ----- PARSE INPUTS -----

p = inputParser;
p.addRequired('sta', @isnumeric);

% optional information to pass (will use if passed)
p.addParamValue('scale_factor', 10, @isnumeric);
p.addParamValue('crop_factor', 2, @isnumeric);

p.parse(sta, varargin{:});

scale_factor = p.Results.scale_factor;
crop_factor = p.Results.crop_factor;

%% -----BODY OF FUNCTION -----

% scale up sta
size_new_sta = size(sta);
size_new_sta(1:2) = size_new_sta(1:2)*scale_factor;
sta_up = zeros(size_new_sta);
for i = 1:size(sta, 3)
    for j = 1:size(sta, 4)
        sta_up(:, :, i, j) = matrix_scaled_up(sta(:, :, i, j), scale_factor);
    end
end

% crop sta
com = rf_com(sta);
if isempty(com)
    warning('no significant stixels found')
    new_sta = [];
else
    % com(2) = size(sta, 2) - com(2);
    com_up = round(com*scale_factor);
    if com_up >= 1/(2*crop_factor)*size_new_sta(1:2) & ...
            com_up <= (1-1/(2*crop_factor))*size_new_sta(1:2)
        a = com_up - 1/(2*crop_factor)*size_new_sta(1:2) + 1;
        b = com_up + 1/(2*crop_factor)*size_new_sta(1:2);
        new_sta = sta_up(a(2):b(2), a(1):b(1), :, :);
    else
        new_sta = [];
        warning('cell not centered enough, try larger crop factor')
    end
end

        
