function [cell_location, varargout] = get_cell_location(datarun, cell_id, varargin)

% SET UP OPTIONAL ARGUMENTS
p = inputParser;

% specify list of optional parameters
p.addParamValue('ei_mode', false, @islogical);
p.addParamValue('foa', -1);
p.addParamValue('com', false, @islogical);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% BODY OF THE FUNCTION

if params.ei_mode
    ei_idx = [];
    rotation = input('EI to STA rotation (in deg):\n');
end

cell_location = zeros(length(cell_id), 2);
for cc = 1:length(cell_id)
    if ismember(cell_id(cc), datarun.cell_ids)
        figure
        if params.ei_mode
            set(gcf, 'Position', [1,1,1000, 500])
            subplot(1,2,2)
            plot_ei(datarun, cell_id(cc), 'rotation', rotation);
                        subplot(1,2,1)
            plot_rf(datarun, cell_id(cc), 'foa', params.foa, 'com', params.com);

            not_done = 1;
            while(not_done)
                mode = input('STA or EI? STA:0  EI:1\n');
                if ismember(mode, [0 1])
                    not_done = 0;
                else
                    warning('please choose from 0 and 1')
                end
            end
            ei_idx = [ei_idx; mode];
            cell_location(cc, :) = ginput();
            % coordinate transformation from EI to STA
            % ...
            
        else
            plot_rf(datarun, cell_id(cc), 'foa', params.foa, 'com', params.com);
            cell_location(cc, :) = ginput;
        end 
        close
    else
        warning('cell ''%.0f'' does not exist in datarun',cell_id(cc))
    end
end

if params.ei_mode
    varargout{1} = ei_idx;
end

end
