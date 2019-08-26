function nnnds = get_nnnds_xy(com, radius)
%
% usage:  nnnds = get_nnnds(datarun, cell_spec, varargin)
%
%
% outputs:     nnnds - nnnd between  cells
%
% optional params, their default values, and what they specify:
%
% fits_to_use           matlab          'vision' or 'matlab' fits can be specified
%
% 2018-04 XY

num_rgcs = size(com, 1);
all_dists = zeros(num_rgcs);

for rgc1 = 1:num_rgcs
    for rgc2 = 1:num_rgcs
        d = nnnd_gdf(com(rgc1,:), com(rgc2,:), 0, 0, radius, radius);
        all_dists(rgc1, rgc2) = d;
    end
end

nnnds = zeros(num_rgcs, 1);
for rgc = 1:num_rgcs
    temp_dists = sort(unique(all_dists(rgc,:)), 'ascend');
    nnnds(rgc) = min(temp_dists(2:end));
end


      
     




    
    
        
        
        
        
        
        
        