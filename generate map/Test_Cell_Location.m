function [id_new, location_new] = Test_Cell_Location(id_old, location_old)

% [id_new, location_new] = Test_Cell_Location(id_old, location_old)
% Delete cells that are not found in WN datarun

map_fail_idx = find(sum(location_old, 2) == 0);
id_old(map_fail_idx) = [];
id_new = id_old;
location_old(map_fail_idx, :) = [];
location_new = location_old;
end
