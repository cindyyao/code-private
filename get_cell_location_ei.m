function elec_num = get_cell_location_ei(datarun, cell_ids)
% return the electrode number with the strongest electrical activity
% (votage change)

% xyao
% 2016-01-25

elec_num = zeros(length(cell_ids), 1);
cell_idx = get_cell_indices(datarun, cell_ids);
for cc = 1:length(cell_idx)
    ei = datarun.ei.eis{cell_idx(cc)};
    ei_max = max(ei, [], 2);
    ei_min = min(ei, [], 2);
    ei_amp = ei_max - ei_min;
    [~, i] = max(ei_amp);
    elec_num(cc) = i;
end

end