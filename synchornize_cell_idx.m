cell_idx = cell(7, 1);
for i = 1:7
cell_idx{i}(:, 1) = get_cell_indices(datarun{1}, cell_id{i}(:, 1));
cell_idx{i}(:, 2) = get_cell_indices(datarun{2}, cell_id{i}(:, 2));
end

cell_id_matched = [];
cell_idx_matched = [];
for i = 1:7
cell_id_matched = [cell_id_matched; cell_id{i}];
cell_idx_matched = [cell_idx_matched; cell_idx{i}];
end

save('RFS131107.mat', 'cell_id', 'cell_idx', 'cell_id_matched', 'cell_idx_matched', '-append')