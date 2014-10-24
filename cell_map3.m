function [cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map3(datarun0, datarun1, datarun2, cell_type)
%[cell_id, cell_idx, cell_id_matched, cell_idx_matched] = cell_map3(datarun0, datarun1, datarun2, cell_type)

[cell_list_map20, ~] = map_ei(datarun2, datarun0);
[cell_list_map21, ~] = map_ei(datarun2, datarun1);

n = length(cell_type);
cell_id = cell(n, 1);
cell_idx = cell(n, 1);
cell_id_matched = [];
cell_idx_matched = [];

for i = 1:n
    [cell_idx2, ~, ~] = get_cell_indices(datarun2, cell_type{i});
    cell_id1 = cell_list_map21(cell_idx2);
    ept = 1-cellfun(@isempty, cell_id1);
    cell_idx2 = cell_idx2.*ept;
    cell_idx2(cell_idx2 == 0) = [];
    cell_id1 = cell2mat(cell_id1);
    cell_idx1 = get_cell_indices(datarun1, cell_id1);
    cell_id2 = datarun2.cell_ids(cell_idx2);
    
    cell_id0 = cell_list_map20(cell_idx2);
    ept = 1-cellfun(@isempty, cell_id0);
    cell_idx2 = cell_idx2.*ept;
    cell_idx1 = cell_idx1.*ept;
    cell_idx2(cell_idx2 == 0) = [];
    cell_idx1(cell_idx1 == 0) = [];
    cell_id0 = cell2mat(cell_id0);
    cell_idx0 = get_cell_indices(datarun0, cell_id0);
    cell_id2 = datarun2.cell_ids(cell_idx2);
    cell_id1 = datarun1.cell_ids(cell_idx1);
    
    id = [cell_id0' cell_id1' cell_id2'];
    idx = [cell_idx0' cell_idx1' cell_idx2'];
    cell_id{i} = id;
    cell_idx{i} = idx;
    cell_id_matched = [cell_id_matched; id];
    cell_idx_matched = [cell_idx_matched; idx];
end
end


