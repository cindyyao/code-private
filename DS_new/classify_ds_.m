function [ds_id, nonds_id] = classify_ds_(datarun, ds_struct, params_idx, ctr_id)

figure
plot(ds_struct.mag{params_idx(1), 1, ctr_id(1)}, ds_struct.mag{params_idx(2), 1, ctr_id(2)}, 'o')
hold on
[x, y] = ginput;
plot(x, y);
IN = inpolygon(ds_struct.mag{params_idx(1), 1, ctr_id(1)}, ds_struct.mag{params_idx(2), 1, ctr_id(2)}, x, y);
[~, I] = find(IN == 1);
id_init = datarun.cell_ids(I);

[~, ~, ib] = intersect(id_init, datarun.cell_ids);
vc = ones(length(datarun.cell_ids),1);
vc(ib) = 2; %initializing ds cells to cluster 2, everything else cluster 1

X(:,1) = log(ds_struct.mag{params_idx(1),1, ctr_id(1)})';
X(:,2) = log(ds_struct.mag{params_idx(2),1, ctr_id(2)})';
[idx, ~, ~] = clustering_analysis_plots(X, 0,1, 2, 0, 1, 0, 0, 0,0, vc);

ds_id = datarun.cell_ids(idx==2);
nonds_id = datarun.cell_ids(idx==1);
% ds_id = id_init;
% nonds_id = setdiff(datarun.cell_ids, ds_id);
