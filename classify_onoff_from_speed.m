function [id_sub, idx_sub] = classify_onoff_from_speed(datarun, ds_id, pc1, pc2, manual)

[NumSpikesCell, ~, StimComb] = get_spikescellstim(datarun,ds_id,0,1);
DG = sort_direction(dscellanalysis(NumSpikesCell, StimComb,datarun));
MAG_all_norm_dg = normalize_MAG(DG);
[id_sub, idx_sub] = classify_onoff(MAG_all_norm_dg', ds_id, pc1, pc2, manual); % select the right cluster!!!

end