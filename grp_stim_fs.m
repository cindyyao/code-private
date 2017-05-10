function index = grp_stim_fs(datarun)

list = datarun.stimulus.trial_list;
repeat = datarun.stimulus.repetitions;

for i = 1:max(list)
    index_temp(i, :) = find(list == i);
end
index(:,:,1) = index_temp * 2 - 1;
index(:,:,2) = index_temp * 2;
end
    
