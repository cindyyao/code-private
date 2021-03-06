function index = grp_stim(datarun)

% index = grp_stim(datarun)
% group stimulus trials by stimulus combination, return the index
% index: M spatial period x N temporal period x R directions x S repeats
%
% xyao
% 2013-12-16

list = datarun.stimulus.trial_list;
repeat = datarun.stimulus.repetitions;
stim_comb = length(datarun.stimulus.combinations);

sp = datarun.stimulus.params.SPATIAL_PERIOD;
tp = datarun.stimulus.params.TEMPORAL_PERIOD;
dr = datarun.stimulus.params.DIRECTION;
cl = datarun.stimulus.params.RGB;



index_temp = zeros(stim_comb, repeat);

for i = 1:stim_comb
    index_temp(i, :) = find(list == i);
end

[sp_temp, tp_temp, dr_temp] = deal(zeros(1, stim_comb));
cl_temp = cell(1, stim_comb);

for i = 1:stim_comb
    sp_temp(i) = datarun.stimulus.trials(i).SPATIAL_PERIOD;
    tp_temp(i) = datarun.stimulus.trials(i).TEMPORAL_PERIOD;
    dr_temp(i) = datarun.stimulus.trials(i).DIRECTION;
    cl_temp(i) = {datarun.stimulus.trials(i).RGB};
end

i_stim = zeros(length(sp), length(tp), length(dr), length(cl));
index = zeros(length(sp), length(tp), length(dr), length(cl), repeat);

for s = 1:length(sp)
    for t = 1:length(tp)
        for d = 1:length(dr)
            for c = 1:length(cl)
                [~, i_s] = find(sp_temp == sp(s));
                [~, i_t] = find(tp_temp == tp(t));
                [~, i_d] = find(dr_temp == dr(d));
                [~, i_c] = find(cellfun(@isequal, cl_temp, repmat(cl(c), size(cl_temp))));
                i_st = intersect(i_s, i_t);
                i_std = intersect(i_st, i_d);
                i = intersect(i_std, i_c);
                i_stim(s, t, d, c) = i;
                index(s, t, d, c, :) = index_temp(i, :);
            end
        end
    end
end


    

