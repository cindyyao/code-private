%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% load data
datarun{1} = load_data('/Analysis/xyao/2014-10-14-0/data000/data000', opt);
datarun{1}.names.stimulus_path = '/Analysis/xyao/2014-10-14-0/stimuli/s00.mat';
datarun{1} = load_stim_matlab(datarun{1});

datarun{2} = load_data('/Analysis/xyao/2014-10-14-0/data003/data003', opt);
datarun{2}.names.stimulus_path = '/Analysis/xyao/2014-10-14-0/stimuli/s03.mat';
datarun{2} = load_stim_matlab(datarun{2});


[NumSpikesCell, StimComb] = get_spikescellstim(datarun{2},datarun{2}.cell_ids,0);
ds_struct = dscellanalysis(NumSpikesCell, StimComb);

% pull out DS cells

figure
plot(ds_struct.mag{3, 1}, ds_struct.mag{4, 1}, 'o')
title('data002 after mapping')
xlabel('TP 24')
ylabel('TP 60')
hold on
[x, y] = ginput;
plot(x, y);

IN = inpolygon(ds_struct.mag{3, 1}, ds_struct.mag{4, 1}, x, y);
[~, I] = find(IN == 1);
id = datarun{2}.cell_ids(I);
idx = get_cell_indices(datarun{2}, id);
I = zeros(length(id), 4);
for i = 1:4
    I(:, i) = ismember(id, datarun{i}.cell_ids);
end
I = sum(I');
id4 = id(I == 4); % DS cells that can be found in all light levels
idx4 = arrayfun(@(x) find(id == x), id4);

