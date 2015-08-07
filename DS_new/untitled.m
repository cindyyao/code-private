ndf = 1;
cc = 1293; %cell_id = 4923
bw = 1;
dt = 1;
cl = 1;

for i = 1:size(temp,1)
    for j = 1:size(temp,2)
        trigger_1 = datamb{ndf}.stimulus.triggers(temp(i,j));
        trigger_2 = datamb{ndf}.stimulus.triggers(temp(i,j)+1);
        spikes{i, j} = datamb{ndf}.spikes{cc}(datamb{ndf}.spikes{cc} > trigger_1 & datamb{ndf}.spikes{cc} <= trigger_2);
    end
end


%%

inbw = find(cell2mat(StimComb(:, 1)) == 120);
indt = find(cell2mat(StimComb(:, 2)) == 2);
a = cell2mat(StimComb(:, 4));
incl = find(a(:,1) == 0.5);

mintersect(inbw,indt,incl)

[NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb(datamb{ndf},4923,duration,1);
