function [datarun_splitted] = split_datarun(datarun, time_points)

time_points = [time_points datarun.duration];
n = length(time_points);
datarun_splitted = cell(1, n);
for i = 1:n
    triggers = datarun.triggers(datarun.triggers < time_points(i));
    datarun.triggers(datarun.triggers < time_points(i)) = [];
    for cc = 1:length(datarun.cell_ids)
        spikes{cc, 1} = datarun.spikes{cc, 1}(datarun.spikes{cc, 1} < time_points(i));
        datarun.spikes{cc, 1}(datarun.spikes{cc, 1} < time_points(i)) = [];
    end
    datarun_splitted{i} = datarun;
    if i >= 2
        datarun_splitted{i}.duration = time_points(i) - time_points(i-1);
        datarun_splitted{i}.triggers = triggers - time_points(i-1);
        for cc = 1:length(datarun.cell_ids)
            spikes{cc, 1} = spikes{cc, 1} - time_points(i-1);
        end
    else
        datarun_splitted{i}.duration = time_points(i);
        datarun_splitted{i}.triggers = triggers;
    end
    datarun_splitted{i}.spikes = spikes;
end
