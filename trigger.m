frame_time = zeros(1, (length(datarun.triggers) - 1)*50);
for j = 1:length(datarun.triggers)-1
    frame_time(50*(j-1)+1:50*j+1) = linspace(datarun.triggers(i), datarun.triggers(i+1), 51);
end

