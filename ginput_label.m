function coordinates = ginput_label(color)
done = 0;
coordinates = [];

while ~done
    c = ginput(1);
    hold on
    if isempty(c)
        done = 1;
    else
        plot(c(1), c(2), '.', 'color', color, 'MarkerSize', 10)
        coordinates = [coordinates; c];
    end
end
    