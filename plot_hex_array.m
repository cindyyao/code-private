function plot_hex_array(datarun)

corner_i = [4 126 195 264 386 455 4];
corner_position = datarun.ei.position(corner_i, :);
plot(corner_position(:, 1), corner_position(:, 2));
end


