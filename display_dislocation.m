n = size(cell_idx_matched, 1);
com_x = zeros(n, 2);
com_y = zeros(n, 2);
for cc = 1:50
    for i = 1:2
        com_x(cc, i) = datarun{i}.stas.rf_coms{cell_idx_matched(cc, i)}(1);
        com_y(cc, i) = datarun{i}.stas.rf_coms{cell_idx_matched(cc, i)}(2);
    end
end
com_x(:, 2) = com_x(:, 2)/3*2;
com_y(:, 2) = com_y(:, 2)/3*2;

diff_com_x = com_x(:, 1) - com_x(:, 2);
diff_com_y = com_y(:, 1) - com_y(:, 2);

diff_x = mean(diff_com_x)
diff_y = mean(diff_com_y)
