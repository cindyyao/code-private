ratio3_mean{1,1} = squeeze(ratio_onoff_mean(2, 1, 1:7))';
ratio3_mean{1,2} = squeeze(ratio_onoff_mean(2, 2, 1:7))';
ratio3_mean{2,1} = squeeze(ratio_onoff_mean(5, 1, 1:7))';
ratio3_mean{2,2} = squeeze(ratio_onoff_mean(5, 2, 1:7))';

ratio3_ste{1,1} = squeeze(ratio_onoff_ste(2, 1, 1:7))';
ratio3_ste{1,2} = squeeze(ratio_onoff_ste(2, 2, 1:7))';
ratio3_ste{2,1} = squeeze(ratio_onoff_ste(5, 1, 1:7))';
ratio3_ste{2,2} = squeeze(ratio_onoff_ste(5, 2, 1:7))';

for ct = 1:4
ratio2_dir_mean{1,1}{ct} = squeeze(ratio_dir_mean(3, ct, 1:7))';
ratio2_dir_mean{2,1}{ct} = squeeze(ratio_dir_mean(6, ct, 1:7))';

ratio2_dir_ste{1,1}{ct} = squeeze(ratio_dir_ste(3, ct, 1:7))';
ratio2_dir_ste{2,1}{ct} = squeeze(ratio_dir_ste(6, ct, 1:7))';
end

for ct = 1:2
ratio2_dir_mean{1,2}{ct} = squeeze(ratio_dir_mean_on(3, ct, 1:7))';
ratio2_dir_mean{2,2}{ct} = squeeze(ratio_dir_mean_on(6, ct, 1:7))';

ratio2_dir_ste{1,2}{ct} = squeeze(ratio_dir_ste_on(3, ct, 1:7))';
ratio2_dir_ste{2,2}{ct} = squeeze(ratio_dir_ste_on(6, ct, 1:7))';
end

ratio2_dir_mean = ratio2_dir_mean';
ratio2_dir_ste = ratio2_dir_ste';

for ct = 1:4
ratio3_dir_mean{1,1}{ct} = squeeze(ratio_dir_mean(2, ct, 1:7))';
ratio3_dir_mean{1,2}{ct} = squeeze(ratio_dir_mean(5, ct, 1:7))';

ratio3_dir_ste{1,1}{ct} = squeeze(ratio_dir_ste(2, ct, 1:7))';
ratio3_dir_ste{1,2}{ct} = squeeze(ratio_dir_ste(5, ct, 1:7))';
end

for ct = 1:3
ratio3_dir_mean{2,1}{ct} = squeeze(ratio_dir_mean_on(2, ct, 1:7))';
ratio3_dir_mean{2,2}{ct} = squeeze(ratio_dir_mean_on(5, ct, 1:7))';

ratio3_dir_ste{2,1}{ct} = squeeze(ratio_dir_ste_on(2, ct, 1:7))';
ratio3_dir_ste{2,2}{ct} = squeeze(ratio_dir_ste_on(5, ct, 1:7))';
end