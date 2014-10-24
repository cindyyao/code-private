surround_mean = zeros(length(surround_ndf0_12.surround), 1);
surround_stev = zeros(length(surround_ndf0_12.surround), 1);

for i = 1:length(surround_ndf0_12.surround)
    surround_mean(i) = mean(surround_ndf0_12.surround{i});
    surround_stev(i) = std(surround_ndf0_12.surround{i})/length(surround_ndf0_12.surround{i});
end
surround_ndf0_12.surround_mean = surround_mean;
surround_ndf0_12.surround_stev = surround_stev;

%%

for j = 1:3
surround_mean = zeros(length(surround_strengths_cone{j}.surround), 1);
surround_stev = zeros(length(surround_strengths_cone{j}.surround), 1);

for i = 1:length(surround_strengths_cone{j}.surround)
    surround_mean(i) = mean(surround_strengths_cone{j}.surround{i});
    surround_stev(i) = std(surround_strengths_cone{j}.surround{i})/length(surround_strengths_cone{j}.surround{i});
end
surround_strengths_cone{j}.surround_mean = surround_mean;
surround_strengths_cone{j}.surround_stev = surround_stev;
end
