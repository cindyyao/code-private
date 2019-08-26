function [] = get_anova1(data)
    % data: A cell array of 1-D vectors
    % Author: xyao 2019-02-27
    group = [];
    data_new = [];
    for i = 1:length(data)
        group = [group ones(1, length(data{i}))*i];
        data_new = [data_new reshape(data{i}, 1, length(data{i}))];
    end
    anova1(data_new, group)
end