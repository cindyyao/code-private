% model_series = zeros(5, 2);
% model_error = zeros(5, 2);

figure

h = bar(model_series);
set(h,'BarWidth',1); % The bars will now touch each other

set(gca,'XTicklabel',cell_type)
ylabel('Surround Strengths')
legend('melatonin','non-melatonin');
title('contrast 0.24')
hold on;

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end