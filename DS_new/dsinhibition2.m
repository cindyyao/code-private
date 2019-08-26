%% load data
load('DS180413.mat', 'response_s', 'response_others')
response_s_0413 = response_s;
response_others_0413 = response_others;
load('DS170629.mat', 'response_s', 'response_others')
response_s_0629 = response_s;
response_others_0629 = response_others;
response_s_0629(3) = [];
response_others_0629(3) = [];

clear response_s response_others
for drug = 1:3
    response_s{drug} = [response_s_0629{drug}; response_s_0413{drug}];
    response_others{drug} = [response_others_0629{drug}; response_others_0413{drug}];
end

%%
ctr_x = [5 10 20 40 80 150 300];
color = 'brgkc';
figure
set(gcf, 'Position', [1 1 600 500])
errorbar(ctr_x, mean(response_s{1}), nanstd(response_s{1})/sqrt(size(response_s{1}, 1)), 'color', color(1));
hold on
errorbar(ctr_x, mean(response_s_ko{1}), nanstd(response_s_ko{1})/sqrt(size(response_s_ko{1}, 1)), 'color', color(2));
errorbar(ctr_x, mean(response_s{2}), nanstd(response_s{2})/sqrt(size(response_s{2}, 1)), 'color', color(3));
errorbar(ctr_x, mean(response_s_ko{2}), nanstd(response_s_ko{2})/sqrt(size(response_s_ko{2}, 1)), 'color', color(4));
legend('WT ctr', 'KO ctr', 'WT SR', 'KO SR', 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike rate')
xlim([3 400])
title('superior')

figure
set(gcf, 'Position', [1 1 600 500])
errorbar(ctr_x, mean(response_others{1}), nanstd(response_others{1})/sqrt(size(response_others{1}, 1)), 'color', color(1));
hold on
errorbar(ctr_x, mean(response_others_ko{1}), nanstd(response_others_ko{1})/sqrt(size(response_others_ko{1}, 1)), 'color', color(2));
errorbar(ctr_x, mean(response_others{2}), nanstd(response_others{2})/sqrt(size(response_others{2}, 1)), 'color', color(3));
errorbar(ctr_x, mean(response_others_ko{2}), nanstd(response_others_ko{2})/sqrt(size(response_others_ko{2}, 1)), 'color', color(4));
legend('WT ctr', 'KO ctr', 'WT SR', 'KO SR', 'location', 'northwest')
set(gca, 'Xscale', 'log')
xlabel('% contrast')
ylabel('spike rate')
xlim([3 400])
title('others')

%%
y = [response_s_ko{2}(:, 1:5); response_others_ko{2}(:, 1:5)]';
g1 = repmat(ctr_x(1:5)', 1, size(y, 2));
g2 = [repmat({'superior'}, 5, size(response_s_ko{2}, 1)) repmat({'others'}, 5, size(response_others_ko{2}, 1))];
p = anovan(y(:), {g1(:), g2(:)})
p = anovan(y(:), {g1(:), g2(:)}, 'model','interaction','varnames',{'g1','g2'})


y = [response_s{2}(:, 1:5); response_others{2}(:, 1:5)]';
g1 = repmat(ctr_x(1:5)', 1, size(y, 2));
g2 = [repmat({'superior'}, 5, size(response_s{2}, 1)) repmat({'others'}, 5, size(response_others{2}, 1))];
p = anovan(y(:), {g1(:), g2(:)})
p = anovan(y(:), {g1(:), g2(:)}, 'model','interaction','varnames',{'g1','g2'})



a = response_s_ko{2}(:, 1:5);
a = a(:);
b = response_others_ko{2}(:, 1:5);
b = b(:);
[~, p] = ttest2(a, b)



figure
for i = 1:size(response_s_ko{2}, 1)
    h1 = plot(ctr_x(1:5), response_s_ko{2}(i, 1:5), 'b');
    hold on
end
for i = 1:size(response_others_ko{2}, 1)
    h2 = plot(ctr_x(1:5), response_others_ko{2}(i, 1:5), 'r');
    hold on
end
set(gca, 'xscale', 'log')
xlim([4 100])
legend([h1, h2], 'superior', 'others')
xlabel('contrast')
ylabel('spike rate')
title('FACx')

figure
for i = 1:size(response_s{2}, 1)
    h1 = plot(ctr_x(1:5), response_s{2}(i, 1:5), 'b');
    hold on
end
for i = 1:size(response_others{2}, 1)
    h2 = plot(ctr_x(1:5), response_others{2}(i, 1:5), 'r');
    hold on
end
set(gca, 'xscale', 'log')
xlim([4 100])
legend([h1, h2], 'superior', 'others')
xlabel('contrast')
ylabel('spike rate')
title('C57')

%%
ya = response_s_ko{2}(:, 1:5);
xa = log10(repmat(ctr_x(1:5), size(ya, 1), 1));
yb = response_others_ko{2}(:, 1:5);
xb = log10(repmat(ctr_x(1:5), size(yb, 1), 1));
yall = [ya; yb];
xall = [xa; xb];

[fall, gall] = fit_nr(xall(:)', yall(:)', 'Upper', [100, 10000, max(xall(:)), Inf])
[fa, ga] = fit_nr(xa(:)', ya(:)', 'Upper', [100, 10000, max(xa(:)), Inf])
[fb, gb] = fit_nr(xb(:)', yb(:)', 'Upper', [100, 10000, max(xb(:)), Inf])

F = ((gall.sse - (ga.sse + gb.sse))/(ga.sse + gb.sse))/((gall.dfe - (ga.dfe + gb.dfe))/(ga.dfe + gb.dfe));
p = fcdf(F, gall.dfe - (ga.dfe + gb.dfe), ga.dfe + gb.dfe)



ya = response_s{2}(:, 1:5);
xa = log10(repmat(ctr_x(1:5), size(ya, 1), 1));
yb = response_others{2}(:, 1:5);
xb = log10(repmat(ctr_x(1:5), size(yb, 1), 1));
yall = [ya; yb];
xall = [xa; xb];

[fall, gall] = fit_nr(xall(:)', yall(:)', 'Upper', [100, 10000, max(xall(:)), Inf])
[fa, ga] = fit_nr(xa(:)', ya(:)', 'Upper', [100, 10000, max(xa(:)), Inf])
[fb, gb] = fit_nr(xb(:)', yb(:)', 'Upper', [100, 10000, max(xb(:)), Inf])

F = ((gall.sse - (ga.sse + gb.sse))/(ga.sse + gb.sse))/((gall.dfe - (ga.dfe + gb.dfe))/(ga.dfe + gb.dfe));
p = fcdf(F, gall.dfe - (ga.dfe + gb.dfe), ga.dfe + gb.dfe)