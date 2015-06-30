[d,h,wf] = my_abfload('/Users/xyao/Documents/data/2015_06_16_0005.abf');
sampRate = h.sampRate;
theta = 4.5;
trigger_raw = d(:, 3);
trigger_idx = find(trigger_raw > theta);
duplicate_idx = diff(trigger_idx) < 20;
duplicate_idx = logical([0; duplicate_idx]);
trigger_idx(duplicate_idx) = [];
trigger_time = trigger_idx/sampRate;
trigger_interval = diff(trigger_time);

a = sort(trigger_interval);
a(end-1)-a(2)

figure
X = [trigger_time; trigger_time];
Y = [zeros(1, length(trigger_time)); ones(1, length(trigger_time))];
line(X, Y);
