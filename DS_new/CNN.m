load('~/Desktop/pred.mat')
act_eval1 = act_eval;
predict_eval1 = predict_eval;
load('~/Desktop/pred1517457856.mat')
act_eval2 = act_eval;
predict_eval2 = predict_eval;


% 74 94 105 124
figure
for i = 100:447
    plot([act_eval(101:200, i) predict_eval(101:200, i)])
    legend('act', 'pred')
    title(['cell #: ' num2str(i)])
    pause
end

for i = 1:447
%     r(i) = corr(act_eval(1:end-1, i), predict_eval(2:end, i));
    r(i) = corr(act_eval(:, i), predict_eval(:, i));
end
r(isnan(r)) = 0;
mean(r)

%% select cells with high firing rate
load('~/Dropbox/CNN_Retina_Project/RetinaCNNCode/BWDb1.mat')
[~, i] = sort(sum(Psth, 2), 'descend');
Psth = Psth(i(1:200), :);
save('/Volumes/lab/Experiments/Array/Shared/xyao/CNN/CNN_Retina_Project/RetinaCNNCode/BWDb3.mat', 'Psth', 'mov', '-v7.3')


%% smooth psth
load('~/Dropbox/CNN_Retina_Project/RetinaCNNCode/BWDb1.mat')
x = [-3:3];
filter = normpdf(x,0,1);
filter = filter/sum(filter);

Psth_temp = conv2(Psth, filter, 'same');
Psth_temp = Psth_temp/max(Psth_temp(:));

i=3;
figure
plot(Psth(i, 1:100));
hold on
plot(Psth_temp(i, 1:100))

Psth = Psth_temp;
save('/Volumes/lab/Experiments/Array/Shared/xyao/CNN/CNN_Retina_Project/RetinaCNNCode/BWDb4.mat', 'Psth', 'mov', '-v7.3')


