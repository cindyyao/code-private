load('dsmodel.mat', 'ERROR_rod','ERROR_cone')
load('dsmodel170504.mat', 'Precision_ratio_rod', 'Precision_ratio_cone')
ymax = 5; alpha1 = 7; c50 = 1.8;
ctr = log10([5 10 20 40 80 150 300]);
input1 = ymax.* ctr.^alpha1./(ctr.^alpha1 + c50.^alpha1);

ymax = 25; alpha1 = 7; c50 = 1.8;
ctr = log10([5 10 20 40 80 150 300]);
input2 = ymax.* ctr.^alpha1./(ctr.^alpha1 + c50.^alpha1);

[input, i] = sort([input1 input2]);
ERROR = [ERROR_rod; ERROR_cone];
ERROR = ERROR(i, :);

Precision_ratio = [Precision_ratio_rod Precision_ratio_cone];
Precision_ratio = Precision_ratio(:, i);
Precision_ratio = Precision_ratio';

theta = unique(sort([linspace(10,50,4)/180*pi linspace(50,80,10)/180*pi linspace(80,170,8)/180*pi]));
figure
for i = 1:size(ERROR,1)
    plot(theta/pi*360, ERROR(i, :))
    hold on
end

[~, optimal] = min(ERROR');
optimalWidth = zeros(size(ERROR,1),1);
for i = 1:size(ERROR,1)
    optimalWidth(i) = ERROR(i,optimal(i));
end
plot(theta(optimal)/pi*180*2, optimalWidth, 'r*')
xlabel('tuning width / degree')
ylabel('prediction error / degree')
ylim([0 100])
xlim([0 400])


figure
for i = 1:size(ERROR,1)
    plot(theta/pi*360, Precision_ratio(i, :))
    hold on
end
for i = 1:size(ERROR,1)
    optimaPrecision(i) = Precision_ratio(i,optimal(i));
end
plot(theta(optimal)/pi*180*2, optimaPrecision, 'r*')
xlabel('tuning width / degree')
ylabel('prediction error / degree')
% ylim([0 100])
xlim([0 400])

ERROR_norm = ERROR/180;
for i = 1:size(ERROR, 1)
    opt = optimal(i);
    while opt < size(ERROR, 2)
        if ERROR_norm(i, opt + 1) - ERROR_norm(i, opt) > Precision_ratio(i, opt + 1) - Precision_ratio(i, opt)
            break
        else
            opt = opt + 1;
        end
    end
    optimal_new(i) = opt;
end

figure
for i = 1:size(ERROR,1)
    plot(theta/pi*360, Precision_ratio(i, :))
    hold on
end

for i = 1:size(ERROR,1)
    optimaPrecision(i) = Precision_ratio(i,optimal_new(i));
end
plot(theta(optimal_new)/pi*180*2, optimaPrecision, 'r*')
xlabel('tuning width / degree')
ylabel('prediction error / degree')
% ylim([0 100])
xlim([0 400])

%% optimal tuning width
clear error theta_e
Ymax = linspace(0.3, 6, 10);
Ymax = Ymax.^2;
theta = unique(sort([linspace(10,50,4)/180*pi linspace(50,80,10)/180*pi linspace(80,170,8)/180*pi]));
Alpha2 = log(0.5)./log(0.5+0.5*cos(theta));


for i = 1:length(ctr)
    for j = 1:length(Alpha2)
        % NDF 4
        % CRF 
        c50 = [1.8 1.8 1.8 1.8];
        alpha1 = [7 7 7 7];
        % tuning curve
        phi = [0 1/2 1 3/2]*pi; % preferred direction
        alpha2 = [1 1 1 1]*Alpha2(j); % tuning width 
        r_base = [1 1 1 1]; % background firing
        ymax = [1 1 1 1]*25; % ymax = 5 for rod light level
        % run_MLdecoder
        [error(i,j), theta_e{i,j}] = run_poiss_MLdecoder(c50,alpha1,phi,alpha2,r_base,ymax,ctr(i));
        [i,j]
    end
end 

%% RUN THIS
contrast = [10 20 40 80 150 300]; % contrast
theta = unique(sort([linspace(10,50,4)/180*pi linspace(50,80,10)/180*pi linspace(80,170,8)/180*pi]));
Alpha2 = log(0.5)./log(0.5+0.5*cos(theta));


c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
r_base = [1 1 1 1]; % background firing
repeat_n = 200;

[Precision_rod, Precision_cone, ERROR_rod, ERROR_cone] = deal(cell(1, 10));
for repeat = 2:10
    clear Precision_ratio
    if repeat <= 5
        ymax = [1 1 1 1]*5; % ymax = 5 for rod light level
    else
        ymax = [1 1 1 1]*30; % ymax = 30 for cone light level
    end
    for width = 1:length(Alpha2)
        alpha2 = [1 1 1 1]*Alpha2(width); % tuning width 
        for ctr = 1:length(contrast)
            Precision_ratio(width, ctr) = run_poiss_Detection(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n);
            [error(width,ctr), theta_e{width,ctr}] = run_poiss_MLdecoder(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n);
    %         [width ctr]
            disp(['repeat ' num2str(repeat) ' detection task: tuning width ' num2str(width) ' contrast ' num2str(ctr) ' finished.'])
        end
    end
    if repeat <= 5
        Precision_rod{repeat} = Precision_ratio;
        ERROR_rod{repeat} = error;
        theta_e_rod{repeat} = theta_e;
    else
        Precision_cone{repeat-5} = Precision_ratio;
        ERROR_cone{repeat-5} = error;
        theta_e_cone{repeat} = theta_e;
    end
    save('dsmodel_direction.mat', 'Precision_rod', 'Precision_cone', 'ERROR_rod', 'ERROR_cone', 'theta_e_rod', 'theta_e_cone')
end

Precision_ratio = Precision_ratio_rod;
theta = unique(sort([linspace(10,50,4)/180*pi linspace(50,80,10)/180*pi linspace(80,170,8)/180*pi]));
figure
for i = 1:size(Precision_ratio, 2)
    plot(theta*360/pi, Precision_ratio(:, i))
    hold on
end

%%
load('dsmodel-1.mat')
load('dsmodel.mat', 'ERROR_cone')
load('dsmodel170504.mat', 'Precision_ratio_cone')
ERROR_cone_mean = ERROR_cone';
Precision_cone_mean = Precision_ratio_cone;

ymax = 5; alpha1 = 7; c50 = 1.8;
ctr = log10([10 20 40 80 150 300]);
input1 = ymax.* ctr.^alpha1./(ctr.^alpha1 + c50.^alpha1);

ymax = 25; alpha1 = 7; c50 = 1.8;
ctr = log10([5 10 20 40 80 150 300]);
input2 = ymax.* ctr.^alpha1./(ctr.^alpha1 + c50.^alpha1);

[input, i] = sort([input1 input2]);
ERROR = [ERROR_rod_mean ERROR_cone'];
ERROR = ERROR(:, i);
ERROR = ERROR(:, 4:end);

Precision_ratio = [Precision_rod_mean Precision_ratio_cone];
Precision_ratio = Precision_ratio(:, i);
Precision_ratio = Precision_ratio(:, 4:end);

theta = unique(sort([linspace(10,50,4)/180*pi linspace(50,80,10)/180*pi linspace(80,170,8)/180*pi]));
figure
for i = 1:size(ERROR,2)
    plot(theta/pi*360, ERROR(:, i), 'k')
    hold on
end

[~, optimal] = min(ERROR);
optimalWidth = zeros(size(ERROR,2),1);
for i = 1:size(ERROR,2)
    optimalWidth(i) = ERROR(optimal(i), i);
end
plot(theta(optimal)/pi*180*2, optimalWidth, 'r*')
xlabel('tuning width / degree')
ylabel('prediction error / degree')
ylim([0 100])
xlim([0 400])



figure
for i = 1:size(ERROR,2)
    plot(theta/pi*360, Precision_ratio(:, i), 'k')
    hold on
end
for i = 1:size(ERROR,2)
    optimaPrecision(i) = Precision_ratio(optimal(i), i);
end
plot(theta(optimal)/pi*180*2, optimaPrecision, 'r*')
xlabel('tuning width / degree')
ylabel('prediction error / degree')
% ylim([0 100])
xlim([0 400])

ERROR_norm = ERROR/180;
for i = 1:size(ERROR, 2)
    opt = optimal(i);
    while opt < size(ERROR, 1)
        if ERROR_norm(opt + 1, i) - ERROR_norm(opt, i) > Precision_ratio(opt + 1, i) - Precision_ratio(opt, i)
            break
        else
            opt = opt + 1;
        end
    end
    optimal_new(i) = opt;
end

figure
for i = 1:size(ERROR,2)
    plot(theta/pi*360, ERROR(:, i), 'k')
    hold on
end

for i = 1:size(ERROR,2)
    optimaERROR(i) = ERROR(optimal_new(i), i);
end
plot(theta(optimal_new)/pi*180*2, optimaERROR, 'r*')
xlabel('tuning width / degree')
ylabel('prediction error / degree')
% ylim([0 100])
xlim([0 400])

%%
for width = 1:20
    for ctr = 1:6
        error_max_rod(width, ctr) = max(mean(theta_e_rod{1}{width, ctr}));
    end
end

for rep = 1:5
for width = 1:20
    for ctr = 1:6
        error_max_cone(width, ctr, rep) = max(mean(theta_e_cone{rep}{width, ctr}));
    end
end
end
[~, i] = min(mean(error_max_cone, 3))