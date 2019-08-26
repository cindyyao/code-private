repeat = 10;
rmax = 25;

%% discrimination and detection task (precision likelihood ratio)
c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
a = 1.8; b = 0.8;
Alpha2 = [a a a a ; a a b a]; % tuning width 
r_base = [1 1 1 1]*0.02; % background firing
ymax = [1 1 1 1]*rmax;
beta = [1 1 1 1]*0.95;
contrast = [40];
repeat_n = 10000;
% [Precision_dir, theta_e_dir] = deal(cell(1, repeat));
for rp = 3:repeat
    for width = 1:size(Alpha2, 1)
        alpha2 = Alpha2(width, :);
        for ctr = 1:length(contrast)
            [~, Precision(width, ctr, :)] = run_poiss_Detection_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,beta);
%             [error(width, ctr), theta_e{width, ctr}] = run_poiss_MLdecoder_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,beta);
            
            disp(['repeat ' num2str(rp) ' detection task: tuning width ' num2str(width) ' contrast ' num2str(ctr) ' finished.'])
        end
    end
    Precision_dir{rp} = Precision;
%     theta_e_dir{rp} = theta_e;
%     save('dsmodel_temp5.mat', 'Precision_ratio_cone', 'error_cone')
    % dsmodel_temp1.mat: contrast = [40 50 60], r_base = 0.1, a = 2; b = 0.8;
    % dsmodel_temp2.mat: contrast = [40 50 60], r_base = 0.02, a = 2; b = 0.8; 
    % dsmodel_temp3.mat: contrast = [40 50 60], r_base = 1, a = 2; b = 0.8;
    % dsmodel_temp4.mat: contrast = [40 50], r_base = 0.02, a = 1.8; b = 0.8;
    % dsmodel_temp5.mat: contrast = [40 50], r_base = 0.02, a = 1.9; b = 0.8;
end

%% beta only

contrast = [40]; % contrast
alpha2 = [1.8 1.8 1.8 1.8]; % tuning width 


c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
r_base = [1 1 1 1]*0.02; % background firing
repeat_n = 10000;
a = 0.95; b = 0.6;
Beta = [a a a a; a a b a];
ymax = [1 1 1 1]*rmax;
% [Precision_beta, theta_e_beta] = deal(cell(1, repeat));
for rp = 3:repeat
    clear Precision_ratio
    for k = 1:size(Beta, 1)
        for ctr = 1:length(contrast)
            [~, Precision(k, ctr, :)] = run_poiss_Detection_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,Beta(k,:));
%             [error(k,ctr), theta_e{k,ctr}] = run_poiss_MLdecoder_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,Beta(k,:));
    %         [width ctr]
            disp(['repeat ' num2str(rp) ' detection task: beta ' num2str(k) ' contrast ' num2str(ctr) ' finished.'])
        end
    end
    Precision_beta{rp} = Precision;
%     ERROR_rod{repeat} = error;
%     theta_e_beta{rp} = theta_e;

%     save('dsmodel_direction.mat', 'Precision_rod', 'Precision_cone', 'ERROR_rod', 'ERROR_cone', 'theta_e_rod', 'theta_e_cone')
end


%% firing rate

contrast = [40]; % contrast
alpha2 = [1.8 1.8 1.8 1.8]; % tuning width 


c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
r_base = [1 1 1 1]*0.02; % background firing
repeat_n = 10000;
a = 0.95; b = 0.6;
Beta = [a a a a; a a a a];
Ymax = [1 1 1 1; 1 1 2 1]*rmax;
% [Precision_rmax, theta_e_rmax] = deal(cell(1, repeat));
for rp = 3:repeat
    clear Precision_ratio
    for k = 1:size(Beta, 1)
        ymax = Ymax(k, :); % ymax = 5 for rod light level
        for ctr = 1:length(contrast)
            [~, Precision(k, ctr, :)] = run_poiss_Detection_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,Beta(k,:));
%             [error(k,ctr), theta_e{k,ctr}] = run_poiss_MLdecoder_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,Beta(k,:));
    %         [width ctr]
            disp(['repeat ' num2str(rp) ' detection task: beta ' num2str(k) ' contrast ' num2str(ctr) ' finished.'])
        end
    end
    Precision_rmax{rp} = Precision;
%     ERROR_rod{repeat} = error;
%     theta_e_rmax{rp} = theta_e;

%     save('dsmodel_direction.mat', 'Precision_rod', 'Precision_cone', 'ERROR_rod', 'ERROR_cone', 'theta_e_rod', 'theta_e_cone')
end
%% all
contrast = [40]; % contrast
Alpha2 = [1.8 1.8 1.8 1.8; 1.8 1.8 0.8 1.8]; % tuning width 


c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
r_base = [1 1 1 1]*0.02; % background firing
repeat_n = 10000;
a = 0.95; b = 0.6;
Beta = [a a a a; a a b a];
Ymax = [1 1 1 1; 1 1 2 1]*rmax;
% [Precision_all, theta_e_all] = deal(cell(1, repeat));
for rp = 11:repeat
    clear Precision_ratio
    for k = 1:size(Beta, 1)
        ymax = Ymax(k, :); % ymax = 5 for rod light level
        alpha2 = Alpha2(k, :);
        for ctr = 1:length(contrast)
            [~, Precision(k, ctr, :)] = run_poiss_Detection_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,Beta(k,:));
%             [error(k,ctr), theta_e{k,ctr}] = run_poiss_MLdecoder_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,Beta(k,:));
    %         [width ctr]
            disp(['repeat ' num2str(rp) ' detection task: beta ' num2str(k) ' contrast ' num2str(ctr) ' finished.'])
        end
    end
    Precision_all{rp} = Precision;
%     ERROR_rod{repeat} = error;
%     theta_e_all{rp} = theta_e;

%     save('dsmodel_direction.mat', 'Precision_rod', 'Precision_cone', 'ERROR_rod', 'ERROR_cone', 'theta_e_rod', 'theta_e_cone')
end

% save('width_beta.mat', 'Precision_dir', 'Precision_beta', 'Precision_both', 'theta_e_dir', 'theta_e_both', 'theta_e_beta')


% %% firing rate & beta
% 
% contrast = [40]; % contrast
% alpha2 = [1.8 1.8 1.8 1.8]; % tuning width 
% 
% 
% c50 = [1.8 1.8 1.8 1.8];
% alpha1 = [7 7 7 7];
% % tuning curve
% phi = [0 1/2 1 3/2]*pi; % preferred direction
% r_base = [1 1 1 1]*0.02; % background firing
% repeat_n = 10000;
% a = 0.95; b = 0.6;
% Beta = [a a a a; a a b a];
% Ymax = [1 1 1 1; 1 1 2 1]*rmax;
% [Precision_fb, theta_e_fb] = deal(cell(1, repeat));
% for rp = 1:repeat
%     clear Precision_ratio
%     for k = 1:size(Beta, 1)
%         ymax = Ymax(k, :); % ymax = 5 for rod light level
%         for ctr = 1:length(contrast)
%             [~, Precision(k, ctr, :)] = run_poiss_Detection_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,Beta(k,:));
% %             [error(k,ctr), theta_e{k,ctr}] = run_poiss_MLdecoder_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,Beta(k,:));
%     %         [width ctr]
%             disp(['repeat ' num2str(rp) ' detection task: beta ' num2str(k) ' contrast ' num2str(ctr) ' finished.'])
%         end
%     end
%     Precision_fb{rp} = Precision;
% %     ERROR_rod{repeat} = error;
% %     theta_e_fb{rp} = theta_e;
% 
% %     save('dsmodel_direction.mat', 'Precision_rod', 'Precision_cone', 'ERROR_rod', 'ERROR_cone', 'theta_e_rod', 'theta_e_cone')
% end

%% figure
idx = 1:repeat;
Precision_dir_all = mean(cat(4,Precision_dir{idx}), 4);
theta_temp = cellfun(@cell2mat, theta_e_dir, 'UniformOutput', 0);
theta_e_dir_all = mean(cat(3, theta_temp{idx}), 3);
THETA = linspace(0, 2*pi, size(Precision_dir_all, 3));
figure
subplot(2,2,2)
plot(THETA, squeeze(Precision_dir_all(2,1,:) - Precision_dir_all(1,1,:)), 'k')
title('detection rate difference')
xlabel('direction')
ylabel('detection rate difference')
xlim([0 2*pi])

subplot(2,2,1)
plot(THETA, squeeze(Precision_dir_all(1,1,:)))
hold on
plot(THETA, squeeze(Precision_dir_all(2,1,:)))
legend('all narrow', 'one broad')
title('detection rate')
xlabel('direction')
ylabel('detection rate')
xlim([0 2*pi])

subplot(2,2,4)
plot(THETA, squeeze(theta_e_dir_all(2,:) - theta_e_dir_all(1,:)), 'k')
title('error difference')
xlabel('direction')
ylabel('error difference')
xlim([0 2*pi])

subplot(2,2,3)
plot(THETA, squeeze(theta_e_dir_all(1,:)))
hold on
plot(THETA, squeeze(theta_e_dir_all(2,:)))
legend('all narrow', 'one broad')
title('errror')
xlabel('direction')
ylabel('error')
xlim([0 2*pi])


% beta
Precision_beta_all = mean(cat(4,Precision_beta{idx}), 4);
theta_temp = cellfun(@cell2mat, theta_e_beta, 'UniformOutput', 0);
theta_e_beta_all = mean(cat(3, theta_temp{idx}), 3);
figure
subplot(2,2,2)
plot(THETA, squeeze(Precision_beta_all(2,1,:) - Precision_beta_all(1,1,:)), 'k')
title('detection rate difference')
xlabel('direction')
ylabel('detection rate difference')
xlim([0 2*pi])

subplot(2,2,1)
plot(THETA, squeeze(Precision_beta_all(1,1,:)))
hold on
plot(THETA, squeeze(Precision_beta_all(2,1,:)))
legend('all narrow', 'one broad')
title('detection rate')
xlabel('direction')
ylabel('detection rate')
xlim([0 2*pi])

subplot(2,2,3)
plot(THETA, squeeze(theta_e_beta_all(1,:)))
hold on
plot(THETA, squeeze(theta_e_beta_all(2,:)))
legend('all narrow', 'one broad')
title('error')
xlabel('direction')
ylabel('error')
xlim([0 2*pi])

subplot(2,2,4)
plot(THETA, squeeze(theta_e_beta_all(2,:) - theta_e_beta_all(1,:)), 'k')
title('error difference')
xlabel('direction')
ylabel('error difference')
xlim([0 2*pi])

% beta only
Precision_rmax_all = mean(cat(4,Precision_rmax{idx}), 4);
theta_temp = cellfun(@cell2mat, theta_e_rmax, 'UniformOutput', 0);
theta_e_rmax_all = mean(cat(3, theta_temp{idx}), 3);
figure
subplot(2,2,2)
plot(THETA, squeeze(Precision_rmax_all(2,1,:) - Precision_rmax_all(1,1,:)), 'k')
title('detection rate difference')
xlabel('direction')
ylabel('detection rate difference')
xlim([0 2*pi])

subplot(2,2,1)
plot(THETA, squeeze(Precision_rmax_all(1,1,:)))
hold on
plot(THETA, squeeze(Precision_rmax_all(2,1,:)))
legend('all narrow', 'one broad')
title('detection rate')
xlabel('direction')
ylabel('detection rate')
xlim([0 2*pi])

subplot(2,2,3)
plot(THETA, squeeze(theta_e_rmax_all(1,:)))
hold on
plot(THETA, squeeze(theta_e_rmax_all(2,:)))
legend('all narrow', 'one broad')
title('error')
xlabel('direction')
ylabel('error')
xlim([0 2*pi])

subplot(2,2,4)
plot(THETA, squeeze(theta_e_rmax_all(2,:) - theta_e_rmax_all(1,:)), 'k')
title('error difference')
xlabel('direction')
ylabel('error difference')
xlim([0 2*pi])

% all 

Precision_all_all = mean(cat(4,Precision_all{idx}), 4);
theta_temp = cellfun(@cell2mat, theta_e_all, 'UniformOutput', 0);
theta_e_all_all = mean(cat(3, theta_temp{idx}), 3);
figure
subplot(2,2,2)
plot(THETA, squeeze(Precision_all_all(2,1,:) - Precision_all_all(1,1,:)), 'k')
title('detection rate difference')
xlabel('direction')
ylabel('detection rate difference')
xlim([0 2*pi])

subplot(2,2,1)
plot(THETA, squeeze(Precision_all_all(1,1,:)))
hold on
plot(THETA, squeeze(Precision_all_all(2,1,:)))
legend('all narrow', 'one broad')
title('detection rate')
xlabel('direction')
ylabel('detection rate')
xlim([0 2*pi])

subplot(2,2,3)
plot(THETA, squeeze(theta_e_all_all(1,:)))
hold on
plot(THETA, squeeze(theta_e_all_all(2,:)))
legend('all narrow', 'one broad')
title('error')
xlabel('direction')
ylabel('error')
xlim([0 2*pi])

subplot(2,2,4)
plot(THETA, squeeze(theta_e_all_all(2,:) - theta_e_all_all(1,:)), 'k')
title('error difference')
title('error difference')
xlabel('direction')
ylabel('error difference')
xlim([0 2*pi])

% % fb
% 
% Precision_fb_all = mean(cat(4,Precision_fb{idx}), 4);
% theta_temp = cellfun(@cell2mat, theta_e_all, 'UniformOutput', 0);
% theta_e_fb_all = mean(cat(3, theta_temp{idx}), 3);
% figure
% subplot(2,2,2)
% plot(THETA, squeeze(Precision_fb_all(2,1,:) - Precision_fb_all(1,1,:)), 'k')
% title('detection rate difference')
% xlabel('direction')
% ylabel('detection rate difference')
% xlim([0 2*pi])
% 
% subplot(2,2,1)
% plot(THETA, squeeze(Precision_fb_all(1,1,:)))
% hold on
% plot(THETA, squeeze(Precision_fb_all(2,1,:)))
% legend('all narrow', 'one broad')
% title('detection rate')
% xlabel('direction')
% ylabel('detection rate')
% xlim([0 2*pi])
% 
% subplot(2,2,3)
% plot(THETA, squeeze(theta_e_fb_all(1,:)))
% hold on
% plot(THETA, squeeze(theta_e_fb_all(2,:)))
% legend('all narrow', 'one broad')
% title('error')
% xlabel('direction')
% ylabel('error')
% xlim([0 2*pi])
% 
% subplot(2,2,4)
% plot(THETA, squeeze(theta_e_fb_all(2,:) - theta_e_fb_all(1,:)), 'k')
% title('error difference')
% title('error difference')
% xlabel('direction')
% ylabel('error difference')
% xlim([0 2*pi])
%% combine figures
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')


repeat = 5;
idx = 1:repeat;
Precision_dir_all = mean(cat(4,Precision_dir{idx}), 4);
theta_temp = cellfun(@cell2mat, theta_e_dir, 'UniformOutput', 0);
theta_e_dir_all = mean(cat(3, theta_temp{idx}), 3);
THETA = linspace(-180, 180, size(Precision_dir_all, 3)+1);
figure
yyaxis left
subplot(2,4,1)
plot(THETA, squeeze(Precision_dir_all(1,1,[1:end 1])), 'bo-', 'MarkerSize', 8)
hold on
plot(THETA, squeeze(Precision_dir_all(2,1,[1:end 1])), 'b*-', 'MarkerSize', 8)
ylabel('detection rate')
xlim([-180 180])


yyaxis right
plot(THETA, squeeze(Precision_dir_all(2,1,[1:end 1]) - Precision_dir_all(1,1,[1:end 1])), 'k')
ylabel('detection rate difference')
% legend('all narrow', 'one broad', 'difference')
title('detection')
xlabel('direction')
ylim([0, inf])

subplot(2,4,5)
yyaxis left
plot(THETA, squeeze(theta_e_dir_all(1,[1:end 1])), 'ro-', 'MarkerSize', 8)
hold on
plot(THETA, squeeze(theta_e_dir_all(2,[1:end 1])), 'r*-', 'MarkerSize', 8)
ylabel('error (degree)')

yyaxis right
plot(THETA, squeeze(theta_e_dir_all(2,[1:end 1]) - theta_e_dir_all(1,[1:end 1])), 'k')
title('discrimination')
xlabel('direction')
ylabel('error difference (degree)')
% legend('all narrow', 'one broad', 'difference')
xlim([-180 180])

% beta
Precision_beta_all = mean(cat(4,Precision_beta{idx}), 4);
theta_temp = cellfun(@cell2mat, theta_e_beta, 'UniformOutput', 0);
theta_e_beta_all = mean(cat(3, theta_temp{idx}), 3);
yyaxis left
subplot(2,4,2)
plot(THETA, squeeze(Precision_beta_all(1,1,[1:end 1])), 'bo-', 'MarkerSize', 8)
hold on
plot(THETA, squeeze(Precision_beta_all(2,1,[1:end 1])), 'b*-', 'MarkerSize', 8)
ylabel('detection rate')

yyaxis right
plot(THETA, squeeze(Precision_beta_all(2,1,[1:end 1]) - Precision_beta_all(1,1,[1:end 1])), 'k')
ylabel('detection rate difference')
% legend('all narrow', 'one broad', 'difference')
title('detection')
xlabel('direction')
xlim([-180 180])

subplot(2,4,6)
yyaxis left
plot(THETA, squeeze(theta_e_beta_all(1,[1:end 1])), 'ro-', 'MarkerSize', 8)
hold on
plot(THETA, squeeze(theta_e_beta_all(2,[1:end 1])), 'r*-', 'MarkerSize', 8)
ylabel('error (degree)')

yyaxis right
plot(THETA, squeeze(theta_e_beta_all(2,[1:end 1]) - theta_e_beta_all(1,[1:end 1])), 'k')
title('discrimination')
xlabel('direction')
ylabel('error difference (degree)')
% legend('all narrow', 'one broad', 'difference')
xlim([-180 180])

% rmax
Precision_rmax_all = mean(cat(4,Precision_rmax{idx}), 4);
theta_temp = cellfun(@cell2mat, theta_e_rmax, 'UniformOutput', 0);
theta_e_rmax_all = mean(cat(3, theta_temp{idx}), 3);
yyaxis left
subplot(2,4,3)
plot(THETA, squeeze(Precision_rmax_all(1,1,[1:end 1])), 'bo-', 'MarkerSize', 8)
hold on
plot(THETA, squeeze(Precision_rmax_all(2,1,[1:end 1])), 'b*-', 'MarkerSize', 8)
ylabel('detection rate')

yyaxis right
plot(THETA, squeeze(Precision_rmax_all(2,1,[1:end 1]) - Precision_rmax_all(1,1,[1:end 1])), 'k')
ylabel('detection rate difference')
% legend('all narrow', 'one broad', 'difference')
title('detection')
xlabel('direction')
xlim([-180 180])

subplot(2,4,7)
yyaxis left
plot(THETA, squeeze(theta_e_rmax_all(1,[1:end 1])), 'ro-', 'MarkerSize', 8)
hold on
plot(THETA, squeeze(theta_e_rmax_all(2,[1:end 1])), 'r*-', 'MarkerSize', 8)
ylabel('error (degree)')

yyaxis right
plot(THETA, squeeze(theta_e_rmax_all(2,[1:end 1]) - theta_e_rmax_all(1,[1:end 1])), 'k')
title('discrimination')
xlabel('direction')
ylabel('error difference (degree)')
% legend('all narrow', 'one broad', 'difference')
xlim([-180 180])

% all
Precision_all_all = mean(cat(4,Precision_all{idx}), 4);
theta_temp = cellfun(@cell2mat, theta_e_all, 'UniformOutput', 0);
theta_e_all_all = mean(cat(3, theta_temp{idx}), 3);
yyaxis left
subplot(2,4,4)
plot(THETA, squeeze(Precision_all_all(1,1,[1:end 1])), 'bo-', 'MarkerSize', 8)
hold on
plot(THETA, squeeze(Precision_all_all(2,1,[1:end 1])), 'b*-', 'MarkerSize', 8)
ylabel('detection rate')

yyaxis right
plot(THETA, squeeze(Precision_all_all(2,1,[1:end 1]) - Precision_all_all(1,1,[1:end 1])), 'k')
ylabel('detection rate difference')
% legend('all narrow', 'one broad', 'difference')
title('detection')
xlabel('direction')
xlim([-180 180])
ylim([0, 0.04])

subplot(2,4,8)
yyaxis left
plot(THETA, squeeze(theta_e_all_all(1,[1:end 1])), 'ro-', 'MarkerSize', 8)
hold on
plot(THETA, squeeze(theta_e_all_all(2,[1:end 1])), 'r*-', 'MarkerSize', 8)
ylabel('error (degree)')

yyaxis right
plot(THETA, squeeze(theta_e_all_all(2,[1:end 1]) - theta_e_all_all(1,[1:end 1])), 'k')
title('discrimination')
xlabel('direction')
ylabel('error difference (degree)')
% legend('all narrow', 'one broad', 'difference')
xlim([-180 180])


%% compare different number of broad tuning

contrast = [40]; % contrast
a = 1.8; b = 0.8;
Alpha2 = [a a a a; a a b a; a b b a; b b b a; b b b b]; % tuning width 


c50 = [1.8 1.8 1.8 1.8];
alpha1 = [7 7 7 7];
% tuning curve
phi = [0 1/2 1 3/2]*pi; % preferred direction
r_base = [1 1 1 1]*0.02; % background firing
repeat_n = 10000;
a = 0.95; b = 0.6;
Beta = [a a a a; a a b a; a b b a; b b b a; b b b b];
a = 1; b = 2;
Ymax = [a a a a; a a b a; a b b a; b b b a; b b b b]*rmax;
[Precision_all, theta_e_all] = deal(cell(1, repeat));
for rp = 1:repeat
    clear Precision_ratio
    for k = 1:size(Beta, 1)
        ymax = Ymax(k, :); % ymax = 5 for rod light level
        alpha2 = Alpha2(k, :);
        for ctr = 1:length(contrast)
            [~, Precision(k, ctr, :)] = run_poiss_Detection_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,Beta(k,:));
            [error(k,ctr), theta_e{k,ctr}] = run_poiss_MLdecoder_beta(c50,alpha1,phi,alpha2,r_base,ymax,contrast(ctr),repeat_n,Beta(k,:));
    %         [width ctr]
            disp(['repeat ' num2str(rp) ' detection task: beta ' num2str(k) ' contrast ' num2str(ctr) ' finished.'])
        end
    end
    Precision_all{rp} = Precision;
    ERROR_rod{repeat} = error;
    theta_e_all{rp} = theta_e;

%     save('dsmodel_direction.mat', 'Precision_rod', 'Precision_cone', 'ERROR_rod', 'ERROR_cone', 'theta_e_rod', 'theta_e_cone')
end

%%
theta = linspace(0, 360, 17);
theta = theta(1:end-1);

tuning1 = (0.5 + 0.5 * cos(theta/180*pi + pi)).^1.8*0.95+(1-0.95);
tuning2 = (0.5 + 0.5 * cos(theta/180*pi + pi)).^0.8*0.95+(1-0.95);
tuning3 = (0.5 + 0.5 * cos(theta/180*pi + pi)).^1.8*0.6+(1-0.6);
tuning4 = (0.5 + 0.5 * cos(theta/180*pi + pi)).^0.8*0.6+(1-0.6);


figure(5)
subplot(5,1,1)
plot(theta, tuning1, 'o-')
hold on 
plot(theta, circshift(tuning1, 4, 2), 'o-')
plot(theta, circshift(tuning1, 8, 2), 'o-')
plot(theta, circshift(tuning1, -4, 2), 'o-')
xlim([0, 360])

subplot(5,1,2)
% plot(theta, tuning1, '*-')
plot(theta, circshift(tuning1, 4, 2), '*-')
hold on 
plot(theta, circshift(tuning1, 8, 2), '*-')
plot(theta, circshift(tuning1, -4, 2), '*-')
xlim([0, 360])


subplot(5,1,3)
plot(theta, tuning2, '*-')
xlim([0, 360])

subplot(5,1,4)
plot(theta, tuning3, '*-')
xlim([0, 360])

subplot(5,1,5)
plot(theta, tuning4, '*-')
xlim([0, 360])
