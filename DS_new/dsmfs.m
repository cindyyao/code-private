%% load data
dataset{1} = load_tostruct('DS160519.mat');
dataset{2} = load_tostruct('DS160829.mat');
dataset{3} = load_tostruct('DS161019.mat');
dataset{4} = load_tostruct('DS161115.mat');
dataset{5} = load_tostruct('DS161122.mat');

%% Superior cells
ONOFF = {'ON', 'OFF'};
figure
for onoff = 1:2
    subplot(1, 2, onoff)
    clear mean_temp ste_temp
%     for ll = 2:3
%         mean_temp(ll-1) = mean(dataset{2}.rf_area_clean{ll}{1}{onoff});
%         ste_temp(ll-1) = std(dataset{2}.rf_area_clean{ll}{1}{onoff})/sqrt(length(dataset{2}.rf_area_clean{ll}{1}{onoff}));
%     end
%     errorbar([2 4], mean_temp, ste_temp, 'color', 'b')
%     hold on
    
    clear mean_temp ste_temp
    for ll = 1:4
        mean_temp(ll) = mean(dataset{3}.rf_area_clean{ll}{1}{onoff});
        ste_temp(ll) = std(dataset{3}.rf_area_clean{ll}{1}{onoff})/sqrt(length(dataset{3}.rf_area_clean{ll}{1}{onoff}));
    end
    errorbar([0 1 2 4], mean_temp, ste_temp, 'color', 'r')
    hold on
    
    clear mean_temp ste_temp
    for ll = 1:5
        mean_temp(ll) = mean(dataset{5}.rf_area_clean{ll}{1}{onoff});
        ste_temp(ll) = std(dataset{5}.rf_area_clean{ll}{1}{onoff})/sqrt(length(dataset{5}.rf_area_clean{ll}{1}{onoff}));
    end
    errorbar(0:4, mean_temp, ste_temp, 'color', 'k')
    legend('stixelsize20', 'stixelsize40', 'stixelsize30')
    xlabel('log(background intensity)')
    ylabel('RF area (mm^2)')
    title(ONOFF{onoff})
    xlim([-1 5])
    ylim([0 0.4])
end

%% 'other' cells
ONOFF = {'ON', 'OFF'};
figure
for onoff = 1:2
    subplot(1, 2, onoff)
    clear mean_temp ste_temp
    for ll = 1:3
        mean_temp(ll) = mean(dataset{1}.rf_area_clean_all{ll}{onoff});
        ste_temp(ll) = std(dataset{1}.rf_area_clean_all{ll}{onoff})/sqrt(length(dataset{1}.rf_area_clean_all{ll}{onoff}));
    end
    errorbar([0 2 4], mean_temp, ste_temp, 'color', 'g')
    hold on

    clear mean_temp ste_temp
    for ll = 2:3
        mean_temp(ll-1) = mean(dataset{2}.rf_area_clean_all{ll}{onoff});
        ste_temp(ll-1) = std(dataset{2}.rf_area_clean_all{ll}{onoff})/sqrt(length(dataset{2}.rf_area_clean_all{ll}{onoff}));
    end
    errorbar([2 4], mean_temp, ste_temp, 'color', 'b')
    
%     clear mean_temp ste_temp
%     for ll = 1:4
%         mean_temp(ll) = mean(dataset{3}.rf_area_clean_all{ll}{onoff});
%         ste_temp(ll) = std(dataset{3}.rf_area_clean_all{ll}{onoff})/sqrt(length(dataset{3}.rf_area_clean_all{ll}{onoff}));
%     end
%     errorbar([0 1 2 4], mean_temp, ste_temp, 'color', 'r')
    
    clear mean_temp ste_temp
    for ll = 1:5
        mean_temp(ll) = mean(dataset{4}.rf_area_clean_all{ll}{onoff});
        ste_temp(ll) = std(dataset{4}.rf_area_clean_all{ll}{onoff})/sqrt(length(dataset{4}.rf_area_clean_all{ll}{onoff}));
    end
    errorbar(0:4, mean_temp, ste_temp, 'color', 'm')

    clear mean_temp ste_temp
    for ll = 1:5
        mean_temp(ll) = mean(dataset{5}.rf_area_clean_all{ll}{onoff});
        ste_temp(ll) = std(dataset{5}.rf_area_clean_all{ll}{onoff})/sqrt(length(dataset{5}.rf_area_clean_all{ll}{onoff}));
    end
    errorbar(0:4, mean_temp, ste_temp, 'color', 'k')
    legend('stixelsize15', 'stixelsize20', 'stixelsize15', 'stixelsize30')
    xlabel('log(background intensity)')
    ylabel('RF area (mm^2)')
    title(ONOFF{onoff})
    xlim([-1 5])
    ylim([0 0.04])
end
