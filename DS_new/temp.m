%% max window
% on-off DSGC
clear Max_i
trial_dur = mean(diff(datamb{1}.stimulus.triggers));
bin_size = 0.01;
xx = bin_size/2:bin_size:trial_dur-bin_size/2;
dscell_type = {'superior', 'anterior', 'inferior', 'posterior'};
condition = {'control', 'AP5', 'AP5+Hex', 'wash'};
step_size = 5;
for drug = 1:4
    for dir = 1:4
        CC = 1;
        for cc = 1:length(idx_dir{dir})
            if ~isempty(raster_p_sum{drug}{idx_dir{dir}(cc)})
                for repeat = 1:10
                    for ctr = 7:-1:1
                        a = raster_p_sum_all{drug}{idx_dir{dir}(cc)}{:, :, ctr, repeat};
                        hist_temp = hist(a, xx);
%                         if drug == 1 && ctr == 7
                            [max_p, max_i] = max(conv(hist_temp, ones(1,step_size), 'valid'));
%                             Max_i{dir}(CC, repeat) = max_i;
%                         else
%                             max_p = conv(hist_temp, ones(1,step_size), 'valid');
%                             max_p = max_p(Max_i{dir}(CC, repeat));
%                         end
                        response_pmax{drug}{dir}(CC, ctr, repeat) = max_p - bgnd_firing{dir}(drug, CC)*bin_size*step_size;
                    end
                end
                CC = CC + 1;
            end
        end
        response_pmax{drug}{dir} = mean(response_pmax{drug}{dir}, 3);
        response_pmax_norm{drug}{dir} = response_pmax{drug}{dir}./repmat(max(response_pmax{1}{dir}, [], 2), 1, size(response_pmax{drug}{dir},2));
%         response_pmax_norm{drug}{dir} = response_pmax{drug}{dir}./repmat(max(response_pmax{drug}{dir}, [], 2), 1, size(response_pmax{drug}{dir},2));
    end
end


ctr_x = [5 10 20 40 80 150 300];
color = 'brgkc';
figure
set(gcf, 'Position', [1 1 900 800])

for dir = 1:4
    subplot(2,2,dir)
    for drug = 1:4
        errorbar(ctr_x, mean(response_pmax_norm{drug}{dir}), std(response_pmax_norm{drug}{dir})/sqrt(size(response_pmax_norm{drug}{dir}, 1)), 'color', color(drug));
        hold on
    end
    legend(condition, 'location', 'northwest')
    set(gca, 'Xscale', 'log')
    xlabel('% contrast')
    ylabel('spike rate')
    title(dscell_type{dir})
%     xlim([3 400])
end
