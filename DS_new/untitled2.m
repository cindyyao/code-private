a = ceil((length(trigger_set_i)+1)/2);
for cc =10:10%length(nds_flash_raster);
    ts = 1;
    b = 1;
    trial_n = length(nds_dark_raster{1});
    FigHandle = figure;
    set(FigHandle, 'Position', [0, 0, 1920, 1080]);
    while ts <= length(trigger_set_i)+1
        if  ts == a+1
            b = 3;
        end
        subplot(a, 6, b)
        if ts > 1
            trial_n = length(nds_flash_raster{cc}{ts-1});
        end
        for j = 1:trial_n
            if ts == 1
                SpikeTime = nds_dark_raster{cc}{j};
            else
                SpikeTime = nds_flash_raster{cc}{ts-1}{j};
            end
            SpikeTime = SpikeTime';
            X = [SpikeTime; SpikeTime];
            Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
            line(X, Y, 'color', 'b');
            axis([0, 3, 0, trial_n]);
            hold on
        end
        if b == 1
            title(num2str(ds_id_flash(cc)))
        end
        b = b+1;
        subplot(a, 6, b)
        if ts > 1
            bar(XX, nds_flash_hist_mean{cc}{ts-1}, 1)
        else
            bar(XX, nds_dark_hist_mean{cc}, 1)
        end
        xlim([0 3])

        b = b+5;
        ts = ts+1;
    end

%     print_close(1, [24 12], num2str(cc))
    subplot(a, 6, [5 6 11 12])
    plot(log10(Irel), nds_Pc(10, :))
    xlabel('log(R*/rod)')
    ylabel('probability')

end


flash = nds_flash_raster{10}{11};
dark = nds_dark_raster{10};
seq = [];
for i = 1:60
    seq = [seq; [randperm(2)-1]];
end

for i = 1:60
    hist_flash(i, :) = hist(flash{i}, XX);
    hist_dark(i, :) = hist(dark{i}, XX);
end
hist_flash_sum = sum(hist_flash);
hist_dark_sum = sum(hist_dark);
for i = 1:60
    template_flash = hist_flash_sum - hist_flash(i, :);
    template_flash = template_flash/norm(template_flash);
    template_flash(isnan(template_flash)) = 0;
    template_dark = hist_dark_sum - hist_dark(i, :);
    template_dark = template_dark/norm(template_dark);
    template_dark(isnan(template_dark)) = 0;
    DV = template_flash - template_dark;
    corr_flash(i) = hist_flash(i, :) * DV';
    corr_dark(i) = hist_dark(i, :) * DV';
end
pc = (sum(corr_flash > corr_dark) + sum(corr_flash == corr_dark)/2)/60;
% pc = (sum(corr_flash > 0) + sum(corr_dark <= 0))/(60*2)
figure
i = 1; j = 1;
while j <= 120
    subplot(12, 10, seq(i, 1)*10+j)
    plot([1:150]/50, hist_flash(i, :))
    subplot(12, 10, seq(i, 2)*10+j)
    plot([1:150]/50, hist_dark(i, :))
    if mod(i, 10) == 0
        j = j + 10;
    end
    i = i + 1;
    j = j + 1;
end



%%
cc = 10;
figure
flash_trials = ds_flash_hist_trial_cd{2}{cc}{11};
dark_trials = ds_dark_hist_trial_cd{2}{11};
for i = 1:60
    subplot(11, 15, floor((i-1)/15)*45+mod(i-1, 15)+1)
    plot(flash_trials{i})
    axis off
    subplot(11, 15, floor((i-1)/15)*45+mod(i-1, 15)+16)
    plot(dark_trials{i})
    axis off
end