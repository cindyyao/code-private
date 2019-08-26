for cc = 1:length(nds_flash_raster)
    raster_1cell = nds_flash_raster{cc};
    for ts = 1:length(trigger_set_i)
        for t = 1:length(raster_1cell{ts})
            hist_1cell{ts}{t} = hist(raster_1cell{ts}{t}, XX);
        end
        raster_1cell_all{ts} = sort(cell2mat(raster_1cell{ts}));
        hist_1cell_all{ts} = hist(raster_1cell_all{ts}, XX)/length(trigger_set_i{ts});
    end
    nds_flash_hist_trial{cc} = hist_1cell;
    nds_flash_hist_mean{cc} = hist_1cell_all;
end

%% get raster and psth for dark trials
% use 400-580 second (roughly) as dark trials

for cc = 1:length(nds_dark_raster)
    for t = 1:length(nds_dark_raster{cc})
        nds_dark_hist_trial{cc}{t} = hist(nds_dark_raster{cc}{t}, XX);
    end
    nds_dark_raster_all{cc} = sort(cell2mat(nds_dark_raster{cc}));
    nds_dark_hist_mean{cc} = hist(nds_dark_raster_all{cc}, XX)/length(trigger);
end
window = 1.5;
bin_n = window/bin_size;

nds_Pc = zeros(length(nds_dark_raster), length(trigger_set_i));
for cc = 1:length(nds_dark_raster)
    for ts = 1:length(trigger_set_i)
        trial_n = min(length(nds_flash_raster{1}{ts}), 60);
%         trial_n = 60;
        corr_flash = zeros(trial_n, 1);
        corr_dark = zeros(trial_n, 1);
        temp = nds_dark_hist_trial{cc}(1:trial_n);
        nds_dark_hist_sum = sum(cell2mat(temp'));
        nds_flash_hist_sum = sum(cell2mat(nds_flash_hist_trial{cc}{ts}'));
        for t = 1:trial_n
            template_flash = nds_flash_hist_sum - nds_flash_hist_trial{cc}{ts}{t};
            template_flash = template_flash(1:bin_n);
            template_flash = template_flash/norm(template_flash);
            template_flash(isnan(template_flash)) = 0;
            template_dark = nds_dark_hist_sum - nds_dark_hist_trial{cc}{t};
            template_dark = template_dark(1:bin_n);
            template_dark = template_dark/norm(template_dark);
            template_dark(isnan(template_dark)) = 0;
            DV = template_flash - template_dark;
            corr_flash(t) = nds_flash_hist_trial{cc}{ts}{t}(1:bin_n) * DV(1:bin_n)';
            corr_dark(t) = nds_dark_hist_trial{cc}{t}(1:bin_n) * DV(1:bin_n)';
        end
        nds_Pc(cc, ts) = (sum(corr_flash > corr_dark) + sum(corr_flash == corr_dark)/2)/trial_n;
%         nds_Pc(cc, ts) = (sum(corr_flash > 0) + sum(corr_dark <= 0))/(trial_n*2);
    end
end
