function [rf, rf_std] = get_fs_rf(fs_spike, field_width, field_height,sub_region)

rf = cell(length(fs_spike), 1);
rf_std = cell(length(fs_spike), 1);
for cc = 1:length(fs_spike)
    if ~isempty(fs_spike{cc})
        repeats = size(fs_spike{cc}, 2);
        fs_spike_all = reshape(fs_spike{cc}, field_width, field_height, repeats, 2);
        fs_spike_all = permute(fs_spike_all,[2 1 3 4]);
        fs_spike_mean = squeeze(mean(fs_spike_all, 3));
        fs_spike_std = squeeze(std(fs_spike_all, [], 3));
        if sub_region
             fs_spike_mean = repmat(fs_spike_mean, 2, 2);
             fs_spike_std = repmat(fs_spike_std, 2, 2);
        end
        rf{cc} = fs_spike_mean;
        rf_std{cc} = fs_spike_std;
    end
end