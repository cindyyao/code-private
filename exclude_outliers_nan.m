function out = exclude_outliers_nan(in, std_n)

sigma = nanstd(in);
miu = nanmean(in);
in(in > miu + std_n*sigma) = nan;
in(in < miu - std_n*sigma) = nan;
out = in;

end