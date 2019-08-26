function out = exclude_outliers(in, std_n)

sigma = std(in);
miu = mean(in);
in(in > miu + std_n*sigma) = [];
in(in < miu - std_n*sigma) = [];
out = in;

end

