function samples = sample_dist(pdf, N)

cdf = cumsum(pdf);
cdf = cdf/cdf(end);
samples = zeros(1, N);
for i = 1:N
    samples(i) = sum(rand > cdf);
end
end

