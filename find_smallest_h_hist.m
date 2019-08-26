function [h, filteredA] = find_smallest_h_hist(samples, max_lag)
A = hist(samples, [0.5:1:max_lag*2+1]);

max_lag = (length(A) - 1)/2;
epsilon = 0.001;
h_low = 0.001;
h_high = max_lag;
sz = max_lag;  % length of gaussFilter vector
x = linspace(-sz, sz, sz*2+1);

while(h_high - h_low > epsilon)
    h_mid = (h_low + h_high)/2;
    gaussFilter = exp(-x .^ 2 / (2 * h_mid ^ 2));
    gaussFilter = gaussFilter / sum (gaussFilter); % normalize
    filteredA = conv(A, gaussFilter, 'same');
    if isUnimodel(filteredA)
        h_high = h_mid;
    else
        h_low = h_mid;
    end
end
h = h_low;

end
