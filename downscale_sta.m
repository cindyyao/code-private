function newm = downscale_sta(m, scale)
        
filter = ones(scale)/(scale^2);
for i = 1:size(m, 3)
    for j = 1:size(m, 4)
        mtemp = filter2(filter, squeeze(m(:, :, i, j)), 'same');
        mtemp = mtemp(1:scale:end, 1:scale:end);
        newm(:, :, i, j) = mtemp;
    end
end

