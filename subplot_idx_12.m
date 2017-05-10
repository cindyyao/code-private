function [idx, xx, yy] = subplot_idx_12(x, y)


xx = 5*x-1; yy = 5*y-1;
idx = zeros(x*y, 16);
for i = 1:y
    idx(i, :) = [yy+4 4 3 2 1 yy+1 2*yy+1 3*yy+1 3*yy+2 3*yy+3 3*yy+4 2*yy+4 yy+2 yy+3 2*yy+2 2*yy+3] + 5*(i-1);
end

if x > 1
    for i = 2:x
        idx(y*(i-1)+1:y*i, :) = idx(1:y, :) + yy*5*(i-1);
    end
end

end

