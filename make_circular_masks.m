function masks = make_circular_masks(center, radius, width, height, scale)
% masks = make_circular_masks(center, radius, width, height)
% size(center) = [N, 2];
% length(radius) = N;

if nargin < 5
    scale = 1;
end

if size(center,1) ~= length(radius) || size(center,2) ~= 2
    error('dimension mismatch')
end

N = length(radius);
masks = cell(N,1);
for i = 1:N
    masks{i} = make_circular_mask(center(i,:), radius(i), width, height, scale);
end

end


