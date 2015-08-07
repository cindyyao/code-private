function[magmax] = max_allspeeds(mag)

%Function returns maximum of magnitude vectors over all temporal periods

% Sneha Ravi 
% Last revision: 12-18-2012

for c = 1:size(mag,2)
for i = 1:size(mag,1)
    magmax_temp(i,1:length(mag{i,1})) = mag{i,c};
end
magmax(c,:) = max(magmax_temp);
end
end
