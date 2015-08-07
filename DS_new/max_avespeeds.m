function[magave] = max_avespeeds(mag)

%Function returns average of magnitude vectors calculated across all temporal periods

% Sneha Ravi 
% Last revision: 12-18-2012
for c = 1:size(mag,2)
for i = 1:size(mag,1)
    magave_temp(i,1:length(mag{i,1})) = mag{i,c};
end
magave(c,:) = mean(magave_temp);
end
end