function[null] = nullmin(rho)

%Function returns spike number of direction with least firing

% Sneha Ravi 
% Last revision: 12-18-2012

null = cell(size(rho));%
for i = 1:size(rho,1)
    for j = 1:size(rho,2)
        null{i,j} = min(rho{i,j},[], 2)';
    end
end