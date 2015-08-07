function[X Y] = get_XY(rho, theta)

%Function converts magnitude of spike numbers in each direction to X and Y polar coordinates


%Sneha Ravi 
  %Last revision: 12-18-2012

X = cell(size(rho));
Y = cell(size(rho));

for i = 1:size(rho,1)
    for j = 1:size(rho,2)
        [X{i,j} Y{i,j}] = pol2cart(theta{i,j}, rho{i,j}); 

    end
end
end