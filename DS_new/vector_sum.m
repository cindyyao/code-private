function[U V angle mag] = vector_sum(X, Y)

%Function that sums spike numbers in all directions     

%Input: X and Y magnitudes of all cells for all temporal periods in all directions

%Output: Cartesian (U,V) and Polar (angle, mag) forms of calculated vector sum

% Sneha Ravi 
% Last revision: 12-18-2012


U = cell(size(X));
V = cell(size(X));
angle = cell(size(X)); 
mag = cell(size(X)); 

for i = 1:size(X,1)
    for j = 1:size(X,2)
        U{i,j} = sum(X{i,j}');
        V{i,j} = sum(Y{i,j}');
        [angle{i,j},mag{i,j}] = cart2pol(U{i,j}, V{i,j}); 
    end
end
end