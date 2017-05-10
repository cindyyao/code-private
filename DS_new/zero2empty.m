function matrix = zero2empty(matrix)
matrix(sum(matrix, 2) == 0,:) = [];
matrix(sum(matrix, 2) == 0,:) = [];
end