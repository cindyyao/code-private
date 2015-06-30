function matrix = nan2empty(matrix)
matrix(any(isnan(matrix),2),:) = [];
matrix(any(isinf(matrix),2),:) = [];
end
