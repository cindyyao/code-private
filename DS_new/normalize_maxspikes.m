function[rho theta] = normalize_maxspikes(rho, theta)

%Function that normalizes by the direction with the maximum spike rate
%Inputs are rho (the values) and theta(the angles)
%Function returns the normalized radius and angle matrix after removing the rows that divide
%by zeros and have Nans and Inf

%Sneha Ravi 
  %Last revision: 12-18-2012
  
for j = 1:size(rho,1)
    for c = 1:size(rho,2)
        norm = [];
        norm = max(rho{j,c}')';
        norm = repmat(norm, 1, size(rho{j,c},2));
        rho{j,c} = rho{j,c}./norm;
    %     [rho{j,1} theta{j,1}] = exciseRows(rho{j,1}, theta{j,1});
        rho{j,c} = exciseRows(rho{j,c});
    end
end
end
