function[dsindex] = DS_index_one(mag, null)

% Function calculates and returns direction selective index, dsi = Vector-Null/Vector+Null

% Sneha Ravi 
% Last revision: 12-18-2012

dsindex = cell(size(mag));
for i = 1:size(mag,1)
    for j = 1:size(mag,2)
        dsindex{i,j} = (mag{i,j}-null{i,j})./(mag{i,j}+null{i,j});
    %     [dsindex{i,1} ntimp] = exciseRows(dsindex{i,1}', zeros(size(dsindex{i,1},2),1));
        dsindex{i,j} = exciseRows(dsindex{i,j}');
    end
end
end