function [dsindex] = DS_index_vector(MAG, RHO);

dsindex = cell(size(MAG));
for i = 1:size(MAG, 1)
    for j = 1:size(MAG, 2)
        dsindex{i, j} = MAG{i, j}./sum(RHO{i, j}');
    %     [dsindex{i,1} ntimp] = exciseRows(dsindex{i,1}', zeros(size(dsindex{i,1},2),1));
        dsindex{i, j} = exciseRows(dsindex{i, j}');
    end
end

end