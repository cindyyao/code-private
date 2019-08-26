function rho = bgnd_subtraction(rho1, bgnd)

rho = cell(length(rho1), 1);
for j = 1:length(rho1)
    rho{j} = rho1{j} - repmat(bgnd', 1, size(rho1{j}, 2))*8; % 8 second
    rho{j} = max(rho{j}, 0);
end

end