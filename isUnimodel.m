function unimodel = isUnimodel(A)
Diff = diff(A);
modality = sum(Diff(1:end-1).*Diff(2:end) < 0);
if modality < 2
    unimodel = true;
else
    unimodel = false;
end

end