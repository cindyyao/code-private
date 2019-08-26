function charge = overageData(datausage, bill, line_n)
% charge = overageData(datausage, bill)
% datausage: vector of data usage of every member
% overaged bill (15/GB)

overage = max(datausage - 15/line_n, 0);
charge = overage/sum(overage)*bill;

end