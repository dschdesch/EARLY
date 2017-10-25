function [COD] = Coeff_Determination(M,P)

ssd = sum((P-M).^2);
ss = sum((M-mean(M)).^2);
COD = 1-(ssd/ss);