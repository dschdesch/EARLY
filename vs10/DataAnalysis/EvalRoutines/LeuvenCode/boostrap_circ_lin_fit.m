function [clfit_intersect,clfit_slope] = boostrap_circ_lin_fit(X,Y,W,nboot)
%X is usually frequency
%Y is phase in cycles
% nboot = 1;
N = length(Y);
for i = 1:nboot
    indx = randi(N,N,1);
    [ph(i),D] = fit_circ_lin(X(indx),Y(indx),W,0,1);
    delay(i) = -D;
end

clfit_intersect = prctile(ph,[5 50 95]);
clfit_slope = prctile(delay,[5 50 95]);