function [I, Yr, rdev] = sparsify(X,Y, rtol);
% sparsify - approximation of a function 
%   I = sparsify(X,Y,dtol) finds index array I such that Y(I) approximates
%   Y in such a way that the reconstruction Yr = interp1(X(I),Y(I),X) does
%   not deviate from Y by more than dtol*abs(Y). That is,
%       abs(Yr-Y)<=dtol*abs(Y).
%   X must be uniformly increasing array. X and Y must have the same size.
%
%   [I, Yr, Rdev] = sparsify(X,Y,dtol) also returns the reconstruction Yr
%   described above, and the relative deviation Rdev = abs((Yr-Y)./Y).

maxDev = rtol*abs(Y); % max absolute deviation

I = 1:numel(Y);
while 1,
    N0 = numel(I);
    itry = [I(1:2:end-1), I(end)];
    iout = setdiff(I,itry);
    Yr = interp1(X(itry), Y(itry), X);
    idev = intersect(iout, find(abs(Yr-Y)>maxDev*0.8));
%<<<<<<< HEAD
    I = sort(union(itry,idev)); % was I = sort(union(itry,idev))' and caused errors; hsin-wei Aug.1, 2017 
%=======
%    I = sort(union(itry,idev))';% no! who did this??? Github?
%>>>>>>> parent of f377d81... Change  sparsify transpose
    if numel(I)>0.85*N0, break; end
end

Yr = interp1(X(I), Y(I), X);
rdev = abs((Yr-Y)./Y);
rdev(Y==0) = 0;







