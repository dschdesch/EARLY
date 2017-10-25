function [meanlowess, stdlowess] = mylowessbootstrap(varargin)
%M. Sayles. KU Leuven. 18/11/14
%Function to compute bootstrapped lowess smoothing fits to data. Returns
%the mean lowess fit, and the standard deviation of the estimate. Assumes
%Gaussian distribution of the bootstrapped estimates.
%Usage: [Ylowessfit,YlowessSTD] =
%mylowessbootstrap(Xdata,Ydata,span,'lowess',nboot);
%Typical values for span ~0.2, and for nboot ~200 works well in most cases.
if length(varargin)==5
    X = varargin{1};
    Y = varargin{2};
    span = varargin{3};
    smoothmethod = varargin{4};
    nboot = varargin{5};
elseif length(varargin)==4
    X = varargin{1};
    Y = varargin{2};
    span = varargin{3};
    smoothmethod = varargin{4};
    nboot = 200;
elseif length(varargin)==3
    X = varargin{1};
    Y = varargin{2};
    span = varargin{3};
    smoothmethod = 'lowess';
    nboot = 200;
end

xy = [X Y];
f = @(xy) mylowess(xy,X,span,smoothmethod);
yboot2 = bootstrp(nboot,f,[X,Y])';
meanlowess = mean(yboot2,2);
stdlowess = std(yboot2,0,2);
return