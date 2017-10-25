function [CP,CP_uw,CD,R,N] = fit_circ_lin(f,theta,weights,varargin)

[m,n] = size(f);
if m>1
    f = f';
end
[m,n] = size(theta);
if n>1
    theta = theta';
end
[m,n] = size(weights);
if m>1
    weights = weights';
end
%% First stage: grid search to find reasonable first guess
if isempty(varargin)
    cdgrid = -6:0.1:6;
    cpgrid = -0.5:0.01:0.5;
    [cdgrid,cpgrid] = meshgrid(cdgrid,cpgrid);
    cdgrid = cdgrid(:);
    cpgrid = cpgrid(:);
    nrgrid = length(cdgrid);
    dgrid = nan(1,nrgrid);
    theta_in = theta;
    for i = 1:nrgrid
        dgrid(i) = objective([cpgrid(i); cdgrid(i)]);
    end
    %Find the index of the cp/cd values with the lowest circular error
    [dummy,indx] = min(dgrid);
    %Make this the start point for the gradient descent in the second stage
    p0 = [cpgrid(indx);cdgrid(indx)];
else
    p0 = [varargin{1};varargin{2}];
end

%% Second stage: gradient descent from this first guess

%Set up the optimization, using gradient descent
options = optimset('HessUpdate','steepdesc','largescale','off','MaxFunEvals',1e12,'MaxIter',1e12,'display','off','TolX',1e-15);
N = sum(weights>0);

%Do the fitting: minimize the circular error between data and model. Return
%CP and CD in P(1) and P(2), respectively, and the value of the objective
%function (the circular error term) in D.
theta_in = theta;
[P,D] = fminunc(@objective,p0,options);

%Get CP and CD: restrict CP to [-0.5 0.5]
CP_uw = -P(1);
CP = wrapToPi(CP_uw*2*pi)/(2*pi);
CD = -P(2);

%Calculate residual error (to estimate goodness of fit)
R = sqrt(D/N);

    function theta_hat = model(p)
        theta_hat = p(1)+(p(2)*f);
    end
    function d = objective(p)   
        theta_hat = model(p);
        d = sum(weights'.*(0.5*(1-cos(2*pi*(theta_in - theta_hat')))));
    end
end













































