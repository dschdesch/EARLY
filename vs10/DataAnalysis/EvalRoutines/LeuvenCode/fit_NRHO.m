function [out] = fit_NRHO(xvals,yvals)


if yvals(1)<yvals(end)
    fpower=@(p,x)p(1)+p(2).*((1+x)/2).^p(3);
else
    fpower=@(p,x)p(1)+p(2).*((1-x)/2).^p(3);
end
options = optimset('display','off');

out.beta = lsqcurvefit(fpower,[1 1 1],xvals,yvals,[-Inf 0 0],[Inf Inf Inf],options);
out.fittedyvals = fpower(out.beta,xvals);

[out.COD] = Coeff_Determination(yvals,out.fittedyvals);