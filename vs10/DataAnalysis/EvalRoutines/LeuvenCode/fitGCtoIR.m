function [xout,resnormout,residualout,exitflagout,outputout] = fitGCtoIR(xvals,yvals)

F = @(x,xdata)x(1).*(xdata.^3).*exp(-xdata./x(2)).*cos(2*pi*((x(3)*xdata)+(0.5*x(4)*(xdata.^2))) + x(5));

% F = @(x,xdata)x(1).*exp(-(xdata.^2)/x(2)).*cos(2*pi*((x(3)*xdata)+(0.5*x(4)*(xdata.^2))) + x(5));

opts = optimoptions('lsqcurvefit','MaxFunEvals',10e4,'MaxIter',10e4,'TolFun',1e-6,'TolX',1e-6);

x0 = [-max(abs(yvals)) 1 1 1 1];
[x{1},resnorm{1},residual{1},exitflag{1},output{1}] = lsqcurvefit(F,x0,xvals',yvals,[],[],opts);
x0 = [max(abs(yvals)) 1 1 1 1];
[x{2},resnorm{2},residual{2},exitflag{2},output{2}] = lsqcurvefit(F,x0,xvals',yvals,[],[],opts);

if resnorm{1}>resnorm{2}
    xout = x{2};
    resnormout = resnorm{2};
    residualout = residual{2};
    exitflagout = exitflag{2};
    outputout = output{2};
else
    xout = x{1};
    resnormout = resnorm{1};
    residualout = residual{1};
    exitflagout = exitflag{1};
    outputout = output{1};
end
return;