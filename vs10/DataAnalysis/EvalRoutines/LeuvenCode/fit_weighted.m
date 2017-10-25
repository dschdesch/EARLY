function [beta_out,gof_out]=fit_weighted(varargin)
if nargin ==3
    x_in = varargin{1};
    y_in = varargin{2};
    weight_in = varargin{3};
    weights = weight_in;
    weights = weights./max(weights);
    s = fitoptions('Method','LinearLeastSquares','Lower',[-Inf -Inf],'Upper',[Inf Inf],'Weights',weights);
elseif nargin ==2
    x_in = varargin{1};
    y_in = varargin{2};
    s = fitoptions('Method','LinearLeastSquares','Lower',[-Inf -Inf],'Upper',[Inf Inf]);
end
f = fittype({'x','1'},'coefficients',{'a','b'},'options',s);
[beta_out,gof_out] = fit(x_in(:),y_in,f);

end