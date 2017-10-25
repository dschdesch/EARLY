function [bayesian_model] = Make_Bayesian_Model(spikes,wv,kernels,dt_s)

%2 ears
Nchan = 2;

%spike times as sample numbers
ST_samps = floor(spikes/(dt_s*1000));

%Filter the input signals
s = zeros(size(wv));
for n = 1:Nchan
    %Now s is scaled as Pascals
    temp = conv(wv(:,n),kernels.h1c(:,n),'full');
    s(:,n) = temp(1:length(wv));
    %useful for the equal-variance scale
    s(:,n) = s(:,n)/std(s(:,n));
    %Get the envelope
    senv(:,n) = abs(hilbert(s(:,n)));
end
%Get likelihood function
ss = s(ST_samps,:);
sigmalim = 4;
[dummy,L,X,Y] = kde2d(ss,2^5,[-sigmalim -sigmalim],[sigmalim sigmalim]);
L = L/max(L(:));

ssenv = senv(ST_samps,:);
sigmalim = 4;
[dummy,Lenv] = kde2d(ssenv,2^5,[0 0],[sigmalim sigmalim]);
Lenv = Lenv/max(Lenv(:));


%Get the prior
Xpdf = normpdf(X(1,:),0,1);
Ypdf = normpdf(Y(:,1),0,1);
[XX,YY] = meshgrid(Xpdf,Ypdf);
P = XX.*YY;
P = P/max(P(:));

[dummy,Penv] = kde2d(senv,2^5,[0 0],[sigmalim sigmalim]);
Penv = Penv/max(Penv(:));

%Get the measured IO non-linearity
IONL = L./P;
ScaleFactor = max(IONL(:));
IONL = IONL/ScaleFactor;

%Make a nice figure showing the IO non-linearity. Plot the sqrt of the IONL
%and scale the color axis as the square of the values - this makes the
%color scale quadratic, so the data display better visually.
figure;
pcolor(X,Y,sqrt(IONL));
colormap (flipud(hot));
axis square;
shading flat;
view(0,90);
set(gca,'xlim',[-(sigmalim+1) sigmalim+1],'ylim',[-(sigmalim+1) sigmalim+1],'clim',[0 1],'fontsize',14);
xlabel ('s_{ipsi}','interpreter','tex');
ylabel ('s_{contra}','interpreter','tex');
hcb = colorbar;
set(hcb,'fontsize',14,'linewidth',1,'ytick',linspace(0,1,6),'yticklabel',linspace(0,1,6),'ylim',[0 1]);
ylabel (hcb,'P(s_{ipsi},s_{contra}|spike)/P(s_{ipsi},s_{contra})','interpreter','tex');
set(hcb,'yticklabel',round((str2num(get(hcb,'yticklabel')).^2)*100)/100); %#ok<ST2NM>

bayesian_model.IONL = IONL;
bayesian_model.ScaleFactor = ScaleFactor;
bayesian_model.Likelihood = L;
bayesian_model.Prior = P;
bayesian_model.sigmalim = sigmalim;
bayesian_model.X = X;
bayesian_model.Y = Y;

return;