function [NDF,BM] = NDF_prediction(AnNum,contrachan,NDF,dt_s,kernels,BM)

%Returns fitted noise-delay functions based on the two impulse responses:
%achieves this in two different ways...

%Method #1: quick and simple:
%Take the h1 kernels expressed as a modulation factor of the mean firing
%rate. Cross correlate them. Then multiply by h0 (the mean firing rate), to
%get the predicted difcor.

%-------------------------------

%Method #2: more complicated:
%Based on the real noise delay functions scale the instantaneous spike
%probability function to give the best possible fit. Uses a constrained
%non-linear optimization approach.

%-------------------------------

%% Implement Method #1
signalRMS = 10^(kernels.SPL/20)*20e-6;
signalPOWER  = signalRMS^2;
ITD = NDF.xvals.ITDx;
ditd = abs(diff(ITD([1,2])));
N = length(ITD);
NDC = [NDF.negative.rate.mean'; NDF.positive.rate.mean'];
[c,lags] = xcorr(kernels.h1(:,2)*(sqrt(signalPOWER)/kernels.h0),kernels.h1(:,1)*(sqrt(signalPOWER)/kernels.h0));
xcorrITDs = lags*dt_s*1000;
c = c(xcorrITDs>=min(ITD)-ditd & xcorrITDs<=max(ITD)+ditd);
NDF.fit.simple.ITDx = xcorrITDs(xcorrITDs>=min(ITD)-ditd & xcorrITDs<=max(ITD)+ditd);
NDF.fit.simple.difcorrate = c*kernels.h0;
try
    [NDF.fit.simple.COD] = Coeff_Determination(NDF.difcor.rate.mean,interp1(NDF.fit.simple.ITDx,NDF.fit.simple.difcorrate,ITD));
catch
    [NDF.fit.simple.COD] = Coeff_Determination(NDF.difcor.rate.mean',interp1(NDF.fit.simple.ITDx,NDF.fit.simple.difcorrate,ITD));
end
try
    NDF.fit.simple.rho = corr(NDF.difcor.rate.mean,interp1(NDF.fit.simple.ITDx,NDF.fit.simple.difcorrate,ITD));
catch
    NDF.fit.simple.rho = corr(NDF.difcor.rate.mean,interp1(NDF.fit.simple.ITDx,NDF.fit.simple.difcorrate,ITD)');
end
figure;set(gcf,'paperpositionmode','auto');
try
    maxrate = max(abs([NDF.fit.simple.difcorrate;NDF.difcor.rate.spline.mean]))+10;
catch
    maxrate = max(abs([NDF.fit.simple.difcorrate;NDF.difcor.rate.spline.mean']))+10;
end
plot(NDF.xvals.ITDxspline,NDF.difcor.rate.spline.mean,'g-','linewidth',2);hold on
plot(ITD,NDF.difcor.rate.mean,'go','markersize',10,'linewidth',2,'markerfacecolor','w');hold on;
plot(NDF.fit.simple.ITDx,NDF.fit.simple.difcorrate,'k-','linewidth',2);
box off;axis square;
set(gca,'xlim',[min(ITD)-1 max(ITD)+1],'ylim',[-maxrate maxrate],'fontsize',14,...
    'linewidth',1,'ticklength',[0.02 0.02],'tickdir','out','layer','top');
ylabel ('\Delta FIRING RATE (spikes/s)','interpreter','tex');
xlabel 'ITD (ms)';
title 'DIFCOR: simple cross-correlation approach';
text(min(ITD)-0.5,-maxrate+10,sprintf('R^2 = %2.2f',NDF.fit.simple.COD),'interpreter','tex','fontsize',10);
pos = get(gca,'outerposition');
set(gca,'outerposition',min(max(pos,0),1));


%% Implement Method #2
SPL = kernels.SPL;
rho = [-1 1];
Nrho = length(rho);
noisedur = 1000;
%Predict noise delay curves
for i = 1:Nrho
    for itd = 1:N
        wv = MakeGaussNoiseBand(10, 20000, noisedur, 1/dt_s, 2, rho(i), ITD(itd), 0, 999,contrachan);
        wv = rescalewv(wv,SPL,'Pa');
        [dummy,dummy,wv] = flipchannels(AnNum,contrachan,wv);%Make sure the order is [ipsi contra]
        s(:,:,i,itd) = zeros(size(wv));
        for j = 1:2
            s(:,j,i,itd) = conv(wv(:,j),kernels.h1wf(:,j),'same');
            s(:,j,i,itd) = s(:,j,i,itd)/std(s(:,j,i,itd));
        end
    end
end
s = min(s,4);s = max(s,-4);%Make sure no very rare samples beyond 4*sigma mess things up.

%Set up the optimization
options = optimset('Algorithm','interior-point');

p0 = [1 0.00000001 0.1];
[BM.P,NDF.fit.MSE] = fmincon(@objective,p0,[],[],[],[],[0.0005 -0.5 1e-10],[3 0.5 1],[],options);

[BM.IOprob,NDF.fit.rate] = model(BM.P);
NDF.fit.difcorrate = NDF.fit.rate(2,:)-NDF.fit.rate(1,:);

[NDF.fit.COD.positive] = Coeff_Determination(NDF.positive.rate.mean,NDF.fit.rate(2,:)');
[NDF.fit.COD.negative] = Coeff_Determination(NDF.negative.rate.mean,NDF.fit.rate(1,:)');
[NDF.fit.COD.difcor] = Coeff_Determination(NDF.difcor.rate.mean,NDF.fit.difcorrate');
[NDF.fit.COD.total] = Coeff_Determination([NDF.negative.rate.mean;NDF.positive.rate.mean], [NDF.fit.rate(1,:) NDF.fit.rate(2,:)]');

[NDF.fit.rho.positive] = corr(NDF.positive.rate.mean,NDF.fit.rate(2,:)');
[NDF.fit.rho.negative] = corr(NDF.negative.rate.mean,NDF.fit.rate(1,:)');
[NDF.fit.rho.difcor] = corr(NDF.difcor.rate.mean,NDF.fit.difcorrate');
[NDF.fit.rho.total] = corr([NDF.negative.rate.mean; NDF.positive.rate.mean], [NDF.fit.rate(1,:) NDF.fit.rate(2,:)]');

%display the real probability function
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.25 0.12 0.45 0.7]);
surf(BM.X,BM.Y,BM.IOprob);
shading flat;
colormap (flipud(hot));
set(gca,'xlim',[-5 5],'ylim',[-5 5],'fontsize',14,...
    'linewidth',1,'tickdir','out','ticklength',[0.02 0.02],'layer','top',...
    'xtick',-4:2:4,'ytick',-4:2:4,'zgrid','off');
axis square;view(0,90);
ylabel ('S_{contra}','interpreter','tex');
xlabel ('S_{ipsi}','interpreter','tex');
title 'Instantaneous Firing Probability';
hcb = colorbar;
set(hcb,'fontsize',14,'linewidth',1,'ytick',linspace(0,round(max(BM.IOprob(:))*100)/100,5),'ylim',[0 max(BM.IOprob(:))]);
ylabel (hcb,'FIRING PROBABILITY','interpreter','tex');

%display the NDF estimates
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.25 0.12 0.45 0.7]);
maxrate = max(abs([NDF.positive.rate.mean; NDF.negative.rate.mean]))+10;
subplot(2,2,1);
plot(NDF.xvals.ITDxspline,NDF.positive.rate.spline.mean,'c-','linewidth',2);hold on
plot(ITD,NDF.positive.rate.mean,'co','markersize',10,'linewidth',2,'markerfacecolor','w');hold on;
plot(ITD,NDF.fit.rate(2,:),'k-','linewidth',2);
box off;axis square;
set(gca,'xlim',[min(ITD)-1 max(ITD)+1],'ylim',[-10 maxrate],'fontsize',14,...
    'linewidth',1,'ticklength',[0.02 0.02],'tickdir','out','layer','top');
ylabel 'FIRING RATE (spikes/s)';
xlabel 'ITD (ms)';
title ('\rho = +1','interpreter','tex');
text(min(ITD)-0.5,0,sprintf('R^2 = %2.2f',NDF.fit.COD.positive),'interpreter','tex','fontsize',10);

subplot(2,2,2);
plot(NDF.xvals.ITDxspline,NDF.negative.rate.spline.mean,'m-','linewidth',2);hold on
plot(ITD,NDF.negative.rate.mean,'mo','markersize',10,'linewidth',2,'markerfacecolor','w');hold on;
plot(ITD,NDF.fit.rate(1,:),'k-','linewidth',2);
box off;axis square;
set(gca,'xlim',[min(ITD)-1 max(ITD)+1],'ylim',[-10 maxrate],'fontsize',14,...
    'linewidth',1,'ticklength',[0.02 0.02],'tickdir','out','layer','top');
ylabel 'FIRING RATE (spikes/s)';
xlabel 'ITD (ms)';
title ('\rho = -1','interpreter','tex');
text(min(ITD)-0.5,0,sprintf('R^2 = %2.2f',NDF.fit.COD.negative),'interpreter','tex','fontsize',10);

subplot(2,2,3);
plot(NDF.xvals.ITDxspline,NDF.difcor.rate.spline.mean,'g-','linewidth',2);hold on
plot(ITD,NDF.difcor.rate.mean,'go','markersize',10,'linewidth',2,'markerfacecolor','w');hold on;
plot(ITD,NDF.fit.difcorrate,'k-','linewidth',2);
box off;axis square;
set(gca,'xlim',[min(ITD)-1 max(ITD)+1],'ylim',[-maxrate maxrate],'fontsize',14,...
    'linewidth',1,'ticklength',[0.02 0.02],'tickdir','out','layer','top');
ylabel ('\Delta FIRING RATE (spikes/s)','interpreter','tex');
xlabel 'ITD (ms)';
title 'DIFCOR';
text(min(ITD)-0.5,-maxrate+10,sprintf('R^2 = %2.2f',NDF.fit.COD.difcor),'interpreter','tex','fontsize',10);

subplot(2,2,4);
plot([NDF.negative.rate.mean; NDF.positive.rate.mean]', [NDF.fit.rate(1,:) NDF.fit.rate(2,:)],'ko','linewidth',2,'markersize',10,'markerfacecolor','w');hold on;
plot([-10 maxrate],[-10 maxrate],'--','linewidth',2,'color',[.5 .5 .5]);
box off;axis square;
set(gca,'xlim',[-10 maxrate],'ylim',[-10 maxrate],'fontsize',14,...
    'linewidth',1,'ticklength',[0.02 0.02],'tickdir','out','layer','top');
ylabel 'PREDICTED FIRING RATE (spikes/s)';
xlabel 'MEASURED FIRING RATE (spikes/s)';
title 'MODEL vs. DATA';
text(maxrate/2,-5,sprintf('R^2 = %2.2f',NDF.fit.COD.total),'interpreter','tex','fontsize',10);
    
    function d = objective(p)
        [dummy,rate_hat] = model(p);
        r = NDC-rate_hat;
        d = mean(r(:).^2);
        
        %Show an output so you can view the minimization in progress
        figure(111);
        cla;
        plot(ITD,NDC','o','markersize',10);hold on;
        ylabel 'Spikes/s';
        xlabel 'ITD (ms)';
        if ~exist('ploth2','var')
            ploth2 = plot(ITD,rate_hat');hold on;
        else
            set(ploth2(1),'YData',rate_hat(1,:));
            set(ploth2(2),'YData',rate_hat(2,:));
        end
    end

    function [IOprob,meanspikerate] = model(p)
        IOprob = ((BM.IONL.^p(1))*p(3))+p(2);
        IOprob = max(0,IOprob);%Don't allow negative probabilities
        meanspikerate = zeros(2,N);
        for ii = 1:2
            for jj = 1:N
                spikeprob = interp2(BM.X,BM.Y,IOprob,s(:,1,ii,jj),s(:,2,ii,jj));
                meanspikerate(ii,jj) = sum(spikeprob)/(noisedur/1000);
            end
        end
    end
end