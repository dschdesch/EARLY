function [] = plot_Binaural_Wiener(data)

colorordermx = [0 0 1;1 0 0;0 0.5 0;0 0.75 0.75;0.75 0 0.75;0.75 0.75 0;0.25 0.25 0.25];
set(0,'defaultAxesFontName','Arial','defaultTextFontName','Arial');
fh1 = figure;
set(fh1,'paperpositionmode','auto','units','normalized',...
    'position',[0.2 0.05 0.55 0.85],'DefaultAxesColorOrder',colorordermx,'color','none');
a = linspace(1,0,64)';
b = ones(64,1);
map  = [a a b; b a a; a a a];
colormap(map);
if ~isempty(data.kernels)
    ah = subplot(3,3,1);
    hold on;
    set(ah,'visible','off','xlim',[0 1],'ylim',[0 1]);
    %noise samples
    N1 = randn(1000,1)/3*0.05;
    N2 = randn(1000,1)/3*0.05;
    T = linspace(0,0.2,1000);
    plot(T,N1+0.45,'r-');
    plot(T,N2+0.65,'b-');
    drawArrow([0.2 0.45],[0.4 0.5],'k');
    drawArrow([0.2 0.65],[0.4 0.6],'k');
    outh = drawArrow([0.7 0.55],[0.9 0.55],'g');
    set(outh,'edgecolor','g');
    rectangle('position',[0.4 0.4 0.3 0.3],'curvature',[1 1],'facecolor',[.7 .7 .7],'linewidth',2)
    text(0.72,0.56,'spikes','fontsize',10,'horizontalalignment','left','verticalalignment','bottom');
    plot([0.1 0.1],[0.35 0.25],'--','linewidth',2,'color','b');
    plot([0.1 0.5],[0.25 0.25],'--','linewidth',2,'color','b');
    plot([0.09 0.09],[0.35 0.05],'--','linewidth',2,'color','r');
    plot([0.09 0.5],[0.05 0.05],'--','linewidth',2,'color','r');
    plot([0.85 0.85],[0.5 0.05],'--','linewidth',2,'color','g');
    plot([0.85 0.6],[0.05 0.05],'--','linewidth',2,'color','g');
    plot([0.85 0.6],[0.25 0.25],'--','linewidth',2,'color','g');
    rectangle('position',[0.5 0.2 0.1 0.1],'curvature',[1 1],'facecolor','w','linewidth',1);
    rectangle('position',[0.5 0 0.1 0.1],'curvature',[1 1],'facecolor','w','linewidth',1);
    plot(0.55,0.25,'kp','markerfacecolor','k','markersize',10);
    plot(0.55,0.05,'kp','markerfacecolor','k','markersize',10);
    text(0.1,0.75,'INPUT','verticalalignment','bottom','horizontalalignment','center');
    text(0.8,0.75,'OUTPUT','verticalalignment','bottom','horizontalalignment','center');
    axis(ah,'square');
    
    %revcors - time domain
    maxx = min(max(data.kernels.h1envif.endlatency)+0.5,25);
    minx = max(min(data.kernels.h1envif.frontlatency)-0.5,0);
    maxy = max(abs(data.kernels.h1wf(:)))*1.2;
    miny = -maxy;
    ah = subplot(3,3,2);
    plot(data.kernels.Time,[squeeze(data.kernels.nullh1wf(:,1,:)) ...
        squeeze(data.kernels.nullh1wf(:,2,:))],'-','color',[.6 .6 .6],'linewidth',1);hold on;
    plot(data.kernels.Time,data.kernels.h1wf,'-','linewidth',2);
    set(ah,'fontsize',12,'tickdir','out','linewidth',1,'layer','top','xlim',[minx maxx],'ylim',[miny maxy],...
        'ticklength',[0.03 0.03]);
    ylabel (ah,'spikes/s/Pa');
    xlabel (ah,'TIME (ms)');
    box(ah,'off');
    axis(ah,'square');
    text(0.9,0.9,'\bfB','units','normalized','fontsize',16,'interpreter','tex','horizontalalignment','center','verticalalignment','middle');
    text(0.98,0.2,'Ipsi.','units','normalized','fontsize',10,'horizontalalignment','right','verticalalignment','middle','color','b');
    text(0.98,0.1,'Contra.','units','normalized','fontsize',10,'horizontalalignment','right','verticalalignment','middle','color','r');
    text(0.05,1,sprintf('%2.0f spikes',data.kernels.TotalSpikes),'units','normalized','fontsize',10,'horizontalalignment','left','verticalalignment','top','color','k');
    
    %power spectra
    maxx = 4;
    minx = 0.04;
    maxy = 5;
    miny = -30;
    psah = subplot(3,3,3);
    hold on;
    scatter3(data.kernels.h1ffax',data.kernels.h1magmask(:,1),data.kernels.h1coeffs.Monaural.weights(:,1),...
        30*ones(size(data.kernels.h1magmask(:,1))),(data.kernels.h1coeffs.Monaural.weights(:,1)/3)-0.01,'fill','o','markeredgecolor','none');
    scatter3(data.kernels.h1ffax',data.kernels.h1magmask(:,2),data.kernels.h1coeffs.Monaural.weights(:,2),...
        30*ones(size(data.kernels.h1magmask(:,2))),(data.kernels.h1coeffs.Monaural.weights(:,2)/3)+(1/3)-0.01,'fill','o','markeredgecolor','none');
    set(psah,'fontsize',12,'tickdir','out','linewidth',1,'layer','top','ylim',[miny maxy],...
        'xlim',[minx maxx],'xtick',[0.1 1],'xticklabel',[0.1 1],'ticklength',[0.03 0.03],'clim',[0 1],'xscale','log');
    ylabel (psah,'GAIN (dB re max)');
    xlabel (psah,'FREQUENCY (kHz)');
    axis(psah,'square');box(psah,'off');
    text(0.9,0.9,'\bfC','units','normalized','fontsize',16,'interpreter','tex','horizontalalignment','center','verticalalignment','middle');
    text(0.98,0.2,'Ipsi.','units','normalized','fontsize',10,'horizontalalignment','right','verticalalignment','middle','color','b');
    text(0.98,0.1,'Contra.','units','normalized','fontsize',10,'horizontalalignment','right','verticalalignment','middle','color','r');
    %Put the difcor spectrum on top of these - need a new axis to display
    %properly - to do with handling of patch objects and line objects in matlabs eps
    %export
    psah2 = axes('position',get(psah,'position'),'color','none');
    Y = data.NDF.difcor.power;
    Y(Y<-20) = nan;
    plot(psah2,data.NDF.difcor.freq/1000,Y,'k-','linewidth',2);
    axis(psah2,'square');
    box(psah2,'off');
    set(psah2,'fontsize',12,'tickdir','out','linewidth',1,'layer','top','ylim',[miny maxy],...
        'xlim',[minx maxx],'xtick',[0.1 1],'xticklabel',[0.1 1],'ticklength',[0.03 0.03],'clim',[0 1],'xscale','log','outerposition',get(psah,'outerposition'),'position',get(psah,'position'),'visible','off');
    ylabel (psah2,'GAIN (dB re max)');
    xlabel (psah2,'FREQUENCY (kHz)');
    
    %unwrapped phase for each ear
    ah = subplot(3,3,4);
    scatter3(data.kernels.h1ffax',data.kernels.h1phase_uw(:,1),data.kernels.h1coeffs.Monaural.weights(:,1),...
        30*ones(size(data.kernels.h1phase_uw(:,1))),(data.kernels.h1coeffs.Monaural.weights(:,1)/3)-0.01,'fill','o','markeredgecolor','none');
    hold on;
    scatter3(data.kernels.h1ffax',data.kernels.h1phase_uw(:,2),data.kernels.h1coeffs.Monaural.weights(:,2),...
        30*ones(size(data.kernels.h1phase_uw(:,2))),(data.kernels.h1coeffs.Monaural.weights(:,2)/3)+(1/3)-0.01,'fill','o','markeredgecolor','none');
    set(ah,'fontsize',12,'tickdir','out','linewidth',1,'layer','top',...
        'ylim',[floor(nanmin(reshape(data.kernels.h1phase_uw(:,[1 2]),numel(data.kernels.h1phase_uw(:,[1 2])),1))) 1],...
        'xlim',[0 2],'ticklength',[0.03 0.03],'clim',[0 1]);
    ylabel (ah,'PHASE (cycles)');
    xlabel (ah,'FREQUENCY (kHz)');
    box(ah,'off');
    axis(ah,'square');
    grid(ah,'off');
    view(ah,[0,90]);
    colormap(map);
    text(0.9,0.9,'\bfD','units','normalized','fontsize',16,'interpreter','tex','horizontalalignment','center','verticalalignment','middle');
    text(0.98,0.2,'Ipsi.','units','normalized','fontsize',10,'horizontalalignment','right','verticalalignment','middle','color','b');
    text(0.98,0.1,'Contra.','units','normalized','fontsize',10,'horizontalalignment','right','verticalalignment','middle','color','r');
    text(0.02,0.2,sprintf('GD: %2.2f ms',data.kernels.h1coeffs.Monaural.Delay(1)),'units','normalized','fontsize',10,'horizontalalignment','left','verticalalignment','middle','color','b');
    text(0.02,0.1,sprintf('GD: %2.2f ms',data.kernels.h1coeffs.Monaural.Delay(2)),'units','normalized','fontsize',10,'horizontalalignment','left','verticalalignment','middle','color','r');
    
    %inter-aural phase
    ah = subplot(3,3,5);
    scatter3(data.kernels.h1ffax,data.kernels.h1phase_uw(:,3),data.kernels.h1coeffs.Binaural.weights,...
        30*ones(size(data.kernels.h1phase_uw(:,3))),(data.kernels.h1coeffs.Binaural.weights/3)+(2/3),'fill','o','markeredgecolor','none');
    set(ah,'fontsize',12,'tickdir','out','linewidth',1,'layer','top',...
        'ylim',[floor(nanmin(reshape(data.kernels.h1phase_uw(:,[1 2]),numel(data.kernels.h1phase_uw(:,[1 2])),1))) 1],...
        'xlim',[0 2],'ticklength',[0.03 0.03],'clim',[0 1],'ylim',[-0.5 0.5]);
    ylabel (ah,'PHASE (cycles)');
    xlabel (ah,'FREQUENCY (kHz)');
    box(ah,'off');
    axis(ah,'square');
    view(ah,[0,90]);
    grid(ah,'off');
    text(0.9,0.9,'\bfE','units','normalized','fontsize',16,'interpreter','tex','horizontalalignment','center','verticalalignment','middle');
    text(0.02,0.2,sprintf('CP: %2.2f cycles',data.kernels.h1coeffs.Binaural.CP),'units','normalized','fontsize',10,'interpreter','tex','horizontalalignment','left','verticalalignment','middle');
    text(0.02,0.1,sprintf('CD: %2.2f ms',data.kernels.h1coeffs.Binaural.CD),'units','normalized','fontsize',10,'interpreter','tex','horizontalalignment','left','verticalalignment','middle');
    
    %instantaneous frequency glides
    ah = subplot(3,3,6);
    Z = data.kernels.h1envif.z(2:end,:);
    Z = bsxfun(@rdivide,Z,max(Z));
    Z = max(0.05,Z);
    Y = data.kernels.h1envif.IF/1000;
    Y(Y<0)=nan;
    minx = min(data.kernels.h1envif.frontlatency);
    maxx = max(data.kernels.h1envif.endlatency);
    scatter3(data.kernels.Time(2:end)',Y(:,1),Z(:,1),...
        30*ones(size(Y(:,1))),(Z(:,1)/3)-0.01,'fill','o','markeredgecolor','none');
    hold on;
    scatter3(data.kernels.Time(2:end)',Y(:,2),Z(:,2),...
        30*ones(size(Y(:,2))),(Z(:,2)/3)+(1/3)-0.01,'fill','o','markeredgecolor','none');
    set(ah,'fontsize',12,'tickdir','out','linewidth',1,'layer','top','ylim',[0 2],...
        'ticklength',[0.03 0.03],'clim',[0 1],'xlim',[minx maxx]);
    ylabel (ah,'FREQUENCY (kHz)');
    xlabel (ah,'TIME (ms)');
    box(ah,'off');
    axis(ah,'square');
    grid(ah,'off');
    view(ah,[0,90]);
    text(0.9,0.9,'\bfF','units','normalized','fontsize',16,'interpreter','tex','horizontalalignment','center','verticalalignment','middle');
    text(0.98,0.2,'Ipsi.','units','normalized','fontsize',10,'horizontalalignment','right','verticalalignment','middle','color','b');
    text(0.98,0.1,'Contra.','units','normalized','fontsize',10,'horizontalalignment','right','verticalalignment','middle','color','r');
    
    %cross correlation of spectra to estimate disparity
    ah = subplot(3,3,7);
    hold on;
    plot([0 0],[0 1.2],'k-','linewidth',1);
    plot(data.kernels.tuning_diff.freq_oct,data.kernels.tuning_diff.normcipsi,'b-','linewidth',2);
    plot(data.kernels.tuning_diff.freq_oct,data.kernels.tuning_diff.normccontra,'r-','linewidth',2);
    set(ah,'fontsize',12,'tickdir','out','linewidth',1,'layer','top','ylim',[0 1.2],...
        'ticklength',[0.03 0.03],'xlim',[-1 1]);
    ylabel (ah,'CORRELATION');
    xlabel (ah,'\Delta FREQUENCY (oct.)','interpreter','tex');
    box(ah,'off');
    axis(ah,'square');
    grid(ah,'off');
    view(ah,[0,90]);
    text(0.9,0.9,'\bfG','units','normalized','fontsize',16,'interpreter','tex','horizontalalignment','center','verticalalignment','middle');
    text(0.98,0.2,'Ipsi.','units','normalized','fontsize',10,'horizontalalignment','right','verticalalignment','middle','color','b');
    text(0.98,0.1,'Contra.','units','normalized','fontsize',10,'horizontalalignment','right','verticalalignment','middle','color','r');
    plot([data.kernels.tuning_diff.ipsivsbin_oct data.kernels.tuning_diff.ipsivsbin_oct],[0 max(data.kernels.tuning_diff.normcipsi)],'b--','linewidth',1);
    plot([data.kernels.tuning_diff.contravsbin_oct data.kernels.tuning_diff.contravsbin_oct],[0 max(data.kernels.tuning_diff.normccontra)],'r--','linewidth',1);
    plot([data.kernels.tuning_diff.ipsivsbin_oct data.kernels.tuning_diff.contravsbin_oct],[0.1 0.1],'k-','linewidth',4);
    text(0.01,1,sprintf('Disparity: \n%2.2f oct\n%2.2f mm',data.kernels.tuning_diff.ipsivscontra_oct,data.kernels.tuning_diff.ipsivscontra_mm),...
        'units','normalized','fontsize',10,'horizontalalignment','left','verticalalignment','top','color','k');
    
    %positive and negative noise delay curves
    ah = subplot(3,3,8);
    hold on;
    maxy = max([data.NDF.positive.ratespline data.NDF.negative.ratespline])+10;
    miny = 0;
    plot([0 0],[miny maxy],'k--','linewidth',1);
    plot(data.NDF.xvals.ITDxspline,data.NDF.positive.ratespline,'-','color',[1 165/255 0],'linewidth',2);
    plot(data.NDF.xvals.ITDx,data.NDF.positive.rate,'o','markerfacecolor',[1 165/255 0],'markeredgecolor',[1 165/255 0],'markersize',5,'linewidth',2);
    plot(data.NDF.xvals.ITDxspline,data.NDF.negative.ratespline,'g-','linewidth',2);
    plot(data.NDF.xvals.ITDx,data.NDF.negative.rate,'go','markerfacecolor','g','markersize',5,'linewidth',2);
    set(ah,'fontsize',12,'tickdir','out','linewidth',1,'layer','top',...
        'ticklength',[0.03 0.03],'xlim',[min(data.NDF.xvals.ITDx) max(data.NDF.xvals.ITDx)],'ylim',[miny maxy]);
    ylabel (ah,'FIRING RATE (spikes/s)');
    xlabel (ah,'ITD (ms)');
    box(ah,'off');
    axis(ah,'square');
    grid(ah,'off');
    view(ah,[0,90]);
    text(0.9,0.9,'\bfH','units','normalized','fontsize',16,'interpreter','tex','horizontalalignment','center','verticalalignment','middle');
    text(0.02,0.9,'\rho = +1','units','normalized','fontsize',10,'horizontalalignment','left','verticalalignment','middle','color',[1 165/255 0],'interpreter','tex');
    text(0.02,0.8,'\rho = -1','units','normalized','fontsize',10,'horizontalalignment','left','verticalalignment','middle','color','g','interpreter','tex');
    
    
    %Difcor and prediction
    ah = subplot(3,3,9);
    hold on;
    maxy = max(abs(data.NDF.difcor.ratespline))+10;
    miny = -maxy;
    plot([0 0],[miny maxy],'k--','linewidth',1);
    plot(data.NDF.xvals.ITDxspline,data.NDF.difcor.ratespline,'k-','linewidth',2);
    plot(data.NDF.xvals.ITDx,data.NDF.difcor.rate,'ko','markerfacecolor','k','markersize',5,'linewidth',2);
    plot(data.NDF.fit.simple.ITDx,data.NDF.fit.simple.difcorrate,'m-','linewidth',2);
    set(ah,'fontsize',12,'tickdir','out','linewidth',1,'layer','top',...
        'ticklength',[0.03 0.03],'xlim',[min(data.NDF.xvals.ITDx) max(data.NDF.xvals.ITDx)],'ylim',[miny maxy]);
    ylabel (ah,'\Delta FIRING RATE (spikes/s)','interpreter','tex');
    xlabel (ah,'ITD (ms)');
    box(ah,'off');
    axis(ah,'square');
    grid(ah,'off');
    view(ah,[0,90]);
    text(0.9,0.9,'\bfI','units','normalized','fontsize',16,'interpreter','tex','horizontalalignment','center','verticalalignment','middle');
    text(0.02,0.9,sprintf('R^2 = %2.2f',data.NDF.fit.simple.COD),'units','normalized','fontsize',10,'horizontalalignment','left','verticalalignment','middle','color','m','interpreter','tex');

    %Make sure that the second axis in the power spectrum plot is correctly
    %alligned with the main axis - tends to jump around during calls to
    %subplot
    drawnow;
    set(psah2,'outerposition',(get(psah,'outerposition')),'position',(get(psah,'position')),'visible','off');
end