function [] = plot_Monaural_Wiener(ParamsIn)

%Plot a raster, rate level function and synchrony level function


if ~isempty(ParamsIn.RLV)
    figure;set(gcf,'units','normalized','position',[0.2069    0.2189    0.3019    0.4533]);
    axes('ActivePositionProperty','position','fontsize',16);
    [ax_h,rate_h,sync_h] = plotyy(ParamsIn.RLV.Level,ParamsIn.RLV.Rate,ParamsIn.RLV.Level,ParamsIn.RLV.SyncMag);
    set(ax_h(1),'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'xlim',[min(RLV.Level)-5 max(RLV.Level)+5],'ycolor','r');
    set(ax_h(2),'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'xlim',[min(RLV.Level)-5 max(RLV.Level)+5],'ylim',[0 1],'ycolor','b','ytick',[0:0.2:1]);
    set(rate_h,'color','r','markersize',10,'marker','o','markerfacecolor','r','markeredgecolor','r','linewidth',1);
    set(sync_h,'color','b','markersize',10,'marker','o','markerfacecolor','none','markeredgecolor','b','linewidth',1);
    hold(ax_h(2),'on');
    plot(ax_h(2),RLV.Level(RLV.SyncPval<0.01),RLV.SyncMag(RLV.SyncPval<0.01),'bo','markersize',10,'markerfacecolor','b');hold on;
    box(ax_h(1),'off');
    box(ax_h(2),'off');
    axis(ax_h(1),'square');
    axis(ax_h(2),'square');
    xlabel (ax_h(1),'SOUND LEVEL (dB SPL)');
    ylabel (ax_h(1),'FIRING RATE (spikes s^{-1})','interpreter','tex');
    ylabel (ax_h(2),'VECTOR STRENGTH','interpreter','tex');
    
    %plot a quick raster
    figure;set(gcf,'units','normalized','position',[0.1406    0.1078    0.3300    0.7833]);
    axes('ActivePositionProperty','position');
    nrspls = length(RLV.Level);
    nrreps = size(RLV.SpikeTimes,2);
    linecols = {'b','r'};
    
    
    [rlvspls,splindx] = sort(RLV.Level);
    count = 1;
    for i = 1:nrspls
        for j = 1:nrreps
            plot([RLV.SpikeTimes{splindx(i),j}; RLV.SpikeTimes{splindx(i),j}],...
                [zeros(1,length(RLV.SpikeTimes{splindx(i),j}))+j+((i-1)*nrreps)-0.75; ...
                zeros(1,length(RLV.SpikeTimes{splindx(i),j}))+j+((i-1)*nrreps)-0.25],...
                '-','color',linecols{count});hold on;
        end
        count = count+1;
        if count>2
            count = 1;
        end
    end
    box off;
    set(gca,'fontsize',16,'xlim',[0 100],'ytick',nrreps/2:nrreps:(nrreps*nrspls)-(nrreps/2),'yticklabel',rlvspls);
    xlabel 'TIME (ms)';
    ylabel 'SOUND LEVEL (dB SPL)';
end
if ~isempty(ParamsIn.TC)
    figure;set(gcf,'units','normalized','position',[0.1631    0.1600    0.3469    0.5322]);
    semilogx(ParamsIn.TC.thr.freq/1000,ParamsIn.TC.thr.thr,'b+','markersize',8,'linewidth',1);hold on
    semilogx(ParamsIn.TC.fit.freq/1000,ParamsIn.TC.fit.thr,'k-','linewidth',2);
    plot(ParamsIn.TC.fit.cf/1000,ParamsIn.TC.fit.minthr,'ro','markersize',10,'linewidth',2);
    set(gca,'ylim',[-20 70],'xlim',[0.05 40],'fontsize',16,'linewidth',1,'layer','top','ytick',-20:20:60,'xtick',[0.1 1 10],'xticklabel',[0.1 1 10]);
    xlabel 'FREQUENCY (kHz)';
    ylabel 'THRESHOLD (dB SPL)';
    box off;
    %Display some text
    text(6,-10,sprintf('SR: %2.1f spikes s^{-1} \nCF: %2.2f kHz \nTh: %2.1f dB SPL \nQ_{10}: %2.1f',...
        ParamsIn.TC.thr.sr,ParamsIn.TC.fit.cf/1000,ParamsIn.TC.fit.minthr,ParamsIn.TC.fit.q10),'interpreter','tex','color','r');
end


%% REVCORs
if ParamsIn.TC.fit.cf<500
    maxtime = 30;
else
    maxtime = 10;
end
% for i = 1:nrspls
    figure;set(gcf,'units','normalized','position',[0.0300 0.1122 0.9300 0.334],...
        'paperpositionmode','auto');
    subplot (1,4,1);
    
    maxy = max(max(abs(ParamsIn.h1.h1)));
    tdnullh = plot(ParamsIn.h1.Time,ParamsIn.h1.nullh1filt,'color',[0.6 0.6 0.6]);hold on;
    tdh = plot(ParamsIn.h1.Time,ParamsIn.h1.h1filt,'b-','linewidth',2);
    set(gca,'fontsize',16,'layer','top','linewidth',1,'ylim',[-1.1*maxy 1.1*maxy],...
        'ytick',[-maxy maxy],'yticklabel',[-1 1],'ticklength',[0.03 0.03],'xlim',[0 maxtime]);
    xlabel 'TIME (ms)';
    ylabel 'A.U.';
    box off;
    %Display some text for information about spike rate etc.
    text(maxtime/2,0.9*maxy,sprintf('Spike Count: %1.0f \nMean Rate: %1.1f s^{-1} \n%2.0f dB SPL',ParamsIn.TotalSpikes,ParamsIn.MeanSpikeRate,ParamsIn.h1.noiseSPL),...
        'fontsize',12,'color','k','interpreter','tex');
    
    %
    %% Plot the magnitude spectra
    subplot(1,4,2);
    
    fxs = max(h1magmask(:,i));
    fdmagnullh = semilogx(ffax,nullh1mag{i}-fxs,'color',[0.6 0.6 0.6],'linewidth',1);hold on;
    fdmagh(i) = semilogx(ffax,h1magmask(:,i)-fxs,'b','linewidth',4);
    set(gca,'xlim',[0.05 10],'xtick',[0.1 1 10],'xticklabel',[0.1 1 10],'linewidth',1,'layer','top',...
        'ticklength',[0.03 0.03],'fontsize',16,'ylim',[-30 5]);
    box off;
    xlabel 'FREQUENCY (kHz)';
    ylabel 'GAIN (dB)';
    %Add some text describing the tuning
    text(0.06,5,sprintf('DF: %2.2f kHz, Q_{3dB}: %2.1f, Q_{6dB}: %2.1f',domfreq(i),qvals(1,i),qvals(2,i)),...
        'fontsize',12,'color','b','horizontalalignment','left','verticalalignment','top','interpreter','tex');
    %
    %% Plot the phase spectra
    ah = subplot(1,4,3);
    %all data
    hold on;
    
    
    if ~all(isnan(h1phasemask_uw(:,i)))
        %Phase for the spectral components above the noise floor
        plot(ffax,h1phasemask_uw(:,i)./(2*pi),'b','linewidth',3);
        plot(ffax(ffax<=4),feval(lin_beta{i},ffax(ffax<=4)),'b:','linewidth',1);
        set(gca,'xlim',[0 3],'xtick',[0:1:3],'xticklabel',[0:1:3],'linewidth',1,'layer','top',...
            'ticklength',[0.03 0.03],'fontsize',16,'ylim',[-12 2]);
        box off;
        xlabel 'FREQUENCY (kHz)';
        ylabel 'PHASE (Cycles)';
        %Display some text with the group-delay estimates for both ears
        text(0.07,-10,sprintf('GD: %2.3f ms (%2.3f, %2.3f)',-coeffs(i,1),-confints{i}(2,1),-confints{i}(1,1)),...
            'interpreter','tex','color','b','fontsize',12,'horizontalalignment','left');
    end
    
    
    
    
    ah = subplot(1,4,4);
    %define two separate colormaps
    IFcolormap = [linspace(1,0,64); linspace(1,0,64); ones(1,64)]';
    colormap(IFcolormap);
    sh = surface([IFTime(:), IFTime(:)], [h1if(:,i)/1000, h1if(:,i)/1000], [zcrit_h1(2:end,i), zcrit_h1(2:end,i)],...
        'EdgeColor','none', 'MarkerFaceColor','flat','marker','o','markersize',3,'linewidth',1,'linestyle','none','facecolor','none');
    view (0,90);
    grid off;box off;hold on;
    set(ah,'yscale','linear','ylim',[0 4],...
        'layer','top','linewidth',1,'ticklength',[0.03 0.03],'fontsize',16,'xlim',[0 10]);
    
    ylabel(ah,'I.F. (kHz)');
    xlabel(ah,'TIME (ms)');
    %
    
    %
    %% Add the envelope to the h1 time-domain plot
    subplot(1,4,1);
    plot(Time,h1env(:,i),'b:','linewidth',1);
    plot(Time(maxenvind(i)),maxenv(i),'b*','markersize',10,'linewidth',1);
    plot(Time(frontlatencyind(i)), h1wfilt(frontlatencyind(i),i),'b+','markersize',10,'linewidth',1);
    text(maxtime/2,-0.7*maxy,sprintf('Latency: %2.1f ms',Latency(i)),'fontsize',12,'color','b');
% end

if nrspls>1
    if plotYN
        %         cols = {'b','g','r','c','m','k','y'};
        cols = {[0 0 1],[0 1 0],[1 0 0],[0 1 1],[1 0 1],[0 0 0],[1 1 0]};
        figure;set(gcf,'units','normalized','position',[0.0300 0.1122 0.9300 0.668],...
            'paperpositionmode','auto','renderer','zbuffer');
        %Plot all time-domain revcors
        subplot (3,4,1);
        count = 1;
        for i =1:nrspls
            xdata = get(tdh(i),'xdata');
            ydata = get(tdh(i),'ydata');
            ydata = ydata./max(abs(ydata));
            plot(xdata,ydata+(i-1),'color',cols{count},'linewidth',2);hold on;
            count = count+1;
            if count>7
                count = 1;
            end
        end
        set(gca,'layer','top','linewidth',1,'fontsize',16,'ticklength',[0.02 0.02],'ytick',0:nrspls-1,'yticklabel',uniquespls,'ylim',[-1.1 nrspls+0.1],'xlim',[0 30]);
        box off;
        xlabel 'TIME (ms)';
        ylabel 'LEVEL (dB SPL)';
        
        %Plot all frequency-domain revors
        subplot(3,4,2);
        count = 1;
        for i =1:nrspls
            xdata = get(fdmagh(i),'xdata');
            ydata = get(fdmagh(i),'ydata');
            semilogx(xdata,ydata-max(ydata),'color',cols{count},'linewidth',4);hold on;
            count = count+1;
            if count>7
                count = 1;
            end
        end
        set(gca,'layer','top','linewidth',1,'fontsize',16,'ticklength',[0.02 0.02],'ylim',[-30 5],'xlim',[0.05 4],'xtick',[0.1 1 10],'xticklabel',[0.1 1 10]);
        box off;
        xlabel 'FREQ. (kHz)';
        ylabel 'GAIN (dB)';
        
        %Plot the unwrapped phase frequency as function of SPL
        %get rid of any whole cycle shifts - align everything to the 70-dB
        %SPL condition
        refcond = find(uniquespls==70);
        if ~isempty(refcond)
            for i = 1:nrspls
                corrphase(:,i) = (ParamsIn.h1phasemask_uw(:,i)./(2*pi));
            end
            for i = 1:nrspls
                if abs(nanmean(corrphase(:,i)-corrphase(:,refcond)))>0.5
                    corrphase(:,i) = corrphase(:,i)-round(nanmean(corrphase(:,i)-corrphase(:,refcond)));
                    corrphase(:,i) = Trifilter(corrphase(:,i)',11)';
                end
            end
            
            subplot(3,4,3);
            count = 1;
            for i = 1:nrspls
                plot(ParamsIn.h1ffax,corrphase(:,i),'color',cols{count});hold on;
                count = count+1;
                if count>7
                    count =1;
                end
            end
            box off;
            set(gca,'layer','top','linewidth',1,'fontsize',16,'ticklength',[0.02 0.02],'xlim',[0 2]);
            xlabel 'FREQ. (kHz)';
            ylabel 'PHASE (cycles)';
            
            
            subplot(3,4,11);
            count = 1;
            for i = 1:nrspls
                plot(ParamsIn.h1ffax(2:end),Trifilter((-diff(corrphase(:,i))/(ParamsIn.h1ffax(2)-ParamsIn.h1ffax(1)))',11),'color',cols{count});hold on;
                count = count+1;
                if count>7
                    count =1;
                end
            end
            box off;
            set(gca,'layer','top','linewidth',1,'fontsize',16,'ticklength',[0.02 0.02],'xlim',[0 2]);
            xlabel 'FREQ. (kHz)';
            ylabel 'G.D. (ms)';
            
            %subtract the linear fit to the 70-dB SPL condition from each of
            %the functions
            for i = 1:nrspls
                corrphase(:,i) = corrphase(:,i)-(ParamsIn.h1phasefit_coeffs(refcond,2)+(ParamsIn.h1ffax.*ParamsIn.h1phasefit_coeffs(refcond,1)))';%corrected phase
                corrphase(:,i) = Trifilter(corrphase(:,i)',11)';
            end
            subplot(3,4,7);
            count = 1;
            for i = 1:nrspls
                plot(ParamsIn.h1ffax,corrphase(:,i),'color',cols{count});hold on;
                count = count+1;
                if count>7
                    count =1;
                end
            end
            box off;
            set(gca,'layer','top','linewidth',1,'fontsize',16,'ticklength',[0.02 0.02],'xlim',[0 2]);
            xlabel 'FREQ. (kHz)';
            ylabel 'PHASE (cycles)';
            
            %Plot group delay as a function of SPL
            subplot(3,4,5);
            count = 1;
            for i = 1:nrspls
                if ~isempty(ParamsIn.h1phasefit_coeffs)
                    plot(uniquespls(i),-ParamsIn.h1phasefit_coeffs(i,1),'markeredgecolor',cols{count},'markerfacecolor',cols{count},'marker','o','linestyle','none','markersize',10);hold on;
                end
                count = count+1;
                if count>7
                    count =1;
                end
            end
            box off;
            set(gca,'layer','top','linewidth',1,'fontsize',16,'ticklength',[0.02 0.02],'xlim',[min(allspls)-10 max(allspls)+10]);
            xlabel 'LEVEL (dB SPL)';
            ylabel 'GD (ms)';
        end
        
        %Plot signal front delay as a function of SPL
        subplot(3,4,9);
        count = 1;
        for i = 1:nrspls
            if ~isempty(ParamsIn.h1phasefit_coeffs)
                plot(uniquespls(i),frontlatency(i),'markeredgecolor',cols{count},'markerfacecolor',cols{count},'marker','+','markersize',10);hold on;
            end
            count = count+1;
            if count>7
                count =1;
            end
        end
        box off;
        set(gca,'layer','top','linewidth',1,'fontsize',16,'ticklength',[0.02 0.02],'xlim',[min(allspls)-10 max(allspls)+10]);
        xlabel 'LEVEL (dB SPL)';
        ylabel 'SFD (ms)';
        
        %Plot dominant frequency and Q_6dB as a function of SPL
        subplot(3,4,6);
        count = 1;
        for i = 1:nrspls
            plot(uniquespls(i),ParamsIn.domfreq(i),'markeredgecolor',cols{count},'markerfacecolor',cols{count},'marker','o','linestyle','none','markersize',10);hold on;
            count = count+1;
            if count>7
                count =1;
            end
        end
        box off;
        set(gca,'layer','top','linewidth',1,'fontsize',16,'ticklength',[0.02 0.02],'xlim',[min(allspls)-10 max(allspls)+10],'ylim',[0 3]);
        xlabel 'LEVEL (dB SPL)';
        ylabel 'DOM. FREQ. (kHz)';
        
        %Plot the Q values
        subplot(3,4,10);
        count = 1;
        for i = 1:nrspls
            plot(uniquespls(i),ParamsIn.qvals(2,i),'markeredgecolor',cols{count},'markerfacecolor',cols{count},'marker','o','linestyle','none','markersize',10);hold on;
            count = count+1;
            if count>7
                count =1;
            end
        end
        box off;
        set(gca,'layer','top','linewidth',1,'fontsize',16,'ticklength',[0.02 0.02],'xlim',[min(allspls)-10 max(allspls)+10],'ylim',[0 3]);
        xlabel 'LEVEL (dB SPL)';
        ylabel ('Q_{6dB}','interpreter','tex');
        
        %Plot the instantaneous frequency at each SPL
        subplot(3,4,4);
        count = 1;
        for i = 1:nrspls
            C = bsxfun(@times,1-cols{count},1-(zcrit_h1(2:end,i)./max(zcrit_h1(2:end,i))));
            C(C==0) = 1;
            scatter3(IFTime', h1if(:,i)/1000, zcrit_h1(2:end,i),10*ones(size(IFTime')),C,'filled');hold on;
            count = count+1;
            if count>7
                count =1;
            end
        end
        box off;view (0,90);grid off;
        set(gca,'layer','top','linewidth',1,'fontsize',16,'ticklength',[0.02 0.02],'xlim',[0 12],'ylim',[0 3]);
        xlabel 'TIME (ms)';
        ylabel 'I.F. (kHz)';
        
        subplot(3,4,8);
        count = 1;
        for i =1:nrspls
            C = bsxfun(@times,1-cols{count},1-(zcrit_h1(3:end,i)./max(zcrit_h1(3:end,i))));
            C(C==0) = 1;
            c_kHzms(:,i) = diff(h1if(:,i)/1000)/(IFTime(2)-IFTime(1));
            scatter3(IFTime(2:end)', c_kHzms(:,i), zcrit_h1(3:end,i),10*ones(size(IFTime(2:end)')),C,'filled');hold on;
            count = count+1;
            if count>7
                count =1;
            end
        end
        box off;view (0,90);grid off;
        set(gca,'layer','top','linewidth',1,'fontsize',16,'ticklength',[0.02 0.02],'xlim',[0 12],'ylim',[-10 10]);
        xlabel 'TIME (ms)';
        ylabel 'Glide Slope (kHz/ms)';
    end
end