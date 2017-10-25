function [] = PlotAllBinRev

%Save the figures to this location
binplotdir = 'C:\work\BinRev\BinRev\pdfs';
datadir = 'C:\work\BinRev\BinRev';
cd (datadir);

files = dir('*.mat');
nrfiles = length(files);
for i = 1:nrfiles
    data = load(fullfile(datadir,files(i).name));
    plotBFSresults;
    plotNOISEresults;
end
    function plotBFSresults
        %Still to add this
    end
    function plotNOISEresults
        if isfield(data.ParamsOut,'NDF')
            fh1 = figure(1);set(gcf,'units','normalized','position',[0.0275    0.1078    0.9300    0.7944],...
                'paperorientation','landscape','paperunits','normalized','paperposition',[0 0 1 1]);
            %Plot the revcors
            if isfield(data.ParamsOut,'Revcor')
                subplot(2,4,1);
                plot(data.ParamsOut.Revcor.Time,data.ParamsOut.Revcor.nullh1filt./max(abs(data.ParamsOut.Revcor.h1filt(:))),'color',[0.6 0.6 0.6]);hold on
                rh1 = plot(data.ParamsOut.Revcor.Time,data.ParamsOut.Revcor.h1filt(:,1)./max(abs(data.ParamsOut.Revcor.h1filt(:))),'b-','linewidth',2);hold on
                rh2 = plot(data.ParamsOut.Revcor.Time,data.ParamsOut.Revcor.h1filt(:,2)./max(abs(data.ParamsOut.Revcor.h1filt(:))),'r-','linewidth',2);hold on
                plot(data.ParamsOut.Revcor.Time,data.ParamsOut.Revcor.h1filtenv(:,1)./max(abs(data.ParamsOut.Revcor.h1filt(:))),'b:','linewidth',0.5);hold on
                plot(data.ParamsOut.Revcor.Time,data.ParamsOut.Revcor.h1filtenv(:,2)./max(abs(data.ParamsOut.Revcor.h1filt(:))),'r:','linewidth',0.5);hold on
                box off;
                set(gca,'xlim',[0 20],'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'ylim',[-1.2 1.2]);
                ylabel 'NORM. AMPLITUDE';
                xlabel 'TIME (ms)';
                lh1 = legend([rh1 rh2],'IPSI','CONTRA');
                set(lh1,'fontsize',8,'box','off');
                title (files(i).name(1:end-11));
                
                subplot(2,4,5);
                reflevel = mean(mean(data.ParamsOut.Revcor.nullh1filtmag((data.ParamsOut.Revcor.h1ffax)<3,:),2));
                semilogx(data.ParamsOut.Revcor.h1ffax,data.ParamsOut.Revcor.nullh1filtmag-reflevel,'color',[.6 .6 .6]);hold on;
                semilogx(data.ParamsOut.Revcor.h1ffax,data.ParamsOut.Revcor.h1filtmag(:,1)-reflevel,'b-');hold on;
                semilogx(data.ParamsOut.Revcor.h1ffax,data.ParamsOut.Revcor.h1filtmag(:,2)-reflevel,'r-');
                semilogx(data.ParamsOut.Revcor.h1ffax(data.ParamsOut.Revcor.h1zscore(:,1)>=5),data.ParamsOut.Revcor.h1filtmag(data.ParamsOut.Revcor.h1zscore(:,1)>=5,1)-reflevel,'b-','linewidth',4);
                semilogx(data.ParamsOut.Revcor.h1ffax(data.ParamsOut.Revcor.h1zscore(:,2)>=5),data.ParamsOut.Revcor.h1filtmag(data.ParamsOut.Revcor.h1zscore(:,2)>=5,2)-reflevel,'r-','linewidth',4);
                set(gca,'ylim',[-10 50],'fontsize',16,'xlim',[0.05 10],'xtick',[0.1 1 10],'xticklabel',[0.1 1 10],'layer','top','ticklength',[0.02 0.02]);
                box off;
                ylabel 'MAGNITUDE (dB)';
                xlabel 'FREQUENCY (kHz)';
                
                subplot(2,4,6);
                plot(data.ParamsOut.Revcor.h1ffax,data.ParamsOut.Revcor.h1phasemask_uw(:,1)./(2*pi),'b-','linewidth',2);hold on;
                plot(data.ParamsOut.Revcor.h1ffax,data.ParamsOut.Revcor.h1phasemask_uw(:,2)./(2*pi),'r-','linewidth',2);hold on;
                plot(data.ParamsOut.Revcor.h1ffax,data.ParamsOut.Revcor.h1phasefit_coeffs(1,2)+...
                    (data.ParamsOut.Revcor.h1phasefit_coeffs(1,1)*data.ParamsOut.Revcor.h1ffax),'b:','linewidth',0.5);
                plot(data.ParamsOut.Revcor.h1ffax,data.ParamsOut.Revcor.h1phasefit_coeffs(2,2)+...
                    (data.ParamsOut.Revcor.h1phasefit_coeffs(2,1)*data.ParamsOut.Revcor.h1ffax),'r:','linewidth',0.5);
                set(gca,'xlim',[0 2],'ylim',[-10 2],'linewidth',1,'layer','top','ticklength',[0.02 0.02],'fontsize',16);
                xlabel 'FREQUENCY (kHz)';
                ylabel 'PHASE (cycles)';
                box off;
                
                subplot(2,4,7);
                plot(data.ParamsOut.Revcor.h1magxcordfreqs, data.ParamsOut.Revcor.h1magxcorcor,'k-','linewidth',2);hold on;
                plot([0 0],[0 1],'k--','linewidth',0.5);
                set(gca,'linewidth',1,'fontsize',16,'xlim',[-1 1],'layer','top','ticklength',[0.02 0.02]);
                xlabel ('\Delta FREQUENCY (kHz)','interpreter','tex');
                ylabel 'CORRELATION';
                box off;
                
                
                subplot(2,4,8);
                weights = 1-[data.ParamsOut.Revcor.h1filtenvzscore(:,1)/max(data.ParamsOut.Revcor.h1filtenvzscore(:,1)) ...
                    data.ParamsOut.Revcor.h1filtenvzscore(:,2)/max(data.ParamsOut.Revcor.h1filtenvzscore(:,2))];
                
                for jj = 1:length(data.ParamsOut.Revcor.h1filtif(:,1))
                    plot(data.ParamsOut.Revcor.IFTime(jj),data.ParamsOut.Revcor.h1filtif(jj,1)/1000,'o','markersize',4,'markerfacecolor',[weights(jj,1) weights(jj,1) 1],'markeredgecolor',[weights(jj,1) weights(jj,1) 1]);hold on;
                end
                for jj = 1:length(data.ParamsOut.Revcor.h1filtif(:,2))
                    plot(data.ParamsOut.Revcor.IFTime(jj),data.ParamsOut.Revcor.h1filtif(jj,2)/1000,'o','markersize',4,'markerfacecolor',[1 weights(jj,2) weights(jj,2)],'markeredgecolor',[1 weights(jj,2) weights(jj,2)]);hold on;
                end

                set(gca,'linewidth',1,'fontsize',16,'layer','top','ticklength',[0.02 0.02],'xlim',[3 15],'ylim',[0 1.8]);                
                box off;
                xlabel 'TIME (ms)';
                ylabel 'INSTANT. FREQ. (kHz)';
                
            end
            %Plot the noise delay curves
            if isfield(data.ParamsOut,'NDF')
                subplot(2,4,2);
                plot(data.ParamsOut.NDF.NTDx,data.ParamsOut.NDF.NTDyPos,'c*','markersize',13,'linewidth',2);hold on;
                ndhp = plot(data.ParamsOut.NDF.NTDx_spline,data.ParamsOut.NDF.NTDyPos_spline,'c-','linewidth',2);hold on;
                plot(data.ParamsOut.NDF.NTDx,data.ParamsOut.NDF.NTDyNeg,'m*','markersize',13,'linewidth',2);
                ndhn = plot(data.ParamsOut.NDF.NTDx_spline,data.ParamsOut.NDF.NTDyNeg_spline,'m-','linewidth',2);hold on;
                set(gca,'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'layer','top','xlim',[min(data.ParamsOut.NDF.NTDx) max(data.ParamsOut.NDF.NTDx)]);
                box off;
                xlabel 'ITD (ms)';
                ylabel ('FIRING RATE (spikes s^{-1})','interpreter','tex');
                
                
                subplot(2,4,3);
                plot(data.ParamsOut.NDF.NTDx,data.ParamsOut.NDF.Difcor_NTD,'g*','markersize',13,'linewidth',2);hold on;
                plot(data.ParamsOut.NDF.NTDx_spline,data.ParamsOut.NDF.Difcor_spline,'g-','linewidth',2);hold on;
                if isfield(data.ParamsOut,'Revcor')
                    plot(data.ParamsOut.Revcor.Difcor_pred_x,data.ParamsOut.Revcor.Difcor_pred_y,'k-','linewidth',2);
                end
                set(gca,'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'layer','top','xlim',[min(data.ParamsOut.NDF.NTDx) max(data.ParamsOut.NDF.NTDx)]);
                box off;
                xlabel 'ITD (ms)';
                ylabel ('\Delta FIRING RATE (spikes s^{-1})','interpreter','tex');
                
                subplot(2,4,4);
                if isfield(data.ParamsOut,'Revcor')
                    plot(data.ParamsOut.NDF.Difcor_Mag_Freq,data.ParamsOut.Revcor.Difcor_Mag_Pred,'k-','linewidth',2);hold on;
                end
                plot(data.ParamsOut.NDF.Difcor_Mag_Freq,data.ParamsOut.NDF.Difcor_Mag_NTD,'g-','linewidth',2);
                set(gca,'linewidth',1,'fontsize',16,'layer','top','ticklength',[0.02 0.02]);
                xlabel 'FREQUENCY (kHz)';
                ylabel ('POWER','interpreter','tex');
                box off;
                
                
                %save the figure as a pdf
                cd(binplotdir);
                print(fh1,'-painters','-dpdf','-r600',[files(i).name(1:end-3) 'pdf']);
                close(1);
            end
            
            
        end
    end
end























