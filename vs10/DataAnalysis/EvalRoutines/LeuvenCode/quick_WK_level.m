function [] = quick_WK_level(AnNum,UnNum)
datadir = fullfile('C:\Users\Mark\Dropbox\MonRev',[AnNum '_' UnNum]);
savedir = 'C:\Users\Mark\Dropbox\MonRevBandwidth';
savefilename = fullfile(savedir,[AnNum '_' UnNum '_Bandwidth.mat']);
datafiles = dir(fullfile(datadir,'*MonRev*.mat'));
clear h1 SPL DF M T F U L h0 SPLpHz ERBHz ERBSPL LowerHz UpperHz SymmetryRatio up_slope lo_slope
load(fullfile(datadir,[AnNum '_' UnNum '_TC.mat']));
CF = TC.fit.cf/1000;
Th = TC.fit.minthr;
SR = TC.thr.sr;
Q10 = TC.fit.q10;
if CF<=2.5 && length(datafiles)>2
    for i = 1:length(datafiles)
        load(fullfile(datadir,datafiles(i).name));
        df_Hz = 1000*(kernels.h1ffax(2)-kernels.h1ffax(1));
        h0(i) = kernels.h0;
        h1(:,i) = kernels.h1c_mod;
        SPL(i) = kernels.SPL;
        SPLpHz(i) = SPL(i)-(10*log10(19990));
        ERBHz(i) = abs((nansum(kernels.h1cmag+20)*df_Hz)/20);
        ERBSPL(i) = SPL(i)-(10*log10(ERBHz(i)));
        
        
        signalRMS(i) = 10^(SPL(i)/20);
        DF(i) = kernels.h1cdomfreq;
        M(:,i) = kernels.h1magmask;
        
        QERB(i) = (DF(i)*1000)/ERBHz(i);
        
        %upper bandwidth
        UpperHz(i) = abs(nansum(kernels.h1cmag(kernels.h1ffax>=DF(i))+20)/20);
        %lower bandwidth
        LowerHz(i) = abs(nansum(kernels.h1cmag(kernels.h1ffax<DF(i))+20)/20);
        %Symmetry ratio (Lower/Upper)
        SymmetryRatio(i) = LowerHz(i)/UpperHz(i);
        
        
        cmag = kernels.h1cmag;
        fax = kernels.h1ffax(~isnan(cmag));
        cmag = cmag(~isnan(cmag));
        [~,di] = max(cmag);
        
        %upper and lower cord slopes - measured at 10-dB down point
        up_f_10dB = interp1(cmag(di:end),fax(di:end),-10);
        lo_f_10dB = interp1(cmag(1:di),fax(1:di),-10);
        
        bw_10dB(i) = up_f_10dB-lo_f_10dB;
        q10dB(i) = DF(i)/bw_10dB(i);
        
        %dy/dx
        up_slope(i) = -10/abs(log2(up_f_10dB/DF(i)));
        lo_slope(i) = -10/abs(log2(lo_f_10dB/DF(i)));
    end
    
    save(savefilename,'CF','SR','Th','Q10','SPL','SPLpHz','ERBSPL','ERBHz','QERB','DF','M',...
        'UpperHz','LowerHz','SymmetryRatio','up_slope','lo_slope','bw_10dB','q10dB');
else
    disp('CF too high for these analyses... limit is 2.5 kHz')
end
%Get the temchin et al data
load('C:\Users\Mark\Dropbox\Disparity_MATLAB\TemchinQ10Data.mat');

% T = kernels.Time;
% F = kernels.h1ffax;
% df = F(2)-F(1);
% 
% shiftFx = (0:length(datafiles)-1)*0.5;
% figure;set(gcf,'units','normalized','position',[0.3531    0.0533    0.1269    0.8556],'paperpositionmode','auto');
% plot(T,bsxfun(@plus,h1,shiftFx),'linewidth',1);
% set(gca,'xlim',[0 25],'ylim',[-0.5 length(datafiles)*0.5],'ytick',[0:1:0.5*(length(datafiles)-1)],...
%     'yticklabel',round(ERBSPL(1:2:end)),'activepositionproperty','outerposition','fontsize',14,'linewidth',1);
% xlabel 'TIME (ms)';
% ylabel 'SOUND LEVEL (dB SPL in ERB)';
% 
% MSpS = 20*log10(abs(fft(h1,2^13)));
% MSpS = bsxfun(@minus,MSpS,max(MSpS));
% MSpS = MSpS(1:(2^13)/2,:);
% MSpS(MSpS<-20)=nan;
% 
% TCfreq = TC.fit.freq;
% TCthresh = TC.fit.thr;
% 
% TCfreq = TCfreq(~isnan(TCthresh));
% TCthresh = TCthresh(~isnan(TCthresh));
% 
% 
% 
% figure;set(gcf,'units','normalized','position',[0.3144    0.2778    0.2194    0.3211],...
%     'paperpositionmode','auto','activepositionproperty','outerposition');
% contour(F',ERBSPL,MSpS',-20:2:-2,'linewidth',2);
% shading flat;
% hold on;
% plot(DF,ERBSPL,'k-o','linewidth',2);
% plot(TCfreq/1000,TCthresh,'-','linewidth',3,'color',[.6 .6 .6]);
% set(gca,'xlim',[0 2],'ylim',[-20 50],'clim',[-20 0],'fontsize',14,...
%     'tickdir','out','activepositionproperty','outerposition',...
%     'ticklength',[0.02 0.02],'linewidth',1);
% xlabel 'FREQUENCY (kHz)';
% ylabel 'LEVEL (dB SPL [in ERB])';
% % h = colorbar;
% % set(h,'fontsize',12);
% % ylabel (h,'GAIN (dB)');
% box off;
% axis square;
% 
% figure;set(gcf,'units','normalized','position',[0.3144    0.2778    0.2194    0.3211],...
%     'paperpositionmode','auto','activepositionproperty','outerposition');
% TCERB = interp1(temchin.erb.freq,temchin.erb.bw,CF*1000);
% plot([min(SPLpHz),max(SPLpHz)],[TCERB TCERB],'k-','linewidth',2);hold on;
% plot(SPLpHz,ERBHz,'ro-','linewidth',1,'markersize',10);
% set(gca,'yscale','log','fontsize',14,'tickdir','out','ytick',[10 100 1000],'yticklabel',[10 100 1000],...
%     'activepositionproperty','outerposition',...
%     'ticklength',[0.02 0.02],'ylim',[10 2000],'linewidth',1);
% ylabel 'ERB (Hz)';
% xlabel 'Sound Level (dBSPL/Hz)';
% box off;
% axis square;
% 
% figure;set(gcf,'units','normalized','position',[0.3144    0.2778    0.2194    0.3211],...
%     'paperpositionmode','auto','activepositionproperty','outerposition');
% TC_Q_ERB = (CF*1000)/TCERB;
% plot([min(SPLpHz),max(SPLpHz)],[TC_Q_ERB TC_Q_ERB],'k-','linewidth',2);hold on;
% plot(SPLpHz,QERB,'ro-','linewidth',1,'markersize',10);
% set(gca,'yscale','log','fontsize',14,'tickdir','out','ylim',[0.1 20],...
%     'activepositionproperty','outerposition','ticklength',[0.02 0.02],...
%     'ytick',[.1 1 10],'yticklabel',[0.1 1 10],'linewidth',1);
% ylabel ('Q_{ERB}','interpreter','tex');
% xlabel 'Sound Level (dBSPL/Hz)';
% box off;
% axis square;




