function [] = Ana_BFS(varargin)
An_Num = varargin{1};
UnEx_Num = varargin{2};
if nargin>2
    NITD_pos_Num = varargin{3};
    NITD_neg_Num = varargin{4};
end

if ~iscell(UnEx_Num)
    tempvar{1} = UnEx_Num; clear UnEx_Num;
    UnEx_Num = tempvar;clear tempvar;
end
if numel(UnEx_Num)>1
    foo;%still to do (merge ds)
else
    ds = dataset(An_Num,UnEx_Num{1});
    data = ds2struct(ds); clear ds;
end
if exist('NITD_pos_Num','var')
    ds = dataset(An_Num,NITD_pos_Num);
    NITDdata_pos = ds2struct(ds); clear ds;
    ds = dataset(An_Num,NITD_neg_Num);
    NITDdata_neg = ds2struct(ds); clear ds;
end
ind = regexp(UnEx_Num{1},'-');
Un_Num = UnEx_Num{1}(1:ind-1);
nrreps = size(data.Data.SpikeTimes,2);
indepvals = data.Stimulus.IndepVar.Values(~isnan(data.Stimulus.IndepVar.Values));
nrconds = length(indepvals);
repdur = data.Stimulus.Special.RepDur;
contrachan = data.Stimulus.StimParam.stimcntrl.contrachan;
stimdur = data.Stimulus.Special.BurstDur(1);
beat_f = data.Stimulus.Special.BeatFreq(1);
beat_p_ms = 1e3/beat_f;


if contrachan==1
    ipsichan = 'R';
else
    ipsichan = 'L';
end
if (strcmp(ipsichan,'R') && beat_f>0) || (strcmp(ipsichan,'L') && beat_f<0)
    correctITD = 1;%Flag to apply correction to ITD conventions
else
    correctITD = 0;
end

fh.AnaBFS = figure;
set(fh.AnaBFS,'units','normalized','position',[0.0225    0.0922    0.9481    0.7900],...
    'Numbertitle','off','Name',['Ana_BFS: ' An_Num '\' UnEx_Num{1} ': Ipsi = ' ipsichan],...
    'paperpositionmode','auto');

%Plot data using nested functions and return updated axis handles from each
[spk_ps,hax] = BFS_plotraster;
[hax] = BFS_plotmeanrate;
[hax,VS,Ph,Pval,spks,CD,CP] = BFS_plotIPD;
[hax,Tone_DifCor] = BFS_plotITD;
%If you have the NITD data for this unit - plot that with the inferred NITD
%from the tonal data
if exist('NITD_pos_Num','var')
    [Noise_DifCor,minNTD,maxNTD] = BFS_plotNTD;
    [rho,pval] = corr(Noise_DifCor/max(abs(Noise_DifCor)),Tone_DifCor'/max(abs(Tone_DifCor)));
    text(minNTD*0.9,1,'Tone difcor','color','r','fontsize',14);
    text(minNTD*0.9,0.75,'Noise difcor','color','b','fontsize',14);
    text(minNTD*0.9,0.5,sprintf('\\rho = %1.2f',rho),'interpreter','tex','fontsize',14,'color','k');
end
[hax,TDF,NDF] = BFS_plotDifCorFFT;

%Display a figure title for the whole thing
taxh = axes;
title (taxh,[An_Num '\' UnEx_Num{1} ': BFS'],'fontsize',20);
axis (taxh,'off');

%Set up the data structure
dataout.CD = CD;
dataout.CP = CP;
dataout.TDF = TDF;
dataout.NDF = NDF;

%save the data
dirstr = 'C:\work\BinRev\BFS\';
str = [dirstr An_Num '_' Un_Num '_BFS.mat'];
if exist(str,'file')
    str = [dirstr An_Num '_' Un_Num 'b_BFS.mat'];
end
save (str,'dataout');
%% Nested functions
%Raster
    function [sp_rate,hax_update]=BFS_plotraster
        %Raster
        hax.rast = subplot(2,4,1);
        linecols = get(hax.rast,'colororder');
        plotstep = 0.5/nrreps;
        col = 1;
        sp_rate = nan(1,nrconds);
        for i=1:nrconds
            for j=1:nrreps
                plot(hax.rast,[data.Data.SpikeTimes{i,j}; data.Data.SpikeTimes{i,j}]./1e3,...
                    [ones(1,length(data.Data.SpikeTimes{i,j})).*(i-0.25+((j-1)*plotstep));
                    ones(1,length(data.Data.SpikeTimes{i,j})).*(i-0.25+((j-1)*plotstep)+plotstep)],...
                    'color',linecols(col,:),'linewidth',1);
                hold on;
            end
            col = col+1;
            if col>7
                col = 1;
            end
            spks = ([data.Data.SpikeTimes{i,:}]);
            sp_rate(i) = sum(spks<=stimdur)/(nrreps*repdur/1e3);
        end
        set(hax.rast,'layer','top','ylim',[0 nrconds+1],'ytick',1:4:nrconds,...
            'yticklabel',indepvals(1:4:end)/1000,'fontsize',16,'linewidth',1,...
            'ticklength',[0.03 0.03],'xlim',[0 repdur/1e3]);
        plot([stimdur stimdur]/1000,[0 nrconds+1],'-','color',[0.6 0.6 0.6],'linewidth',2);
        ylabel (hax.rast,'FREQUENCY (kHz)');
        xlabel (hax.rast,'TIME (s)');
        axis square;box off;
        hax_update = hax;
    end
%Mean Firing Rate
    function [hax_update] = BFS_plotmeanrate
        linecols = get(hax.rast,'colororder');
        hax.meanrate = subplot(2,4,2);
        hold on;
        col = 1;
        for i=1:nrconds
            plot(hax.meanrate,i,spk_ps(i),'o','markerfacecolor',linecols(col,:),...
                'markeredgecolor','k','markersize',10,'linestyle','none');
            col = col+1;
            if col>7
                col = 1;
            end
        end
        set(hax.meanrate,'layer','top','xlim',[0 nrconds+1],'xtick',1:4:nrconds,...
            'xticklabel',indepvals(1:4:end)/1000,'ticklength',[0.03 0.03],...
            'fontsize',16,'linewidth',1);
        axis square;box off;
        ylabel (hax.meanrate,'FIRING RATE (Spikes s^-^1)','interpreter','tex');
        hax_update = hax;
        xlabel(hax.meanrate,'FREQUENCY (kHz)');
    end
%IPD functions, mean phase and vector strength
    function [hax_update,VS,phase_spks,Pval,spks,CD,CP] = BFS_plotIPD
        %plot normalized IPD curves in 2D
        hax.IPD = subplot(2,4,5);
        linecols = get(hax.IPD,'colororder');
        hold on;
        nrbins = 20;
        x.IPD = linspace(0,1,nrbins);
        col = 1;
        phase_spks = nan(1,nrconds); VS = nan(1,nrconds); Pval = nan(1,nrconds);
        spks = cell(1,nrconds);
        for i=1:nrconds
            spks{i} = ([data.Data.SpikeTimes{i,:}]);
            spks{i} = spks{i}(spks{i}<=stimdur);%only driven spikes
            beat_p = 1e3/data.Stimulus.Special.BeatFreq(i);%Beat period in millisec
            spks_p = mod(spks{i},beat_p)./beat_p;%Spikes times modulo beat period
            if correctITD
                spks_p=1-spks_p;
            end
            sum_sines = sum(sin((2*pi).*spks_p));%Sum of sines and cosines
            sum_cosines = sum(cos((2*pi).*spks_p));
            phase_spks(i) = atan2(sum_sines,sum_cosines);%Mean phase
            VS(i) = sqrt(sum_sines^2 + sum_cosines^2)/length(spks{i});%Vector strength
            Pval(i) = exp(-length(spks{i})*(VS(i)^2));
            y.IPD(i,:) = hist(spks_p,x.IPD);
%             y.IPD(i,:) = y.IPD(i,:)./max(y.IPD(i,:));%Normalize to max.
            y.IPD(i,[1 end]) = sum(y.IPD(i,[1 end])); %Take care of edge effects in histogram (for display purposes only since they're circular data)
        end
        for i=1:nrconds
            if Pval(i)<=0.001
                plot(hax.IPD,x.IPD,y.IPD(i,:)./max(max(y.IPD)),'color',linecols(col,:),'linewidth',2);
            else
                plot(hax.IPD,x.IPD,y.IPD(i,:)./max(max(y.IPD)),'color',[0.5 0.5 0.5],'linewidth',1);
            end
            col = col+1;
            if col>7
                col = 1;
            end
        end
        set(hax.IPD,'layer','top','ylim',[0 1],'ticklength',[0.03 0.03],...
            'fontsize',16,'linewidth',1);
        ylabel(hax.IPD,'NORMALIZED FIRING RATE');
        xlabel(hax.IPD,'IPD (cycles)');
        axis square; box off;
        
        %plot phase and vector strength
        phase_spks = unwrap(phase_spks);%Unwrap the phase data;
        hax.P = subplot(2,4,4);
        linecols = get(hax.P,'colororder');
        hax.VS = subplot(2,4,3);
        hold(hax.P,'on');
        hold(hax.VS,'on');
        col = 1;
        for i=1:nrconds
            if Pval(i)<=0.001
                plot(hax.VS,i,VS(i),'marker','o','markerfacecolor',linecols(col,:),'markeredgecolor','k','markersize',10,'linestyle','none');
                plot(hax.P,i,phase_spks(i)/(2*pi),'marker','o','markerfacecolor',linecols(col,:),'markeredgecolor','k','markersize',10,'linestyle','none');
            else
                plot(hax.VS,i,VS(i),'marker','o','markerfacecolor','none','markeredgecolor',[0.6 0.6 0.6],'markersize',10,'linestyle','none');
            end
            col = col+1;
            if col>7
                col = 1;
            end
        end
        set(hax.P,'layer','top','ticklength',[0.03 0.03],'xlim',[0 nrconds+1],'xtick',1:4:nrconds,...
            'xticklabel',indepvals(1:4:end)/1000,...
            'fontsize',16,'linewidth',1);
        set(hax.VS,'layer','top','xlim',[0 nrconds+1],'xtick',1:4:nrconds,...
            'xticklabel',indepvals(1:4:end)/1000,'ticklength',[0.03 0.03],...
            'fontsize',16,'linewidth',1,'ylim',[0 1]);
        ylabel(hax.P,'PHASE (cycles)');
        xlabel(hax.P,'FREQUENCY (kHz)');
        ylabel(hax.VS,'VECTOR STRENGTH');
        xlabel(hax.VS,'FREQUENCY (kHz)');
        axis (hax.P,'square'); box(hax.P,'off');
        axis (hax.VS,'square'); box(hax.VS,'off');
        
        %Use Rayleigh statistic to compute a set of weights for the linear
        %regression of phase vs frequency
        weights = Pval(Pval<=0.001);
        weights = -log(weights);%This is the Rayleigh stat
        weights =  weights./max(weights);
        options = fitoptions('Method','linearleastsquares','weights',weights);
        f = fittype({'1','x'},'coefficients',{'beta1','beta2'},'options',options);
        beta = fit(indepvals(Pval<=0.001)/1000,(phase_spks(Pval<=0.001)/(2*pi))',f);
        fity = beta.beta1 + beta.beta2.*(indepvals/1000);
        plot(hax.P,1:nrconds,fity,'--k','linewidth',4);
        %Display some information in text
        str = sprintf('CD = %1.2f ms\nCP = %1.2f cycles',beta.beta2,beta.beta1);
        th = text (indepvals(1)/1000,max(get(hax.P,'ylim'))*0.9,str,'fontsize',14,'verticalalignment','top','interpreter','tex');
        set(th,'parent',hax.P);
        hax_update = hax;
        CD = beta.beta2;
        CP = beta.beta1;
    end

    function [hax_update,Difcor]=BFS_plotITD
        hax.ITD = subplot(2,4,6);
        if exist('NITDdata_neg','var')
            %If you have the NTD data, compute the tone-based NITD
            %simulation at the same ITDs. Otherwise, pick something
            %sensible
            minITD = min(NITDdata_neg.Stimulus.IndepVar.Values)/1e3;
            maxITD = max(NITDdata_neg.Stimulus.IndepVar.Values)/1e3;
            nITD = numel(NITDdata_neg.Stimulus.IndepVar.Values);
        else
            minITD = -5;
            maxITD = 5;
            nITD = 41;
        end
        x.ITD = linspace(minITD,maxITD,nITD);
        ITD = cell(1,nrconds);
        ITD_flip = cell(1,nrconds);
        for i=1:nrconds
            carrp = 1e3/indepvals(i);
            spks_p = mod(spks{i}./beat_p_ms,1);
            spks_c = ceil(spks{i}./beat_p_ms);
            spks_itd = nan(1,length(spks{i}));
            spks_itd(spks_p<0.5) = spks_p(spks_p<0.5)*carrp;
            spks_itd(spks_p>=0.5) = (spks_p(spks_p>=0.5)-1)*carrp;
            if correctITD
                spks_itd = spks_itd*-1;
            end
            spks_itd = [spks_itd+((spks_c-1)*carrp) spks_itd(spks_c>1)-((spks_c(spks_c>1)-1)*carrp)];
            
            ITD{i} = spks_itd;
            ITD_flip{i} = spks_itd+(carrp/2);
            y.ITD(i,:) = hist(ITD{i}(ITD{i}>=minITD & ITD{i}<=maxITD),x.ITD);
            y.ITD_flip(i,:) = hist(ITD_flip{i}(ITD_flip{i}>=minITD & ITD_flip{i}<=maxITD),x.ITD);
        end
        colormap (flipud(gray));
        pcolor(hax.ITD,repmat(x.ITD,nrconds,1),repmat(indepvals/1000,1,numel(x.ITD)),y.ITD);hold on;
        set(hax.ITD,'layer','top','linewidth',1,'fontsize',16,'ticklength',[0.03 0.03],'ylim',[min(indepvals/1000) max(indepvals/1000)],'xlim',[minITD maxITD],'xtick',[minITD 0 maxITD]);
        shading (hax.ITD,'flat'); axis (hax.ITD,'square');box (hax.ITD,'off');view (hax.ITD,[0,90]);grid(hax.ITD,'off');
        %add some lines indicating ethological range of ITDs and the pi
        %limit
        plot([-0.25 0.25;-0.25 0.25],[min(indepvals/1000) min(indepvals/1000); max(indepvals/1000) max(indepvals/1000)],'r--','linewidth',2);
        plot(0.5*(1000./indepvals),indepvals/1000,'r-',-0.5*(1000./indepvals),indepvals/1000,'r-','linewidth',2);
%         plot(0.25*0.5*(1000./indepvals),indepvals/1000,'b--',0.25*-0.5*(1000./indepvals),indepvals/1000,'b--','linewidth',2);
        xlabel(hax.ITD,'ITD (ms)');
        ylabel(hax.ITD,'FREQUENCY (kHz)');
        
        hax.DC = subplot(2,4,7);
        hold(hax.DC,'on');
        Difcor = y.ITD-y.ITD_flip;
        Difcor = sum(Difcor);        
        area(hax.DC,[minITD -0.25],[1.1 1.1],-1.1,'facecolor',[0.7 0.7 0.7],'edgecolor','k','linewidth',1);
        area(hax.DC,[0.25 maxITD],[1.1 1.1],-1.1,'facecolor',[0.7 0.7 0.7],'edgecolor','k','linewidth',1);
%         plot(hax.DC,x.ITD,Difcor/max(Difcor),'k*','linewidth',2,'markersize',10);
        Difcor_spline = spline(x.ITD,Difcor/max(abs(Difcor)),linspace(minITD,maxITD,1000));
        plot(hax.DC,linspace(minITD,maxITD,1000),Difcor_spline,'r-','linewidth',2);

        set(hax.DC,'layer','top','linewidth',1,'fontsize',16,'ticklength',[0.03 0.03],'ylim',[-1.1 1.1],'xlim',[minITD maxITD],'xtick',[minITD 0 maxITD]);
        axis (hax.DC,'square'); box(hax.DC,'off');
        xlabel(hax.DC,'ITD (ms)');
        ylabel(hax.DC,'NORMALIZED FIRING RATE');
        hax_update = hax;
    end
    function [DifCor,minNTD,maxNTD] = BFS_plotNTD
        nrcond = numel(NITDdata_neg.Stimulus.IndepVar.Values);
        nrreps = size(NITDdata_neg.Data.SpikeTimes,2);
        stimlen = NITDdata_neg.Stimulus.Special.BurstDur(1);%ms
        x.NTD = NITDdata_neg.Stimulus.IndepVar.Values/1e3; %NTD x axis in ms
        NTD_neg = nan(nrcond,1);NTD_pos = nan(nrcond,1);
        for i=1:nrcond
            NTD_neg(i) = length(NITDdata_neg.Data.SpikeTimes{i}(NITDdata_neg.Data.SpikeTimes{i}<=stimlen))/(stimlen*1e-3*nrreps);
            NTD_pos(i) = length(NITDdata_pos.Data.SpikeTimes{i}(NITDdata_pos.Data.SpikeTimes{i}<=stimlen))/(stimlen*1e-3*nrreps);
        end
        if correctITD
            NTD_neg = flipud(NTD_neg);
            NTD_pos = flipud(NTD_pos);
        end
        minNTD = min(x.NTD);
        maxNTD = max(x.NTD);
        SumCor = NTD_neg+NTD_pos;
        DifCor = NTD_pos-NTD_neg;
        DifCor_spline = spline(x.NTD,DifCor/max(abs(DifCor)),linspace(minNTD,maxNTD,1000));
%         plot(hax.DC,x.NTD,DifCor/max(DifCor),'b*');
        plot(hax.DC,linspace(minNTD,maxNTD,1000),DifCor_spline,'b-','linewidth',2);
    end
    function [hax_update,TDF,NDF] = BFS_plotDifCorFFT
        if exist('NITDdata_neg','var')
            minITD = min(NITDdata_neg.Stimulus.IndepVar.Values)/1e3;
            maxITD = max(NITDdata_neg.Stimulus.IndepVar.Values)/1e3;
            nITD = numel(NITDdata_neg.Stimulus.IndepVar.Values);
        else
            minITD = -5;
            maxITD = 5;
            nITD = 41;
        end
        dt = (maxITD-minITD)*1e-3/nITD;
        sr = 1/dt;
        freq = (sr/2)*linspace(0,1,floor(nITD/2)+1);
        
        if exist('Noise_DifCor','var')
            s = fft([Tone_DifCor'/max(abs(Tone_DifCor)) Noise_DifCor/max(abs(Noise_DifCor))]);
        else
            s = fft(Tone_DifCor'/max(abs(Tone_DifCor)));
        end
        s = s(1:floor(nITD/2)+1,:);
        s = 20*log10(abs(s));
        hax.FFT = subplot(2,4,8);
        set(hax.FFT,'colororder',[1 0 0; 0 0 1],'nextplot','replacechildren','xscale','log',...
            'xtick',[0.1 1],'xticklabel',[0.1 1],'fontsize',16,'ylim',[-42 2],'xlim',[0.1 3],'ticklength',[0.03 0.03],'layer','top');
        plot(hax.FFT,freq(2:end)/1000,s(2:end,:)-repmat(max(s),length(freq(2:end)),1),'linewidth',2);hold on;
        [dummy,maxind] = max(s);
        if exist('NITDdata_neg','var')
            plot(hax.FFT,freq(maxind(1))/1000,0,'ro','markerfacecolor','r','markersize',10);
            plot(hax.FFT,freq(maxind(2))/1000,0,'bo','markerfacecolor','b','markersize',10);
            TDF = freq(maxind(1))/1000;
            str = sprintf('Tone DF = %1.2f kHz',freq(maxind(1))/1000);
            text (0.15,-30,str,'fontsize',14,'verticalalignment','top','horizontalalignment','left','interpreter','tex','color','r','backgroundcolor','w');
            NDF = freq(maxind(2))/1000;
            str = sprintf('Noise DF = %1.2f kHz',freq(maxind(2))/1000);
            text (0.15,-35,str,'fontsize',14,'verticalalignment','top','horizontalalignment','left','interpreter','tex','color','b','backgroundcolor','w');
        else
            plot(hax.FFT,freq(maxind(1))/1000,0,'ro','markerfacecolor','r','markersize',10);
        end
        xlabel (hax.FFT,'FREQUENCY (kHz)');
        ylabel (hax.FFT,'MAGNITUDE (dB re. max.)');
        axis (hax.FFT,'square'); box (hax.FFT,'off');
        hax_update = hax;
    end
end