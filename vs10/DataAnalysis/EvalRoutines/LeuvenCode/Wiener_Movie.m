function [] = Wiener_Movie(varargin)

%Makes a movie showing a 1-st order WK building up rep by rep


%M. Sayles 2015. Leuven.

%Inputs:
%

%Outputs:
%-------------------------------------------------------------------------
%% Set up the figure window and axes.
fh = figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9],'color','w','paperpositionmode','auto');
tcax = subplot(2,2,1);
irax = subplot(2,2,3);
msax = subplot(2,2,2);
phax = subplot(2,2,4);
%% Get input data in correct format
AnNum = varargin{1};
count = 0;
for i = 2:nargin
    if strcmp(varargin{i-1},'thr')
        thrdata = varargin{i};
    elseif strcmp(varargin{i-1},'spl')
        spldata = varargin{i};
    elseif ~strcmp(varargin{i},'thr') && ~strcmp(varargin{i},'spl')
        count = count+1;
        ndzdata{count} = varargin{i};
    end
end
if exist('ndzdata','var')
    nrdata = length(ndzdata);
    for i = 1:nrdata
        try
            dsz{i} = dataset(AnNum,[ndzdata{i} '-NTD']);
        catch
            dsz{i} = dataset(AnNum,[ndzdata{i} '-NSPL']);
        end
    end
    
    
    %% Set some parameters
    minST = 50;%Minimum spike time in ms
    maxST = dsz{1}.Stimulus.Special.BurstDur(1); %Maximum spike time in ms
    RCwin = 30;%Length of reverse-correlation window in ms
    Nboot = 10;
    [dummy, dt] = StimSam(dsz{1},1);
    dt = 1e-3*dt; %us -> ms
    Nwin = round(RCwin/dt);
    SpikeStep = 20;
end
%% Get the tuning curve data
dsthr = dataset(AnNum,[thrdata '-THR']);
TC = EvalTHR(dsthr);
ParamsOut.TC = TC;
maxSPL = dsthr.stimulus.StimParam.SPL;
CF_frq = TC.fit.cf;
CF_spl = TC.fit.minthr;
Q10 = TC.fit.q10;
if isnan(Q10)
    Q10 = TC.thr.q10;
end

%Plot the TC
tch = semilogx(tcax,TC.fit.freq(TC.fit.thr<maxSPL),TC.fit.thr(TC.fit.thr<maxSPL),'b-','linewidth',5);
hold(tcax, 'on');
semilogx(tcax,TC.thr.freq(TC.thr.thr<maxSPL),TC.thr.thr(TC.thr.thr<maxSPL),'k+','markersize',8,'linewidth',1);

set(tcax,'xlim',[50 10e3],'xtick',[10e1 10e2 10e3],'xticklabel',[0.1 1 10],...
    'fontsize',16,'ticklength',[0.02 0.02],'ylim',[-10 80]);
ylabel(tcax,'THRESHOLD (dB SPL)');
xlabel(tcax,'TONE FREQUENCY (kHz)');
box(tcax,'off');
set(fh,'currentaxes',tcax);
text(70,80,sprintf('CF: %2.2f kHz \nThreshold: %2.1f dB SPL \nQ_{10}: %2.2f',CF_frq/1000,CF_spl,Q10),...
    'fontsize',12,'color','k','interpreter','tex','verticalalignment','top','horizontalalignment','left');
plot(tcax,[CF_frq CF_frq],[min(get(tcax,'ylim')) max(get(tcax,'ylim'))],'--',...
    [min(get(tcax,'xlim')) max(get(tcax,'xlim'))],[CF_spl CF_spl],'--','linewidth',1,'color',[.6 .6 .6]);
plot(tcax,CF_frq,CF_spl,'ro','markersize',10,'markerfacecolor','r');

%% 1. Get revcor cumulatively and plot on update on each rep
if exist('ndzdata','var')
    [wv, dt] = StimSam(dsz{1},1);
    Nchan = size(wv, 2);
    dt = 1e-3*dt; %us -> ms
    Nsamples = length(wv);
    allspikes = dsz{1}.SPT;%get spike times
    allspikes =cat(2,allspikes{4,:});
    allspikes = allspikes(allspikes>=minST & allspikes<=maxST);%limit spike times to only include driven spikes
    SpikeCounts = 1:SpikeStep:length(allspikes);
    wbh = waitbar(0,sprintf('Bootstrapping %2.0f revcors: Please wait...',length(SpikeCounts)));
    for i = 1:length(SpikeCounts)
        spikes = allspikes(1:SpikeCounts(i));%concatenate all spikes from all reps
        [h1{i},Time,SpikeCount(i)] = local_h1_kernel(spikes);
        for j = 1:Nboot
            nullspikes = scramble_spikes(spikes)';
            nullspikes = nullspikes(nullspikes>=minST & nullspikes<=maxST);
            [nullh1{i,j},dummy1,dummy2] = local_h1_kernel(nullspikes);
            waitbar((((i-1)*Nboot)+j)/(Nboot*length(SpikeCounts)));
        end
    end
    close (wbh);
    for i = 1:length(h1)
        %Impulse response (Time domain)
        cla(irax);
        cla(msax);
        cla(phax);
        TotalSpikes = (i-1)*SpikeStep+1;
        sh1_lp = local_low_pass(h1{i});
        scaleF = max(abs(h1{i}));
        sh1 = h1{i}./scaleF;
        sh1_lp = sh1_lp./scaleF;
        for j = 1:Nboot
            nullsh1{i,j} = nullh1{i,j}./scaleF;
            plot(irax,Time,nullsh1{i,j},'color',[.6 .6 .6]);hold(irax,'on');
        end
        plot(irax,Time,sh1,'-','color',[0.5 0.5 1],'linewidth',2);
        plot(irax,Time,sh1_lp,'b-','linewidth',4);
%         set(irax,'ylim',[-1.1 1.1],'fontsize',16,'xlim',[0 30],'ytick',[-1 1],'layer','top','ticklength',[0.02 0.02]);
set(irax,'ylim',[-1.1 1.1],'fontsize',16,'xlim',[0 15],'ytick',[-1 1],'layer','top','ticklength',[0.02 0.02]);
        ylabel(irax,'AMPLITUDE');
        xlabel(irax,'TIME (ms)');
        box(irax,'off');
        hold(irax,'off');
        set(fh,'currentaxes',irax);
        text(10,-1,sprintf('# Spikes: %2.0f',TotalSpikes),'fontsize',12,'color','k',...
            'verticalalignment','bottom','horizontalalignment','right');
        
        %Magnitude spectrum
        [sh1mag,nullsh1mag,ffax,sh1phase,sh1zscore] = Phase_Freq_analysis(sh1,horzcat(nullsh1{i,:}));
        correctionF = mean(mean(nullsh1mag(ffax<4000,:)));
        sh1mag = sh1mag-correctionF;
        sh1mag = Trifilter(sh1mag',13);
        nullsh1mag =nullsh1mag-correctionF;
        semilogx(msax,ffax,nullsh1mag,'color',[0.6 0.6 0.6]);
        hold(msax,'on');
        semilogx(msax,ffax,sh1mag,'-','color',[0.5 0.5 1],'linewidth',2);
        set(msax,'ylim',[-10 40],'xlim',[50 4500],'xtick',[100 1000],...
            'xticklabel',[0.1 1],'fontsize',16,'ticklength',[0.02 0.02],'layer','top');
        box(msax,'off');
        
        ylabel(msax,'MAGNITUDE (dB)');
        xlabel(msax,'FREQUENCY (kHz)');
        set(fh,'currentaxes',msax);
        
        %find the region which is above the noise floor
        [maxval,ind] = max(sh1mag);
        startind = find(sh1zscore(1:ind)<3,1,'last');
        stopind = find(sh1zscore(ind:end)<3,1,'first')+ind;
        semilogx(msax,ffax(startind:stopind),sh1mag(startind:stopind),'b-','linewidth',4);
        %Find the dominant frequency within that region
        [maxval,ind] = max(sh1mag(startind:stopind));
        df = ffax(startind+ind-1);
        plot(df,maxval,'ro','markersize',10,'markerfacecolor','r');
        if df<4e3
            text(70,40,sprintf('# Spikes: %2.0f \nDF: %2.2f kHz',TotalSpikes,df/1000),'fontsize',12,'color','k',...
                'verticalalignment','top','horizontalalignment','left');
        else
            text(70,40,sprintf('# Spikes: %2.0f',TotalSpikes),'fontsize',12,'color','k',...
                'verticalalignment','top','horizontalalignment','left');
        end
        hold(msax,'off');
        
        %Phase spectrum
        sh1phase(1:startind-1) = NaN;
        sh1phase(stopind+1:end) = NaN;
        sh1phase = unwrap(sh1phase)./(2*pi);
        set(fh,'currentaxes',phax);
        if df<4e3
            [beta_out,gof_out]=local_fit_linear(ffax(startind:stopind),sh1phase(startind:stopind));
            semilogx(phax,1:1:4500,feval(beta_out,1:1:4500),'b--');
            text(70,-7,sprintf('# Spikes: %2.0f \nGroup Delay: %2.2f ms',TotalSpikes,-1000*beta_out.a),'fontsize',12,'color','k',...
                'verticalalignment','bottom','horizontalalignment','left');
        end
        
        hold(phax,'on');
        semilogx(phax,ffax(startind:stopind),sh1phase(startind:stopind),'b-','linewidth',4);
        set(phax,'ylim',[-8 2],'xlim',[50 4500],'xtick',[100 1000],...
            'xticklabel',[0.1 1],'fontsize',16,'ticklength',[0.02 0.02],'layer','top','xscale','log');
        box(phax,'off');
        
        ylabel(phax,'PHASE (cycles)');
        xlabel(phax,'FREQUENCY (kHz)');
        
        F(i) = getframe(fh);
    end
    %         movie2avi(F,'revcormovie','compression','FFDS','fps',6,'quality',100);
    VV = VideoWriter('RevcorVideo.mp4','MPEG-4');
    VV.FrameRate = 10;
    VV.Quality = 100;
    open(VV);
    for i = 1:length(h1)
        writeVideo(VV,F(i));
    end
    close(VV);
end

%% Local functions
    function [filtout] = local_low_pass(IN)
        fsamp = 1000/dt; %sample frequency
        fcuts = [4500 5500]; % Transition bands
        mags = [1 0];
        devs = [0.01 10^(-40/20)];% 1% passband ripple and 40-dB attenuation of stopband
        
        [n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);
        b = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale'); %Filter coefficients to match the required parameters
        filtout = filtfilt(b,1,IN); %Zero-phase filtering
    end

    function [h1out,Time,Nspikes] = local_h1_kernel(ST)
        Nspikes = length(ST);
        tx = (0:Nsamples)*dt; %Binning edges in ms ...
        pst = histc(ST,tx);%Make PSTH with same sampling frequency as stimulus
        for ii = 1:Nchan
            h1out(:,ii) = flipud(xcorr(wv(:,ii),pst,Nwin))/Nspikes;%Get the revcor (h1 kernel) - cross correlation of the pst and the stimulus waveform
        end
        h1out = h1out((1:Nwin)+(Nwin+1),:);
        Time = tx(1:Nwin);
    end
    function [scramspikes] = scramble_spikes(spikes)
        if length(spikes)==1
            scramspikes = spikes+(rand*10);
        else
            % shuffle spike train by randomizing the ISI's
            D = diff([0; sort(spikes')]); % first sort to make sure diff gives the ISI's
            NInt = length(D);
            scramspikes = cumsum(D(randperm(NInt)));
        end
    end
    function [h1magout,nullh1magout,freqscale,h1phaseout,zcritout] = Phase_Freq_analysis(y,nully)
        nfft = 2^13;
        h1fft = fft(y,nfft);
        h1magout = 20*log10(abs(h1fft));
        h1magout = h1magout(1:floor(nfft/2),:);
        nullh1magout = 20*log10(abs(fft(nully,nfft)));
        nullh1magout = nullh1magout(1:floor(nfft/2),:);
        nf = 0.5*(1000/dt);
        freqscale = nf*linspace(0,1,floor(nfft/2));
        h1phaseout = angle(h1fft(1:floor(nfft/2),:));
        %get z-scores for magnitude spectrum
        zcritout = zeros(size(h1magout));
        for ii = 1:Nchan
            zcritout(:,ii)= ((10.^(h1magout/20))-mean(10.^(nullh1magout/20),2))./std(10.^(nullh1magout/20),1,2);
        end
        %Apply smoothing (Triangular smoothing function)
        for ii = 1:Nchan
            [zcritout(:,ii)] = Trifilter(zcritout(:,ii)',7)';
        end
    end
    function [beta_out,gof_out]=local_fit_linear(x_in,y_in)
        s = fitoptions('Method','LinearLeastSquares','Lower',[-Inf -Inf],'Upper',[Inf Inf]);
        f = fittype({'x','1'},'coefficients',{'a','b'},'options',s);
        [beta_out,gof_out] = fit(x_in(:),y_in,f);
    end
    function [otvect]=Trifilter(invect,nfw)
        nfwi = 2*floor(nfw/2) + 1;
        filt = zeros(1,nfwi);
        summ = 0;
        nfw2 = floor(nfwi/2);
        for jj=1:nfw2
            filt(jj) = jj;
            filt(nfwi+1-jj) = jj;
            summ = summ + 2*jj;
        end
        nfw3 = nfw2 + 1;
        filt(nfw3) = nfw3;
        summ = summ + nfw3;
        filt = filt./summ;
        svect = size(invect,2) + 2*nfw2;
        vect1 = zeros(1,svect);
        vect1(1:nfw2) = invect(1)*ones(1,nfw2);
        vect1(nfw3:svect-nfw2) = invect;
        vect1(svect-nfw2+1:svect) = invect(size(invect,2))*ones(1,nfw2);
        vect2 = conv(vect1, filt);
        otvect = vect2(2*nfw2+1:svect);
    end
end