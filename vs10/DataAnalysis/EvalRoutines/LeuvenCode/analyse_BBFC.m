function [BFS] = analyse_BBFC(data,params,ITDx)

count = 0;
for ddx = 1:length(data)
    startFreq = min(data{ddx}.stim_param.StartFreq);
    endFreq = min(data{ddx}.stim_param.EndFreq);
    stepFreq = data{ddx}.stim_param.StepFreq;
    Freqs = startFreq:stepFreq:endFreq;
    [nrconds,nrreps] = size(data{ddx}.spikes);
    for idx = 1:nrconds
        for jdx = 1:nrreps
            if ~isempty(data{ddx}.spikes{idx,jdx})
                count = count+1;
                playedfreqs(count) = Freqs(idx);
                spikes{count} = data{ddx}.spikes{idx,jdx};
                spl(count) = data{ddx}.stim_param.SPL(1);
            end
        end
    end
    dataSPL(ddx) = data{ddx}.stim_param.SPL(1);
end
Freqs = playedfreqs;

uSPLs = unique(spl);
nrSPLs = length(uSPLs);

dt = params.dt;
ITDx = sort(ITDx,'ascend');
ditd = diff(ITDx([1,2]));
stimdur = data{1}.stim_param.BurstDur;
beatF = data{1}.stim_param.BeatFreq; %Hz
beatP = 1e3/beatF;%beat period (in milliseconds)

minITD = min(ITDx);
maxITD = max(ITDx);

for SPLidx = 1:nrSPLs
    
    thesedata = find(dataSPL==uSPLs(SPLidx));
    
    SPL = uSPLs(SPLidx);
    
    ufreq = unique(Freqs(spl==SPL));
    Nfreq = length(ufreq);
    
    monauralF = [ufreq'+beatF ufreq'];
    
    ITDpos = zeros(Nfreq,length(ITDx));
    ITDneg = zeros(Nfreq,length(ITDx));
    meanphase = zeros(Nfreq,1);
    VS = zeros(Nfreq,1);
    VSipsi = zeros(Nfreq,1);
    VScontra = zeros(Nfreq,1);
    Pval = zeros(Nfreq,1);
    Pvalipsi = zeros(Nfreq,1);
    Pvalcontra = zeros(Nfreq,1);
    syncrate = zeros(Nfreq,1);
    spikerate = zeros(Nfreq,1);
    
    for i = 1:Nfreq
        clear STpos; clear STneg;
        carrP = 1e3/ufreq(i);
        ind = find(Freqs==ufreq(i) & spl==SPL);
        S = cat(2,spikes{ind});
        S = S(S<=stimdur);
        Nspikes = length(S);
        Nreps = length(ind);
        
        %Spike rate
        spikerate(i) = Nspikes/Nreps/(stimdur/1000);%Spikes per second
        
        SCpos = mod(S/beatP,1);%driven spike times in cycles re beat period
        SCneg = mod((S+(beatP/2))/beatP,1);%driven spike times in cycles re beat period
        
        SCipsi = mod(S/(1e3/monauralF(i,1)),1);%driven spike times in cycles re ipsi period
        SCcontra = mod(S/(1e3/monauralF(i,2)),1);%driven spike times in cycles re contra period

        
        %convert cycles to ITD - simulate 5 beat cycles (so that it is
        %consistent with the noise data)
        STpos = [SCpos(SCpos<=0.5) -1+SCpos(SCpos>0.5)];
        STneg = [SCneg(SCneg<=0.5) -1+SCneg(SCneg>0.5)];
        
        STpos = bsxfun(@plus,STpos,(-5:5)');
        STneg = bsxfun(@plus,STneg,(-5:5)');
        
        STpos = STpos(:)*carrP;
        STneg = STneg(:)*carrP;
        
        if strcmp(params.AnNum,'F15075')%ITD convention wrong for this animal
            STpos = -STpos;
            STneg = -STneg;
        end
        
        ITDpos(i,:) = hist(STpos(STpos>=minITD & STpos<=maxITD),ITDx);
        ITDneg(i,:) = hist(STneg(STneg>=minITD & STneg<=maxITD),ITDx);
        
        ITDpos(i,:) = ITDpos(i,:)/(Nreps*stimdur/1000)/ditd;
        ITDneg(i,:) = ITDneg(i,:)/(Nreps*stimdur/1000)/ditd;
        
        %Mean phase
        sines = sin(2*pi*SCpos);
        cosines = cos(2*pi*SCpos);
        meanphase(i) = atan2(sum(sines),sum(cosines));%Mean phase
        meanphase(i) = meanphase(i)/(2*pi); %radians --> cycles
        
        sinesipsi = sin(2*pi*SCipsi);
        cosinesipsi = cos(2*pi*SCipsi);
        meanphaseipsi(i) = atan2(sum(sinesipsi),sum(cosinesipsi));%Mean phase
        meanphaseipsi(i) = meanphaseipsi(i)/(2*pi); %radians --> cycles
        
        sinescontra = sin(2*pi*SCcontra);
        cosinescontra = cos(2*pi*SCcontra);
        meanphasecontra(i) = atan2(sum(sinescontra),sum(cosinescontra));%Mean phase
        meanphasecontra(i) = meanphasecontra(i)/(2*pi); %radians --> cycles
        
        %Vector strength
        VS(i) = sqrt(sum(sines)^2 + sum(cosines)^2)/Nspikes;%Vector strength
        VSipsi(i) = sqrt(sum(sinesipsi)^2 + sum(cosinesipsi)^2)/Nspikes;%Vector strength
        VScontra(i) = sqrt(sum(sinescontra)^2 + sum(cosinescontra)^2)/Nspikes;%Vector strength
        
        %P-value
        Pval(i) = exp(-Nspikes*(VS(i)^2));
        Pvalipsi(i) = exp(-Nspikes*(VSipsi(i)^2));
        Pvalcontra(i) = exp(-Nspikes*(VScontra(i)^2));
        
        %Synchronized rate
        syncrate(i) = VS(i)*spikerate(i);
        
        
    end
    
    if ~length(VSipsi)==length(VScontra)
        error('Inconsistent VS vectors')
    end
    if ~length(VSipsi)==length(ufreq)
        error('Inconsistent VS vectors')
    end
    
    %Do the circular-linear regression on the freq-phase data to get CP and CD
    weights = syncrate(Pval<0.001);
    weights = weights/max(weights);
    [CP,CD,dummy,dummy] = fit_circ_lin(ufreq(Pval<0.001)/1000,meanphase(Pval<0.001),weights');
    
    %Get the difcor
    difcorrate = nanmean(ITDpos)-nanmean(ITDneg);
    
    %spline fit
    ITDxspline = min(ITDx):dt:max(ITDx);
    Dspline = spline(ITDx,difcorrate,ITDxspline);
    
    %Get the BD
    [dummy,ind] = max(Dspline);
    BD = ITDxspline(ind);
    
    %Get the power spectrum
    NFFT = 2^13;
    nf = 0.5/(dt/1000);
    freqax = linspace(0,1,NFFT/2)*nf;
    vals = Dspline;
    hanwin = hann(length(vals));
    spec = fft((vals'-mean(vals)).*hanwin,NFFT);
    Dphase = angle(spec(1:NFFT/2));
    Dspec = 20*log10(abs(spec(1:NFFT/2)));
    Dspec = Dspec-max(Dspec);
    Dphase(Dspec<-20) = nan;
    [dummy,ind] = max(Dspec);
    DspecPeakHz = freqax(ind);
    
    
    %Make the output structure
    BFS{SPLidx}.beat.CP = CP;
    BFS{SPLidx}.beat.CD = CD;
    BFS{SPLidx}.beat.meanphase = meanphase;
    BFS{SPLidx}.beat.ipsiphase = meanphaseipsi;
    BFS{SPLidx}.beat.contraphase = meanphasecontra;
    BFS{SPLidx}.beat.VS = VS;
    BFS{SPLidx}.beat.VSipsi = VSipsi;
    BFS{SPLidx}.beat.VScontra = VScontra;
    BFS{SPLidx}.beat.Pval = Pval;
    BFS{SPLidx}.beat.Pvalipsi = Pvalipsi;
    BFS{SPLidx}.beat.Pvalcontra = Pvalcontra;
    BFS{SPLidx}.beat.syncrate = syncrate;
    BFS{SPLidx}.beat.ufreq = ufreq;
    BFS{SPLidx}.beat.freq = playedfreqs(spl==SPL);
    BFS{SPLidx}.beat.spiketimes = spikes(spl==SPL);
    BFS{SPLidx}.positive.rate = ITDpos;
    BFS{SPLidx}.negative.rate = ITDneg;
    BFS{SPLidx}.xvals.ITDx = ITDx;
    BFS{SPLidx}.difcor.rate = difcorrate;
    BFS{SPLidx}.xvals.ITDxspline = ITDxspline;
    BFS{SPLidx}.difcor.ratespline = Dspline;
    BFS{SPLidx}.difcor.bd = BD;
    BFS{SPLidx}.difcor.power = Dspec;
    BFS{SPLidx}.difcor.phase = Dphase;
    BFS{SPLidx}.difcor.peakhz = DspecPeakHz;
    BFS{SPLidx}.difcor.freq = freqax;
    BFS{SPLidx}.SPL = SPL;
    BFS{SPLidx} = add_extra_peaks(BFS{SPLidx});
    clear ITDpos ITDneg meanphase VS Pval syncrate spikerate CD CP;
end
    function [data] = add_extra_peaks(data)
        bBD = data.difcor.bd;
        Xspline = data.xvals.ITDxspline;
        yvals = data.difcor.ratespline-min(data.difcor.ratespline);
        [peakval,peakind] = findpeaks(yvals);
        [troughval,troughind] = findpeaks(-yvals);
        troughval = -troughval;
        if min(peakind)<min(troughind)
            troughind = [1 troughind];
            troughval = [yvals(1) troughval];
        end
        if max(peakind)>max(troughind)
            troughind = [troughind length(yvals)];
            troughval = [troughval yvals(end)];
        end
        Dpeaks = Xspline(peakind);
        Dtroughs = Xspline(troughind);
        keeppeak = zeros(size(peakval));
        keeptrough = zeros(size(troughval));
        for ii = 1:length(peakval)
            clo = (peakval(ii)-troughval(find(troughind<peakind(ii),1,'last')))/...
                (peakval(ii)+troughval(find(troughind<peakind(ii),1,'last')));
            chi = (peakval(ii)-troughval(find(troughind>peakind(ii),1,'first')))/...
                (peakval(ii)+troughval(find(troughind>peakind(ii),1,'first')));
            if clo>0.5
                keeppeak(ii) = 1;
                keeptrough(find(troughind<peakind(ii),1,'last'))=1;
            elseif chi>0.5
                keeppeak(ii) = 1;
                keeptrough(find(troughind>peakind(ii),1,'first'))=1;
            end
        end
        Dpeaks = Dpeaks(logical(keeppeak));
        data.difcor.troughs = Dtroughs(logical(keeptrough));
        data.difcor.peaks = Dpeaks(Dpeaks~=bBD);
    end
end