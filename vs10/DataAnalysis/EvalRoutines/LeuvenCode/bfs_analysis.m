function [BFS] = bfs_analysis(dsbfs,params,ITDx)

dt = params.dt;
correctITD = params.correctITD;
ITDx = sort(ITDx,'ascend');
ditd = abs(diff(ITDx([1,2])));
Ndata = length(dsbfs);
stimdur = dsbfs{1}.Stimulus.Special.BurstDur(1);%stimulus length (in millseconds)
beatF = dsbfs{1}.Stimulus.Special.BeatFreq(1);%beat frequency (Hz)
beatP = 1e3/beatF;%beat period (in milliseconds)
minITD = min(ITDx);
maxITD = max(ITDx);
count = 1;
%Gather all spikes from all reps of all conditions in all datasets
for i = 1:Ndata
    spls(i) = dsbfs{i}.Stimulus.StimParam.indiv.stim{1}.spl;
    thesefreq = dsbfs{i}.Stimulus.IndepVar.Values;
    inds = find(~isnan(thesefreq));
    nrreps = size(dsbfs{i}.Data.SpikeTimes,2);
    for j = 1:length(inds)
        for k = 1:nrreps
            spikes{count} = dsbfs{i}.Data.SpikeTimes{j,k};
            freq(count) = thesefreq(inds(j));
            spl(count) = spls(i);
            count = count+1;
        end
    end
end

uSPLs = unique(spls);
NSPLs = length(uSPLs);

for j = 1:NSPLs
    %Re-sort data so that all reps from a single condition are together,
    %regardless of their dataset
    ufreq = unique(freq(spl==uSPLs(j)));
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
        carrP = 1e3/ufreq(i);
        ind = find(freq==ufreq(i) & spl==uSPLs(j));
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

        
        %correct ITD for sgsr conventions
        if correctITD
            SCpos = 1-SCpos;
            SCneg = 1-SCneg;
        end
        
        %convert cycles to ITD - simulate 5 beat cycles (so that it is
        %consistent with the noise data)
        STpos = [SCpos(SCpos<=0.5) -1+SCpos(SCpos>0.5)];
        STneg = [SCneg(SCneg<=0.5) -1+SCneg(SCneg>0.5)];
        
        STpos = bsxfun(@plus,STpos,(-5:5)');
        STneg = bsxfun(@plus,STneg,(-5:5)');
        
        STpos = STpos(:)*carrP;
        STneg = STneg(:)*carrP;
        
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
    
    %Do the circular-linear regression on the freq-phase data to get CP and CD
    weights = syncrate(Pval<0.001);
    weights = weights/max(weights);
    [CP,CD,~,~] = fit_circ_lin(ufreq(Pval<0.001)/1000,meanphase(Pval<0.001),weights');
    
    %Get the difcor
    difcorrate = nanmean(ITDpos)-nanmean(ITDneg);
    
    %spline fit
    ITDxspline = min(ITDx):dt:max(ITDx);
    Dspline = spline(ITDx,difcorrate,ITDxspline);
    
    %Get the BD
    [~,ind] = max(Dspline);
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
    [~,ind] = max(Dspec);
    DspecPeakHz = freqax(ind);
    
    
    %Make the output structure
    BFS{j}.beat.CP = CP;
    BFS{j}.beat.CD = CD;
    BFS{j}.beat.meanphase = meanphase;
    BFS{j}.beat.ipsiphase = meanphaseipsi;
    BFS{j}.beat.contraphase = meanphasecontra;
    BFS{j}.beat.VS = VS;
    BFS{j}.beat.VSipsi = VSipsi;
    BFS{j}.beat.VScontra = VScontra;
    BFS{j}.beat.Pval = Pval;
    BFS{j}.beat.Pvalipsi = Pvalipsi;
    BFS{j}.beat.Pvalcontra = Pvalcontra;
    BFS{j}.beat.syncrate = syncrate;
    BFS{j}.beat.ufreq = ufreq;
    BFS{j}.beat.freq = freq(spl==uSPLs(j));
    BFS{j}.beat.spiketimes = {spikes{spl==uSPLs(j)}};
    BFS{j}.positive.rate = ITDpos;
    BFS{j}.negative.rate = ITDneg;
    BFS{j}.xvals.ITDx = ITDx;
    BFS{j}.difcor.rate = difcorrate;
    BFS{j}.xvals.ITDxspline = ITDxspline;
    BFS{j}.difcor.ratespline = Dspline;
    BFS{j}.difcor.bd = BD;
    BFS{j}.difcor.power = Dspec;
    BFS{j}.difcor.phase = Dphase;
    BFS{j}.difcor.peakhz = DspecPeakHz;
    BFS{j}.difcor.freq = freqax;
    BFS{j}.SPL = uSPLs(j);
    BFS{j} = add_extra_peaks(BFS{j});
    clear ITDpos ITDneg meanphase VS Pval syncrate spikerate...
        meanphaseipsi meanphasecontra VSipsi VScontra Pvalipsi Pvalcontra; 
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