function [] = simulate_CPCD(SPL)

% Simulate cross correlation of fiber pairs and get CP as a function of tuning difference and reference CF

currdir = cd;
datadir = 'C:\work\BinRev\MonRev';
savedir = 'C:\work\BinRev\Analyses';
cd (datadir);
datafiles = dir('*_MonRev.mat');
nrdata = numel(datafiles);
count = 0;
h = waitbar(0,'Please wait...');
spec_crit = 3;
for i =1:nrdata
    h1indx = [];
    data = load (datafiles(i).name);
    ParamsOut = data.ParamsOut;clear data;
    annum(i,:) = str2num(datafiles(i).name(2:6));
    if isfield(ParamsOut,'h1')
        if isfield(ParamsOut,'noiseSPLs')
            h1indx = find(unique(ParamsOut.noiseSPLs)==SPL);
        else
            if SPL==70
                h1indx = 1;
            end
        end
        if isempty(h1indx)
            h1mat(:,i) = nan(1803,1);
            h1phmat(:,i) = nan(4096,1);
            h1magmat(:,i) = nan(4096,1);
            h1latency(i) = nan;
            h1z(:,i) = nan(4096,1);
            df(:,i) = nan(1803,1);
            continue;
        end
        Time = ParamsOut.Time;
        FreqAx = ParamsOut.h1ffax;
        h1mat(:,i) = ParamsOut.h1filt(:,h1indx)./max(abs(ParamsOut.h1filt(:,h1indx)));
        h1magmat(:,i) = ParamsOut.h1mag(:,h1indx);
        h1phmat(:,i) = ParamsOut.h1phase(:,h1indx);
        h1latency(i) = ParamsOut.h1latency(h1indx);
        h1z(:,i) = ParamsOut.h1zscore(:,h1indx);
        df(:,i) = ones(length(ParamsOut.h1filt(:,h1indx)),1).*ParamsOut.domfreq(h1indx);
    else
        h1mat(:,i) = nan(length(h1mat),1);
        h1phmat(:,i) = nan(length(h1phmat),1);
        h1magmat(:,i) = nan(length(h1magmat),1);
        h1latency(i) = nan;
        h1z(:,i) = nan(length(h1z),1);
        df(:,i) = nan(length(df),1);
    end
    cf(:,i) = ones(length(df),1).*ParamsOut.TC.thr.cf;
    waitbar(i/nrdata);
end
close (h);

cd (currdir);
%% Tidy up the CF/DF axis - to take account of "DF" being much lower than CF for fibers around 3-4 kHz.
df = df*1000;
FreqX = df;
FreqX(cf>2000) = cf(cf>2000);

fffind = find(~isnan(FreqX(1,:)) & FreqX(1,:)<=2000);
fff = FreqX(1,fffind);
[fff,sortind] = sort(fff);
fffind =fffind(sortind);
cff = cf(1,fffind);
phase = h1phmat(:,fffind);
zscore = h1z(:,fffind);
mag = h1magmat(:,fffind);
h1time = h1mat(:,fffind);
weights = zeros(length(mag(:,1)),2);
animalnum = annum(fffind,:);
dt = Time(2)-Time(1);
count = 1;
h = waitbar(0,'Please wait...');
for i = 1:length(fff)
    if fff(i)<=2000;
        indx = find(abs(log2(fff/fff(i)))<=0.5);
        indx = indx(indx~=i);
        mph(:,1) = phase(:,i);
        [maxval(1),maxind(1)] = max(mag(:,i));
        DF(1) = FreqAx(maxind(1));
        startind(1) = find(zscore(1:maxind-1,i)>=spec_crit,1,'first');
        stopind(1) = find(zscore(maxind+1:end,i)>=spec_crit,1,'last')+maxind(1);
        weights(startind(1):stopind(1),1) = 10.^(mag(startind(1):stopind(1),i)./20);
        weights(:,1) = weights(:,1)./max(weights(:,1));
        weights(find(zscore(startind(1):stopind(1),i)<spec_crit)+startind(1)-1,1)=0;
        for j = 1:length(indx)
            if animalnum(i,:)==animalnum(indx(j),:)
                comptype(count) = 1;
            else
                comptype(count) = 2;
            end
            mph(:,2) = phase(:,indx(j));
            [maxval(2),maxind(2)] = max(mag(:,indx(j)));
            DF(2) = FreqAx(maxind(2));
            startind(2) = find(zscore(1:maxind(2)-1,indx(j))>=spec_crit,1,'first');
            stopind(2) = find(zscore(maxind(2)+1:end,indx(j))>=spec_crit,1,'last')+maxind(2);
            weights(startind(2):stopind(2),2) = 10.^(mag(startind(2):stopind(2),indx(j))./20);
            weights(:,2) = weights(:,2)./max(weights(:,2));
            weights(find(zscore(startind(2):stopind(2),indx(j))<spec_crit)+startind(2)-1,2)=0;
            
            startat = max(startind);
            stopat = min(stopind);
            BPhase_cycles = (mph(startat:stopat,2)-mph(startat:stopat,1))./(2*pi);
            Bweights = sqrt(prod(weights(startat:stopat,:),2));
            Bweights = Bweights/max(Bweights);
            [dummy,ind]=max(Bweights);
            ind1 = find(Bweights(1:ind-1)==0,1,'last')+1;
            ind2 = find(Bweights(ind+1:end)==0,1,'first')+ind-1;
            if isempty(ind1)
                ind1 = 1;
            end
            if isempty(ind2)
                ind2 = length(Bweights);
            end
            Bweights(1:ind1-1)=0;
            Bweights(ind2+1:end)=0;
            Freq = FreqAx(startat:stopat);
            [CP(count),CD(count),dummy1,dummy2] = fit_circ_lin(Freq,BPhase_cycles,Bweights');
            [BD(count),LPD(count),UPD(count)] = local_get_peaks(h1time(:,i),h1time(:,indx(j)));
            [DomFreq(count)] = local_difcor_spec(xcorr(h1time(:,i),h1time(:,indx(j))));
            CDcycles(count) = CD(count)*DomFreq(count);
            [deltaF(count)] = local_get_deltaF(mag(:,i),zscore(:,i),mag(:,indx(j)),zscore(:,indx(j)));
            count = count+1;
        end
    end
    waitbar(i/length(fff));
end
[CDCPinter,CDCPslope,dum1,dum2] = fit_circ_lin(CP(deltaF>=0)',CDcycles(deltaF>=0),ones(size(CDcycles(deltaF>=0)))');
CDCPinter = -CDCPinter;
CDCPslope = -CDCPslope;
close (h);
save(fullfile(savedir,['AN_CD_CP_' num2str(SPL)]),'DomFreq','deltaF','CD','CP','CDcycles','BD','LPD','UPD','comptype','CDCPinter','CDCPslope');

    function [DFdifcor] = local_difcor_spec(IN)
        diffcorsr = 1000/dt;
        NFFT = 2^21;
        nsam = length(IN);
        DiffCorSpecFreq = diffcorsr*linspace(0,0.5,NFFT/2)/1000; % freq in kHz
        
        % compute & apply hann window, compute complex fft spectrum
        hanwin = hann(nsam);
        
        PredDifSpec = fft(IN.*hanwin,NFFT); % complex spec after windowing
        % magnitude ->power
        PredDifMagSpec = abs(PredDifSpec(1:NFFT/2)).^2;
        [dum,Pind] = max(PredDifMagSpec);
        DFdifcor = DiffCorSpecFreq(Pind);
    end
    function [bd,lpd,upd] = local_get_peaks(ipsi,contra)
        yvals = xcorr(ipsi,contra,'coeff');
        xvals = -sort([(1:floor(numel(yvals)/2))*-dt 0 (1:floor(numel(yvals)/2))*dt]);
        [dummy,bdind] = max(yvals);
        bd = xvals(bdind);
        [dummy,updind] = findpeaks(smooth(yvals(1:bdind),5),'sortstr','descend','minpeakheight',0.2);
        if ~isempty(updind)
            upd = xvals(updind(1));
        else
            upd = nan;
        end
        [dummy,lpdind] = findpeaks(smooth(yvals(bdind:end),5),'sortstr','descend','minpeakheight',0.2);
        if ~isempty(lpdind)
            lpd = xvals(lpdind(1)+bdind-1);
        else
            lpd = nan;
        end
    end

    function [deltaout] = local_get_deltaF(ipsimag,ipsiz,contramag,contraz)
        ipsimag(ipsiz>=spec_crit) = ipsimag(ipsiz>=spec_crit)-min(ipsimag(ipsiz>=spec_crit));
        ipsimag(ipsiz<spec_crit)=nan;ipsimag(isnan(ipsimag))=0;
        scaled(:,1) = ipsimag/max(ipsimag);
        contramag(contraz>=spec_crit) = contramag(contraz>=spec_crit)-min(contramag(contraz>=spec_crit));
        contramag(contraz<spec_crit)=nan;contramag(isnan(contramag))=0;
        scaled(:,2) = contramag/max(contramag);
        
        [c,lags] = xcorr(scaled(:,1),scaled(:,2));
        auto1 = xcorr(scaled(:,1),0);
        auto2 = xcorr(scaled(:,2),0);
        normc = c/sqrt(auto1*auto2); % normalize cross correlation
        [val indc] = max(normc);
        difff = FreqAx(2)-FreqAx(1); % Get increment of freqUp
        xcor.dfreqs = lags*difff; % Convention: positive frequency differences mean ipsi is tuned higher than contra.
        xcor.cor = normc;
        xcor.indPeak = indc;
        xcor.freqDif = xcor.dfreqs(indc);
        
        %reference point for octave comparison is the difcor DF
        if xcor.freqDif>0 % Ipsi > Contra
            deltaout = 2*-log2((DomFreq(count) - xcor.freqDif)/DomFreq(count)); % log(<1) = -
        elseif xcor.freqDif<0% Ipsi < Contra
            deltaout = 2*log2((DomFreq(count) - abs(xcor.freqDif))/DomFreq(count)); % log(>1) = +
        else
            deltaout = 0;
        end
    end
end