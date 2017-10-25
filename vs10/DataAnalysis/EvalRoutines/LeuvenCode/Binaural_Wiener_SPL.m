function [] = Binaural_Wiener_SPL(varargin)

minST = 50;
maxST = 5000;
RCwin = 30;%Length of reverse-correlation window in ms
Nboot = 20;

correctITD = 1;
AnNum = varargin{1};
UnNum = varargin{2};


ndndata = [];
ndpdata = [];
ndzdata=[];
for i = 3:nargin
    if strcmp(varargin{i-1},'ndp')
        ndpdata = varargin{i};
    elseif strcmp(varargin{i-1},'ndn')
        ndndata = varargin{i};
    elseif strcmp(varargin{i-1},'ndz')
        ndzdata = varargin{i};
    end
end
if ~isempty(ndzdata)
    if iscell(ndzdata)
        for i = 1:length(ndzdata)
            dsz{i} = dataset(AnNum,ndzdata{i});
        end
    else
        dsz{1} = dataset(AnNum,ndzdata);
    end
    nrdata = length(dsz);
end
if ~isempty(ndndata)
    if iscell(ndndata)
        for i = 1:length(ndndata)
            dsn{i} = dataset(AnNum,ndndata{i});
        end
    else
        dsn{1} = dataset(AnNum,ndndata);
    end
end
if ~isempty(ndpdata)
    if iscell(ndpdata)
        for i = 1:length(ndpdata)
            dsp{i} = dataset(AnNum,ndpdata{i});
        end
    else
        dsp{1} = dataset(AnNum,ndpdata);
    end
end


[dummy, dt] = StimSam(dsz{1},1);
dt = 1e-3*dt; %us -> ms
Nwin = round(RCwin/dt);

%Get the SPL for each revcor dataset
for i = 1:length(dsz)
    SPL(i) = dsz{i}.Stimulus.StimParam.SPL(1);
end
uSPLs = unique(SPL);
NrSPLs = length(uSPLs);
%Tell the user something
fprintf(1,'Found %1.0f revcor datasets at %1.0f SPLs\n',length(dsz),NrSPLs);

%Get the measured noise delay function
NTD_SPL = dsp{1}.Stimulus.StimParam.SPL(1);
maxST = dsp{1}.Stimulus.Special.BurstDur(1); %Max spike time to consider for noise-delay data
ITDx = dsp{1}.Stimulus.IndepVar.Values./1e3; %ITD axis - us->ms
if correctITD
    ITDx = flipud(ITDx);
end
NrITD = length(ITDx);
NTD_pos = zeros(NrITD,1);NTD_neg = zeros(NrITD,1);NTD_difcor = zeros(NrITD,1);
for i = 1:NrITD
    NTD_pos(i) = numel(dsp{1}.SPT{i}(dsp{1}.SPT{i}<maxST))./(maxST/1e3);
    NTD_neg(i) = numel(dsn{1}.SPT{i}(dsn{1}.SPT{i}<maxST))./(maxST/1e3);
    NTD_difcor(i) = NTD_pos(i) - NTD_neg(i);
end
ITDx_spline = linspace(min(ITDx),max(ITDx),1000);
NTD_pos_spline = spline(ITDx,NTD_pos,ITDx_spline);
NTD_neg_spline = spline(ITDx,NTD_neg,ITDx_spline);
NTD_difcor_spline = NTD_pos_spline - NTD_neg_spline;

%Get the binaural revcors at each unique SPL
for u = 1:NrSPLs
    nullh1 = zeros(Nwin,Nboot*2);
    SPLind = find(SPL==uSPLs(u));
    nrdata = length(SPLind);
    NrReps = 0;
    for i = 1:nrdata
        [wv, dt] = StimSam(dsz{SPLind(i)},1);
        dt = 1e-3*dt; %us -> ms
        Nchan = size(wv, 2);
        Nsamples = length(wv);
        spikes = dsz{SPLind(i)}.SPT;%get spike times
        NrReps = NrReps+size(spikes,2);
        spikes = cat(2,spikes{:});%concatenate all spikes from all reps
        spikes = spikes(spikes>=minST & spikes<=maxST);%limit spike times to only include driven spikes
        fprintf('Calculating %2.0f of %2.0f revcors...\n',i,nrdata);
        [h1{i},Time,SpikeCount(i)] = local_h1_kernel(spikes);
        
        [newh1,recside,channelflipper] = mapchannel2side(dsz{i},h1{i});
        if ~isempty(newh1)
            h1{i} = newh1;
        else
            disp 'Unable to figure out contra/ipsi channel!!!!';
            %First assume recording was on right side and that you connected the left
            %(red) DAC to the left ear.
            if strcmp(AnNum,'L14593')
                recside = 'Left';
                channelflipper = 0;
                correctITD = 0;
            else
                if contrachan == 1 || strcmp(AnNum,'L5004'); %Exception here for this animal
                    recside = 'Right';
                    channelflipper = 0;
                else
                    recside = 'Left';
                    channelflipper = 0;
                end
            end
        end
        
        %Now we want to put the ipsi channel in the left column and the contra
        %in the right column
        if strcmp(recside,'Right')
            isflipped = 1;
            h1{i} = h1{i}(:,[2 1]);%reverse [left(contra) right(ipsi)] to make [right(ipsi) left(contra)]
        else
            isflipped = 0;
            %Leave channels as they are - left (ipsi) is already on the left
        end
        
        fprintf('Done %2.0f of %2.0f revcors at %2.0f dB SPL!\n',i,nrdata,uSPLs(u));
        
        h = waitbar(0,sprintf('Bootstrapping %2.0f of %2.0f revcors: Please wait...',i,nrdata));
        for j = 1:Nboot
            nullspikes = scramble_spikes(spikes)';
            nullspikes = nullspikes(nullspikes>=minST);
            [nullvals,dummy] = local_h1_kernel(nullspikes);
            if isflipped
                nullvals = nullvals(:,[2 1]);
            end
            nullh1(:,(j*2)-1:j*2) = nullh1(:,(j*2)-1:j*2)+nullvals;
            waitbar(j/Nboot);
        end
        close (h);
    end
    TotalSpikes(u) = sum(SpikeCount);
    MeanSpikeRate(u) = TotalSpikes(u)/NrReps/((maxST-minST)/1000);
    h1mat = [h1{:}];%Put all the estimates into a matrix
    clear h1;
    h1spl{u}(:,1) = mean(h1mat(:,1:2:end),2);
    h1spl{u}(:,2) = mean(h1mat(:,2:2:end),2);
    clear h1mat;
    nullh1spl{u} = nullh1/nrdata;
    clear nullh1;
    clear SpikeCount;
end


for i = 1:NrSPLs
    normfactor(i) = max(max(abs(h1spl{i})));
    for j = 1:2
        h1spl{i}(:,j) = local_low_pass(h1spl{i}(:,j));
        h1spl{i}(:,j) = h1spl{i}(:,j)/normfactor(i);
    end
    for j = 1:Nboot*2
        nullh1spl{i}(:,j) = local_low_pass(nullh1spl{i}(:,j));
        nullh1spl{i}(:,j) = nullh1spl{i}(:,j)/normfactor(i);
    end
    [h1mag{i},nullh1mag{i},freqs,h1phase{i},h1zscore{i}] = Phase_Freq_analysis(h1spl{i},nullh1spl{i});
    for j = 1:2
        h1mag{i}(:,j) = Trifilter(h1mag{i}(:,j)',13)';
    end
end
spec_crit = 5; %Z-score(f) criterion for "significance" in the frequency domain
ffax = freqs/1e3; %Hz -> kHz
for i = 1:NrSPLs
    [h1bw{i},qvals{i},domfreq{i},domind{i}] = Bandwidth_analysis(h1mag{i},h1zscore{i},spec_crit,ffax);
end

% %Get the CD and CP from the monaural phase-frequency
% %functions
clear ind;
for i = 1:NrSPLs
    weights{i} = zeros(size(h1phase{i}));
    for j = 1:2
        ind(i,j) = find(ffax==domfreq{i}(j));
        startind(i,j) = find(h1zscore{i}(1:ind(i,j)-1,j)>=5,1,'first')+1;
        stopind(i,j) = find(h1zscore{i}(ind(i,j)+1:end,j)>=5,1,'last')+ind(i,j)-1;
        weights{i}(startind(i,j):stopind(i,j),j) = 10.^(h1mag{i}(startind(i,j):stopind(i,j),j)/20);
        weights{i}(:,j) = weights{i}(:,j)/max(weights{i}(:,j));
        weights{i}(find(h1zscore{i}(startind(i,j):stopind(i,j),j)<5)+startind(i,j)-1,j)=0;
    end
end
for i = 1:NrSPLs
    startat = max(startind(i,:));
    stopat = min(stopind(i,:));
    Phase = [h1phase{i}(startat:stopat,1) h1phase{i}(startat:stopat,2)];
    Freq = ffax(startat:stopat);
    BPhase_cycles = (Phase(:,2)-Phase(:,1))./(2*pi);
    Bweights = prod([weights{i}(startat:stopat,1) weights{i}(startat:stopat,2)],2);
    Bweights = Bweights/max(Bweights);
    [dummy,indx]=max(Bweights);
    ind1 = find(Bweights(1:indx-1)==0,1,'last')+1;
    ind2 = find(Bweights(indx+1:end)==0,1,'first')+indx-1;
    if isempty(ind1)
        ind1 = 1;
    end
    if isempty(ind2)
        ind2 = length(Bweights);
    end
    Bweights(1:ind1-1)=0;
    Bweights(ind2+1:end)=0;
    [CP(i),CD(i),dummy1,dummy2] = fit_circ_lin(Freq,BPhase_cycles,Bweights');
end


%Get the xcorr at each SPL
for i = 1:NrSPLs
    difcor(:,i) = xcorr(h1spl{i}(:,1),h1spl{i}(:,2),'coeff');
    difcor_df(i) = local_difcor_spec(difcor(:,i));
end
itd = -[fliplr(-Time) Time(2:end)];

%scale the difcors for plotting
for j = 1:NrSPLs
    difcor(:,j) = ((2*j)-2)+difcor(:,j)./max(abs(difcor(:,j)));
    [bd(j).y,bd(j).ind] = max(difcor(:,j));
    bd(j).x = itd(bd(j).ind);
end
BD = [bd(:).x];


%Plotting
cols = {[0 0 1],[0 1 0],[1 0 0],[0 1 1],[1 0 1],[0 0 0],[.8 .8 0],[0 0.5 0],[1 0.5 0],[0.5 0 0.5]};
figure;
subplot(2,4,[1 5]);
for i = 1:NrSPLs
    plot(Time,h1spl{i}(:,1)+(i*2)-2,'b','linewidth',2);hold on;
    plot(Time,h1spl{i}(:,2)+(i*2)-2,'r','linewidth',2);
end
set(gca,'xlim',[0 20],'ylim',[-1.5 2*(NrSPLs-1)+1.5],'fontsize',16,'ytick',[0:2:2*(NrSPLs-1)+1],'yticklabel',uSPLs,'linewidth',1,'layer','top');
box off;
xlabel 'TIME (ms)';
ylabel 'SOUND LEVEL (dB SPL)';

subplot(2,4,[2 6]);
plot([0 0],[-1.5 2*(NrSPLs-1)+1.5],'--','color',[0.6 0.6 0.6]);hold on;
plot(itd,difcor,'k','linewidth',2);hold on;
n = 1;
for i = 1:NrSPLs
    plot(bd(i).x, bd(i).y,'marker','o','markersize',10,'markerfacecolor',cols{n},'markeredgecolor', cols{n});
    n = n+1;
    if n >10
        n = 1;
    end
end
set(gca,'xlim',[-5 5],'ylim',[-1.5 2*(NrSPLs-1)+1.5],'fontsize',16,'ytick',[0:2:2*(NrSPLs-1)+1],'yticklabel',uSPLs,'linewidth',1,'layer','top');
box off;
xlabel 'ITD (ms)';
ylabel 'SOUND LEVEL (dB SPL)';

%Dominant frequency of monaural revcors and binaural difcor as a
%fuction of SPL
subplot(2,4,3);
h1df = vertcat(domfreq{:});
plot(uSPLs,h1df(:,1),'b^','markerfacecolor','none','markersize',10,'linewidth',1);hold on;
plot(uSPLs,h1df(:,2),'rv','markerfacecolor','none','markersize',10,'linewidth',1);hold on;
n = 1;
for i = 1:NrSPLs
    plot(uSPLs(i),difcor_df(i),'marker','o','markerfacecolor',cols{n},'markeredgecolor',cols{n},'markersize',10);
    n = n+1;
    if n >10
        n = 1;
    end
end
set(gca,'xlim',[min(uSPLs)-5 max(uSPLs)+5],'ylim',[0.2 0.4],'fontsize',16,'xtick',uSPLs(1:2:end),'xticklabel',uSPLs(1:2:end),'linewidth',1,'layer','top');
box off;
ylabel 'DOMINANT FREQUENCY (kHz)';
xlabel 'SOUND LEVEL (dB SPL)';

%BD and BD*DF as a function of SPL
subplot(2,4,7);
n = 1;
for i = 1:NrSPLs
    plot(uSPLs(i),bd(i).x,'marker','o','markerfacecolor',cols{n},'markeredgecolor',cols{n},'markersize',10);hold on;
    n = n+1;
    if n >10
        n = 1;
    end
end
set(gca,'xlim',[min(uSPLs)-5 max(uSPLs)+5],'ylim',[0 1],'fontsize',16,'xtick',uSPLs(1:2:end),'xticklabel',uSPLs(1:2:end),'linewidth',1,'layer','top');
box off;
ylabel 'BEST DELAY (ms)';
xlabel 'SOUND LEVEL (dB SPL)';

%Characteristic Phase as a function of SPL
subplot(2,4,4);
n = 1;
for i = 1:NrSPLs
    plot(uSPLs(i),CP(i),'marker','o','markerfacecolor',cols{n},'markeredgecolor',cols{n},'markersize',10);hold on;
    n = n+1;
    if n >10
        n = 1;
    end
end
set(gca,'xlim',[min(uSPLs)-5 max(uSPLs)+5],'ylim',[-0.5 0.5],'fontsize',16,'xtick',uSPLs(1:2:end),'xticklabel',uSPLs(1:2:end),'linewidth',1,'layer','top');
box off;
ylabel 'CHARACTERISTIC PHASE (cycles)';
xlabel 'SOUND LEVEL (dB SPL)';

%Characteristic Delay as a function of SPL
subplot(2,4,8);
n = 1;
for i = 1:NrSPLs
    plot(uSPLs(i),CD(i),'marker','o','markerfacecolor',cols{n},'markeredgecolor',cols{n},'markersize',10);hold on;
    n = n+1;
    if n >10
        n = 1;
    end
end
set(gca,'xlim',[min(uSPLs)-5 max(uSPLs)+5],'ylim',[0 1],'fontsize',16,'xtick',uSPLs(1:2:end),'xticklabel',uSPLs(1:2:end),'linewidth',1,'layer','top');
box off;
ylabel 'CHARACTERISTIC DELAY (ms)';
xlabel 'SOUND LEVEL (dB SPL)';

%Save the output
save (fullfile('C:\work\BinRev\Analyses',['MSO_SPL_' AnNum '_' UnNum]),'CD','CP','difcor_df','h1df','h1spl','uSPLs','BD','difcor');


    function [h1bw,qvals,domfreq,domind] = Bandwidth_analysis(magIN,zIN,critIN,fIN)
        %Returns the 3 and 6-dB bandwidths of the h1 kernels (in kHz), the dominant
        %frequency (in kHz) and the index of that value
        h1bw = nan(2,Nchan);
        qvals = nan(2,Nchan);
        domfreq = nan(1,Nchan);
        domind = nan(1,Nchan);
        bwlevels = [3 6];
        for ii=1:Nchan
            [maxval,maxind] = max(magIN(:,ii));
            domfreq(ii) = fIN(maxind);
            domind(ii) = maxind;
            lovals = magIN(1:maxind,ii);
            hivals = magIN(maxind:end,ii);
            flo = fIN(1:maxind);
            fhi = fIN(maxind:end);
            zlo = zIN(1:maxind,ii);
            zhi = zIN(maxind:end,ii);
            for jj = 1:length(bwlevels)
                ind1 = find(lovals<maxval-bwlevels(jj),1,'last');
                ind2 = find(lovals>maxval-bwlevels(jj),1,'first');
                if zlo(ind1)>=critIN && zlo(ind2)>=critIN
                    lointercept = interp1(lovals([ind1 ind2]),flo([ind1 ind2]),maxval-bwlevels(jj));
                else
                    lointercept = nan;
                end
                ind1 = find(hivals<maxval-bwlevels(jj),1,'first');
                ind2 = find(hivals>maxval-bwlevels(jj),1,'last');
                if zhi(ind1)>=critIN && zhi(ind2)>=critIN
                    hiintercept = interp1(hivals([ind1 ind2]),fhi([ind1 ind2]),maxval-bwlevels(jj));
                else
                    hiintercept = nan;
                end
                h1bw(jj,ii) = hiintercept-lointercept;
            end
            qvals(:,ii) = domfreq(ii)./h1bw(:,ii);
        end
    end
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
        % shuffle spike train by randomizing the ISI's
        %         RandStream.setDefaultStream (RandStream('mt19937ar','seed',sum(100*clock)));
        D = diff([0; sort(spikes')]); % first sort to make sure diff gives the ISI's
        NInt = length(D);
        scramspikes = cumsum(D(randperm(NInt)));
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
        %         h1phaseout = unwrap(angle(h1fft(1:floor(nfft/2),:)));
        h1phaseout = angle(h1fft(1:floor(nfft/2),:));
        %get z-scores for magnitude spectrum
        zcritout = zeros(size(h1magout));
        for ii = 1:Nchan
            zcritout(:,ii)= ((10.^(h1magout(:,ii)/20))-mean(10.^(nullh1magout(:,ii:2:end)/20),2))./std(10.^(nullh1magout(:,ii:2:end)/20),1,2);
        end
        %Apply smoothing (Triangular smoothing function)
        for ii = 1:Nchan
            [zcritout(:,ii)] = Trifilter(zcritout(:,ii)',7)';
        end
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







