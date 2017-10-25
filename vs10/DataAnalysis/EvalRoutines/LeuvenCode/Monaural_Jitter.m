function [] = Monaural_Jitter(Frange,thr)

%Inputs:
%Frange: bandwidth of search (octaves)

Nchan = 2;
RCwin = 30;%Length of reverse-correlation window in ms
Nboot = 20;
cols = {[0 0 1],[0 1 0],[1 0 0],[0 1 1],[1 0 1],[0 0 0],[.8 .8 0],[0 0.5 0],[1 0.5 0],[0.5 0 0.5]};
spec_crit = 3;

ANdata = load('C:\LeuvenDataAnalysis\MonRevPopData.mat');
nrdata = length(ANdata.AN);
for i = 1:nrdata
    ind = find(unique(ANdata.AN(i).noiseSPLs)==70);
    if isempty(ind)
        ANdfs(i) = NaN;
    else
        ANdfs(i) = ANdata.AN(i).df(ind);
    end
    ANcfs(i) = ANdata.AN(i).cf;
end
ANfreqs = ANdfs;
ANfreqs(ANcfs>2000) = ANcfs(ANcfs>2000)/1000;

ind = find(ANfreqs<=2);
for i = 1:length(ind)
    if ismember(111,ANdata.AN(ind(i)).noiseseed)
        if length(unique(ANdata.AN(ind(i)).noiseseed))>1
            type(i) = 1;%We have both 111 and some other noise seed
        else
            type(i) = 2;%We only have 111
        end
    else
        type(i) = 3;%We don't have 111
    end
end

ind = ind(type~=3); %Those fibers that can be used as "ipsi" inputs (and some could also contribute to "contra" if they are type 1)

targetCFs = 0.15*2.^([0:0.1:3.7]); %Roughly 150 Hz to 2000 Hz in 0.1 octave steps

for i = 1:length(targetCFs)
        clear h1;clear nullh1;clear temph1;clear tempnullh1;
        lo = targetCFs(i)*2^(-Frange/2);
        hi = targetCFs(i)*2^(Frange/2);
        CFind = ind(ANdfs(ind)<=hi & ANdfs(ind)>=lo);
        if numel(CFind)>10
            CFind = randsample(CFind,10);
        end
        if numel(CFind)<4
            continue;
        end
        
        [cispikes,wv,dt,h1contra{i},h1ipsi{i}] = Mon_Jitter_get_coincident_spikes(CFind);
        
        if ~isempty(cispikes)
            
            Nwin = round(RCwin/dt);
            
            Nsamples = length(wv);
            for j = 1:length(cispikes)
                cispikes{j} = cispikes{j}(cispikes{j}<2500 & cispikes{j}>50);
                [h1{i,j},tvect,dummy] = local_h1_kernel(cispikes{j});
                h1{i,j} = local_low_pass(h1{i,j});
            end
            %     Bootstrap to get the noise floor for each estimate of the
            %     binaural impulse responses
            for j = 1:length(cispikes)
                nullh1{1} = zeros(Nwin,Nboot);
                nullh1{2} = zeros(Nwin,Nboot);
                for k = 1:Nboot
                    nullspikes = scramble_spikes(cispikes{j});
                    [nullvals,dummy,dummy] = local_h1_kernel(nullspikes);
                    nullh1{1}(:,k) = nullvals(:,1);
                    nullh1{2}(:,k) = nullvals(:,2);
                end
                nullh1{1} = local_low_pass(nullh1{1});
                nullh1{2} = local_low_pass(nullh1{2});
                
                h1time{i,j}(:,1) = h1{i,j}(:,1);
                h1time{i,j}(:,2) = h1{i,j}(:,2);
                
                tempnullh1{1} = nullh1{1};
                tempnullh1{2} = nullh1{2};
                clear nullh1;
                nullh1 = zeros(Nwin,2*Nboot);
                nullh1(:,1:2:(2*Nboot)) = tempnullh1{1};
                nullh1(:,2:2:(2*Nboot)) = tempnullh1{2};
                
                [h1mag{i,j},freqscale,h1phase{i,j},h1zscore{i,j},h1wf{i,j},dummy] = Phase_Freq_analysis(h1time{i,j},nullh1);
                
                for k = 1:2
                    h1mag{i,j}(:,k) = Trifilter(h1mag{i,j}(:,k)',13)';
                end
                
                [h1bw{i,j},qvals{i,j},domfreq{i,j},dummy] = Bandwidth_analysis(h1mag{i,j},h1zscore{i,j},spec_crit,freqscale);
                
                clear nullh1;
                
                %Get the predicted difcor ("bincor") for
                %the pseudo-binaural output
                [bincor_bin{i,j},BD_bin(i,j),ITD,bincor_power{i,j},bincor_domfreq(i,j),bincor_freqax] = local_get_bincor(h1wf{i,j}(:,1),h1wf{i,j}(:,2));
                bincor_bin{i,j} = bincor_bin{i,j}/max(abs(bincor_bin{i,j}));
                
                
                %Get the CD and CP for the pseudo-binaural revcors
                weights = zeros(size(h1mag{i,j}));
                for k = 1:2
                    mph(:,k) = h1phase{i,j}(:,k);
                    [maxval(k),maxind(k)] = max(h1mag{i,j}(:,k));
                    DF(k) = freqscale(maxind(k));
                    if length(find(h1zscore{i,j}(1:maxind-1,k)>=spec_crit))>=3 && length(find(h1zscore{i,j}(maxind+1:end,k)>=spec_crit))>=3
                        startind(k) = find(h1zscore{i,j}(1:maxind-1,k)>=spec_crit,1,'first');
                        stopind(k) = find(h1zscore{i,j}(maxind+1:end,k)>=spec_crit,1,'last')+maxind(1);
                    else
                        startind(k) = NaN;
                        stopind(k) = NaN;
                        break;
                    end
                    weights(startind(k):stopind(k),k) = 10.^(h1mag{i,j}(startind(k):stopind(k),k)./20);
                    weights(:,k) = weights(:,k)./max(weights(:,k));
                    weights(find(h1zscore{i,j}(startind(k):stopind(k),k)<spec_crit)+startind(k)-1,k)=0;
                end
                if any(isnan([startind stopind]))
                    continue;
                end
                startat = max(startind);
                stopat = min(stopind);
                BPhase_cycles = (mph(startat:stopat,2)-mph(startat:stopat,1))./(2*pi);
                Bweights = sqrt(prod(weights(startat:stopat,:),2));
                Bweights = Bweights/max(Bweights);
                [dummy,Bweightind]=max(Bweights);
                ind1 = find(Bweights(1:Bweightind-1)==0,1,'last')+1;
                ind2 = find(Bweights(Bweightind+1:end)==0,1,'first')+Bweightind-1;
                if isempty(ind1)
                    ind1 = 1;
                end
                if isempty(ind2)
                    ind2 = length(Bweights);
                end
                Bweights(1:ind1-1)=0;
                Bweights(ind2+1:end)=0;
                Freq = freqscale(startat:stopat)/1000;%Frequency in kHz!
                [CP(i,j),CD(i,j),dummy1,dummy2] = fit_circ_lin(Freq,BPhase_cycles,Bweights');
                CDcycles(i,j) = CD(i,j)*bincor_domfreq(i,j);
            end
            
            %Draw the results
            figure(1);
            subplot(2,2,1);cla; %Magnitude spectrum - ipsi
            for j = 1:length(h1ipsi{i}.h1ipsi)
                yvals = h1ipsi{i}.h1ipsimag{j};
                yvals(h1ipsi{i}.h1ipsizscore{j}<spec_crit) = NaN;
                yvals = yvals-max(yvals);
                plot(freqscale/1000,yvals,'color',cols{j},'linewidth',1);hold on;
                plot(h1ipsi{i}.h1ipsidomfreq{j},0,'o','markeredgecolor',cols{j},'markerfacecolor',cols{j},'linewidth',1,'markersize',10);
            end
            for j = 1:10
                yvals = h1mag{i,j}(:,1);
                yvals(h1zscore{i,j}(:,1)<spec_crit) = NaN;
                yvals = yvals-max(yvals);
                plot(freqscale/1000,yvals,'k-','linewidth',3);
            end
            %Plot the geometric mean of the dominant frequencies of the
            %contralateral input fibers
            averageipsidomfreq = geomean([h1ipsi{i}.h1ipsidomfreq{:}]);
            plot(averageipsidomfreq,1,'v','color',[0.6 0.6 0.6],'markerfacecolor',[.6 .6 .6],'markersize',12);
            
            %Plot the dominant frequency at the output of the binaural revcor
            %for this ear
            for j = 1:10
                plot(domfreq{i,j}(1)/1000,1,'kv','markerfacecolor','k','markersize',12);
            end
            set(gca,'xscale','log','xlim',[0.05 5],'ylim',[-20 2],'fontsize',16,'linewidth',1,'xtick',[0.1 1],'xticklabel',[0.1 1],'layer','top');
            ylabel 'GAIN (dB)';
            xlabel 'FREQUENCY (kHz)';
            box off;
            
            
            subplot(2,2,2);cla; %Magnitude spectrum - contra
            for j = 1:size(h1contra{i}.h1contra,1)
                yvals = h1contra{i}.h1contramag{j};
                yvals(h1contra{i}.h1contrazscore{j}<spec_crit) = NaN;
                yvals = yvals-max(yvals);
                plot(freqscale/1000,yvals,'color',cols{j},'linewidth',1);hold on;
                plot(h1contra{i}.h1contradomfreq{j},0,'o','markeredgecolor',cols{j},'markerfacecolor',cols{j},'linewidth',1,'markersize',10);
            end
            for j = 1:10
                yvals = h1mag{i,j}(:,2);
                yvals(h1zscore{i,j}(:,2)<spec_crit) = NaN;
                yvals = yvals-max(yvals);
                plot(freqscale/1000,yvals,'k-','linewidth',3);
            end
            %Plot the geometric mean of the dominant frequencies of the
            %contralateral input fibers
            averagecontradomfreq = geomean([h1contra{i}.h1contradomfreq{:}]);
            plot(averagecontradomfreq,1,'v','color',[0.6 0.6 0.6],'markerfacecolor',[.6 .6 .6],'markersize',12);
            
            %Plot the dominant frequency at the output of the binaural revcor
            %for this ear
            for j = 1:10
                plot(domfreq{i,j}(2)/1000,1,'kv','markerfacecolor','k','markersize',12);
            end
            set(gca,'xscale','log','xlim',[0.05 5],'ylim',[-20 2],'fontsize',16,'linewidth',1,'xtick',[0.1 1],'xticklabel',[0.1 1],'layer','top');
            ylabel 'GAIN (dB)';
            xlabel 'FREQUENCY (kHz)';
            box off;
            
            subplot(2,2,3);cla; %Time domain - ipsi
            for j = 1:length(h1ipsi{i}.h1ipsi)
                yvals = h1ipsi{i}.h1ipsi{j};
                yvals = yvals/max(abs(yvals));
                plot(tvect,yvals,'color',cols{j},'linewidth',1);hold on;
            end
            for j = 1:10
                yvals = h1wf{i,j}(:,1);
                yvals = yvals/max(abs(yvals));
                plot(tvect,yvals,'k--','linewidth',2);
            end
            box off
            set(gca,'xlim',[0 15],'ylim',[-1.2 1.2],'fontsize',16,'linewidth',1,'layer','top');
            xlabel 'TIME (ms)';
            ylabel 'AMPLITUDE';
            
            subplot(2,2,4);cla; %Time domain - contra
            for j = 1:size(h1contra{i}.h1contra,1)
                for k = 1:size(h1contra{i}.h1contra,2)
                    yvals = h1contra{i}.h1contra{j,k};
                    yvals = yvals/max(abs(yvals));
                    plot(tvect,yvals,'color',cols{j},'linewidth',1);hold on;

                end
            end
            for j = 1:10
                yvals = h1wf{i,j}(:,2);
                yvals = yvals/max(abs(yvals));
                plot(tvect,yvals,'k--','linewidth',2);
            end
            
            box off
            set(gca,'xlim',[0 15],'ylim',[-1.2 1.2],'fontsize',16,'linewidth',1,'layer','top');
            xlabel 'TIME (ms)';
            ylabel 'AMPLITUDE';
            
            
            figure(2);cla;
            plot([0 0],[-1.2 1.2],'--','color',[0.6 0.6 0.6],'linewidth',3);hold on;
            for j = 1:10
                plot(ITD,bincor_bin{i,j},'k-','linewidth',3);hold on;
            end
            box off
            set(gca,'xlim',[-5 5],'ylim',[-1.2 1.2],'fontsize',16,'linewidth',1,'layer','top');
            xlabel 'ITD (ms)';
            ylabel 'AMPLITUDE';
        end
end
return;
    function [bincorout,BDout,ITDout,BinSpec,BinDomFreq,BinCorSpecFreq] = local_get_bincor(IN1,IN2)
        bincorout = xcorr(IN1,IN2);
        ITDout = -sort([(1:floor(numel(bincorout)/2))*-dt 0 (1:floor(numel(bincorout)/2))*dt]);
        [dum,indout] = max(bincorout);
        BDout = ITDout(indout);
        bincorsr = 1/dt;
        NFFT = 2^15;
        nsam = length(bincorout);
        BinCorSpecFreq = bincorsr*linspace(0,0.5,NFFT/2); % freq in kHz
        
        % compute & apply hann window, compute complex fft spectrum
        hanwin = hann(nsam);
        
        BinSpec = fft((bincorout./max(abs(bincorout))).*hanwin,NFFT); % complex spec after windowing
        % magnitude ->power
        BinSpec = abs(BinSpec(1:NFFT/2)).^2;
        
        [dummy,maxbinspecind] = max(BinSpec);
        BinDomFreq = BinCorSpecFreq(maxbinspecind);
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
        D = diff([0; sort(spikes)']); % first sort to make sure diff gives the ISI's
        NInt = length(D);
        scramspikes = cumsum(D(randperm(NInt)));
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
    function [h1magout,freqscale,h1phaseout,zcritout,wf,nullwf] = Phase_Freq_analysis(y,nully)
        wf = zeros(size(y));
        nullwf = zeros(size(nully));
        nfft = 2^13;
        h1fft = fft(y,nfft);
        h1magout = 20*log10(abs(h1fft));
        h1magout = h1magout(1:floor(nfft/2),:);
        nullh1fft = fft(nully,nfft);
        nullh1magout = 20*log10(abs(nullh1fft));
        nullh1magout = nullh1magout(1:floor(nfft/2),:);
        nf = 0.5*(1000/dt);
        freqscale = nf*linspace(0,1,floor(nfft/2));
        h1phaseout = angle(h1fft(1:floor(nfft/2),:));
        %get z-scores for magnitude spectrum
        zcritout = zeros(size(h1magout));
        for ii = 1:Nchan
            zcritout(:,ii) = ((10.^(h1magout(:,ii)/20))-mean(10.^(nullh1magout(ii:2:(2*Nboot))/20),2))./std(10.^(nullh1magout(ii:2:(2*Nboot))/20),1,2);
        end
        %
        %Do the Wiener filtering
        for ii = 1:Nchan
            zz = [zcritout(:,ii); flipud(zcritout(:,ii))];
            h1fft(zz<spec_crit,ii) = 0; %Set noise components to zeros
            vals = ifft(h1fft(:,ii),nfft);
            vals = real(vals);
            wf(:,ii) = vals(1:length(y));
            for jj = ii:2:size(nully,2)
                nullh1fft(zz<spec_crit,jj) = 0; %Set noise components to zeros
                vals = ifft(nullh1fft(:,jj),nfft);
                vals = real(vals);
                nullwf(:,jj) = vals(1:length(y));
            end
        end
        
        %Apply smoothing (Triangular smoothing function)
        for ii = 1:Nchan
            [zcritout(:,ii)] = Trifilter(zcritout(:,ii)',7)';
        end
    end
    function [h1bw,qvals,domfreq,domind] = Bandwidth_analysis(magIN,zIN,critIN,fIN)
        %Returns the 3 and 6-dB bandwidths of the h1 kernels (in kHz), the dominant
        %frequency (in kHz) and the index of that value
        Nchan = 2;
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
            hivals = magIN(maxind:length(magIN),ii);
            flo = fIN(1:maxind);
            fhi = fIN(maxind:length(magIN));
            zlo = zIN(1:maxind,ii);
            zhi = zIN(maxind:length(magIN),ii);
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













