function [] = Monaural_Convergence(Frange,thr,plotYN)

%Inputs:
%Frange: bandwidth of search (octaves)
%Ctype: type of comparison; within animal (1), or across animal (2)

%Monaural coincidence analysis for ANFs in response to noise. Finds the
%revcor from just coincident spikes from fibers with neighbouring CFs in
%the region of CFin (+/- Frange);

Nchan = 2;
RCwin = 30;%Length of reverse-correlation window in ms
Nboot = 20;
cols = {[0 0 1],[0 1 0],[1 0 0],[0 1 1],[1 0 1],[0 0 0],[.8 .8 0],[0 0.5 0],[1 0.5 0],[0.5 0 0.5]};
spec_crit = 3;


ANdata = load('C:\LeuvenDataAnalysis\MonRevPopData.mat');
nrdata = length(ANdata.AN);
for i = 1:nrdata
    if ~isempty(ANdata.AN(i).noiseSPLs)
        ind = find(unique(ANdata.AN(i).noiseSPLs)==70);
        if isempty(ind)
            ANdfs(i) = NaN;
        else
            ANdfs(i) = ANdata.AN(i).df(ind);
        end
    else
        if ~isempty(ANdata.AN(i).df)
            ANdfs(i) = ANdata.AN(i).df;
        else
            ANdfs(i) = NaN;
        end
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

indA = ind(type~=2); %Those fibers that can be used as "ipsi" inputs (and some could also contribute to "contra" if they are type 1)
indB = ind(type==2 | type==1); %Those fibers that can only be used as "contra" inputs

for i = 1:length(indA)
    try
        clear h1;clear nullh1;clear temph1;clear tempnullh1;
        dfA = ANdfs(indA(i));
        lo = dfA*2^(-Frange/2);
        hi = dfA*2^(Frange/2);
        setB = indB(ANdfs(indB)<=hi & ANdfs(indB)>=lo);
        if numel(setB)>10
            setB = randsample(setB,10);
        end
        if numel(setB)<2
            continue;
        end
        [cispikes,dummy,wv,dt,h1ipsi{i},h1contra{i}] = Mon_Conv_get_coincident_spikes(indA(i),setB,thr);
        if ~isempty(cispikes)
            NrEst = length(cispikes);
            
            Nwin = round(RCwin/dt);
            for j = 1:NrEst
                Nsamples = length(wv{j});
                cispikes{j} = cispikes{j}(cispikes{j}<5000 & cispikes{j}>50);
                [h1{j},tvect,dummy] = local_h1_kernel(cispikes{j});
            end
            %         Bspikes = Bspikes(Bspikes<5000 & Bspikes>50);
            %         [Bh1{i},dummy,dummy] = local_h1_kernel(Bspikes);
            
            %     Bootstrap to get the noise floor for the summed spike train
            for j = 1:NrEst
                nullh1{j,1} = zeros(Nwin,Nboot);
                nullh1{j,2} = zeros(Nwin,Nboot);
                for k = 1:Nboot
                    nullspikes = scramble_spikes(cispikes{j});
                    [nullvals,dummy,dummy] = local_h1_kernel(nullspikes);
                    nullh1{j,1}(:,k) = nullvals(:,1);
                    nullh1{j,2}(:,k) = nullvals(:,2);
                end
            end
            %         nullBh1 = zeros(Nwin,Nboot);
            %         for k = 1:Nboot
            %             nullspikesB = scramble_spikes(Bspikes);
            %             [nullvals,dummy,dummy] = local_h1_kernel(nullspikesB);
            %             nullBh1(:,k) = nullvals(:,2);
            %         end
            
            temph1 = ([h1{:}]);
            clear h1;
            h1time{i}(:,1) = mean(temph1(:,1:2:end),2);
            h1time{i}(:,2) = mean(temph1(:,2:2:end),2);
            
            tempnullh1{1} = zeros(Nwin,Nboot);
            tempnullh1{2} = zeros(Nwin,Nboot);
            for j = 1:NrEst
                tempnullh1{1} = tempnullh1{1}+nullh1{j,1};
                tempnullh1{2} = tempnullh1{2}+nullh1{j,2};
            end
            clear nullh1;
            tempnullh1{1} = tempnullh1{1}/NrEst;
            tempnullh1{2} = tempnullh1{2}/NrEst;
            
            nullh1 = zeros(Nwin,2*Nboot);
            nullh1(:,1:2:end) = tempnullh1{1};
            nullh1(:,2:2:end) = tempnullh1{2};
            
            h1time{i} = local_low_pass(h1time{i});
            
            [h1mag{i},freqscale,h1phase{i},h1zscore{i},h1wfilt{i}] = Phase_Freq_analysis(h1time{i},nullh1);
            %         [h1Bmag{i},dummy,h1Bphase{i},h1Bzscore{i},h1Bwfilt{i}] = Phase_Freq_analysis(Bh1{i},nullBh1);
            %         [h1mag{i},freqscale,h1phase{i},h1zscore{i}] = Phase_Freq_analysis(h1time{i},nullh1);
            for j = 1:2
                h1mag{i}(:,j) = Trifilter(h1mag{i}(:,j)',13)';
            end
            
            [h1bw{i},qvals{i},domfreq{i},dummy] = Bandwidth_analysis(h1mag{i},h1zscore{i},spec_crit,freqscale);
            
            %Get the predicted difcor ("bincor") for the input fibers and for
            %the pseudo-binaural output
            [bincor_bin{i},BD_bin{i},ITD] = local_get_bincor(h1wfilt{i}(:,1),h1wfilt{i}(:,2));
            bincor_bin{i} = bincor_bin{i}/max(abs(bincor_bin{i}));
            for j = 1:length(h1contra{i}.h1contra)
                [bincor_mon{i}(:,j),BD_mon{i}(j),dum] = local_get_bincor(h1ipsi{i}.h1ipsi{1},h1contra{i}.h1contra{j});
            end
            
            if plotYN
                %Draw the results
                figure(1);
                subplot(2,2,1);cla; %Magnitude spectrum - ipsi
                
                yvals = h1ipsi{i}.h1ipsimag{1};
                yvals(h1ipsi{i}.h1ipsizscore{1}<spec_crit) = NaN;
                yvals = yvals-max(yvals);
                plot(freqscale/1000,yvals,'color',cols{1},'linewidth',1);hold on;
                plot(h1ipsi{i}.h1ipsidomfreq{1},0,'o','markeredgecolor',cols{1},'markerfacecolor',cols{1},'linewidth',1,'markersize',10);
                
                yvals = h1mag{i}(:,1);
                yvals(h1zscore{i}(:,1)<spec_crit) = NaN;
                yvals = yvals-max(yvals);
                plot(freqscale/1000,yvals,'k-','linewidth',3);
                
                set(gca,'xscale','log','xlim',[0.05 5],'ylim',[-20 2],'fontsize',16,'linewidth',1,'xtick',[0.1 1],'xticklabel',[0.1 1],'layer','top');
                ylabel 'GAIN (dB)';
                xlabel 'FREQUENCY (kHz)';
                box off;
                
                
                subplot(2,2,2);cla; %Magnitude spectrum - contra
                for j = 1:length(h1contra{i}.h1contra)
                    yvals = h1contra{i}.h1contramag{j};
                    yvals(h1contra{i}.h1contrazscore{j}<spec_crit) = NaN;
                    yvals = yvals-max(yvals);
                    plot(freqscale/1000,yvals,'color',cols{j},'linewidth',1);hold on;
                    plot(h1contra{i}.h1contradomfreq{j},0,'o','markeredgecolor',cols{j},'markerfacecolor',cols{j},'linewidth',1,'markersize',10);
                    
                    yvals = h1mag{i}(:,2);
                    yvals(h1zscore{i}(:,2)<spec_crit) = NaN;
                    yvals = yvals-max(yvals);
                    plot(freqscale/1000,yvals,'k-','linewidth',3);
                end
                
                %Plot the geometric mean of the dominant frequencies of the
                %contralateral input fibers
                averagecontradomfreq = geomean([h1contra{i}.h1contradomfreq{:}]);
                plot(averagecontradomfreq,1,'v','color',[0.6 0.6 0.6],'markerfacecolor',[.6 .6 .6],'markersize',12);
                
                %Plot the dominant frequency at the output of the binaural revcor
                %for this ear
                plot(domfreq{i}(2)/1000,1,'kv','markerfacecolor','k','markersize',12);
                
                set(gca,'xscale','log','xlim',[0.05 5],'ylim',[-20 2],'fontsize',16,'linewidth',1,'xtick',[0.1 1],'xticklabel',[0.1 1],'layer','top');
                ylabel 'GAIN (dB)';
                xlabel 'FREQUENCY (kHz)';
                box off;
                
                subplot(2,2,3);cla; %Time domain - ipsi
                yvals = h1ipsi{i}.h1ipsi{1};
                yvals = yvals/max(abs(yvals));
                plot(tvect,yvals,'b','linewidth',1);hold on;
                
                yvals = h1wfilt{i}(:,1);
                yvals = yvals/max(abs(yvals));
                plot(tvect,yvals,'k--','linewidth',2);
                box off
                set(gca,'xlim',[0 15],'ylim',[-1.2 1.2],'fontsize',16,'linewidth',1,'layer','top');
                xlabel 'TIME (ms)';
                ylabel 'AMPLITUDE';
                
                subplot(2,2,4);cla; %Time domain - contra
                for j = 1:length(h1contra{i}.h1contra)
                    yvals = h1contra{i}.h1contra{j};
                    yvals = yvals/max(abs(yvals));
                    plot(tvect,yvals,'color',cols{j},'linewidth',1);hold on;
                end
                yvals = h1wfilt{i}(:,2);
                yvals = yvals/max(abs(yvals));
                plot(tvect,yvals,'k--','linewidth',2);
                box off
                set(gca,'xlim',[0 15],'ylim',[-1.2 1.2],'fontsize',16,'linewidth',1,'layer','top');
                xlabel 'TIME (ms)';
                ylabel 'AMPLITUDE';
                
                
                figure(2);cla;
                plot([0 0],[-1.2 1.2],'--','color',[0.6 0.6 0.6],'linewidth',3);hold on;
                for j = 1:length(h1contra{i}.h1contra)
                    yvals = bincor_mon{i}(:,j)/max(abs(bincor_mon{i}(:,j)));
                    plot(ITD,yvals,'color',cols{j},'linewidth',1);hold on;
                    plot(BD_mon{i}(j),yvals(ITD==BD_mon{i}(j)),'o','color',cols{j},'markerfacecolor',cols{j},'linewidth',1);
                end
                
                plot(ITD,bincor_bin{i},'k-','linewidth',3);hold on;
                box off
                set(gca,'xlim',[-5 5],'ylim',[-1.2 1.2],'fontsize',16,'linewidth',1,'layer','top');
                xlabel 'ITD (ms)';
                ylabel 'AMPLITUDE';
            end
            
        end
    catch
        fprintf('failed to analyze unit %2.0f \n',i);
    end
end
return;
    function [bincorout,BDout,ITDout] = local_get_bincor(IN1,IN2)
        bincorout = xcorr(IN1,IN2);
        ITDout = -sort([(1:floor(numel(bincorout)/2))*-dt 0 (1:floor(numel(bincorout)/2))*dt]);
        [dum,ind] = max(bincorout);
        BDout = ITDout(ind);
    end
    function [h1out,Time,Nspikes] = local_h1_kernel(ST)
        Nspikes = length(ST);
        tx = (0:Nsamples)*dt; %Binning edges in ms ...
        pst = histc(ST,tx);%Make PSTH with same sampling frequency as stimulus
        for ii = 1:Nchan
            h1out(:,ii) = flipud(xcorr(wv{j}(:,ii),pst,Nwin))/Nspikes;%Get the revcor (h1 kernel) - cross correlation of the pst and the stimulus waveform
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
    function [h1magout,freqscale,h1phaseout,zcritout,wf] = Phase_Freq_analysis(y,nully)
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
            zcritout(:,ii)= ((10.^(h1magout(:,ii)/20))-mean(10.^(nullh1magout(:,ii:2:end)/20),2))./std(10.^(nullh1magout(:,ii:2:end)/20),1,2);
        end
        
        
        %Do the Wiener filtering
        for ii = 1:Nchan
            zz = [zcritout(:,ii); flipud(zcritout(:,ii))];
            h1fft(zz<spec_crit,ii) = 0; %Set noise components to zeros
            vals = ifft(h1fft(:,ii),nfft);
            vals = real(vals);
            wf(:,ii) = vals(1:length(y));
            for jj = 1:size(nully,2)
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