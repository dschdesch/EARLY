function [] = Monaural_Coincidence_Pairs(Fstart,Ctype)

%Inputs:
%Fstart: approximate CF of reference fiber
%Ctype: type of comparison; within animal (1), or across animal (2)

%Monaural coincidence analysis for ANFs in response to noise. Finds the
%revcor from just coincident spikes from fibers with neighbouring CFs in
%the region of CFin (+/- Frange);

RevCorDir = ('C:\work\BinRev\MonRev\');

ANdata = load('C:\LeuvenDataAnalysis\MonRevPopData.mat');
nrdata = length(ANdata.AN);
for i = 1:nrdata
    if ismember(111,ANdata.AN(i).noiseseed)
        if ~isempty(ANdata.AN(i).noiseSPLs)
            ind = find(unique(ANdata.AN(i).noiseSPLs)==70);
            ANdfs(i) = ANdata.AN(i).df(ind);
            if isempty(ind)
                ANdfs(i) = NaN;
            end
        else
            if ~isempty(ANdata.AN(i).df)
                ANdfs(i) = ANdata.AN(i).df;
            else
                ANdfs(i) = NaN;
            end
        end
        ANcfs(i) = ANdata.AN(i).cf;
    else
        ANcfs(i) = ANdata.AN(i).cf;
        ANdfs(i) = NaN;
    end
end
ANfreqs = ANdfs;
ANfreqs(ANcfs>2000) = ANcfs(ANcfs>2000)/1000;

allind = find(ANfreqs<=1.5);
[dummy,CFind] = min(abs(ANfreqs-Fstart));
if ~ismember(CFind,allind)
    error ('This CF is not available');
end
allpairs = [allind' repmat(CFind,length(allind),1)];%all possible fiber pairs
allpairsdfs = ANfreqs(allpairs);
octdiff = log2(allpairsdfs(:,1)./allpairsdfs(:,2));

ind = [];
indpairs = allpairs(abs(octdiff)<=0.5,:);
for i = 1:size(indpairs,1)
    isok = 0;
    switch Ctype
        case 1
            if  strcmp(ANdata.AN(indpairs(i,1)).animal,ANdata.AN(indpairs(i,2)).animal)%Same animal or different animal
                isok = 1;
            end
        case 2
            if  ~strcmp(ANdata.AN(indpairs(i,1)).animal,ANdata.AN(indpairs(i,2)).animal)%Same animal or different animal
                isok = 1;
            end
    end
    if isok
        if ismember(111,ANdata.AN(indpairs(i,1)).noiseseed) && ismember(111,ANdata.AN(indpairs(i,2)).noiseseed)%Same noise seed
%             foo;
            if length(ANdata.AN(indpairs(i,1)).noiseseed)==length(ANdata.AN(indpairs(i,1)).noiseSPLs)
                if ismember(70,ANdata.AN(indpairs(i,1)).noiseSPLs(ANdata.AN(indpairs(i,1)).noiseseed==111))
                    pass1 = 1;
                else
                    pass1 = 0;
                end
            else
                if ANdata.AN(indpairs(i,1)).noiseSPLs == 70
                    pass1 = 1;
                else
                    pass1 = 0;
                end
            end
            if length(ANdata.AN(indpairs(i,2)).noiseseed)==length(ANdata.AN(indpairs(i,2)).noiseSPLs)
                if ismember(70,ANdata.AN(indpairs(i,2)).noiseSPLs(ANdata.AN(indpairs(i,2)).noiseseed==111))
                    pass2 = 1;
                else
                    pass2 = 0;
                end
            else
                if ANdata.AN(indpairs(i,2)).noiseSPLs == 70
                    pass2 = 1;
                else
                    pass2 = 0;
                end
            end
            if pass1 && pass2
                ind = [ind i];
            end
        end
    end
end
indpairs = indpairs(ind,:);
clear ind;

% Initialize some variables and allocate memory
Nboot = 20;
nrpairs = size(indpairs,1);
% hpairs = waitbar(0,'Pairwise ANF coincidence analysis: Please wait...');
d = 0;
for i = 1:nrpairs;
    clear dataout;
%     dataout.h1ci = zeros(1803,7);dataout.h1cifilt = zeros(1803,7);
%     dataout.h1cifiltmag = zeros(4096,7);
%     dataout.h1cifiltphase = zeros(4096,7);
%     dataout.zcritcifilt = zeros(4096,7);
%     dataout.SpikeCountci = zeros(1,7);
%     dataout.domfreqA = 0; dataout.domfreqB = 0; dataout.cidomfreq = zeros(1,7);
%     dataout.meandomfreq = 0;
%     dataout.domfreqAdiff_oct = 0; dataout.domfreqBdiff_oct = 0;
%     dataout.cidomfreqdiff_oct = zeros(1,7);
%     dataout.h1cibw = zeros(2,7);
%     dataout.qvalsci = zeros(2,7);
%     dataout.domindci = zeros(2,7);
%     splind = zeros(1,2);
%     spiketimes = cell(1,2);
%     data = cell(1,2);
    ind = cell(2,1);
    for j = 1:2
        data{j} = load(fullfile(RevCorDir,ANdata.AN(indpairs(i,j)).filename));
        spiketimes{j} = [];
        for k = 1:length(data{j}.ParamsOut.StimulusInfo)
            if data{j}.ParamsOut.StimulusInfo{k}.StimParam.RandomSeed==111
                if data{j}.ParamsOut.StimulusInfo{k}.StimParam.SPL(1)==70
                    ind{j} = [ind{j} k];
                end
            end
        end
    end
    if ~isempty(ind{1}) && ~isempty(ind{2})
        for j = 1:2
            %Get the spike times for these datasets
            for k = 1:length(ind{j})
                if data{j}.ParamsOut.StimulusInfo{ind{j}(k)}.StimParam.SPL(1)==70
                    spikedata = dataset(data{j}.ParamsOut.DatasetIDs.Animal,[char(data{j}.ParamsOut.DatasetIDs.noise(ind{j}(k))) '-NTD']);
                    s = spikedata.SPT;
                    spikes = [s{:}];
                    spiketimes{j} = [spiketimes{j} spikes];
                end
            end
            SPLs = zeros(1,length(data{j}.ParamsOut.StimulusInfo));
            for k = 1:length(data{j}.ParamsOut.StimulusInfo)
                SPLs(k) = data{j}.ParamsOut.StimulusInfo{k}.StimParam.SPL(1);
            end
            spl_list = unique(SPLs);
            splind(j) = find(spl_list==70);
            clear SPLs;
        end
        
        if ~exist('wv','var')
            [wv, dt] = StimSam(spikedata,1);
            Nsamples = length(wv);
            RCwin = 30;%Length of reverse-correlation window in ms
            dt = 1e-3*dt; %us -> ms
            Nwin = round(RCwin/dt);
        end
        
        dataout.domfreqA = data{1}.ParamsOut.domfreq(splind(1));
        dataout.domfreqB = data{2}.ParamsOut.domfreq(splind(2));
        dataout.meandomfreq = mean([dataout.domfreqA dataout.domfreqB]);
        dataout.dataset{1} = ANdata.AN(indpairs(i,1)).filename;
        dataout.dataset{2} = ANdata.AN(indpairs(i,2)).filename;
        
        for k = 1:length(d)
            cispikes = findcoincidentspikes(spiketimes{1},spiketimes{2},0.1,d(k));
            [dataout.h1ci(:,k),Time,dataout.SpikeCountci(k)] = local_h1_kernel(cispikes);
            
            h = waitbar(0,'Bootstrapping coincidence revcor: Please wait...');
            for j = 1:Nboot
                nullspikes = scramble_spikes(cispikes)';
                [nullvals,dummy] = local_h1_kernel(nullspikes);
                nullh1ci(:,j) = nullvals;
                nullh1cifilt(:,j) = local_low_pass(nullh1ci(:,j));
                waitbar(j/Nboot,h);
            end
            close (h);
            dataout.h1cifilt(:,k) = local_low_pass(dataout.h1ci(:,k));
            [dataout.h1cifiltmag(:,k),freqscale,dataout.h1cifiltphase(:,k),dataout.zcritcifilt(:,k)] = Phase_Freq_analysis(dataout.h1cifilt(:,k),nullh1cifilt);
            
            dataout.h1cifiltmag(:,k) = Trifilter(dataout.h1cifiltmag(:,k)',13)';
            [dataout.h1cibw(:,k),dataout.qvalsci(:,k),dataout.cidomfreq(k),dataout.domindci(k)] = Bandwidth_analysis(dataout.h1cifiltmag(:,k),dataout.zcritcifilt(:,k),5,freqscale/1000);
            
            dataout.cidomfreqdiff_oct(k) = log2(dataout.cidomfreq(k)/dataout.meandomfreq);
        end
        
        dataout.domfreqAdiff_oct = log2(dataout.domfreqA/dataout.meandomfreq);
        dataout.domfreqBdiff_oct = log2(dataout.domfreqB/dataout.meandomfreq);
        
        [ITDx,BDci,WDci,difcor_ci] = local_predict_difcor(data{2}.ParamsOut.h1filt,dataout.h1cifilt);
                [ITDx,BD,WD,difcor] = local_predict_difcor(data{2}.ParamsOut.h1filt,data{1}.ParamsOut.h1filt);
%         switch Ctype
%             case 1
%                 save(fullfile(CIpath,sprintf('cirevcordata_within_%04.0f.mat',i)),'dataout');
%             case 2
%                 save(fullfile(CIpath,sprintf('cirevcordata_across_%04.0f.mat',i)),'dataout');
%         end
    end
    clear ind;
%     waitbar(i/nrpairs,hpairs);
end
% close(hpairs);


    function [outspikes] = findcoincidentspikes(Aspikes,Bspikes,w,dd)
        if dataout.domfreqA<dataout.domfreqB
            if length(Bspikes)<length(Aspikes)
                outspikes = Bspikes(ismemberf(Bspikes,Aspikes+dd,'tol',w));
            else
                outspikes = Aspikes(ismemberf(Aspikes+dd,Bspikes,'tol',w));
            end
        else
            if length(Bspikes)<length(Aspikes)
                outspikes = Bspikes(ismemberf(Bspikes+dd,Aspikes,'tol',w));
            else
                outspikes = Aspikes(ismemberf(Aspikes,Bspikes+dd,'tol',w));
            end
        end
    end
    function [ITDout,BDout,WDout,preddifcor_scaled_out] = local_predict_difcor(IN1,IN2)
        preddifcor = xcorr(IN1,IN2);
        ITDout = -sort([(1:floor(numel(preddifcor)/2))*-dt 0 (1:floor(numel(preddifcor)/2))*dt]);
        [dum,ind] = max(preddifcor);
        BDout = ITDout(ind);
        [dum,ind] = min(preddifcor);
        WDout = ITDout(ind);
        preddifcor_scaled_out = preddifcor./max(abs(preddifcor));
    end
    function [h1out,Time,Nspikes] = local_h1_kernel(ST)
        Nspikes = length(ST);
        tx = (0:Nsamples)*dt; %Binning edges in ms ...
        pst = histc(ST,tx);%Make PSTH with same sampling frequency as stimulus
        h1out = flipud(xcorr(wv,pst,Nwin))/Nspikes;%Get the revcor (h1 kernel) - cross correlation of the pst and the stimulus waveform
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
    function [h1magout,freqscale,h1phaseout,zcritout] = Phase_Freq_analysis(y,nully)
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
        zcritout= ((10.^(h1magout/20))-mean(10.^(nullh1magout/20),2))./std(10.^(nullh1magout/20),1,2);
        
        %Apply smoothing (Triangular smoothing function)
        [zcritout] = Trifilter(zcritout',7)';
        
    end
    function [h1bw,qvals,domfreq,domind] = Bandwidth_analysis(magIN,zIN,critIN,fIN)
        %Returns the 3 and 6-dB bandwidths of the h1 kernels (in kHz), the dominant
        %frequency (in kHz) and the index of that value
        Nchan = 1;
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













