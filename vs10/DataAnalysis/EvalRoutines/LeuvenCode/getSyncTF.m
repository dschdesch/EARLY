function [stf] = getSyncTF(DataIN)
Nfft = 2^13;
nf = 0.5*(1000/DataIN.dt);
freq = nf*linspace(0,1,floor(Nfft/2));
for i = 1:length(DataIN.dataID)
    spikes = DataIN.dsz{DataIN.dataID(i)}.SPT;
    spikes = cat(2,spikes{DataIN.Subseq(i),:});
    spikes = spikes(spikes>=DataIN.minST & spikes<=DataIN.maxST);
    [wv,dt] = StimSam(DataIN.dsz{DataIN.dataID(i)},1);
    dt = dt/1000; % us --> ms
    Nwin = round((DataIN.h1.h1endlatency-DataIN.h1.h1frontlatency)/dt);
    dZERO=round((spikes-DataIN.h1.h1frontlatency)/dt);%first sample to take for each spike
    dMAX=dZERO-Nwin;%corresponding max time lag samples
    yy=[dMAX' dZERO']; %Compute a matrix of sample numbers to get
    yi=flipud(round(interp1([1;Nwin],yy',(1:Nwin)')));
    %spike triggered fft
    stfft = fft(wv(yi),Nfft);
    stfftmag = abs(stfft(1:Nfft/2,:));
    stfftmag = stfftmag(freq<=4000,:);
    stfftphase = angle(stfft(1:Nfft/2,:));
    stfftphase = stfftphase(freq<=4000,:);
    clear stfft;
    %random fft
    dZERO = Nwin+1:Nwin:round(DataIN.maxST/dt);
    dMAX=dZERO-Nwin;
    yy=[dMAX' dZERO']; %Compute a matrix of sample numbers to get
    yi=flipud(round(interp1([1;Nwin],yy',(1:Nwin)')));
    rnfft = fft(wv(yi),Nfft);
    rnfftmag = abs(rnfft(1:Nfft/2,:));
    rnfftmag = rnfftmag(freq<=4000,:);
    crit_level = mean(rnfftmag,2);%+3*std(rnfft,0,2);%look for everything above the mean of the noise
    clear rnfft;
    [indi, dummy] = ind2sub(size(stfftmag),find(bsxfun(@gt,stfftmag,crit_level)));clear dummy;
    frvals = freq(indi);%Frequency vector
    phvals = stfftphase(bsxfun(@gt,stfftmag,crit_level));%Phase vector (radians)
    mvals = stfftmag(bsxfun(@gt,stfftmag,crit_level));%Magnitude vector
    ufreqs = unique(frvals);%Unique frequencies
    for j = 1:length(ufreqs)
        pp = phvals(frvals==ufreqs(j));
        sum_sines = sum(sin(pp));%Sum of sines and cosines
        sum_cosines = sum(cos(pp));
        nspikes(i,j) = length(pp);
        VS(i,j) = sqrt(sum_sines^2 + sum_cosines^2)/nspikes(i,j);%Vector strength       
        mm = mvals(frvals==ufreqs(j));
        mag(i,j) = sum(mm);%/nspikes(i,j);
    end
    mag = mag/sum(nspikes);
    NrST = sum(nspikes,1);
    stf.vs = mean(VS,1);
    stf.rs = 2*(stf.vs.^2).*NrST;
    stf.pval = exp(-NrST.*(stf.vs.^2));
    stf.freq = ufreqs;
end
            