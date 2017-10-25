function [] = updateMonRev
%Applies an update to the data structure, with h1filt being now a
%weiner filtered version of the raw h1.
cd 'C:\work\BinRev\MonRev';
datafiles = dir('*.mat');
nrfiles = length(datafiles);
spec_crit = 3;
for i = 1:nrfiles
    clear ParamsOut;
    saveYN = 0;
    load(datafiles(i).name);
    if isfield(ParamsOut,'h1');
        Time = ParamsOut.Time;
        dt = Time(2)-Time(1);
        if ~isfield(ParamsOut,'noiseSPLs')
            if length(ParamsOut.StimulusInfo)==1
                ParamsOut.noiseSPLs = ParamsOut.StimulusInfo{1}.StimParam.SPL(1);
            else
                for j = 1:length(ParamsOut.StimulusInfo)
                    spls(j) = ParamsOut.StimulusInfo{j}.StimParam.SPL(1);
                end
                ParamsOut.noiseSPLs = unique(spls);
                clear spls;
            end
            saveYN =1;
        end
        if isfield(ParamsOut,'noiseSPLs')
            ParamsOut.noiseSPLs = unique(ParamsOut.noiseSPLs(~isnan(ParamsOut.noiseSPLs)));
            saveYN = 1;
        end
        nrSPLs = length(ParamsOut.noiseSPLs);
        clear 'h1mag' 'h1wf' 'nullh1wf' 'h1env' 'nullh1env' 'h1if' 'zcrit_h1' 'ampz_h1' 'maxenv' 'maxenvind' 'frontlatencyind' 'frontlatency' 'Latency';
        for j = 1:nrSPLs
            if ~iscell(ParamsOut.nullh1)
                temp = ParamsOut.nullh1;
                ParamsOut = rmfield(ParamsOut,'nullh1');
                ParamsOut.nullh1{1} = temp; clear temp;
            end
            ParamsOut.h1 = local_low_pass(ParamsOut.h1);
            [h1mag(:,j),dummy,dummy,dummy,dummy,h1wf(:,j),nullh1wf{j}] = Phase_Freq_analysis(ParamsOut.h1(:,j),ParamsOut.nullh1{j});
            h1mag(:,j) = Trifilter(h1mag(:,j)',13)';
            [h1env(:,j),nullh1env{j},h1if(:,j),zcrit_h1(:,j),ampz_h1(:,j)] = ENV_IF_analysis(h1wf(:,j),nullh1wf{j});
            ParamsOut.nullh1{j} = local_low_pass(ParamsOut.nullh1{j});
            
            %set all z-scores <spec_crit to 0.
            zcrit_h1(zcrit_h1(:,j)<spec_crit,j)=0;
            ampz_h1(ampz_h1(:,j)<spec_crit,j)=0;
            
            h1if(zcrit_h1(2:end,j)==0,j)=NaN;
            
            [maxenv(j),maxenvind(j)] = max(h1env(:,j));
            frontlatencyind(j) = min([find(ampz_h1(1:maxenvind(j),j)==0,1,'last')+1, numel(Time)]);
            if ~isempty(frontlatencyind)
                frontlatency(j) = Time(frontlatencyind(j));
            else
                frontlatency(j) = NaN;
            end
            
            Latency(j) = Time(maxenvind(j));
        end
        ParamsOut.h1env = h1env(:,1:nrSPLs);
        ParamsOut.h1mag = h1mag(:,1:nrSPLs);
        ParamsOut.nullh1env = {nullh1env{1:nrSPLs}};
        ParamsOut.h1if = h1if(:,1:nrSPLs);
        ParamsOut.h1envzscore = zcrit_h1(:,1:nrSPLs);
        ParamsOut.h1latency = Latency(1:nrSPLs);
        ParamsOut.h1frontlatency = frontlatency(1:nrSPLs);
        ParamsOut.h1filt = h1wf(:,1:nrSPLs);
        ParamsOut.nullh1filt = {nullh1wf{1:nrSPLs}};
        try
            ParamsOut = rmfield(ParamsOut,{'h1filtmag','nullh1filtmag','h1filtphase','h1filtzscore','h1filtenv','nullh1filtenv','h1filtif','h1filtenvzscore'});
        catch
            %Do nothing;
        end
        saveYN = 1;
    end
    if saveYN
        save(datafiles(i).name,'ParamsOut');
    end
end

function [envout,nullenvout,ifout,zcritout,ampz] = ENV_IF_analysis(y,nully)
        z = hilbert(y);
        ifout = ((1000/dt)/(2*pi)*diff(unwrap(angle(z))));%instantaneous frequency (Hz)
        envout = abs(z);%Hilbert envelope of the revcors
        znull = hilbert(nully);
        nullenvout = abs(znull);
        %convert to z scores
        zcritout = (envout-mean(nullenvout,2))./std(nullenvout,1,2);
        %Apply smoothing (Triangular smoothing function)
        [zcritout] = Trifilter(zcritout',15)';
        
        ampz = (abs(y)-mean(abs(nully),2))./std(abs(nully),1,2);
        [ampz] = Trifilter(ampz',15)';
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
    function [h1magout,nullh1magout,freqscale,h1phaseout,zcritout,wf,nullwf] = Phase_Freq_analysis(y,nully)
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
        zcritout = ((10.^(h1magout/20))-mean(10.^(nullh1magout/20),2))./std(10.^(nullh1magout/20),1,2);
        
        %Do the Wiener filtering
        
        zz = [zcritout; flipud(zcritout)];
        h1fft(zz<spec_crit) = 0; %Set noise components to zeros
        vals = ifft(h1fft,nfft);
        vals = real(vals);
        wf = vals(1:length(y));
        for jj = 1:size(nully,2)
            nullh1fft(zz<spec_crit,jj) = 0; %Set noise components to zeros
            vals = ifft(nullh1fft(:,jj),nfft);
            vals = real(vals);
            nullwf(:,jj) = vals(1:length(y));
        end
        
        
        %Apply smoothing to the output z-score estimate (Triangular
        %smoothing function)
        [zcritout] = Trifilter(zcritout',7)';
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