function [] = efferent_Wiener(varargin)

%Analyses for conditions where you have monaural and binaural noise data
%from ANFs - looking for a contralateral efferent effect on the monaural
%filtering.

%Parse the inputs
AnNum = varargin{1};
for i = 2:nargin
    ndzdata{i-1} = varargin{i};
end

%set some useful parameters
minST = 200;%Minimum spike time in ms
maxST = 5000; %Maximum spike time in ms
RCwin = 30;%Length of reverse-correlation window in ms
Nboot = 20;%Number of boostrap samples

%Get the data and figure out which are monaural and which are binaural
if length(ndzdata)<2
    error('Not enough input data for these analyses');
end
for i = 1:length(ndzdata)
    ds(i) = dataset(AnNum,ndzdata{i});
    S(i) = struct(ds(i));
    stimtype{i} = S(i).ID.StimType;
    switch stimtype{i}
        case 'NTD'
            datatype(i) = 2;
        otherwise
            if (S(i).Stimulus.StimParam.startSPL(1)==S(i).Stimulus.StimParam.startSPL(2))...
                    && (S(i).Stimulus.StimParam.endSPL(1)==S(i).Stimulus.StimParam.endSPL(2))
                %This is a monaural stimulus - the "reference"
                disp 'Found reference data :)';
                datatype(i) = 1;
            else
                %This is a test "binaural" stimulus
                %Which ear is varied?
                
                
                if (S(i).Stimulus.StimParam.startSPL(1)~=S(i).Stimulus.StimParam.endSPL(1))...
                        && (S(i).Stimulus.StimParam.startSPL(2)==S(i).Stimulus.StimParam.endSPL(2))
                    %The ipsi ear is varied and the contra ear is fixed
                    disp 'Found contra efferent data :)';
                    datatype(i) = 2;
                elseif (S(i).Stimulus.StimParam.startSPL(1)==S(i).Stimulus.StimParam.endSPL(1))...
                        && (S(i).Stimulus.StimParam.startSPL(2)~=S(i).Stimulus.StimParam.endSPL(2))
                    %The contra ear is varied and the ipsi ear is fixed
                    disp 'Found contra efferent growth function data :)';
                    datatype(i) = 3;
                else
                    error('Unknown data type!');
                end
            end
    end
    [wv,dt] = StimSam(ds(i),1);
    wv = wv(:,1);
    Nsamples = length(wv);
    dt = 1e-3*dt; %us -> ms
    Nwin = round(RCwin/dt);
    SPLs{i} = S(i).Stimulus.IndepVar.Values;
    SPLs{i} = SPLs{i}(~isnan(SPLs{i}));
    for j = 1:length(SPLs{i})
        nullh1{j,i} = zeros(Nwin,Nboot);
        spikes = cat(2,S(i).Data.SpikeTimes{j,:});
        spikes = spikes(spikes>=minST & spikes<=maxST);%limit spike times to only include driven spikes
        %Get the revcor
        [h1{j,i},Time,SpikeCount(j,i)] = local_h1_kernel(spikes);
        h1{j,i} = local_low_pass(h1{j,i});
        
        %Bootstrap to define the noise floor
        h = waitbar(0,'Bootstrapping revcor: Please wait...');
        for k = 1:Nboot
            nullspikes = scramble_spikes(spikes)';
            nullspikes = nullspikes(nullspikes>=minST);
            [nullvals,dummy] = local_h1_kernel(nullspikes);
            nullh1{j,i}(:,k) = nullh1{j,i}(:,k)+nullvals;
            waitbar(k/Nboot);
        end
        close (h);
        
        %Get the magnitude and phase response
        [h1mag{j,i},nullh1mag{j,i},ffax,h1phase{j,i},h1zscore{j,i}] = Phase_Freq_analysis(h1{j,i},nullh1{j,i});
        h1mag{j,i} = Trifilter(h1mag{j,i}',41)';
        
    end
end


%% Plotting

if ismember(1,datatype) && ismember(2,datatype)
    %You have both the control revcors and the same revcors with the contra
    %ear at a constant SPL
    
    %Which SPLs do you have in common between the two datasets
    [commonSPLs,IA,IB] = intersect(SPLs{1},SPLs{2});
    nrSPLs = length(commonSPLs);
    
    figure;
    for i  = 1:nrSPLs
        plot(Time,h1{IA(i),1}./max(abs(h1{IA(i),1}))+(i-1),'k',Time,h1{IB(i),2}./max(abs(h1{IB(i),2}))+(i-1),'r','linewidth',2);
    end
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
        h1out = flipud(xcorr(wv,pst,Nwin))/Nspikes;%Get the revcor (h1 kernel) - cross correlation of the pst and the stimulus waveform
        h1out = h1out((1:Nwin)+(Nwin+1),:);
        Time = tx(1:Nwin);
    end
    function [scramspikes] = scramble_spikes(spikes)
        % shuffle spike train by randomizing the ISI's
        D = diff([0; sort(spikes')]); % first sort to make sure diff gives the ISI's
        NInt = length(D);
        scramspikes = cumsum(D(randperm(NInt)));
    end
    function [h1magout,nullh1magout,freqscale,h1phaseout,zcritout] = Phase_Freq_analysis(y,nully)
        nfft = 2^15;
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
        zcritout = ((10.^(h1magout/20))-mean(10.^(nullh1magout/20),2))./std(10.^(nullh1magout/20),1,2);
        %Apply smoothing (Triangular smoothing function)
        zcritout = Trifilter(zcritout',7)';
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