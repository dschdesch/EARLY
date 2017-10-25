function [] = AN_X_h2(varargin)

AnNum = varargin{1};
dsid = varargin{2};
dsz = dataset(AnNum,dsid);
ds = struct(dsz);
maxST = 5000;
minST = 30;
SpikeTimes = cat(2,ds.Data.SpikeTimes{1,:});
SpikeTimes = SpikeTimes(SpikeTimes>minST & SpikeTimes<=maxST);
ST{1} = SpikeTimes(SpikeTimes<maxST/2);
ST{2} = SpikeTimes(SpikeTimes>=maxST/2);
ST{2} = ST{2}-maxST/2;
ST{2} = ST{2}(ST{2}>minST);
[wv,dt] = StimSam(dsz,1);
dt = dt/1000;
Nsamples = floor(maxST/dt);
Nwin = round(30/dt);
wv = wv(1:Nsamples);
temp = wv;clear wv;
wv(:,1) = temp(1:Nsamples/2);
wv(:,2) = temp(Nsamples/2+1:end);clear temp;
delaysamps = round(2/dt);
wv(:,2) = [wv(delaysamps+1:end,2); zeros(delaysamps,1)];
h2wv = local_low_pass(wv,3000,3500);
[h1,tvect,spikecount] = local_h1_kernel(ST);
[h2,h2X] = local_h2_kernel(ST);


    function [h1out,Time,Nspikes] = local_h1_kernel(ST)
        for ii = 1:2
            Nspikes(ii) = length(ST{ii});
            tx = (0:Nsamples)*dt; %Binning edges in ms ...
            pst = histc(ST{ii},tx);%Make PSTH with same sampling frequency as stimulus
            temph1 = flipud(xcorr(wv(:,ii),pst,Nwin))/Nspikes(ii);%Get the revcor (h1 kernel) - cross correlation of the pst and the stimulus waveform
            h1out(:,ii) = temph1((1:Nwin)+(Nwin+1));
        end
        Time = tx(1:Nwin);
    end
    function [h2out,h2Xout] = local_h2_kernel(Spikes)
        for ii = 1:2
            Spikes = ST{ii}';
            Nspikes(ii) = length(Spikes);
            dZERO=round(Spikes/dt);%sample for each spike in response to this stimulus
            dMAX=dZERO-Nwin+1;%corresponding max time lag samples
            yi = fliplr(cell2mat(arrayfun(@colon,dMAX,dZERO,'uniform',false)));
            hh2wv = h2wv(:,ii);
            wvh2{ii} = hh2wv(yi);
            autocorr = xcorr(hh2wv,Nwin-1,'unbiased');
            autocorrmatrix = toeplitz(autocorr(Nwin:end));%Toeplitz matrix
            h2out(:,:,ii)=(wvh2{ii}'*wvh2{ii})-(Nspikes(ii)*autocorrmatrix);%R2 is the outer product of the two matrices corrected for the signal correlation
        end
        Nlines = min(Nspikes);
        h2Xout = wvh2{1}(1:Nlines,:)'*wvh2{2}(1:Nlines,:);
    end
    function [filtout] = local_low_pass(IN,varargin)
        if isempty(varargin)
            fcuts = [4500 5500]; % Transition bands
        else
            fcuts = [varargin{1} varargin{2}];
        end
        fsamp = 1000/dt; %sample frequency
        
        mags = [1 0];
        devs = [0.01 10^(-40/20)];% 1% passband ripple and 40-dB attenuation of stopband
        
        [n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);
        b = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale'); %Filter coefficients to match the required parameters
        filtout = filtfilt(b,1,IN); %Zero-phase filtering
    end
end