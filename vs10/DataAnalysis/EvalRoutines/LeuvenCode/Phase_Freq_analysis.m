function [h1magout,nullh1magout,freqscale,h1phaseout,zcritout,h1phasemask,h1magmask,wf,nullwf] = Phase_Freq_analysis(h1,nullh1,params)

dt = params.dt;
spec_crit = params.spec_crit;
Nchan = size(h1,2);
wf = zeros(size(h1));
nullwf = zeros(size(nullh1));
nfft = 2^13;
h1fft = fft(h1,nfft);
h1magout = 20*log10(abs(h1fft));
h1magout = h1magout(1:floor(nfft/2),:);
ScaleFactor = max(h1magout);
h1magout = bsxfun(@minus,h1magout,ScaleFactor);%dB gain re. peak
nullh1fft = fft(nullh1,nfft);
nullh1magout = 20*log10(abs(nullh1fft));
nullh1magout = nullh1magout(1:floor(nfft/2),:,:);
nullh1magout = bsxfun(@minus,nullh1magout,ScaleFactor);%dB gain re. peak
nf = 0.5*(1000/dt);
freqscale = nf*linspace(0,1,floor(nfft/2))/1000;
h1phaseout = angle(h1fft(1:floor(nfft/2),:));
%get z-scores for magnitude spectrum
zcritout = zeros(size(h1magout));
zcritout = ((10.^(h1magout/20))-mean(10.^(nullh1magout/20),3))./std(10.^(nullh1magout/20),1,3);

%Do the Wiener filtering
for ii = 1:Nchan
    zz = [zcritout(:,ii); flipud(zcritout(:,ii))];
    h1fft(zz<spec_crit,ii) = 0; %Set noise components to zeros
    vals = ifft(h1fft(:,ii),nfft);
    vals = real(vals);
    wf(:,ii) = vals(1:length(h1));
    for jj = 1:size(nullh1,3)%NBoot
        nullh1fft(zz<spec_crit,ii,jj) = 0; %Set noise components to zeros
        vals = ifft(nullh1fft(:,ii,jj),nfft);
        vals = real(vals);
        nullwf(:,ii,jj) = vals(1:length(h1));
    end
end

%Apply smoothing to the output z-score estimate (Triangular
% smoothing function)
[zcritout] = Trifilter(zcritout',7)';
h1magout = Trifilter(h1magout',13)';%Some smoothing

%Make sure the gain is zero at the peak
ScaleFactor = max(h1magout);
h1magout = bsxfun(@minus,h1magout,ScaleFactor);
nullh1magout = bsxfun(@minus,nullh1magout,ScaleFactor);

h1phasemask = h1phaseout;
h1magmask = h1magout;
for i=1:Nchan
    h1phasemask(zcritout(:,i)<spec_crit,i) = NaN;
    h1magmask(zcritout(:,i)<spec_crit,i) = NaN;
    h1magmask(freqscale>6,i) = NaN;
end
end

