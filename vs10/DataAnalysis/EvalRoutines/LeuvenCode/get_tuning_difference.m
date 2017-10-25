% function [tuning_diff] = get_tuning_difference(freqs,h1mag,NDF)

function [kernels] = get_tuning_difference(kernels,NDF)

freqs = kernels.h1ffax;
h1mag = kernels.h1cmag;
binmag = NDF.difcor.power;
binmag(binmag<-20) = NaN;
binmag = 10.^(binmag/20);
binmag(isnan(binmag))=0;

ipsimag = h1mag(:,1);
ipsimag = 10.^(ipsimag/20);
ipsimag(isnan(ipsimag))=0;

contramag = h1mag(:,2);
contramag = 10.^(contramag/20);
contramag(isnan(contramag))=0;

%cross correlate each ear with the binaural tuning
[cipsi,dummy] = xcorr(ipsimag,binmag);
[ccontra,lags] = xcorr(contramag,binmag);
autoipsi = xcorr(ipsimag,0);
autocontra = xcorr(contramag,0);
autobin = xcorr(binmag,0);
normcipsi = cipsi/sqrt(autoipsi*autobin); % normalize cross correlation
normccontra = ccontra/sqrt(autocontra*autobin); % normalize cross correlation
[dummy, indcipsi] = max(normcipsi);
[dummy, indccontra] = max(normccontra);
difff = freqs(2)-freqs(1); % Get increment of freqUp
xcor.dfreqskHz = lags*difff; % Convention: positive frequency differences mean the test ear is tuned higher than the binaural response.
difcorpeak_kHz = NDF.difcor.peakhz/1000;
N = length(lags);
dfreqs_kHz = xcor.dfreqskHz(ceil(N/2):end)+difcorpeak_kHz;
dfreqs_oct = log2(dfreqs_kHz/difcorpeak_kHz);
xcor.dfreqsoct = [-fliplr(dfreqs_oct(2:end)) dfreqs_oct];

%Tuning differences
tuning_diff.ipsivsbin_kHz = xcor.dfreqskHz(indcipsi);
tuning_diff.contravsbin_kHz = xcor.dfreqskHz(indccontra);
tuning_diff.ipsivsbin_oct = xcor.dfreqsoct(indcipsi);
tuning_diff.contravsbin_oct = xcor.dfreqsoct(indccontra);
tuning_diff.ipsivscontra_kHz = tuning_diff.ipsivsbin_kHz - tuning_diff.contravsbin_kHz;
tuning_diff.ipsivscontra_oct = tuning_diff.ipsivsbin_oct - tuning_diff.contravsbin_oct;
tuning_diff.ipsivscontra_mm = tuning_diff.ipsivscontra_oct*2.55;%Muller et al. 2010. Hear Res.
tuning_diff.normcipsi = normcipsi;
tuning_diff.normccontra = normccontra;
tuning_diff.normcipsi = normcipsi;
tuning_diff.freq_kHz = xcor.dfreqskHz;
tuning_diff.freq_oct = xcor.dfreqsoct;

%output
kernels.tuning_diff = tuning_diff;
return;