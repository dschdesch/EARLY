function [NDFpred] = difcor_prediction_h1(h1,dt,signalPOWER,NDF)

%predict the difcor based on the contra and ipsi h1 kernels scaled as
%spikes per second per Pascal.
[c, lags] = xcorr(h1(:,2),h1(:,1),'none');
NDFpred.difcor.rate = c*(2*signalPOWER^2);
NDFpred.xvals.ITDx = lags*dt;

%ITD vals to get
ITDs = NDF.xvals.ITDx;
NDFpred.xvals.ITDxsampled = ITDs;
%predicted rate at those ITDs
NDFpred.difcor.ratesampled = interp1(NDFpred.xvals.ITDx,NDFpred.difcor.rate,ITDs);
%coefficient of determination for this linear fit
NDFpred.difcor.COD = Coeff_Determination(NDF.difcor.rate',NDFpred.difcor.ratesampled);

%get peaks
vals = spline(ITDs,NDFpred.difcor.ratesampled,NDF.xvals.ITDxspline);
[~,ind] = max(vals);%BD
NDFpred.bd.bd = NDF.xvals.ITDxspline(ind);
%Side peaks
vals = vals-min(vals);
[peakval,peakind] = findpeaks(vals);
Dpeaks = NDF.xvals.ITDxspline(peakind);
[troughval,troughind] = findpeaks(-vals);
troughval = -troughval;
Dtroughs = NDF.xvals.ITDxspline(troughind);
if min(peakind)<min(troughind)
    troughind = [1 troughind];
    troughval = [vals(1) troughval];
end
if max(peakind)>max(troughind)
    troughind = [troughind length(vals)];
    troughval = [troughval vals(length(vals))];
end
keeppeak = zeros(size(peakval));
keeptrough = zeros(size(troughval));
for i = 1:length(peakval)
    clo = (peakval(i)-troughval(find(troughind<peakind(i),1,'last')))/...
        (peakval(i)+troughval(find(troughind<peakind(i),1,'last')));
    chi = (peakval(i)-troughval(find(troughind>peakind(i),1,'first')))/...
        (peakval(i)+troughval(find(troughind>peakind(i),1,'first')));
    if clo>0.5
        keeppeak(i) = 1;
        keeptrough(find(troughind<peakind(i),1,'last'))=1;
    elseif chi>0.5
        keeppeak(i) = 1;
        keeptrough(find(troughind>peakind(i),1,'first'))=1;
    end
end
Dpeaks = Dpeaks(logical(keeppeak));
Dtroughs = Dtroughs(logical(keeptrough));
Dpeaks = Dpeaks(Dpeaks~=NDFpred.bd.bd);
NDFpred.bd.dpeaks = Dpeaks;
NDFpred.bd.dtroughs = Dtroughs;

%Frequency domain signal
NFFT = 2^13;
nf = 0.5/(dt/1000);
freq = linspace(0,1,NFFT/2)*nf;
vals = spline(ITDs,NDFpred.difcor.ratesampled,NDF.xvals.ITDxspline);
hanwin = hann(length(vals));
spec = fft((vals'-mean(vals)).*hanwin,NFFT);
Dphase = angle(spec(1:NFFT/2));
Dspec = 20*log10(abs(spec(1:NFFT/2)));
Dspec = Dspec-max(Dspec);
Dphase(Dspec<-20) = nan;
[dummy,ind] = max(Dspec);
DspecPeakHz = freq(ind);

NDFpred.difcor.mag = Dspec;
NDFpred.difcor.phase = Dphase;
NDFpred.difcor.peakhz = DspecPeakHz;
NDFpred.difcor.freq = freq;