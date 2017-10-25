function [kernels] = Phase_Freq_Fit(kernels)

%Returns the monaural phase-frequency fits and the broad-band inter-aural
%CP and CD.
s = fft(kernels.h1c,2^13);
s = s(1:4096,:);
phase = angle(s);
mag = 20*log10(abs(s));
mag = bsxfun(@minus,mag,max(mag));
mag(mag<-20) = nan;
freq = kernels.h1ffax;
for i = 1:size(mag,2)
    [dummy,ind(i)] = max(mag(:,i));
end
Nchan = length(ind);
%Weights for the linear fit, based on the amplitude re. noise floor
weights = 10.^(mag/20);
weights(isnan(weights)) = 0;
cycles = phase./(2*pi);
for i = 1:Nchan
    %Do the regression based on the wrapped data using circular-linear
    %regression, a la victor.
    [PhaseDelay(i),GroupDelay(i),ResidError(i),N(i)] = fit_circ_lin(freq(weights(:,i)>0),cycles(weights(:,i)>0,i),weights(weights(:,i)>0,i));
    inds = find(weights(:,i)>0);
    %Calculate the unwrapped phase - but don't use this for the fitting -
    %it's just returned so that it can be used for plotting purposes
    phase_uw(:,i) = nan(size(phase(:,i)));
    phase_uw(inds,i) = unwrap(phase(inds,i))/(2*pi);
end
if Nchan==2
    IAP = diff(cycles,1,2);%inter-aural phase
    Bweights = prod(weights,2);
    Bweights = Bweights./max(Bweights);
    [dummy,indx]=max(Bweights);
    [CP,CD,BinError,BinN] = fit_circ_lin(freq(Bweights>0),IAP(Bweights>0),Bweights(Bweights>0));
    %Now that you have the real CP you can go back and correct h1phasemask
    %to give the real relationship between the two phase-frequency
    %functions
    %find the inter-aural phase relation at the binaural dominant frequency
    startphase = (freq(indx)*CD)+CP;
    %equalize the phase at that frequency
    phasenow = diff(phase_uw(indx,:));
    phase_uw(:,2) = phase_uw(:,2)-phasenow;
    %correct to the right phase
    phase_uw(:,2) = phase_uw(:,2)+startphase;
    phase_uw = phase_uw-ceil(max(phase_uw(:)));
    %add the unwrapped inter-aural phase;
    phase_uw(:,3) = phase_uw(:,2)-phase_uw(:,1);
end
%Make an output structure
Monaural = struct('Phase',PhaseDelay,'Delay',GroupDelay,'Error',ResidError,'N',N,'weights',weights);
if Nchan==2
    Binaural = struct('CP',CP,'CD',CD,'Error',BinError,'N',BinN,'weights',Bweights);
    coeffs = struct('Monaural',Monaural,'Binaural',Binaural);
else
    coeffs = struct('Monaural',Monaural);
end
kernels.h1coeffs = coeffs;
kernels.h1phase_uw = phase_uw;
kernels.h1cmag = mag;
kernels.h1cphase = phase;