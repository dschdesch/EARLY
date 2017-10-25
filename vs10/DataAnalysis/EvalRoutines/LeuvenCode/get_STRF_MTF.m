function [MTFmag,tMTFdB,MTFyScale] = get_STRF_MTF(LSTRF,df,NFFT)


MTF=fft2(LSTRF)/size(LSTRF,1)/size(LSTRF,2);

% sumLSTRF = sum(LSTRF);
% figure;plot(abs(fft(sumLSTRF)));

MTFamp=(MTF(2:end/2+1,2:end/2+1)) + ...
    (MTF(end:-1:end/2+1,2:end/2+1)) + ...
    (MTF(2:end/2+1,end:-1:end/2+1)) + ...
    (MTF(end:-1:end/2+1,end:-1:end/2+1)); %This collapses the 4 quadrants, max tM frequency is 25 kHz. Throw out 0-Hz component.
MTFamp=MTFamp(:,1:end/2); %max tM frequency is 12.5 kHz
MTFmag=20*log10(abs(MTFamp));
MTFyScale=linspace(0,1/df/2,NFFT/8+1)*1e3;
MTFyScale = MTFyScale(2:end);

tMTF=abs(sum(MTFamp));
tMTFdB=20*log10(tMTF);

%Note to self:
% This calculates the tMTF by collapsing across all spectral modulation
% frequencies - that is probably not ideal since you end up averaging in a lot
% of noise - therefore might be better to calculate the tMTF based on a
% band around CF.