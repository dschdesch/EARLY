function [filtout] = zero_phase_filt(IN,fsamp)

fcuts = [1 30]; % Transition bands
mags = [0 1];
devs = [10^(-40/20) 0.01];% 1% passband ripple and 40-dB attenuation of stopband

[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);
b = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale'); %Filter coefficients to match the required parameters
filtout = filtfilt(b,1,IN); %Zero-phase filtering