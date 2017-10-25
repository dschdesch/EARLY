function [waveform] = MovingNoise(fLow,fHigh,Fs,Dur,Rho,Seed,StartITD,Speed)

%Get Basic Noise Waveforms
W = MakeGaussNoiseBand(fLow, fHigh, Dur, Fs, 2, Rho, StartITD, 0, Seed,1);

FsN = (Fs*(1-(Speed/1000000)));        %calculate new sample rate [samples/s]

FsN = (max(1,FsN));
[FsN_temp,Fs_temp] = rat( FsN/Fs, 2.5e-5 );   %reduce up/down sample factors
FsN_temp = max(1,FsN_temp);
w = resample(W(:,2),FsN_temp,Fs_temp);          %resample waveform and avoid fsr_temp=0
Nw = numel(w);                              %#elements in wfr
FsN = Fs * (FsN_temp/Fs_temp);                %recalculate new sample rate



