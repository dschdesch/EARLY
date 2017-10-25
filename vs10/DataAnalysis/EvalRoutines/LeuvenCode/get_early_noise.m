function [waveform,dt,signalPOWER] = get_early_noise(nchan,corrchan,burstdur,noiseseed,cutoffs,fs,SPL)


[noiseseed, Nmax] = SetRandState(noiseseed);
noiseseed2 = rem(noiseseed+31^2, Nmax);
for idx = 1:nchan
    if idx==corrchan
        thisseed = noiseseed;
    else
        thisseed = noiseseed2;
    end
    [waveform(:,idx),dt] = local_NoiseSpec(fs, burstdur, thisseed, cutoffs);
end
[waveform,signalPOWER] = rescalewv(waveform,SPL,'Pa');

    function [wv,dt] = local_NoiseSpec(fs, burstdur, thisseed, cutoffs)
        df = 1e3/burstdur; % freq spacing in Hz
        Nsam = round(burstdur*fs*1e-3);
        dt = 1e3/fs;
        SetRandState(thisseed);
        Buf = Nsam*randn(Nsam,2)*[1; 1i]; % the factor Nsam anticipates real(ifft(Buf)) as time waveform
        Buf(1) = Buf(1)/2; % compensate for doubling the DC
        Freq = Xaxis(Buf,df);
        Buf(Freq<cutoffs(1) | Freq>cutoffs(2)) = 0;
        wv = real(ifft(Buf));
    end
end