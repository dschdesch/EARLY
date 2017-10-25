function [S] = Noise_IAP(dt,dur,IPD,ITD)

sr = 1/dt;
nf = sr/2;
Nsamps = ceil(dur/dt);

NFFT = 2^nextpow2(Nsamps);
f = linspace(0,nf,NFFT/2)';

S(:,1) = randn(Nsamps,1);

spec = fft(S,NFFT);
p = (2*pi)*(IPD + (f*(ITD/1000)));

spec(1:NFFT/2) = spec(1:NFFT/2).*(exp(sqrt(-1)*p));%real part
spec(NFFT/2+1:NFFT) = spec(NFFT/2+1:NFFT).*(exp(sqrt(-1)*fliplr(p)));%imaginary part

temp = real(ifft(spec,NFFT));
S(:,2) = temp(1:Nsamps);
[c,lags] = xcorr(S(:,1),S(:,2));
figure(1);plot(lags*dt,c);

return;