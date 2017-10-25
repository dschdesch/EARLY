function [varargout] = bootstrap_h1_h2(spikes,wv,params,Dur,SPL)

Nboot = params.Nboot;
Nwin = params.Nwin;
Nchan = size(wv,2);
maxorder = params.MaxOrder;
nullh1 = zeros(Nwin,Nchan,Nboot);
if maxorder==2
    nullh2 = zeros(Nwin,Nwin,Nchan,Nboot);
    if Nchan==2
        nullh2X = zeros(Nwin,Nwin,Nboot);
    end
    h = waitbar(0,sprintf('Bootstrapping first- and second-order kernels at %2.0f dB SPL: Please wait...',SPL));
else
    h = waitbar(0,sprintf('Bootstrapping first-order kernel at %2.0f dB SPL: Please wait...',SPL));
end

spec = fft(wv);
for j = 1:Nboot
    %For each bootstrap re-sampling - make a new noise with the same
    %magnitude spectrum but a randomized phase spectrum
    spec = spec.*exp(1i*2*pi*randn(size(wv)));
    rwv = real(ifft(spec));
    [rwv,signalPOWER] = rescalewv(rwv,SPL,'Pa');
    [dummy,nullvalsh1,dummy,dummy] = get_h0_h1_kernel(spikes,rwv,params,signalPOWER,Dur,'SpSpPa');
    nullh1(:,:,j) = nullh1(:,:,j)+nullvalsh1;
    if maxorder==2
        [nullvalsh2,nullvalsh2X] = get_h2_kernel(spikes,rwv,params,signalPOWER,Dur,'SpSpPa2');
        nullh2(:,:,:,j) = nullh2(:,:,:,j)+nullvalsh2;
        if Nchan==2
            nullh2X(:,:,j) = nullh2X(:,:,j)+nullvalsh2X;
        end
    end
    clear nullvalsh1 nullvalsh2 nullvalsh2X;
    waitbar(j/Nboot);
end
varargout{1} = nullh1;clear nullh1;
if maxorder==2
    varargout{2} = nullh2;clear nullh2;
    if Nchan==2
        varargout{3} = nullh2X;clear nullh2X;
    end
end
close (h);