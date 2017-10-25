function [h0,h1,Time,Nspikes] = get_h0_h1_kernel(ST,wv,params,signalPOWER,Dur,outscale)
dt = params.dt;
Nwin = params.Nwin;
Nchan = size(wv,2);
Nsamples = size(wv,1);
Nspikes = length(ST);
tx = (0:Nsamples)*dt; %Binning edges in ms ...
pst = histc(ST,tx);%Make PSTH with same sampling frequency as stimulus
h1 = zeros(Nwin,Nchan);
for ii = 1:Nchan
    temp = flipud(xcorr(wv(:,ii),pst,Nwin))/Nspikes;%Get the revcor (h1 kernel) - cross correlation of the pst and the stimulus waveform
    h1(:,ii) = temp((1:Nwin)+(Nwin+1),:);clear temp;
end
Time = tx(1:Nwin);
%Scale the kernel appropriately
switch outscale
    case 'SpSpPa'%Spikes per Second per Pascal
        h0 = Nspikes/(Dur/1000);%spikes per second
        h1 = h1*(h0/signalPOWER);%spike per second per pascal
    case 'SpSpV'%Spikes per Second per Volt
    otherwise
        error('Check input arguments')
end
return;