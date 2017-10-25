function [h2out,h2Xout] = get_h2_kernel(ST,wv,params,signalPOWER,Dur,outscale)
%2-nd order Wiener kernel from spike times and GWN input
%------------------
%Inputs: 
%ST: spike times (single column vector containing all spikes in response to all repetitions of stimulus wv);
%wv: stimulus waveform (Gaussian white noise)
%dt: stimulus sample period
%Nwin: number of samples to consider in the reverse correlation (i.e., Nwin = round(TimeWindow/dt), where TimeWindow is typically 20-30 ms);
%signalPOWER: Pa^2
%Dur: Total duration of the signal used (in milliseconds). This is summed across all reps of all noise tokens.
%outscale: a string, 'SpSpPa2' (Spikes per second per pascal^2)
%-----------------------
%Output:
%h2out: the scaled second order kernel
%h2Xout: the scaled cross-second order kernel (If there are two inputs to the system)
%---------------------------
% M Sayles. 2016.
%--------------------------
persistent yi;

Nwin = params.Nwin;
dt = params.dt;
Nchan = size(wv,2);
Nspikes = length(ST);%column vector of spike times
mfilelist = dbstack;
if ~strcmp(mfilelist(2).name,'bootstrap_h1_h2')
    dZERO = round(ST/dt);%sample for each spike in response to this stimulus
    dMAX = dZERO-Nwin+1;%corresponding max time lag samples
    yi = fliplr(cell2mat(arrayfun(@colon,dMAX,dZERO,'uniform',false)));%matrix of stimulus samples to get
end
h2out = double(zeros(Nwin,Nwin,Nchan));
wvh2 = cell(1,Nchan);
for ii = 1:Nchan
    wvwv = double(wv(:,ii));
    wvh2{ii} = wvwv(yi);%Get those samples
    autocorr = xcorr(wvwv,Nwin-1,'unbiased');%stimulus autocorrelation
    autocorrmatrix = toeplitz(autocorr(Nwin:end)); %Toeplitz matrix
    h2out(:,:,ii) = double(((wvh2{ii}'*wvh2{ii})-(Nspikes*autocorrmatrix))/Nspikes);%R2 is the matrix product corrected for the input signal correlation and number of spikes
end
if Nchan==2
    h2Xout = double((wvh2{1}'*wvh2{2})/Nspikes);%No need to correct for signal correlation, since theoretically they are independent
else
    h2Xout = [];
end

%Scale the kernel(s) appropriately
switch outscale
    case 'SpSpPa2'
        MeanSpikeRate = Nspikes/(Dur/1000);
        h2out = h2out*(MeanSpikeRate/(2*signalPOWER^2));
        if ~isempty(h2Xout)
            h2Xout = h2Xout*(MeanSpikeRate/(2*signalPOWER^2));
        end
    case 'SpSpV2'
        %To do
    otherwise
        error('Check input arguments')
end
return;