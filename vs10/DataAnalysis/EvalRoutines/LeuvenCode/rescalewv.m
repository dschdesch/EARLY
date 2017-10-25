function [wvOUT,signalPOWER] = rescalewv(wvIN,SPL,outscale)
%--------------------------------------------------------------------------
%Inputs:
%wvIN: noise waveform samples
%SPL: overall dB SPL
%outscale: a string, either 'Pa' (Pascals), or 'V' (Volts - not done yet)

%Output:
%wvOUT: scaled waveform
%signalPOWER
%--------------------------------------------------------------------------

signalRMS = 10^(SPL/20)*20e-6;
signalPOWER  = signalRMS^2;
if isa(wvIN,'cell')
wvOUT = zeros(size(wvIN{1,1},1),size(wvIN,2));
Nchan = size(wvIN,2);
switch outscale
    case 'Pa'%Pascals
        for i = 1:Nchan
            wvOUT(:,i) = wvIN{:,i}./rms(wvIN{:,i});
            wvOUT(:,i) = wvOUT(:,i)*signalRMS;
        end
    case 'V'%Volts
        %To do
end
else 
wvOUT = zeros(size(wvIN,1),size(wvIN,2));
Nchan = size(wvIN,2);
switch outscale
    case 'Pa'%Pascals
        for i = 1:Nchan
            wvOUT(:,i) = wvIN(:,i)./rms(wvIN(:,i));
            wvOUT(:,i) = wvOUT(:,i)*signalRMS;
        end
    case 'V'%Volts
        %To do
end
end