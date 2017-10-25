function [IRate,MRate] = predictRate(varargin)

%Predict firing rate of a unit in response to noise input based on h1, h2
%and h2X. Returns the instantaneous rate vector, and a the mean rate as a
%scalar.

wv = varargin{1};
h0 = varargin{2};
h1 = varargin{3};
h2 = varargin{4};
if nargin==5
    h2X = varargin{5};
elseif nargin==4
    h2X = [];
else
    error('Check input arguments');
end



end

figure;plot(0.2*(MeanSpikeRate+conv(wv,ParamsOut.h1.h1,'same')+conv(max(0,wv).^2,abs(U(:,1)),'same')+conv(max(0,wv).^2,-abs(U(:,2)),'same')))