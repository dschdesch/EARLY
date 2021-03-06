function [filtout] = my_low_pass(IN,dt,flag,varargin)
if isempty(varargin)
    fcuts = [4500 5500]; % Transition bands
else
    fcuts = [varargin{1} varargin{2}];
end
fsamp = 1000/dt; %sample frequency
mags = [1 0];
devs = [0.05 0.05];% 5% ripple
[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);
b = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale'); %Filter coefficients to match the required parameters
if flag==1
    filtout = filtfilt(b,1,IN); %Zero-phase filtering
elseif flag==2
    L = length(b);
    M = L/2;
    xx = linspace(-M,M,L);
    [x,y] = meshgrid(xx);
    r = sqrt( x.^2 + y.^2 );
    b2D = zeros(L);
    b2D(r<=M) = interp1(xx,b,r(r<=M));
    b2D = b2D./sum(b2D(:));
    filtout = conv_fft2(IN,b2D,'same');
end