function [h1] = rescaleh1(h1,outscale)

n = size(h1,2);
mag = zeros(size(h1));
phs = zeros(size(h1));
for chan = 1:n
    s = fft(h1(:,chan));
    mag(:,chan) = abs(s);
    phs(:,chan) = angle(s);
end

switch outscale
    case 'att'
        mag = mag/max(mag(:));
        h1 = real(ifft(mag.*exp(1i*phs)));
    case 'SpSpPa'
end
end