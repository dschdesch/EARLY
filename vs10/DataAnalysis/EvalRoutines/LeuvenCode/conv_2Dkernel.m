function [h2pst,h2Xpst] = conv_2Dkernel(stimulus,h2,h2X)

% stimulus = stimulus';
stim_pts = length(stimulus);
h2_pts = length(h2);

%% Perform convolution
% Adapted from "Kernel 2D Convolution of signal in MATLAB" posted by sebas
% at http://stackoverflow.com/questions/21686781/kernel-2d-convolution-of-signal-in-matlab

nbins = h2_pts+stim_pts-1;      % Number of bins for the full convolution
h2pst = zeros(nbins,2);
h2Xpst = zeros(nbins,1);
k = 1:nbins;
imin = k+1-min(k,h2_pts);
imax = min(k,stim_pts);
jmin = max(1,k-stim_pts+1);
jmax = min(k,h2_pts);

for i = 1:nbins
    ivals = imax(i):-1:imin(i);
    jvals = jmin(i):jmax(i);
    for k = 1:2
        h2pst(i,k) = stimulus(ivals,k)'*(h2(jvals,jvals,k)*stimulus(ivals,k));
    end
    h2Xpst(i) = stimulus(ivals,1)'*(h2X(jvals,jvals)*stimulus(ivals,2));
end

%% Extract the central portion of the output, of length stim_pts
i_start = round(h2_pts/2);
h2pst = h2pst(i_start:i_start+stim_pts-1,:);
h2Xpst = h2Xpst(i_start:i_start+stim_pts-1);
return;