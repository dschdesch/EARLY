function [scramspikes] = scramble_spikes(spikes)
% shuffle spike train by randomizing the ISI's
D = diff([0; sort(spikes')]); % first sort to make sure diff gives the ISI's
NInt = length(D);
scramspikes = cumsum(D(randperm(NInt)));
end