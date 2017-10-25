function [spikes,spikerate,meanspikerate] = getPoissonSpikes(spikeprob,dt,N)

n = length(spikeprob);
T = 0:dt:dt*(n-1); % Time axis
dur = dt*n;
refper = 1e-3;
spikes = cell(1,N);
spikerate = zeros(1,N);
for i = 1:N
    spikes = rand(length(spikeprob),1)<spikeprob;%Poisson spikes (0s and 1s for each sample)
    spikes = T(logical(spikes));%Spike times
    if ~isempty(spikes)
        intervals = [spikes(1); diff(spikes)'];%interspike intervals in ms
        %calculate spike rate after applying refractoriness to the
        %spike times
        spikerate(i) = sum(intervals>refper)/dur; %--> ms to s
    else
        spikerate(i) = 0;
    end
end
meanspikerate = mean(spikerate);