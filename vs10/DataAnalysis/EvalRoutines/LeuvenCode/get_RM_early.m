function [] = get_RM_early(AnNum,DataID,UnNum)


clc;

datapath = 'C:\Users\Mark\Dropbox\SPON\RMs';
nrdata = length(DataID);


D = read(dataset,AnNum,DataID);
data.Stimulus = D.stimparam;
data.SpikeTimes = D.spiketimes;

SPL = data.Stimulus.SPL;
freqs = data.Stimulus.Fcar;

burstdur = data.Stimulus.BurstDur;
repdur = data.Stimulus.ISI;

Nreps = data.Stimulus.Nrep;

nrconds = length(SPL);

for j = 1:nrconds
    Spikes = cat(2,data.SpikeTimes{j,:});
    OnSpikes = Spikes(Spikes<burstdur);
    OffSpikes = Spikes(Spikes>=burstdur);
    OnRate(j) = length(OnSpikes)/(burstdur/1000);
    OnRate(j) = OnRate(j)/Nreps;
    OffRate(j) = length(OffSpikes)/((repdur-burstdur)/1000);
    OffRate(j) = OffRate(j)/Nreps;
end

nrSPLs = length(unique(SPL));
nrFreqs = length(unique(freqs));


OnRate = reshape(OnRate,nrFreqs,nrSPLs);
OffRate = reshape(OffRate,nrFreqs,nrSPLs);
levels = reshape(SPL,nrFreqs,nrSPLs);
freqs = reshape(freqs,nrFreqs,nrSPLs);


if ~exist(fullfile(datapath,[AnNum '_' UnNum]),'dir');
    mkdir(fullfile(datapath,[AnNum '_' UnNum]));
end

filename2save = [AnNum '_' UnNum '_RM.mat'];

save(fullfile(datapath,[AnNum '_' UnNum],filename2save)...
    ,'OnRate','OffRate','levels','freqs');
disp done
clear;



