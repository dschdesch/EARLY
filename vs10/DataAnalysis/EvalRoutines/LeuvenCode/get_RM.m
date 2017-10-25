function [] = get_RM(AnNum,DataID,UnNum)


clc;

datapath = 'C:\Users\Mark\Dropbox\SPON\RMs';
nrdata = length(DataID);
for i = 1:nrdata
    dataobj(i) = dataset(AnNum,DataID{i});
    data(i) = struct(dataobj(i));
    SPL(i) = data(i).Stimulus.StimParam.indiv.stim{1}.spl;
    burstdur(i) = data(i).Stimulus.Special.BurstDur;
    repdur(i) = data(i).Stimulus.Special.RepDur;
    freqs(i,:) = data(i).Stimulus.IndepVar.Values;
    nrconds = length(freqs(i,:));
    levels(i,:) = repmat(SPL(i),1,nrconds);
    for j = 1:nrconds
        Spikes = cat(2,data(i).Data.SpikeTimes{j,:});
        OnSpikes = Spikes(Spikes<burstdur(i));
        OffSpikes = Spikes(Spikes>=burstdur(i));
        OnRate(i,j) = length(OnSpikes)/(burstdur(i)/1000);
        OffRate(i,j) = length(OffSpikes)/((repdur(i)-burstdur(i))/1000);
    end
end

[~,ix] = sort(levels(:,1));
OnRate = OnRate(ix,:);
OffRate = OffRate(ix,:);
levels = levels(ix,:);

if ~exist(fullfile(datapath,[AnNum '_' UnNum]),'dir');
    mkdir(fullfile(datapath,[AnNum '_' UnNum]));
end

filename2save = [AnNum '_' UnNum '_RM.mat'];

save(fullfile(datapath,[AnNum '_' UnNum],filename2save)...
    ,'OnRate','OffRate','levels','freqs');
disp done
clear;


















