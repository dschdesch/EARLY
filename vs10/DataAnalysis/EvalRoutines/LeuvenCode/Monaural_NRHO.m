function [] = Monaural_NRHO(AnID,SeqID)

%Load the monaural NRHO data described by the animal ID (AnID) and sequence
%ID (SeqID) and does a coincidence analysis between this and another
%fiber's response to the same noise stimulus.

RevCorDir = ('C:\work\BinRev\MonRev\');

ANdata = load('C:\LeuvenDataAnalysis\MonRevPopData.mat');
ANdfs = vertcat(ANdata.AN(:).df);
ANcfs = vertcat(ANdata.AN(:).cf);
ANfreqs = ANdfs;
ANfreqs(ANcfs>2000) = ANcfs(ANcfs>2000)/1000;

if ~iscell(SeqID)
    temp = SeqID;clear SeqID;
    SeqID = cell(1,1);
    SeqID{1} = temp;clear temp;
end
nrdata = length(SeqID);
for i = 1:nrdata
    nrhodata{i} = dataset(AnID,[SeqID{i} '-NRHO']);
end
nrhoseed = nrhodata{1}.Stimulus.StimParam.Rseed;

%load this unit's revcor
ind = regexp(SeqID,'-');
Revcorfilename = [AnID '_' SeqID{1}(1:ind{1}-1) '_MonRev.mat'];
revcordata = load(fullfile(RevCorDir,Revcorfilename));

dif = abs(log2(ANfreqs./revcordata.ParamsOut.domfreq));
indx = find(dif<=0.25);%potential matches

%Now find those with the same random seed
counti = 1;
for i =1:length(indx)
    foundit = 0;
    countj = 1;
    testdata = load(fullfile(RevCorDir,ANdata.AN(indx(i)).filename));
    for j = 1:length(testdata.ParamsOut.StimulusInfo)
        revcorseed = testdata.ParamsOut.StimulusInfo{j}.StimParam.RandomSeed;
        if revcorseed==nrhoseed
            foundit = 1;
            founddata{counti}.animalID = testdata.ParamsOut.DatasetIDs.Animal;
            founddata{counti}.seqID{countj} = testdata.ParamsOut.DatasetIDs.noise{j};
            founddata{counti}.revcor = testdata.ParamsOut.h1filt;
            countj = countj+1;
        end
    end
    if foundit
        counti = counti+1;
    end
end
dt = revcordata.ParamsOut.Time(2)-revcordata.ParamsOut.Time(1);
%Now do the analysis between this unit and the found units with the same
%noise seed and within the half octave range of DFs
%NRHO data will always be "ipsi" and static noise data "contra"
clear ind;
for i = 1:length(founddata)
    difcor{i} = xcorr(revcordata.ParamsOut.h1filt,founddata{i}.revcor);
    [dum,ind(i)] = max(difcor{i});
end
ITDx = -sort([(1:floor(numel(difcor{1})/2))*-dt 0 (1:floor(numel(difcor{1})/2))*dt]);
BD_ms = ITDx(ind);
BD_samps = numel(revcordata.ParamsOut.h1filt)-ind;

for i = 1:length(founddata)
    nCI{i} = zeros(1,21);
    for j = 1:length(founddata{i}.seqID)
         ntddata = dataset(founddata{i}.animalID,[founddata{i}.seqID{j} '-NTD']);
         ntdSPT = ntddata.SPT;
         ntdSPT = cat(2,ntdSPT{:});
         ntdSPT = ntdSPT(ntdSPT<=5000);
         ntdSPT = ntdSPT-BD_ms(i);
         for k = 1:length(nrhodata)
             for ii = 1:length(nrhodata{k}.SPT)
                 nrhoSPT = nrhodata{k}.SPT{ii};
                 nrhoSPT = nrhoSPT(nrhoSPT<=5000);
                 for jj =1:length(nrhoSPT)
                     dif = abs(nrhoSPT(jj)-ntdSPT);
                     nCI{i}(ii) = nCI{i}(ii)+numel(dif(dif<0.05));
                 end
             end
         end
    end
    %Now scale the spike rate appropriately as spikes per second
    nCIsps{i} = nCI{i}./(length(founddata{i}.seqID)*length(nrhodata)*5);
    %Normalize each predicted NRHO function to the max spike rate
    nCInorm{i} = nCIsps{i}./max(nCIsps{i});
end
foo;






































