function [spikesout,Bspikesout,wv,dt,h1ipsi,h1contra] = Mon_Conv_get_coincident_spikes(indA,setB,thr)

ANdata = load('C:\LeuvenDataAnalysis\MonRevPopData.mat');

spikesout=[];
wv=[];
dt=[];
h1ipsi=[];
h1contra =[];
Bspikesout = [];
%% Get the data and error check
dataA = load(fullfile('C:\work\BinRev\MonRev',ANdata.AN(indA).filename));
% if isfield(dataA.ParamsOut,'noiseSPLs')
SPLind = find(unique(dataA.ParamsOut.noiseSPLs)==70);
if isempty(SPLind)
    return;
end
% else
%     if isfield(dataA.ParamsOut.StimulusInfo{1}.StimParam,'SPL')
%         if dataA.ParamsOut.StimulusInfo{1}.StimParam.SPL(1)==70
%             SPLind =1;
%         else
%             disp('Wrong SPL! These data were not included in the analysis');
%             return;
%         end
%     else
%         error('Something went wrong here!');
%     end
% end
h1ipsi.h1ipsi{1} = dataA.ParamsOut.h1filt(:,SPLind);
h1ipsi.h1ipsimag{1} = dataA.ParamsOut.h1mag(:,SPLind);
h1ipsi.h1ipsiphase{1} = dataA.ParamsOut.h1phase(:,SPLind);
h1ipsi.h1ipsizscore{1} = dataA.ParamsOut.h1zscore(:,SPLind);
h1ipsi.h1ipsibw{1} = dataA.ParamsOut.h1bw(:,SPLind);
h1ipsi.h1ipsiqvals{1} = dataA.ParamsOut.qvals(:,SPLind);
h1ipsi.h1ipsidomfreq{1} = dataA.ParamsOut.domfreq(SPLind);


animalA = dataA.ParamsOut.DatasetIDs.Animal;
seqindxA = find(ANdata.AN(indA).noiseseed~=111);
%Check that these are all at the same SPL
for i = 1:length(seqindxA)
    if dataA.ParamsOut.StimulusInfo{seqindxA(i)}.StimParam.SPL(1) ==70
        checkSPLA(i) = 1;
    else
        checkSPLA(i) = 0;
        disp 'Wrong SPL found!: excluded from analysis';
    end
end
if ~any(checkSPLA)
    return;
end
seqindxA = seqindxA(logical(checkSPLA));
seqidA = dataA.ParamsOut.DatasetIDs.noise(seqindxA);
for i = 1:length(seqidA)
    try
        dsA{i} = dataset(animalA,[seqidA{i} '-NTD']);
    catch
        dsA{i} = dataset(animalA,[seqidA{i} '-NSPL']);
    end
    spikedataA{i} = struct(dsA{i});
    [wvA{i},dt] = StimSam(dsA{i},1);
end


%Lots of error checking here....
for i = 1:length(setB)
    dataB{i} = load(fullfile('C:\work\BinRev\MonRev',ANdata.AN(setB(i)).filename));
    
    SPLind = find(unique(dataB{i}.ParamsOut.noiseSPLs)==70);
    if isempty(SPLind)
        disp ('This dataset doesn''t contain a 70-dB condition');
        continue;
    end
    
    h1contra.h1contra{i} = dataB{i}.ParamsOut.h1filt(:,SPLind);
    h1contra.h1contramag{i} = dataB{i}.ParamsOut.h1mag(:,SPLind);
    h1contra.h1contraphase{i} = dataB{i}.ParamsOut.h1phase(:,SPLind);
    h1contra.h1contrazscore{i} = dataB{i}.ParamsOut.h1zscore(:,SPLind);
    h1contra.h1contrabw{i} = dataB{i}.ParamsOut.h1bw(:,SPLind);
    h1contra.h1contraqvals{i} = dataB{i}.ParamsOut.qvals(:,SPLind);
    h1contra.h1contradomfreq{i} = dataB{i}.ParamsOut.domfreq(SPLind);
    
    animalB{i} = dataB{i}.ParamsOut.DatasetIDs.Animal;
    seqidB{i} = dataB{i}.ParamsOut.DatasetIDs.noise;
    ok = 0;
    for j = 1:length(seqidB{i})
        try
            dsB{j} = dataset(animalB{i},[seqidB{i}{j} '-NTD']);
        catch
            dsB{j} = dataset(animalB{i},[seqidB{i}{j} '-NSPL']);
        end
        switch dsB{j}.stimtype
            case 'NTD'
                if dsB{j}.Stimulus.StimParam.SPL(1)==70
                    ok = 1;
                else
                    ok = 0;
                end
            case 'NSPL'
                if ismember(70,dsB{j}.Stimulus.IndepVar.Values)
                    ok = 1;
                else
                    error ('Something is wrong here! Looks like the NSPL data doesn''t contain 70-dB as a condition');
                end
            otherwise
                error ('Unknown stimulus type!')
        end
        if ok==1
            spikedataB{i} = struct(dsB{j});
            if strcmp(dsB{j}.stimtype,'NSPL')
                indspl = find(spikedataB{i}.Stimulus.IndepVar.Values==70);
                tempdata = spikedataB{i}.Data.SpikeTimes;
                spikedataB{i}.Data.SpikeTimes = {tempdata{indspl,1:end}};
            end
            break;
        end
    end
end
[wvB,dt] = StimSam(dsB{1},1);

%% The meat of the function

%Figure out how many reps we have to work with, usually 25 or 20
for j = 1:length(setB)
    tempnrreps(j) = size(spikedataB{j}.Data.SpikeTimes,2);
end
nrrepsB = min(tempnrreps);clear tempnrreps;


for j = 1:length(dsA)
    tempnrreps(j) = size(spikedataA{j}.Data.SpikeTimes,2);
end
nrrepsA = min(tempnrreps);clear tempnrreps;
% 
% rand('state',sum(100.*clock));
% delayB = 3*rand(length(setB),1);
% delayA = 3*rand;
for u = 1:length(dsA)%For each "ipsi" noise response (i.e., different non-111 noise seeds)
    for i = 1:nrrepsA
        spikesA{i} = spikedataA{u}.Data.SpikeTimes{i};%+delayA;
        for j = 1:nrrepsB
            ciind{u,i,j} = zeros(size(spikesA{i}));
            for k = 1:length(setB)
                spikesB = spikedataB{k}.Data.SpikeTimes{j};%+delayB(k);
                [dummy logind] = findcoincidentspikes(spikesA{i},spikesB,0.2);
                ciind{u,i,j} = ciind{u,i,j}+logind;
            end
            cispikes{u,i,j} = spikesA{i}(ciind{u,i,j}>=thr);
        end
    end
end

%Outputs
dt = dt*1e-3;%micro to milli-seconds
for u = 1:length(dsA)
    spikesout{u} = cat(2,cispikes{u,:,:});%all spikes from all reps
    wv{u} = [wvA{u},wvB];%arranged ["ipsi","contra"]
end
%added 8/2/16 - also return a vector of all the spikes from the "B set"....
%to make sure the low pass in the revcors is actually coming from the
%coincidence.
Bspikesout = [];
for k = 1:length(setB)
    Bspikesout = [Bspikesout [spikedataB{k}.Data.SpikeTimes{1:nrrepsB}]];
end

return;