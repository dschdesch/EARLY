function [spikesout,wv,dt,h1contra,h1ipsi] = Mon_Jitter_get_coincident_spikes(CFind)

ANdata = load('C:\LeuvenDataAnalysis\MonRevPopData.mat');

spikesout=[];
wv=[];
dt=[];
h1ipsi=[];
h1contra =[];

%% Lots of error checking here....
for i = 1:length(CFind)
    data{i} = load(fullfile('C:\work\BinRev\MonRev',ANdata.AN(CFind(i)).filename));
    SPLind = find(unique(data{i}.ParamsOut.noiseSPLs)==70);
    if isempty(SPLind)
        disp ('This dataset doesn''t contain a 70-dB condition');
        continue;
    end
    
    h1ipsi.h1ipsi{i} = data{i}.ParamsOut.h1filt(:,SPLind);
    h1ipsi.h1ipsiphase{i} = data{i}.ParamsOut.h1phase(:,SPLind);
    h1ipsi.h1ipsimag{i} = data{i}.ParamsOut.h1mag(:,SPLind);
    h1ipsi.h1ipsizscore{i} = data{i}.ParamsOut.h1zscore(:,SPLind);
    h1ipsi.h1ipsibw{i} = data{i}.ParamsOut.h1bw(:,SPLind);
    h1ipsi.h1ipsiqvals{i} = data{i}.ParamsOut.qvals(:,SPLind);
    h1ipsi.h1ipsidomfreq{i} = data{i}.ParamsOut.domfreq(SPLind);
    
    h1contra.h1contramag{i} = data{i}.ParamsOut.h1mag(:,SPLind);
    h1contra.h1contrazscore{i} = data{i}.ParamsOut.h1zscore(:,SPLind);
    h1contra.h1contrabw{i} = data{i}.ParamsOut.h1bw(:,SPLind);
    h1contra.h1contraqvals{i} = data{i}.ParamsOut.qvals(:,SPLind);
    h1contra.h1contradomfreq{i} = data{i}.ParamsOut.domfreq(SPLind);
    
    animal{i} = data{i}.ParamsOut.DatasetIDs.Animal;
    seqid{i} = data{i}.ParamsOut.DatasetIDs.noise;
    ok = 0;
    for j = 1:length(seqid{i})
        try
            ds{j} = dataset(animal{i},[seqid{i}{j} '-NTD']);
        catch
            ds{j} = dataset(animal{i},[seqid{i}{j} '-NSPL']);
        end
        switch ds{j}.stimtype
            case 'NTD'
                if ds{j}.Stimulus.StimParam.SPL(1)==70 && ds{j}.Stimulus.StimParam.RandomSeed==111
                    ok = 1;
                else
                    ok = 0;
                end
            case 'NSPL'
                if ismember(70,ds{j}.Stimulus.IndepVar.Values) && ds{j}.Stimulus.StimParam.RandomSeed==111
                    ok = 1;
                else
                    error ('Something is wrong here! Looks like the NSPL data doesn''t contain 70-dB as a condition');
                end
            otherwise
                error ('Unknown stimulus type!')
        end
        if ok==1
            spikedata{i} = struct(ds{j});
            if strcmp(ds{j}.stimtype,'NSPL')
                indspl = find(spikedata{i}.Stimulus.IndepVar.Values==70);
                tempdata = spikedata{i}.Data.SpikeTimes;
                spikedata{i}.Data.SpikeTimes = {tempdata{indspl,1:end}};
            end
            noiseds{i} = ds{j};
            break;
        end
    end
    clear ds;
end
[wv,dt] = StimSam(noiseds{1},1);
dt = dt/1000;
nsamps = floor(5000/dt);
wv = wv(1:nsamps,:);
tempwv(:,1) = wv(1:floor(nsamps/2),1);
tempwv(:,2) = wv(floor(nsamps/2)+1:end,1);
clear wv;
wv = tempwv;

%% Split everything in two (first half and second half of each trial are considered independent noises)
maxST = 5000;
minST = 50;
for i = 1:length(spikedata)
    spikedataA{i}.Data.SpikeTimes = spikedata{i}.Data.SpikeTimes;
    for j = 1:length(spikedataA{i}.Data.SpikeTimes)
        spikedataA{i}.Data.SpikeTimes{j} = spikedataA{i}.Data.SpikeTimes{j}(spikedataA{i}.Data.SpikeTimes{j}<=maxST/2 & spikedataA{i}.Data.SpikeTimes{j}>minST);
    end
    spikedataB{i}.Data.SpikeTimes = spikedata{i}.Data.SpikeTimes;
    for j = 1:length(spikedataB{i}.Data.SpikeTimes)
        spikedataB{i}.Data.SpikeTimes{j} = spikedataB{i}.Data.SpikeTimes{j}(spikedataB{i}.Data.SpikeTimes{j}>maxST/2 & spikedataB{i}.Data.SpikeTimes{j}<=maxST)-(maxST/2);
        spikedataB{i}.Data.SpikeTimes{j} = spikedataB{i}.Data.SpikeTimes{j}(spikedataB{i}.Data.SpikeTimes{j}>minST);
    end
end
clear spikedata;


% %% Debugging purposes only
% tempA = spikedataA;
% tempB = spikedataB;
% spikedataA = tempB;
% spikedataB = tempA;
% clear tempA;clear tempB;
% wv = wv(:,[2 1]);

%% The meat of the function

%Figure out how many reps we have to work with, usually 25 or 20
for j = 1:length(spikedataB)
    tempnrreps(j) = size(spikedataB{j}.Data.SpikeTimes,2);
end
nrreps = min(tempnrreps);clear tempnrreps;

% Apply random delays to the "B" spike trains
nrinputs = length(spikedataB);
for ii = 1:10
    %     randn('state',sum(100*clock));
    %     delays = (0.5*randn(nrinputs,1))+0.5;
    
        rand('state',sum(100*clock));
        delays = -0.5 + 2.*rand(nrinputs,1);
        delays = delays(randperm(nrinputs));
    
%     delays = zeros(nrinputs,1);
    %Make the delays an integer number of samples
    delays = round(delays/dt)*dt;
    for i = 1:nrinputs
        for j = 1:nrreps
            delayedspikes{i,j} = spikedataB{i}.Data.SpikeTimes{j}+delays(i);
        end
    end
    
    %Monaural coincidence stage
%     count = 1;
    for i = 1:nrinputs-1
        AA = [];
        DD = [];
        for j = i+1:nrinputs
            AA = [AA [spikedataA{j}.Data.SpikeTimes{1:nrreps}]];
            DD = [DD [delayedspikes{j,1:nrreps}]];
        end
        [mcispikesA{i}] = findcoincidentspikes([spikedataA{i}.Data.SpikeTimes{1:nrreps}],AA,0.2);
        [mcispikesB{i}] = findcoincidentspikes([delayedspikes{i,1:nrreps}],DD,0.2);
        
        %         count = count+1;
        
    end
    
    %Binaural coincidence stage
    for i = 1:length(mcispikesA)
        
            [bcispikes{i}] = findcoincidentspikes(mcispikesA{i},mcispikesB{i},0.2);
        
    end
    
    spikesout{ii} = [bcispikes{:,:}];
    clear bcispikes;
    clear delayedspikes;
    
    %Adjust the phase and the timing (same thing) for the "contra" side
    %based on the random delays per input
    ffax = data{1}.ParamsOut.h1ffax;
    for i = 1:nrinputs
        if delays(i)>0
            delaysamps = delays(i)/dt;
            h1contra.h1contra{i,ii} = [zeros(round(delaysamps),1); h1ipsi.h1ipsi{i}(1:end-round(delaysamps))];
        elseif delays(i)<0
            delaysamps = abs(delays(i)/dt);
            h1contra.h1contra{i,ii} = [h1ipsi.h1ipsi{i}(round(delaysamps)+1:end);zeros(round(delaysamps),1)];
        else
            h1contra.h1contra{i,ii} = h1ipsi.h1ipsi{i};
        end
        h1contra.h1contraphase{i,ii} = h1ipsi.h1ipsiphase{i}+(2*pi.*ffax*delays(i))';
    end
end



return;