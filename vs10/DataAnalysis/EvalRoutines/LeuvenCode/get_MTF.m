function [] = get_MTF(AnNum,DataID,UnNum)


% AnNum  = 'L15018';
% DataID = '2-10';
% UnNum = '2';

clc;

datapath = 'C:\Users\Mark\Dropbox\SPON\MTFs';

dataobj = dataset(AnNum,DataID);
data = struct(dataobj);

StimLen = data.Stimulus.Special.BurstDur(1);
ModDepth = data.Stimulus.StimParam.indiv.stim{1}.modpercent;
SPL = data.Stimulus.StimParam.indiv.stim{1}.spl;

fm = data.Stimulus.IndepVar.Values;
fm = fm(~isnan(fm));
nrcond = length(fm);
nrep = data.Sizes.Nrep;

fcar = data.Stimulus.Special.CarFreq;
for i = 1:nrcond
    ModP = 1/fm(i);%Modulation period in seconds
    
    Spikes = cat(2,data.Data.SpikeTimes{i,:});
    Spikes = Spikes(Spikes<StimLen);
    Spikes = Spikes/1000; %ms --> s
    
    raster_spikes = data.Data.SpikeTimes(1:nrcond,:);
    
    Nspikes = length(Spikes);
    
    spikerate(i) = Nspikes/(StimLen/1000)/nrep;
    
    if Nspikes
        Spikes_cycles = mod(Spikes,ModP)/ModP;
        
        %Mean phase
        sines = sin(2*pi*Spikes_cycles);
        cosines = cos(2*pi*Spikes_cycles);
        meanphase(i) = atan2(sum(sines),sum(cosines));%Mean phase
        meanphase(i) = meanphase(i)/(2*pi); %radians --> cycles
        
        %Vector strength
        VS(i) = sqrt(sum(sines)^2 + sum(cosines)^2)/Nspikes;%Vector strength
        
        %P-value
        Pval(i) = exp(-Nspikes*(VS(i)^2));
        
        %synchronized rate
        syncrate(i) = spikerate(i)*VS(i);
    else
        meanphase(i) = nan;
        VS(i) = nan;
        Pval(i) = nan;
    end
end

if ~exist(fullfile(datapath,[AnNum '_' UnNum]),'dir');
    mkdir(fullfile(datapath,[AnNum '_' UnNum]));
end

filename2save = [AnNum '_' UnNum '_MTF_' sprintf('%2.0f_%2.0fdB',ModDepth,SPL) '.mat'];

save(fullfile(datapath,[AnNum '_' UnNum],filename2save)...
    ,'spikerate','meanphase','VS','Pval','syncrate','fm','fcar','nrep','nrcond','raster_spikes');
disp done
clear


















