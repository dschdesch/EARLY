function [NRHO] = nrho_analysis(contrachan,dsnrho)

ITD = dsnrho{1}.Stimulus.StimParam.delay;%delay was probably always applied to the left ear, but check it
ITDear = find(ITD~=0.0);
if isempty(ITDear)
    ITD = 0;%easy
else
    ITD = ITD(ITDear);
    if (contrachan==1 && ITDear==1) || (contrachan==2 && ITDear==2)
        ITD = -1*ITD;
    end
end

nrrho = 0;
nrreps = 0;
for i = 1:length(dsnrho)
    rho = dsnrho{i}.Stimulus.IndepVar.Values;
    rho = rho(~isnan(rho));
    nrreps = max(nrreps,size(dsnrho{i}.Data.SpikeTimes,2));
    nrrho = max(nrrho,length(rho));
end
spikerate = nan(length(dsnrho),nrrho,nrreps);
rhoval = nan(length(dsnrho),nrrho,nrreps);
splval = nan(length(dsnrho),nrrho,nrreps);

stimlen = dsnrho{1}.Stimulus.StimParam.burstDur;
for i = 1:length(dsnrho)
    rho = dsnrho{i}.Stimulus.IndepVar.Values;
    ind = find(~isnan(rho));
    rho = rho(ind);
    thisSPL = dsnrho{i}.Stimulus.StimParam.SPL;
    for j = 1:length(ind)
        nrreps = size(dsnrho{i}.Data.SpikeTimes,2);
        for k = 1:nrreps
            spikes = dsnrho{i}.Data.SpikeTimes{ind(j),k};
            spikes = spikes(spikes<=stimlen);
            spikerate(i,j,k) = length(spikes)/(stimlen/1000);
            rhoval(i,j,k) = rho(ind(j));
            splval(i,j,k) = thisSPL;
        end
    end
end
rhoval = rhoval(:);
rhoval = rhoval(~isnan(rhoval));
spikerate = spikerate(:);
spikerate = spikerate(~isnan(spikerate));
splval = splval(:);
splval = splval(~isnan(splval));

uspl = unique(splval);
for j = 1:length(uspl)
    NRHO(j).SPL = uspl(j);
    urho = unique(rhoval(splval==uspl(j)));
    NRHO(j).data.rate = nan(length(urho),1);
    NRHO(j).data.rho = nan(length(urho),1);
    for i = 1:length(urho)
        NRHO(j).data.rate(i) = mean(spikerate(rhoval==urho(i) & splval==uspl(j)));
        NRHO(j).data.rho(i) = urho(i);
    end
    NRHO(j).data.powerfit = fit_NRHO(NRHO(j).data.rho,NRHO(j).data.rate);
    NRHO(j).data.itd = ITD;
end
return;