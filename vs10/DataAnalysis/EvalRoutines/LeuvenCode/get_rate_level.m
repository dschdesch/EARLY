function [RLVout] = get_rate_level(dsspl)
splspikes = dsspl{1}.Data.SpikeTimes;
[nrcond,nrreps] = size(splspikes);
STmax = dsspl{1}.Stimulus.Special.BurstDur;
Freq = dsspl{1}.Stimulus.Special.CarFreq;
for ii = 1:nrcond
    SPL(ii) = dsspl{1}.Stimulus.IndepVar.Values(ii);
    if isnan(SPL(ii))
        SPL = SPL(1:ii-1);
        break;
    end
    ST = [splspikes{ii,:}];
    ST = ST(ST<=STmax);
    SPS(ii) = numel(ST)/nrreps/(STmax/1000);
    ST_c = mod(ST,1000/Freq)./(1000/Freq); %Spike times modulo stimulus frequency
    sum_sines = sum(sin((2*pi).*ST_c));%Sum of sines and cosines
    sum_cosines = sum(cos((2*pi).*ST_c));
    Phase(ii) = atan2(sum_sines,sum_cosines);%Mean phase
    Mag(ii) = sqrt(sum_sines^2 + sum_cosines^2)/numel(ST);%Vector strength
    Pval(ii)=exp(-length(ST)*(Mag(ii)^2));
end
RLVout = struct('Freq',Freq,'Level',SPL,'Rate',SPS,'SyncMag',Mag,'SyncPhase',Phase,'SyncPval',Pval,'SpikeTimes',{splspikes});
end