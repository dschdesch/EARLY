function [] = update_bbCPCD

BinRevDir = 'C:\work\BinRev\BinRev';

cd (BinRevDir);

datafiles = dir('*.mat');
nrdata = length(datafiles);

for i = 1:nrdata
    cd (BinRevDir);
    S = load(datafiles(i).name);
    if isfield(S.ParamsOut,'Revcor')
%         if ~isfield(S.ParamsOut.Revcor,'BPhaseFits')
            ParamsOut = S.ParamsOut;clear S;
            weights = zeros(size(ParamsOut.Revcor.h1zscore));
            for j = 1:2
                ind(j) = find(ParamsOut.Revcor.h1ffax==ParamsOut.Revcor.domfreq(j));
                startind(j) = find(ParamsOut.Revcor.h1zscore(1:ind(j)-1,j)>=5,1,'first')+1;
                stopind(j) = find(ParamsOut.Revcor.h1zscore(ind(j)+1:end,j)>=5,1,'last')+ind(j)-1;
                weights(startind(j):stopind(j),j) = 10.^(ParamsOut.Revcor.h1mag(startind(j):stopind(j),j)/20);
                weights(:,j) = weights(:,j)/max(weights(:,j));
                weights(find(ParamsOut.Revcor.h1zscore(startind(j):stopind(j),j)<5)+startind(j)-1,j)=0;
%                 weights(ParamsOut.Revcor.h1zscore(:,j)<5,j)=0;
            end
            
            startat = max(startind);
            stopat = min(stopind);
            Phase = ParamsOut.Revcor.h1phase(startat:stopat,:);
            Freq = ParamsOut.Revcor.h1ffax(startat:stopat);
            BPhase_cycles = (Phase(:,2)-Phase(:,1))./(2*pi);
            Bweights = prod(weights(startat:stopat,:),2);
            Bweights = Bweights/max(Bweights);
            [dummy,indx]=max(Bweights);
            ind1 = find(Bweights(1:indx-1)==0,1,'last')+1;
            ind2 = find(Bweights(indx+1:end)==0,1,'first')+indx-1;
            if isempty(ind1)
                ind1 = 1;
            end
            if isempty(ind2)
                ind2 = length(Bweights);
            end
            Bweights(1:ind1-1)=0;
            Bweights(ind2+1:end)=0;
            cd C:\LeuvenDataAnalysis;
            [CP,CD,R,N] = fit_circ_lin(Freq,BPhase_cycles,Bweights');
            ParamsOut.Revcor.BPhaseFits.bbCP = CP;
            ParamsOut.Revcor.BPhaseFits.bbCD = CD;
            ParamsOut.Revcor.BPhaseFits.Resid = R;
            ParamsOut.Revcor.BPhaseFits.N = N;
            ParamsOut.Revcor.BPhaseFits.Bweights = Bweights;
            ParamsOut.Revcor.BPhaseFits.BPhase_cycles = BPhase_cycles;
            ParamsOut.Revcor.BPhaseFits.Freq = Freq;
            
            cd (BinRevDir);
            save (datafiles(i).name,'ParamsOut');
%         end
    end
end