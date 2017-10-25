function [h1bw,qvals,domfreq,domind] = Bandwidth_analysis(magIN,zIN,critIN,fIN)
%Returns the 3 and 6-dB bandwidths of the h1 kernels (in kHz), the dominant
%frequency (in kHz) and the index of that value
bwlevels = 1:10;
Nchan = size(magIN,2);
h1bw = nan(length(bwlevels),Nchan);
qvals = nan(length(bwlevels),Nchan);
domfreq = nan(1,Nchan);
domind = nan(1,Nchan);
for ii=1:Nchan
    [maxval,maxind] = max(magIN(:,ii));
    domfreq(ii) = fIN(maxind);
    domind(ii) = maxind;
    lovals = magIN(1:maxind,ii);
    hivals = magIN(maxind:end,ii);
    flo = fIN(1:maxind);
    fhi = fIN(maxind:end);
    zlo = zIN(1:maxind,ii);
    zhi = zIN(maxind:end,ii);
    for jj = 1:length(bwlevels)
        ind1 = find(lovals<maxval-bwlevels(jj),1,'last');
        ind2 = find(lovals>maxval-bwlevels(jj),1,'first');
        if ~isempty(ind1) && ~isempty(ind2)
            if zlo(ind1)>=critIN && zlo(ind2)>=critIN
                lointercept = interp1(lovals([ind1 ind2]),flo([ind1 ind2]),maxval-bwlevels(jj));
            else
                lointercept = nan;
            end
        else
            lointercept = nan;
        end
        ind1 = find(hivals<maxval-bwlevels(jj),1,'first');
        ind2 = find(hivals>maxval-bwlevels(jj),1,'last');
        if ~isempty(ind1) && ~isempty(ind2)
            if zhi(ind1)>=critIN && zhi(ind2)>=critIN
                hiintercept = interp1(hivals([ind1 ind2]),fhi([ind1 ind2]),maxval-bwlevels(jj));
            else
                hiintercept = nan;
            end
        else
            hiintercept = nan;
        end
        h1bw(jj,ii) = hiintercept-lointercept;
    end
    qvals(:,ii) = domfreq(ii)./h1bw(:,ii);
end
end