function [updatedBD,NDF] = correctNDF_BD(kernels,NDF)

%Get the two biggest peaks in the NDF spline function
[pks,locs] = findpeaks(NDF.difcor.rate.spline.mean,'npeaks',2,'sortstr','descend');
possBDs = NDF.xvals.ITDxspline(locs);
% ndvar = (NDF.difcor.rate.spline.std(locs)).^2;
% dprime = abs((diff(pks))/sqrt(0.5*(sum(ndvar))));

% if dprime<5
oldBD_ndf = NDF.bd.bd;
% [~,loc] = max(NDF.fit.simple.difcorrate);
kernelBD = kernels.h1coeffs.Binaural.CD+(kernels.h1coeffs.Binaural.CP/(NDF.difcor.peakhz/1000));
oldBDind = find(possBDs==oldBD_ndf);
if oldBDind==1
    if abs(diff([kernelBD,possBDs(2)]))<abs(diff([kernelBD,possBDs(1)]))
        newBD_ndf = possBDs(2);
        updatedBD = 1;
    else
        updatedBD = 0;
    end
elseif oldBDind==2
    if abs(diff([kernelBD,possBDs(1)]))<abs(diff([kernelBD,possBDs(2)]))
        newBD_ndf = possBDs(1);
        updatedBD = 1;
    else
        updatedBD = 0;
    end
end
% else
%     updatedBD = 0;
% end

if updatedBD
    %update the BD
    NDF.bd.bd = newBD_ndf;
    %Update the side peaks
    vals = NDF.difcor.rate.spline.mean-min(NDF.difcor.rate.spline.mean);
    [peakval,peakind] = findpeaks(vals);
    [troughval,troughind] = findpeaks(-vals);
    troughval = -troughval;
    if min(peakind)<min(troughind)
        troughind = [1;troughind];
        troughval = [vals(1);troughval];
    end
    if max(peakind)>max(troughind)
        troughind = [troughind;length(NDF.difcor.rate.spline.mean)];
        troughval = [troughval;vals(length(NDF.difcor.rate.spline.mean))];
    end
    Dpeaks = NDF.xvals.ITDxspline(peakind);
    Dtroughs = NDF.xvals.ITDxspline(troughind);
    keeppeak = zeros(size(peakval));
    keeptrough = zeros(size(troughval));
    for i = 1:length(peakval)
        clo = (peakval(i)-troughval(find(troughind<peakind(i),1,'last')))/...
            (peakval(i)+troughval(find(troughind<peakind(i),1,'last')));
        chi = (peakval(i)-troughval(find(troughind>peakind(i),1,'first')))/...
            (peakval(i)+troughval(find(troughind>peakind(i),1,'first')));
        if clo>0.5
            keeppeak(i) = 1;
            keeptrough(find(troughind<peakind(i),1,'last'))=1;
        elseif chi>0.5
            keeppeak(i) = 1;
            keeptrough(find(troughind>peakind(i),1,'first'))=1;
        end
    end
    Dpeaks = Dpeaks(logical(keeppeak));
    Dtroughs = Dtroughs(logical(keeptrough));
    Dpeaks = Dpeaks(Dpeaks~=newBD_ndf);
    NDF.bd.dpeaks = Dpeaks;
    NDF.bd.dtroughs = Dtroughs;
end