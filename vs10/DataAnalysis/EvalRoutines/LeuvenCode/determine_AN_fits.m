function [foundmatch,delaysamps,ANfilename] = determine_AN_fits(h1,DFs)
%Gets two binaural revcors (ipsi and contra) as input and tries to find AN
%fibers in the data with matching revcors by minimizing the summed squared
%error between MSO and AN data.

%Get the list of AN dominant frequencies so you can narrow your search
%space: only consider fibers with DFs within 0.25 octaves of the MSO DF
ANdata = load('C:\LeuvenDataAnalysis\MonRevPopData.mat');
ANdfs = zeros(length(ANdata.AN),1);
for i = 1:length(ANdata.AN)
    ind = find(unique(ANdata.AN(i).noiseSPLs)==70,1);
    if ~isempty(ind)
        ANdfs(i) = ANdata.AN(i).df(ind);
    else
        ANdfs(i) = NaN;
    end
end
ANcfs = vertcat(ANdata.AN(:).cf);
ANfreqs = ANdfs;
ANfreqs(ANcfs>2000) = ANcfs(ANcfs>2000)/1000;
for j = 1:2
    dif = abs(log2(ANfreqs./DFs(j)));
    indx{j} = find(dif<=0.25);
    %For each potential match, load the AN data and cross correlate the
    %impulse response from the MSO fiber with that from the AN fiber
    for k = 1:length(indx{j})
        load(fullfile('C:\work\BinRev\MonRev',ANdata.AN(indx{j}(k)).filename));
        if isfield(ParamsOut,'noiseSPLs')
            spls = unique(ParamsOut.noiseSPLs);
            splind = find(spls==70);
        else
            splind = 1;
        end
        if ~isempty(splind)
            ANh1 = ParamsOut.h1(:,splind);
            XC{j,k} = xcorr(ANh1,h1(:,j),'coeff');
            XCM{j,k} = max(XC{j,k});
        else
            XC{j,k} = zeros(length(ParamsOut.h1)*2-1,1);
            XCM{j,k} = 0;
        end
    end
    %Take the one with the best correlation
    [dum,ind]= max([XCM{j,:}]);
    foundmatch(j) = indx{j}(ind);
    [dum,ind] = max(XC{j,ind});
    delaysamps(j) = length(h1(:,j))-ind;
    ANfilename{j} = ANdata.AN(foundmatch(j)).filename;
end
return;