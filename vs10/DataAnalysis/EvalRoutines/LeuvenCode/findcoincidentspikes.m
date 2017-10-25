function [outspikes] = findcoincidentspikes(Aspikes,Bspikes,w)
outspikes = [];
for i = 1:length(Aspikes)
    ciB = Bspikes(Bspikes<=Aspikes(i)+(w/2) & Bspikes>=Aspikes(i)-(w/2));
    outspikes = [outspikes max([repmat(Aspikes(i),1,length(ciB));ciB])];
end


% inds = find(abs(bsxfun(@minus,Aspikes,Bspikes'))<=w);
% [Binds,Ainds] = ind2sub([length(Bspikes),length(Aspikes)],inds);
% outspikes = max([Aspikes(Ainds);Bspikes(Binds)]);

return;