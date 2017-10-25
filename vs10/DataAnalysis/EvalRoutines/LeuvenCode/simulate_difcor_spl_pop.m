function [] = simulate_difcor_spl_pop

%Get the data list;
load('C:\LeuvenDataAnalysis\MonRevPopData.mat');
indx = [];
count = 1;
for i = 1:length(AN)
    if length(unique(AN(i).noiseSPLs(~isnan(AN(i).noiseSPLs))))>1
        if mean(AN(i).df)<3 && AN(i).cf<=3000
            indx = [indx i];
            mdf(count) = mean(AN(i).df);
            count = count+1;
        end
    end
end
iindx = 1:length(indx);
allpairs = nchoosek(iindx,2);
alldf = mdf(allpairs);
alldeltaf = log2(alldf(:,2)./alldf(:,1));
newpairs = allpairs(alldeltaf>0 & alldeltaf<0.5,:);
fiberpairs = indx(newpairs);

for i = 1:length(fiberpairs)
    SPLs1 = AN(fiberpairs(i,1)).noiseSPLs;
    SPLs1 = unique(SPLs1(~isnan(SPLs1)));
    SPLs2 = AN(fiberpairs(i,2)).noiseSPLs;
    SPLs2 = unique(SPLs2(~isnan(SPLs2)));
    nrcommonspls = length(intersect(SPLs1,SPLs2));
    if nrcommonspls>1
        keepit(i) = 1;
    end
end
fiberpairs = fiberpairs(logical(keepit),:);

%Now simulate the bincor for each of these pairs of fibers
h = waitbar(0,'Please wait...');
for i = 1:length(fiberpairs)
    for j = 1:2
        ind = regexp(AN(fiberpairs(i,j)).filename,'_');
        unitnum{j} = AN(fiberpairs(i,j)).filename(ind(1)+1:ind(2)-1);
    end
    [CD{i},CP{i},DF{i},BD{i},SPL{i}] = simulate_difcor_spl({AN(fiberpairs(i,1)).animal,AN(fiberpairs(i,2)).animal},{unitnum{1},unitnum{2}},0);
    waitbar(i/length(fiberpairs));
    BDDF{i} = BD{i}.*DF{i};
    if strcmp(AN(fiberpairs(i,1)).animal,AN(fiberpairs(i,2)).animal)
        comptype(i) = 1;
    else
        comptype(i) = 2;
    end
    if any(SPL{i}==70)
        BDDFnorm{i} = BDDF{i}-BDDF{i}(SPL{i}==70);
    else
        BDDFnorm{i} = nan(1,length(BDDF{i}));
    end
    mDF(i) = mean(DF{i});
end
close (h);

save(fullfile('C:\work\BinRev\Analyses','AN_simulate_difcor_spl_pop'),'CD','CP','DF','mDF','SPL','BDDF','BDDFnorm','comptype');

figure;
colormap jet;
for i = 1:length(SPL)
    if mDF(i)<=1.5;
        x = SPL{i};if size(x,1)>size(x,2);x = x';end;
        y = BDDFnorm{i};if size(y,1)>size(y,2);y = y';end;
        z = ones(1,length(SPL{i}))*mDF(i);
        col = z;
        surface([x;x],[y;y],[z;z],[col;col],...
            'facecol','no',...
            'edgecol','flat',...
            'linew',2);hold on;
%         scatter3(x,y,z,100,col,'filled');
    end
end
set(gca,'fontsize',16,'linewidth',1,'layer','top');
xlabel 'SOUND LEVEL (dB SPL)';
ylabel ('\Delta(BD*DF) (cycles)','interpreter','tex');

figure;
colormap jet;
for i = 1:length(SPL)
    if mDF(i)<=1.5;
        x = SPL{i};if size(x,1)>size(x,2);x = x';end;
        y = BDDF{i};if size(y,1)>size(y,2);y = y';end;
        z = ones(1,length(SPL{i}))*mDF(i);
        col = z;
        surface([x;x],[y;y],[z;z],[col;col],...
            'facecol','no',...
            'edgecol','flat',...
            'linew',2);hold on;
%         scatter3(x,y,z,100,col,'filled');
    end
end
set(gca,'fontsize',16,'linewidth',1,'layer','top');
xlabel 'SOUND LEVEL (dB SPL)';
ylabel ('BD*DF (cycles)','interpreter','tex');

figure;
colormap jet;
for i = 1:length(SPL)
    if mDF(i)<=1.5;
        x = SPL{i};if size(x,1)>size(x,2);x = x';end;
        y = BD{i};if size(y,1)>size(y,2);y = y';end;
        z = ones(1,length(SPL{i}))*mDF(i);
        col = z;
        surface([x;x],[y;y],[z;z],[col;col],...
            'facecol','no',...
            'edgecol','flat',...
            'linew',2);hold on;
%         scatter3(x,y,z,100,col,'filled');
    end
end
set(gca,'fontsize',16,'linewidth',1,'layer','top');
xlabel 'SOUND LEVEL (dB SPL)';
ylabel ('BD (ms)','interpreter','tex');




