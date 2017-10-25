function [] = get_AN_CP_CD
clear;
clc;
%load list of AN data
ANdatalist = load('C:\Users\Mark\Dropbox\Disparity_MATLAB\MonRevData.mat');
%load the real MSO tuning difference parameters from file
load('C:\Users\Mark\Dropbox\Disparity_MATLAB\tuning_diff_params.mat');

CF = ANdatalist.CF;
DF = ANdatalist.DF;

inds = find(CF<=1500 & ~isnan(DF));

ANdata.CF = ANdatalist.CF(inds);
ANdata.DF = ANdatalist.DF(inds);
ANdata.unitID = ANdatalist.unitID(inds);

nrdata = length(ANdata.DF);

for i = 1:nrdata
    loF = ANdata.DF(i)*2^-1;
    hiF = ANdata.DF(i)*2^1;
    indx{i} = find(ANdata.DF>=loF & ANdata.DF<=hiF);
    indx{i} = indx{i}(indx{i}~=i);
    deltaDF{i} = log2(ANdata.DF(i)./ANdata.DF(indx{i}));
    DFind{i} = i*ones(1,length(deltaDF{i}));
end

alldeltaDF = [deltaDF{:}];
allindx = [indx{:}];
allDFind = [DFind{:}];

X = -1:0.01:1;
requireddist = normpdf(X,meantuningdiff,stdtuningdiff);
requireddist = requireddist/max(requireddist);

currentdist = hist(alldeltaDF,X);
currentdist  = currentdist/max(currentdist);
currentdist = 1-currentdist;

probdeltaDF = interp1(X,requireddist,alldeltaDF).*interp1(X,currentdist,alldeltaDF);

flag = probdeltaDF>rand(1,length(alldeltaDF));
newdeltaDF = alldeltaDF(flag);
newindx = allindx(flag);
newDFind = allDFind(flag);

clear deltaDF indx DFind;

deltaDF = newdeltaDF;
indx = newindx;
DFind = newDFind;
A = [indx' DFind'];
A = sort(A,2);
[~,vv,~] = unique(A,'rows');
indx = indx(vv);
DFind = DFind(vv);
nrpairs = length(indx);
datapath = 'C:\Users\Mark\Dropbox\MonRev';
savedir = 'C:\Users\Mark\Dropbox\Disparity_MATLAB';
h = waitbar(0,'wait');
CD = nan(1,nrpairs);
CP = nan(1,nrpairs);
TD = nan(1,nrpairs);
difcorDF = nan(1,nrpairs);
for i = 1:nrpairs
    IPSI = load(fullfile(datapath,ANdata.unitID{DFind(i)},[ANdata.unitID{DFind(i)} '_MonRev_70.00dB.mat']));
    CONTRA = load(fullfile(datapath,ANdata.unitID{indx(i)},[ANdata.unitID{indx(i)} '_MonRev_70.00dB.mat']));
    if length(IPSI.kernels.h1c)==length(CONTRA.kernels.h1c)
        kernels.h1c = [IPSI.kernels.h1c CONTRA.kernels.h1c];
    elseif length(IPSI.kernels.h1c)>length(CONTRA.kernels.h1c)
        paddzeros = zeros(length(IPSI.kernels.h1c)-length(CONTRA.kernels.h1c),1);
        kernels.h1c = [IPSI.kernels.h1c [CONTRA.kernels.h1c;paddzeros]];
    elseif length(IPSI.kernels.h1c)<length(CONTRA.kernels.h1c)
        paddzeros = zeros(length(CONTRA.kernels.h1c)-length(IPSI.kernels.h1c),1);
        kernels.h1c = [[IPSI.kernels.h1c;paddzeros] CONTRA.kernels.h1c];
    end
    kernels.h1ffax = IPSI.kernels.h1ffax;
    [kernels] = Phase_Freq_Fit(kernels);
    NDF.difcor.power = 20*log10(kernels.h1coeffs.Binaural.weights);
    [~,loc] = max(NDF.difcor.power);
    NDF.difcor.peakhz = 1000*kernels.h1ffax(loc);
    [kernels] = get_tuning_difference(kernels,NDF);
    ipsiID{i} = ANdata.unitID{DFind(i)};
    contraID{i} = ANdata.unitID{indx(i)};
    if strcmp(ipsiID{i}(1:6),contraID{i}(1:6))
        sameanimal(i) = 1;
    else
        sameanimal(i) = 0;
    end
    CD(i) = kernels.h1coeffs.Binaural.CD;
    CP(i) = kernels.h1coeffs.Binaural.CP;
    TD(i) = kernels.tuning_diff.ipsivscontra_oct;
    ipsiDF(i) = IPSI.kernels.h1domfreq;
    contraDF(i) = CONTRA.kernels.h1domfreq;
    difcorDF(i) = NDF.difcor.peakhz/1000;
    waitbar(i/nrpairs);
end
save(fullfile(savedir,'AN_CP_CD.mat'),'CP','CD','difcorDF','TD','ipsiDF','contraDF','ipsiID','contraID','sameanimal');








































