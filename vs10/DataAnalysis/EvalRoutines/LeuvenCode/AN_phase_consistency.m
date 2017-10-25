clc;clear;
ANdatapath = 'C:\Users\Mark\Dropbox\MonRev';
ANdatalist = load('C:\Users\Mark\Dropbox\Disparity_MATLAB\MonRevData.mat');
CF = ANdatalist.CF;
DF = ANdatalist.DF;
inds = find(CF<=1600 & ~isnan(DF));
unitIDs = ANdatalist.unitID(inds);
DF = ANdatalist.DF(inds);


targetDF = 0.1*2.^(0:0.2:4);
nrF = length(targetDF);

for i = 2:nrF
    loF = targetDF(i)*2^-0.075;
    hiF = targetDF(i)*2^0.075;
    idx = find(DF>=loF & DF<=hiF);
    nrfibers = length(idx);
    h1mat = zeros(4096,nrfibers);
    for j = 1:nrfibers
        load(fullfile(ANdatapath,unitIDs{idx(j)},[unitIDs{idx(j)} '_MonRev_70.00dB.mat']));
        h1mat(:,j) = kernels.h1phase_uw;
    end
    figure;
    plot(kernels.h1ffax,h1mat);
    set(gca,'xlim',[0 3]);
end