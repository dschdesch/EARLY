%Build monaural wiener kernel population database

datapath = 'C:\Users\Mark\Dropbox\MonRev';
cd (datapath);
datafolders = dir;
datafolders = datafolders(4:end-1);
nrunits = length(datafolders);
unitID = cell(1,nrunits);
RC = cell(1,nrunits);
CF = nan(1,nrunits);Th = nan(1,nrunits); 
Q10 = nan(1,nrunits); SR = nan(1,nrunits); 
DF = nan(1,nrunits); GD = nan(1,nrunits);
for i  = 1:nrunits
    TC=[];
    kernels = [];
    unitID{i} = datafolders(i).name;
    cd(unitID{i});
    
    %get the TC
    load([unitID{i} '_TC.mat']);
    CF(i) = TC.fit.cf;
    Th(i) = TC.fit.minthr;
    Q10(i) = TC.fit.q10;
    if isnan(Q10(i))
        Q10(i) = TC.thr.q10;
    end
    SR(i) = TC.thr.sr;
    TCFreqs{i} = TC.fit.freq;
    TCThresh{i} = TC.fit.thr;
    
    %get the 70-dB revcor info if it exists
    if ~isempty(dir('*MonRev_70.00dB.mat'))
        load([unitID{i} '_MonRev_70.00dB.mat']);
        DF(i) = kernels.h1cdomfreq;
        GD(i) = kernels.h1coeffs.Monaural.Delay;
    end
    
    %get the SPL function is it exists
    if ~isempty(dir('*RLV.mat'))
        load([unitID{i} '_RLV.mat']);
        [RC{i}.SPL,ix] = sort(RLV.Level);
        RC{i}.Rate = RLV.Rate(ix);
        RC{i}.SyncMag = RLV.SyncMag(ix);
        RC{i}.SyncPhase = RLV.SyncPhase(ix);
        RC{i}.SyncPval = RLV.SyncPval(ix);
        RC{i}.Freq = RLV.Freq;
    end
    
    cd(datapath);
end

%linear regression of GD on log2(CF)
X = CF(CF<=1300);
Y = GD(CF<=1300);
X = X(~isnan(Y))/1000;
Y = Y(~isnan(Y));
[fo,gof] = fit(log2(X'),Y','poly1');

%Save the data for use in MSO analyses
savedir = 'C:\Users\Mark\Dropbox\Disparity_MATLAB';
save(fullfile(savedir,'MonRevData.mat'),'unitID','CF','SR','Th','Q10','DF','GD','RC','fo','gof','TCFreqs','TCThresh');