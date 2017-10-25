function [] = MonRev_plotpop

currdir = cd;
datadir = 'C:\Users\Mark\Dropbox\MonRev';
cd (datadir);
datafiles = dir('*_MonRev.mat');
nrdata = numel(datafiles);
count = 0;
rscount = 0;
h = waitbar(0,'Please wait...');
for i =1:nrdata
    data = load (datafiles(i).name);
    ParamsOut = data.ParamsOut;clear data;
    annum(i,:) = str2num(datafiles(i).name(2:6));
    if isfield(ParamsOut,'h1')
        h1indx = find(unique(ParamsOut.noiseSPLs)==70);
        nrspls(i) = length(unique(ParamsOut.noiseSPLs));
        if ~isempty(h1indx)
            for j  = 1:numel(ParamsOut.StimulusInfo)
                if isfield(ParamsOut.StimulusInfo{j},'StimParam')
                    rscount = rscount+1;
                    RSeed(rscount) = ParamsOut.StimulusInfo{j}.StimParam.RandomSeed;
                end
            end
            Time = ParamsOut.Time;
            FreqAx = ParamsOut.h1ffax;
            try
            h1mat(:,i) = ParamsOut.h1filt(:,h1indx)./max(abs(ParamsOut.h1filt(:,h1indx)));
            catch
                foo
            end
            h1magmat(:,i) = ParamsOut.h1mag(:,h1indx);
            h1phmat(:,i) = ParamsOut.h1phase(:,h1indx);
            h1latency(i) = ParamsOut.h1latency(h1indx);
            h1z(:,i) = ParamsOut.h1zscore(:,h1indx);
            df(:,i) = ones(length(ParamsOut.h1filt(:,h1indx)),1).*ParamsOut.domfreq(h1indx);
            if ~isempty(ParamsOut.h1phasefit_coeffs)
                try
                    gd(i) = ParamsOut.h1phasefit_coeffs(h1indx,1)*-1;
                catch
                    gd(i) = nan;
                end
            else
                gd(i) = nan;
            end
            frontlatency(i) = ParamsOut.h1frontlatency(h1indx);
            bw_3dB(i) = ParamsOut.h1bw(1,h1indx);
            bw_6dB(i) = ParamsOut.h1bw(2,h1indx);
        else
            df(:,i) = nan(length(ParamsOut.h1filt(:,1)),1);
            frontlatency(i) = nan;
            bw_3dB(i) = nan;
            bw_6dB(i) = nan;
            gd(i) = nan;
        end
    else
        h1mat(:,i) = nan(length(h1mat),1);
        h1phmat(:,i) = nan(length(h1phmat),1);
        h1magmat(:,i) = nan(length(h1magmat),1);
        h1latency(i) = nan;
        h1z(:,i) = nan(length(h1z),1);
        df(:,i) = nan(length(df),1);
        gd(i) = nan;
        frontlatency(i) = nan;
        bw_3dB(i) = nan;
        bw_6dB(i) = nan;
        nrspls(i) = 0;
    end
    cf(:,i) = ones(length(df),1).*ParamsOut.TC.fit.cf;
    cfthr(i) = ParamsOut.TC.fit.minthr;
    sr(i) = ParamsOut.TC.thr.sr;
    TC{i}.freq = ParamsOut.TC.fit.freq;
    TC{i}.thr = ParamsOut.TC.fit.thr;
    q10fit(i) = ParamsOut.TC.fit.q10;
    q10thr(i) = ParamsOut.TC.thr.q10;
    if ~isempty(ParamsOut.RLV)
        nrrlvspls = numel(ParamsOut.RLV.Level);
        for j = 1:nrrlvspls
            ss = [ParamsOut.RLV.SpikeTimes{j,:}];
            nrspikes(j) = length(ss(ss<=25));
        end
        if ~isempty(max(ParamsOut.RLV.SyncMag(ParamsOut.RLV.SyncPval<0.01 & nrspikes>100)))
            sync(i,1) = max(ParamsOut.RLV.SyncMag(ParamsOut.RLV.SyncPval<0.01 & nrspikes>100));
            sync(i,2) = 1;
        elseif ~isempty(max(ParamsOut.RLV.SyncMag(nrspikes>100)))
            sync(i,1) = max(ParamsOut.RLV.SyncMag(nrspikes>100));
            sync(i,2) = 2;
        else
            sync(i,[1 2]) = [NaN NaN];
        end
        clear nrspikes;
    else
        sync(i,[1 2]) = [NaN NaN];
    end
    AN(i).filename = datafiles(i).name;
    AN(i).DatasetIDs = ParamsOut.DatasetIDs;
    if isfield(ParamsOut,'h1')
        AN(i).df = ParamsOut.domfreq;
        if isfield(ParamsOut,'noiseSPLs')
            AN(i).noiseSPLs = ParamsOut.noiseSPLs;
        else
            AN(i).noiseSPLs = 70;
        end
        for j = 1:length(ParamsOut.StimulusInfo)
            AN(i).noiseseed(j) = ParamsOut.StimulusInfo{j}.StimParam.RandomSeed;
        end
        AN(i).animal = ParamsOut.DatasetIDs.Animal;
    else
        AN(i).df = NaN;
    end
    AN(i).cf = ParamsOut.TC.thr.cf;
    waitbar(i/nrdata);
end
close (h);

cd (currdir);
%% Tidy up the CF/DF axis - to take account of "DF" being much lower than CF for fibers around 3-4 kHz.
df = df*1000;
FreqX = df;
FreqX(cf>2000) = cf(cf>2000);

%% 
fffind = find(~isnan(FreqX(1,:)) & FreqX(1,:)<=2000);
fff = FreqX(1,fffind);
[fff,sortind] = sort(fff);
fffind =fffind(sortind);
zzz = h1mat(:,fffind);
xxx = repmat(fff,1803,1);
yyy = repmat(Time',1,226);

%% Plot CF vs DF
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.4056    0.2467    0.3181    0.4244]);
DFCFax = axes('ActivePositionProperty','outerposition');
DFCFs = cf(1,cf(1,:)<=2500 & cf(1,:)>=50);
DFDFs = df(1,cf(1,:)<=2500 & cf(1,:)>=50);
mycmap = [[zeros(32,1);(linspace(0,1,16))';ones(16,1)] [(linspace(0,1,16))';ones(32,1);(linspace(1,0,16))'] [ones(16,1); linspace(1,0,16)'; zeros(32,1)]];
colormap(mycmap);
colormapyvals = mycmap;
minCF = 50;
maxCF = 2500;

xdata = DFCFs;
ydata = DFDFs;
xdata = xdata(~isnan(ydata));
ydata = ydata(~isnan(ydata));
rho = corr(xdata',ydata');
[beta,gof] = local_fit_linear(xdata'/1000,ydata'/1000);
[coeffpval,coefftval] = local_stats(xdata'/1000,ydata'/1000,beta,gof.dfe);

colormapxvals = logspace(log10(minCF),log10(maxCF),64);
plot([0.05 3000],[0.05 3000],'color',[.6 .6 .6],'linewidth',1);hold on;
for i = 1:length(DFCFs)
    currentrgb = [interp1(colormapxvals,colormapyvals(:,1),DFCFs(i)) interp1(colormapxvals,colormapyvals(:,2),DFCFs(i)) interp1(colormapxvals,colormapyvals(:,3),DFCFs(i))];
    plot(DFCFax,DFCFs(i)/1000,DFDFs(i)/1000,'color',currentrgb,'markersize',10,'marker','*');hold on;
end
box off;
plot([0.05 3],feval(beta,[0.05 3]),'k--','linewidth',1);
set(gca,'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'xtick',[0.1 1],'xticklabel',[0.1 1],'yticklabel',[0.1 1],...
    'xscale','log','XLIM',[0.05 3],'yLIM',[0.05 3],'yscale','log');
xlabel 'CHARACTERISTIC FREQUENCY (kHz)';
ylabel ('DOMINANT FREQUENCY (kHz)');
axis square;
set(gca,'position',[0.1822    0.1100    0.5608    0.8150]);
text(0.1,2,sprintf('\\rho: %2.2f \n\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',rho,beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex','color','k');



%% Plot Group Delay vs. CF
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
GDCFax = axes('ActivePositionProperty','outerposition');
latCFs = cf(1,cf(1,:)<=2500 & cf(1,:)>=50);
latLats = gd(cf(1,:)<=2500 & cf(1,:)>=50);
mycmap = [[zeros(32,1);(linspace(0,1,16))';ones(16,1)] [(linspace(0,1,16))';ones(32,1);(linspace(1,0,16))'] [ones(16,1); linspace(1,0,16)'; zeros(32,1)]];
colormap(mycmap);
colormapyvals = mycmap;
minCF = 50;
maxCF = 2500;

colormapxvals = logspace(log10(minCF),log10(maxCF),64);
for i = 1:length(latCFs)
    currentrgb = [interp1(colormapxvals,colormapyvals(:,1),latCFs(i)) interp1(colormapxvals,colormapyvals(:,2),latCFs(i)) interp1(colormapxvals,colormapyvals(:,3),latCFs(i))];
    plot(GDCFax,latCFs(i),latLats(i),'color',currentrgb,'markersize',10,'marker','*');hold on;
end
%Fit the data
[beta,gof]=local_fit_linear(log10(latCFs(~isnan(latLats))'),latLats(~isnan(latLats))');
[coeffpval,coefftval] = local_stats(log10(latCFs(~isnan(latLats))'),latLats(~isnan(latLats))',beta,gof.dfe);
%add the regression line
plot(80:3000,feval(beta,log10(80:3000)),'k--','linewidth',2);
%add some text describing the regression coefficient and stats
text(300,7,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex')

box off;
set(gca,'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'xtick',[100 1000],'xticklabel',[0.1 1],'xscale','log','XLIM',[50 5000],'ylim',[0 8]);
xlabel 'CF (kHz)';
ylabel ('GROUP DELAY (ms)');
axis square;
set(gca,'position',[0.1822    0.1100    0.5608    0.8150]);

%% Plot Signal-Front Delay vs. CF
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
SFDCFax = axes('ActivePositionProperty','outerposition');
SFDCFs = cf(1,cf(1,:)<=2500 & cf(1,:)>=50);
SFDs = frontlatency(cf(1,:)<=2500 & cf(1,:)>=50);
mycmap = [[zeros(32,1);(linspace(0,1,16))';ones(16,1)] [(linspace(0,1,16))';ones(32,1);(linspace(1,0,16))'] [ones(16,1); linspace(1,0,16)'; zeros(32,1)]];
colormap(mycmap);
colormapyvals = mycmap;
minCF = 50;
maxCF = 2500;
colormapxvals = logspace(log10(minCF),log10(maxCF),64);
for i = 1:length(SFDs)
    currentrgb = [interp1(colormapxvals,colormapyvals(:,1),SFDCFs(i)) interp1(colormapxvals,colormapyvals(:,2),SFDCFs(i)) interp1(colormapxvals,colormapyvals(:,3),SFDCFs(i))];
    plot(SFDCFax,SFDCFs(i),SFDs(i),'color',currentrgb,'markersize',10,'marker','*');hold on;
end
%Fit the data
[beta,gof]=local_fit_linear(log10(SFDCFs(~isnan(SFDs))'),SFDs(~isnan(SFDs))');
[coeffpval,coefftval] = local_stats(log10(SFDCFs(~isnan(SFDs))'),SFDs(~isnan(SFDs))',beta,gof.dfe);
%add the regression line
plot(80:3000,feval(beta,log10(80:3000)),'k--','linewidth',2);
%add some text describing the regression coefficient and stats
text(300,7,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex')
box off;
set(gca,'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'ylim',[0 8],'xtick',[100 1000],'xticklabel',[0.1 1],'xscale','log','XLIM',[50 5000]);
xlabel 'CF (kHz)';
ylabel ('SIGNAL-FRONT DELAY (ms)');
axis square;
set(gca,'position',[0.1822    0.1100    0.5608    0.8150]);

%% Plot Filter Delay (i.e., Group delay - signal-front delay) vs. CF
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
FDCFax = axes('ActivePositionProperty','outerposition');
FDCFs = cf(1,cf(1,:)<=2500 & cf(1,:)>=50);
FDs = gd(cf(1,:)<=2500 & cf(1,:)>=50) -  frontlatency(cf(1,:)<=2500 & cf(1,:)>=50);
mycmap = [[zeros(32,1);(linspace(0,1,16))';ones(16,1)] [(linspace(0,1,16))';ones(32,1);(linspace(1,0,16))'] [ones(16,1); linspace(1,0,16)'; zeros(32,1)]];
colormap(mycmap);
colormapyvals = mycmap;
minCF = 50;
maxCF = 2500;

colormapxvals = logspace(log10(minCF),log10(maxCF),64);
for i = 1:length(FDs)
    currentrgb = [interp1(colormapxvals,colormapyvals(:,1),FDCFs(i)) interp1(colormapxvals,colormapyvals(:,2),FDCFs(i)) interp1(colormapxvals,colormapyvals(:,3),FDCFs(i))];
    plot(FDCFax,FDCFs(i),FDs(i),'color',currentrgb,'markersize',10,'marker','*');hold on;
end
%Fit the data
[beta,gof]=local_fit_linear(log10(FDCFs(~isnan(FDs))'),FDs(~isnan(FDs))');
[coeffpval,coefftval] = local_stats(log10(FDCFs(~isnan(FDs))'),FDs(~isnan(FDs))',beta,gof.dfe);

%add the regression line
plot(80:3000,feval(beta,log10(80:3000)),'k--','linewidth',2);
%add some text describing the regression coefficient and stats
text(300,7,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex')
box off;
set(gca,'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'ylim',[0 8],'xtick',[100 1000],'xticklabel',[0.1 1],'xscale','log','XLIM',[50 5000]);
xlabel 'CF (kHz)';
ylabel ('FILTER DELAY (ms)');
axis square;
set(gca,'position',[0.1822    0.1100    0.5608    0.8150]);

%% Plot Group Delay vs. Signal-Front Delay
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
GDFDax = axes('ActivePositionProperty','outerposition');
GDs = gd(cf(1,:)<=2000 & cf(1,:)>=50);
GDCFs = cf(1,cf(1,:)<=2000 & cf(1,:)>=50);
SFDs = frontlatency(cf(1,:)<=2000 & cf(1,:)>=50);
mycmap = [[zeros(32,1);(linspace(0,1,16))';ones(16,1)] [(linspace(0,1,16))';ones(32,1);(linspace(1,0,16))'] [ones(16,1); linspace(1,0,16)'; zeros(32,1)]];
colormap(mycmap);
colormapyvals = mycmap;
minCF = 50;
maxCF = 2000;

colormapxvals = logspace(log10(minCF),log10(maxCF),64);
for i = 1:length(GDs)
    currentrgb = [interp1(colormapxvals,colormapyvals(:,1),GDCFs(i)) interp1(colormapxvals,colormapyvals(:,2),GDCFs(i)) interp1(colormapxvals,colormapyvals(:,3),GDCFs(i))];
    plot(GDFDax,SFDs(i),GDs(i),'color',currentrgb,'markersize',10,'marker','*');hold on;
end
box off;
set(gca,'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'ylim',[0 8],'xlim',[0 8]);
xlabel 'Signal-Front Delay (ms)';
ylabel ('Group Delay (ms)');
axis square;
set(gca,'position',[0.1822    0.1100    0.5608    0.8150]);

%relationship looks almost like FD  = 0.5GD... plot some lines
plot([0 8],[0 8],'--','color',[.6 .6 .6],'linewidth',1);
plot([0 4],[0 8],':','color',[.6 .6 .6],'linewidth',1);

%% Plot Filter Delay vs. 3-DB bandwidth
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
GDFDax = axes('ActivePositionProperty','outerposition');
GDCFs = cf(1,cf(1,:)<=2500 & cf(1,:)>=50);
BW_3dB = bw_3dB(cf(1,:)<=2500 & cf(1,:)>=50);
mycmap = [[zeros(32,1);(linspace(0,1,16))';ones(16,1)] [(linspace(0,1,16))';ones(32,1);(linspace(1,0,16))'] [ones(16,1); linspace(1,0,16)'; zeros(32,1)]];
colormap(mycmap);
colormapyvals = mycmap;
minCF = 50;
maxCF = 2500;

colormapxvals = logspace(log10(minCF),log10(maxCF),64);
for i = 1:length(FDs)
    currentrgb = [interp1(colormapxvals,colormapyvals(:,1),GDCFs(i)) interp1(colormapxvals,colormapyvals(:,2),GDCFs(i)) interp1(colormapxvals,colormapyvals(:,3),GDCFs(i))];
    plot(GDFDax,BW_3dB(i),FDs(i),'color',currentrgb,'markersize',10,'marker','*');hold on;
end

%Fit the data
[beta,gof]=local_fit_linear(log10(BW_3dB(~isnan(FDs) & ~isnan(BW_3dB))'),FDs(~isnan(FDs) & ~isnan(BW_3dB))');
[coeffpval,coefftval] = local_stats(log10(BW_3dB(~isnan(FDs) & ~isnan(BW_3dB))'),FDs(~isnan(FDs) & ~isnan(BW_3dB))',beta,gof.dfe);

%add the regression line
plot(logspace(log10(0.08),log10(1),10),feval(beta,log10(logspace(log10(0.08),log10(1),10))),'k--','linewidth',2);
%add some text describing the regression coefficient and stats
text(0.1,7,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex')
box off;
set(gca,'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'ylim',[0 8],'xscale','log');
xlabel '3-dB BANDWIDTH (kHz)';
ylabel ('FILTER DELAY (ms)');
axis square;
set(gca,'position',[0.1822    0.1100    0.5608    0.8150]);
set(gca,'xlim',[0.05 1.5],'xtick',[0.1 1],'xticklabel',[0.1 1]);

%% Plot Filter Delay vs. 6-DB bandwidth
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
GDFDax = axes('ActivePositionProperty','outerposition');
GDCFs = cf(1,cf(1,:)<=2500 & cf(1,:)>=50);
BW_6dB = bw_6dB(cf(1,:)<=2500 & cf(1,:)>=50);
mycmap = [[zeros(32,1);(linspace(0,1,16))';ones(16,1)] [(linspace(0,1,16))';ones(32,1);(linspace(1,0,16))'] [ones(16,1); linspace(1,0,16)'; zeros(32,1)]];
colormap(mycmap);
colormapyvals = mycmap;
minCF = 50;
maxCF = 2500;

colormapxvals = logspace(log10(minCF),log10(maxCF),64);
for i = 1:length(FDs)
    currentrgb = [interp1(colormapxvals,colormapyvals(:,1),GDCFs(i)) interp1(colormapxvals,colormapyvals(:,2),GDCFs(i)) interp1(colormapxvals,colormapyvals(:,3),GDCFs(i))];
    plot(GDFDax,BW_6dB(i),FDs(i),'color',currentrgb,'markersize',10,'marker','*');hold on;
end

%Fit the data
[beta,gof]=local_fit_linear(log10(BW_6dB(~isnan(FDs) & ~isnan(BW_6dB))'),FDs(~isnan(FDs) & ~isnan(BW_6dB))');
[coeffpval,coefftval] = local_stats(log10(BW_6dB(~isnan(FDs) & ~isnan(BW_6dB))'),FDs(~isnan(FDs) & ~isnan(BW_6dB))',beta,gof.dfe);

%add the regression line
plot(logspace(log10(0.1),log10(1),10),feval(beta,log10(logspace(log10(0.1),log10(1),10))),'k--','linewidth',2);
%add some text describing the regression coefficient and stats
text(0.1,7,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex')


box off;
set(gca,'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'ylim',[0 8],'xscale','log');
xlabel '6-dB BANDWIDTH (kHz)';
ylabel ('FILTER DELAY (ms)');
axis square;
set(gca,'position',[0.1822    0.1100    0.5608    0.8150]);
set(gca,'xlim',[0.05 1.5],'xtick',[0.1 1],'xticklabel',[0.1 1]);


%% Plot a bunch of tuning curves
CFs_to_find = logspace(2,4,10); %10 CFs spaced logarithmically between 0.1 and 10 kHz
for i=1:length(CFs_to_find)
    %Find a fiber matching this CF
    dif = abs(FreqX(1,:)-CFs_to_find(i));
    ind = find(dif==min(dif));
    if CFs_to_find(i)>3000
        if length(ind)>1
            for j = 1:length(ind)
                if cfthr(ind(j))>-10
                    ind = ind(j);
                    break
                end
            end
        end
    else
        if length(ind)>1
            ind = ind(end);
        end
    end
    %what is the exact CF for this fiber
    foundCF(i) = FreqX(1,ind);
    foundTC{i}.freq = TC{ind}.freq;
    foundTC{i}.thr = TC{ind}.thr;
end
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2044    0.4544    0.4719    0.3000]);
TCax = axes('ActivePositionProperty','outerposition');
mycmap = [[zeros(32,1);(linspace(0,1,16))';ones(16,1)] [(linspace(0,1,16))';ones(32,1);(linspace(1,0,16))'] [ones(16,1); linspace(1,0,16)'; zeros(32,1)]];
colormap(mycmap);
colormapyvals = mycmap;
colormapxvals = linspace(1,length(CFs_to_find),64);
for i = 1:length(CFs_to_find)
    currentrgb = [interp1(colormapxvals,colormapyvals(:,1),i) interp1(colormapxvals,colormapyvals(:,2),i) interp1(colormapxvals,colormapyvals(:,3),i)];
    plot(TCax,foundTC{i}.freq,foundTC{i}.thr,'color',currentrgb,'linewidth',3);
    hold on;
end
view(0,90);
set(gca,'xscale','log','xlim',[50 20e3],'fontsize',16,'xtick',[10e1 10e2 10e3],'xticklabel',[0.1 1 10],'ylim',[-20 68],'ticklength',[0.013 0.013]);
xlabel 'FREQUENCY (kHz)';
ylabel 'THRESHOLD (dB SPL)';
box off;

%% Plot all CF threshold data
figure;set(gcf,'paperpositionmode','auto','units','inches','position',[3 2 6.5 6.5]);
set(gcf,'units','normalized');
axes('position',[0.15 0.4 0.5 0.5])
lsrh = semilogx(cf(1,(sr<18)),cfthr(sr<18),'ro','markersize',10,'linewidth',1);%low and medium SR
hold on;
hsrh = semilogx(cf(1,(sr>=18)),cfthr(sr>=18),'k+','markersize',10,'linewidth',1);%high SR
set(gca,'fontsize',16,'xlim',[80 20e3],'xtick',[100 1000 10000],'xticklabel',[0.1 1 10],'ticklength',[0.02 0.02]);
xlabel 'CF (kHz)';
ylabel 'THRESHOLD (dB SPL)';
box off;axis normal;
legend([lsrh hsrh],'L/MSR','HSR');

%Add some histograms to project the data along the two axes

cfx = 50*2.^([0:0.2:9]);
lsrcfy = hist(cf(1,sr<18),cfx);
hsrcfy = hist(cf(1,sr>=18),cfx);
axes('position',[0.15 0.15 0.5 0.15]);
cfbarh = bar(log10(cfx),[lsrcfy; hsrcfy]',1,'stacked');
set(cfbarh(1),'facecolor','r','edgecolor','none');
set(cfbarh(2),'facecolor','k','edgecolor','none');
set(gca,'xtick',[2 3 4],'xticklabel',10.^([2 3 4])./1000,'xlim',log10([80 20e3]),'ylim',[0 60],'fontsize',16);
box off;
ylabel '# FIBERS';
xlabel 'CF (kHz)';

thx = -20:2:80;
axes('position',[0.75 0.4 0.15 0.5]);
lsrthy = hist(cfthr(sr<18),thx);
hsrthy = hist(cfthr(sr>=18),thx);
thbarh = barh(thx,[lsrthy; hsrthy]',1,'stacked');
set(thbarh(1),'facecolor','r','edgecolor','none');
set(thbarh(2),'facecolor','k','edgecolor','none');
set(gca,'xlim',[0 80],'ylim',[-20 80],'fontsize',16,'xdir','reverse','yaxislocation','right');
box off;
xlabel '# FIBERS';
ylabel 'THRESHOLD (dB SPL)';



%And a lowess fit to the data
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);

XdataLSR = cf(1,(sr<18));
YdataLSR = cfthr(sr<18);
[XdataLSR,ind] = sort(XdataLSR);
YdataLSR = YdataLSR(ind);
[YfitLSR,YfitSTDLSR] = mylowessbootstrap(XdataLSR',YdataLSR',0.2,'lowess',200);

XdataHSR = cf(1,(sr>=18));
YdataHSR = cfthr(sr>=18);
[XdataHSR,ind] = sort(XdataHSR);
YdataHSR = YdataHSR(ind);
[YfitHSR,YfitSTDHSR] = mylowessbootstrap(XdataHSR',YdataHSR',0.2,'lowess',200);

lsrh = semilogx(XdataLSR,YfitLSR,'r-','linewidth',4);hold on;
semilogx(XdataLSR,YfitLSR-YfitSTDLSR,'r:',XdataLSR,YfitLSR+YfitSTDLSR,'r:','linewidth',1);hold on;

hsrh = semilogx(XdataHSR,YfitHSR,'k-','linewidth',4);hold on;
semilogx(XdataHSR,YfitHSR-YfitSTDHSR,'k:',XdataHSR,YfitHSR+YfitSTDHSR,'k:','linewidth',1);hold on;

set(gca,'fontsize',16,'xlim',[80 20e3],'xtick',[100 1000 10000],'xticklabel',[0.1 1 10],'ticklength',[0.02 0.02]);
xlabel 'CF (kHz)';
ylabel 'THRESHOLD (dB SPL)';
box off;axis square;
legend([lsrh hsrh],'L/MSR','HSR');
%% Plot all Q10 data from the fit
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
semilogx(cf(1,(sr<18)),q10fit(sr<18),'ro','markersize',10,'linewidth',1);%low and medium SR
hold on;
semilogx(cf(1,(sr>=18)),q10fit(sr>=18),'k+','markersize',10,'linewidth',1);%high SR
set(gca,'fontsize',16,'xlim',[80 20e3],'xtick',[100 1000 10000],'xticklabel',[0.1 1 10],'ticklength',[0.02 0.02],'yscale','log','ylim',[0.5 10],'ytick',[1 10],'yticklabel',[1 10]);
xlabel 'CF (kHz)';
ylabel ('Q_{10}','interpreter','tex');
box off;
%% Plot the Spont Rate distribution
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.0869    0.1100    0.3413    0.4756]);
xx = 0:1:200;
srd = hist(sr,xx);
srd = srd./sum(srd);
bar_hh = bar(xx,srd*100);
baselineh = get(bar_hh,'baseline');
set(baselineh,'linestyle','none');
set(bar_hh,'barwidth',1,'facecolor','k','edgecolor','none');
hold on;
gm = gmdistribution.fit(sr',3);
plot(xx,pdf(gm,xx')*100,'b','linewidth',3);
cols = {'g','r','c'};
for j=1:3
line(xx,gm.PComponents(j)*normpdf(xx,gm.mu(j),sqrt(gm.Sigma(j)))*100,'color',cols{j},'linewidth',3,'linestyle','--')
end
set(gca,'fontsize',16,'xlim',[0 200],'ticklength',[0.02 0.02]);
xlabel 'SPONTANEOUS RATE (sps)';
ylabel 'PERCENT';
box off;
%% Plot synchrony vs. CF
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
VSax = axes('ActivePositionProperty','outerposition');
mycmap = [[zeros(32,1);(linspace(0,1,16))';ones(16,1)] [(linspace(0,1,16))';ones(32,1);(linspace(1,0,16))'] [ones(16,1); linspace(1,0,16)'; zeros(32,1)]];
colormap(mycmap);
colormapyvals = mycmap;
minCF = 100;
maxCF = 5000;
colormapxvals = logspace(log10(minCF),log10(maxCF),64);
VSFreqLSR = FreqX(1,(sync(:,2)==1 & sr'<18));
VSSyncLSR = sync((sync(:,2)==1 & sr'<18),1);
VSFreqHSR = FreqX(1,(sync(:,2)==1 & sr'>=18));
VSSyncHSR = sync((sync(:,2)==1 & sr'>=18),1);
VSFreqNS = FreqX(1,(sync(:,2)==2));
VSSyncNS = sync((sync(:,2)==2),1);
% for i = 1:length(VSSyncLSR)
%     currentrgb = [interp1(colormapxvals,colormapyvals(:,1),VSFreqLSR(i)) interp1(colormapxvals,colormapyvals(:,2),VSFreqLSR(i)) interp1(colormapxvals,colormapyvals(:,3),VSFreqLSR(i))];
    lsrh = semilogx(VSax,VSFreqLSR,VSSyncLSR,'color','r','marker','o','markersize',10,'linestyle','none','linewidth',1); hold on;%Low and Medium SR
% end
% for i = 1:length(VSSyncHSR)
%     currentrgb = [interp1(colormapxvals,colormapyvals(:,1),VSFreqHSR(i)) interp1(colormapxvals,colormapyvals(:,2),VSFreqHSR(i)) interp1(colormapxvals,colormapyvals(:,3),VSFreqHSR(i))];
    hsrh = semilogx(VSax,VSFreqHSR,VSSyncHSR,'color','k','marker','*','markersize',10,'linestyle','none','linewidth',1); hold on;%Low and Medium SR

    semilogx(VSFreqNS,VSSyncNS,'+','color',[.6 .6 .6],'linewidth',1,'markersize',10);

box off;axis square;
set(gca,'fontsize',16,'xtick',[100 1000 10000],'xticklabel',[0.1 1 10],'ticklength',[0.02 0.02],'xlim',[50 20e3],'ylim',[0 1]);
xlabel 'CF (kHz)';
ylabel 'SYNCHRONIZATION INDEX';

%% Plot a nice figure showing impulse reponses from 100 Hz to 2000 Hz

FFs = [200:100:1000]; %CFs to search for
for i=1:length(FFs)
    %Find a fiber matching this CF
    dif = abs(FreqX(1,:)-FFs(i));
    ind = find(dif==min(dif));
    if length(ind)>1
        ind = ind(1);
    end
    %what is the exact CF for this fiber
    foundCF(i) = FreqX(1,ind);
    hh1(:,i) = (h1mat(:,ind)./max(abs(h1mat(:,ind))))+i;
    hh1mag(:,i) = (h1magmat(:,ind)-max(abs(h1magmat(:,ind))));
    %     hh1mag(h1z(:,i)<3,i) = NaN;
end
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.0869    0.1100    0.2800    0.7544]);
IRax = axes('ActivePositionProperty','outerposition');
plot(IRax,Time,hh1,'linewidth',3);
set(IRax,'ytick',[1:length(FFs)],'yticklabel',[FFs/1e3],'fontsize',20,'linewidth',1,'ticklength',[0.02 0.02]);
box off
ylabel 'FREQUENCY (kHz)';
xlabel 'TIME (ms)';

figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.0869    0.1100    0.2800    0.7544]);
MSax = axes('ActivePositionProperty','outerposition');
plot(MSax,hh1mag,FreqAx,'linewidth',3);
hold on;
hhh = plot(MSax,0,foundCF'./1e3,'o','markersize',15,'linewidth',3);
for i=1:length(hhh)
    set(hhh(i),'markerfacecolor',get(hhh(i),'color'));
end
set(MSax,'ytick',[FFs/1e3],'yticklabel',[FFs/1e3],'fontsize',20,'linewidth',1,'ticklength',[0.02 0.02],'xlim',[-20 1],'ylim',[0.1 1.1]);
box off
ylabel 'FREQUENCY (kHz)';
xlabel 'dB re max.';


%% Cross correlate pairs of impulse responses

%%% First simulate the orthodoxy with delay lines
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.1006    0.1067    0.2800    0.5000]);
ACax = axes('ActivePositionProperty','outerposition');
plot(ACax,[0 0],[-3.5 3.5],'k--','linewidth',1);hold on;
dt = Time(2)-Time(1);
%500 Hz fiber cross correlated with delays of [-0.8,-0.4,0,0.4,0.8] ms
delays = [-0.8:0.4:0.8];
delayind = [2:-1:-2];
cols = {'b','r','k','g','m'};
nsamps = round(delays./dt);%Number of delay samples to add
for i = 1:length(delays)
    h1ipsi = (hh1(:,4)-4);
    h1contra = h1ipsi;
    %Get the envelope - Hilbert transform
    h1ipsienv = abs(hilbert(h1ipsi));
    h1contraenv = h1ipsienv;
    if nsamps(i)<0%delay ipsi
        h1ipsi = [zeros(abs(nsamps(i)),1);h1ipsi];
        h1contra = [h1contra;zeros(abs(nsamps(i)),1)];
    elseif nsamps(i)>0
        h1contra = [zeros(abs(nsamps(i)),1);h1contra];
        h1ipsi = [h1ipsi;zeros(abs(nsamps(i)),1)];
    end
    msoout = xcorr(h1contra,h1ipsi);
    msoout = msoout./max(abs(msoout));
    msooutenv = abs(hilbert(msoout));
    msooutenv = msooutenv./max(msooutenv);
    timeout = [-dt*floor(length(msoout)/2):dt:-dt 0:dt:dt*floor(length(msoout)/2)];
    plot(ACax,timeout,msoout+delayind(i),'color',cols{i},'linewidth',3);hold on;
    plot(ACax,timeout,msooutenv+delayind(i),'--','color',[0.6 0.6 0.6],'linewidth',2);hold on;
    plot(ACax,timeout,-msooutenv+delayind(i),'--','color',[0.6 0.6 0.6],'linewidth',2);hold on;
    plot(ACax,delays(i),delayind(i)+1,'o','markersize',15,'color',cols{i},'markerfacecolor',cols{i});
    plot(ACax,delays(i),delayind(i)+1,'*','markersize',15,'color',[.6 .6 .6],'linewidth',2);
end
set(ACax,'xlim',[-5 5],'ytick',[-2:1:2],'yticklabel',fliplr(delays),'fontsize',20,'ticklength',[0.02 0.02],'ylim',[-3.5 3.5]);
xlabel 'ITD (ms)';
ylabel 'AXONAL DELAY (ms)';
box off;


%%% Now simulate the heresy with frequency differences
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.1006    0.1067    0.2800    0.5000]);
XCax = axes('ActivePositionProperty','outerposition');
plot(XCax,[0 0],[-3.5 3.5],'k--','linewidth',1);hold on;
dt = Time(2)-Time(1);
%500 Hz fiber cross correlated with delays of [300,400,500,600,700] Hz
delays = [200:-100:-200];
delayind = [2:-1:-2];
flipdelayind = fliplr(delayind);
cols = {'b','r','k','g','m'};
% cols = {'m','g','k','r','b'};
% nsamps = round(delays./dt);%Number of delay samples to add
for i = 1:length(delays)
    h1ipsi = (hh1(:,4)-4);
    h1contra = hh1(:,4-flipdelayind(i))-(4-flipdelayind(i));
    msoout = xcorr(h1contra,h1ipsi);
    msoout = msoout./max(abs(msoout));
    msooutenv = abs(hilbert(msoout));
    timeout = [-dt*floor(length(msoout)/2):dt:-dt 0:dt:dt*floor(length(msoout)/2)];
    [dum,maxind] = max(msoout);
    maxpeak = timeout(maxind);
    [dum,maxindenv] = max(msooutenv);
    maxpeakenv = timeout(maxindenv);
    plot(XCax,timeout,msoout+delayind(i),'color',cols{i},'linewidth',3);hold on;
    plot(XCax,timeout,msooutenv+delayind(i),'--','color',[.6 .6 .6],'linewidth',2);hold on;
    plot(XCax,timeout,-msooutenv+delayind(i),'--','color',[.6 .6 .6],'linewidth',2);hold on;
    plot(XCax,maxpeak,delayind(i)+1,'o','markersize',15,'color',cols{i},'markerfacecolor',cols{i});
    plot(XCax,maxpeakenv,delayind(i)+1,'*','markersize',15,'color',[.6 .6 .6],'linewidth',2);
end
set(XCax,'xlim',[-5 5],'ytick',[-2:1:2],'yticklabel',delays/1e3,'fontsize',20,'ticklength',[0.02 0.02],'ylim',[-3.5 3.5]);
xlabel 'ITD (ms)';
ylabel ('\Delta FREQUENCY (kHz)','interpreter','tex');
box off;


%% Save an index to the database - a look-up table of dominant frequencies
cd ('C:\LeuvenDataAnalysis');
save('MonRevPopData','AN');

%% Local Functions
    function [DFdifcor] = local_difcor_spec(IN)
        diffcorsr = 1000/dt;
        NFFT = 2^13;
        nsam = length(IN);
        DiffCorSpecFreq = diffcorsr*linspace(0,0.5,NFFT/2)/1000; % freq in kHz
        
        % compute & apply hann window, compute complex fft spectrum
        hanwin = hann(nsam);
        
        PredDifSpec = fft(IN.*hanwin,NFFT); % complex spec after windowing
        % magnitude ->power
        PredDifMagSpec = abs(PredDifSpec(1:NFFT/2)).^2;
        [dum,Pind] = max(PredDifMagSpec);
        DFdifcor = DiffCorSpecFreq(Pind);
    end
    function [otvect]=Trifilter(invect,nfw)
        nfwi = 2*floor(nfw/2) + 1;
        filt = zeros(1,nfwi);
        summ = 0;
        nfw2 = floor(nfwi/2);
        for jj=1:nfw2
            filt(jj) = jj;
            filt(nfwi+1-jj) = jj;
            summ = summ + 2*jj;
        end
        nfw3 = nfw2 + 1;
        filt(nfw3) = nfw3;
        summ = summ + nfw3;
        filt = filt./summ;
        svect = size(invect,2) + 2*nfw2;
        vect1 = zeros(1,svect);
        vect1(1:nfw2) = invect(1)*ones(1,nfw2);
        vect1(nfw3:svect-nfw2) = invect;
        vect1(svect-nfw2+1:svect) = invect(size(invect,2))*ones(1,nfw2);
        vect2 = conv(vect1, filt);
        otvect = vect2(2*nfw2+1:svect);
    end
    function [beta_out,gof_out]=local_fit_linear(varargin)
        if nargin ==3
            x_in = varargin{1};
            y_in = varargin{2};
            weight_in = varargin{3};
            weights = weight_in;
            weights = weights./max(weights);
            s = fitoptions('Method','LinearLeastSquares','Lower',[-Inf -Inf],'Upper',[Inf Inf],'Weights',weights,'robust','on');
        elseif nargin ==2
            x_in = varargin{1};
            y_in = varargin{2};
            s = fitoptions('Method','LinearLeastSquares','Lower',[-Inf -Inf],'Upper',[Inf Inf],'robust','on');
        end
        f = fittype({'x','1'},'coefficients',{'a','b'},'options',s);
        [beta_out,gof_out] = fit(x_in(:),y_in,f);
    end
    function [coeffpval,coefftval] = local_stats(xdata,ydata,beta,dfe)
        r=ydata-feval(beta,xdata);%residuals
        SSres=sum(r.^2);%sum square of residuals
        coeffSE = sqrt(SSres/dfe)/sqrt(sum((xdata-mean(xdata)).^2));%standard error of the regression coefficient
        coefftval = beta.a/coeffSE;
        coeffpval=2*(1-tcdf(abs(coefftval),dfe));
    end
end