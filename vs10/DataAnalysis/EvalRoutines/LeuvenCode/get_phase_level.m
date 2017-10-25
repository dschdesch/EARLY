function [] = get_phase_level(AnNum,UnNum)
datapath = 'C:\Users\Mark\Dropbox\MonRev\';
AnUnStr = [AnNum '_' UnNum];
datafolder = [datapath, AnNum '_' UnNum];
datafiles = dir(fullfile(datafolder,'*MonRev*.mat'));
nrdata = length(datafiles);
phase_uw = nan(2^12,nrdata);
load(fullfile(datafolder,[AnUnStr '_TC.mat']));
for i = 2:nrdata
    data(1) = load(fullfile(datafolder,datafiles(i-1).name));
    data(2) = load(fullfile(datafolder,datafiles(i).name));
    
    for j = 1:2
        spec = fft(data(j).kernels.h1c_mod,2^13);
        temp = angle(spec);
        phase(:,j) = temp(1:2^12);
        temp = 20*log10(abs(spec));
        mag(:,j) = temp(1:2^12);
        mag(:,j) = mag(:,j)-max(mag(:,j));
        mag(mag(:,j)<-20,j) = nan;
    end
    phase(isnan(mag))=nan;
    freq = data(1).kernels.h1ffax;
    SPL(i) = data(2).kernels.SPL;
    SPL(i-1) = data(1).kernels.SPL;
    [~,ind] = max(mag);
    weights = 10.^(mag/20);
    weights(isnan(weights)) = 0;
    cycles = phase./(2*pi);
    for j = 1:2
        try
        ind_start(j) = find(isnan(phase((1:ind(j)),j)),1,'last')+1;
        catch
            ind_start(j) = 1;
        end
        ind_stop(j) = find(isnan(phase((ind(j):end),j)),1,'first')+ind(j)-2;
        phase(1:ind_start(j)-1,j) = NaN;
        phase(ind_stop(j)+1:end,j) = NaN;
        cycles(1:ind_start(j)-1,j) = NaN;
        cycles(ind_stop(j)+1:end,j) = NaN;
    end
    if i==2
        for j = 1:2
            phase_uw(ind_start(j):ind_stop(j),j) = unwrap(phase(ind_start(j):ind_stop(j),j))/(2*pi);
            GD(j) = data(j).kernels.h1coeffs.Monaural.Delay;
        end
    else
        phase_uw(ind_start(2):ind_stop(2),i) = unwrap(phase(ind_start(2):ind_stop(2),2))/(2*pi);
        GD(i) = data(2).kernels.h1coeffs.Monaural.Delay;
    end
    startat = max(ind_start);
    stopat = min(ind_stop);
    IAP = diff(cycles,1,2);%inter-aural phase
    Bweights = zeros(length(weights),1);
    Bweights(startat:stopat) = prod(weights(startat:stopat,:),2);
    Bweights = Bweights./max(Bweights);
    [dummy,indx]=max(Bweights);
    ind1 = find(Bweights(1:indx-1)==0,1,'last')+1;
    ind2 = max([find(Bweights(indx+1:end)==0,1,'first')+indx-1, length(Bweights)]);
    if isempty(ind1)
        ind1 = 1;
    end
    if isempty(ind2)
        ind2 = length(Bweights);
    end
    Bweights(1:ind1-1)=0;
    Bweights(ind2+1:end)=0;
    [CP,CD,BinError,BinN] = fit_circ_lin(freq(startat:stopat),IAP(startat:stopat),Bweights(startat:stopat)');
    %Now that you have the real CP you can go back and correct h1phasemask
    %to give the real relationship between the two phase-frequency
    %functions
    %find the inter-aural phase relation at the binaural dominant frequency
    [~,indx] = min(abs(freq-(TC.fit.cf/1000)));
    startphase = (freq(indx)*CD)+CP;
    %equalize the phase at that frequency
    phasenow = diff(phase_uw(indx,[i-1 i]));
    phase_uw(:,i) = phase_uw(:,i)-phasenow;
    %correct to the right phase
    phase_uw(:,i) = phase_uw(:,i)+startphase;
end
% % 
% phase_uw = phase_uw(:,1:2:end);
% SPL = SPL(1:2:end);

meandelay = freq'*mean(-GD);
CF = TC.fit.cf/1000;
meandelay = meandelay-interp1(freq,meandelay,CF);


figure;
minSPL = min(SPL);
maxSPL = max(SPL);
X = linspace(minSPL,maxSPL,64)';
Y = colormap;
h = plot(freq,phase_uw,'linewidth',2);
for i = 1:length(h)
    rgbvals = [interp1(X,Y(:,1),SPL(i)) interp1(X,Y(:,2),SPL(i)) interp1(X,Y(:,3),SPL(i))];
    set(h(i),'color',rgbvals);
end
set(gca,'xlim',[0 3],'fontsize',14,'activepositionproperty','outerposition','clim',[0 1]);
box off;
axis square;
xlabel 'FREQUENCY (kHz)';
ylabel 'PHASE (cycles)';
h = colorbar;
set(h,'fontsize',14,'ytick',[0:0.25:1],'yticklabel',linspace(minSPL,maxSPL,5));
ylabel(h,'OVERALL LEVEL (dB SPL)');

figure;
minSPL = min(SPL);
maxSPL = max(SPL);
X = linspace(minSPL,maxSPL,64)';
Y = colormap;
h = plot(freq,bsxfun(@minus,phase_uw,meandelay),'linewidth',2);
for i = 1:length(h)
    rgbvals = [interp1(X,Y(:,1),SPL(i)) interp1(X,Y(:,2),SPL(i)) interp1(X,Y(:,3),SPL(i))];
    set(h(i),'color',rgbvals);
end
set(gca,'xlim',[0 3],'fontsize',14,'activepositionproperty','outerposition','clim',[0 1]);
box off;
axis square;
xlabel 'FREQUENCY (kHz)';
ylabel 'PHASE (cycles)';
h = colorbar;
set(h,'fontsize',14,'ytick',[0:0.25:1],'yticklabel',linspace(minSPL,maxSPL,5));
ylabel(h,'OVERALL LEVEL (dB SPL)');

maxphase = max(phase_uw(:));
minphase = min(phase_uw(:));
phaserange = maxphase-minphase;
nrlines = floor(phaserange/0.125);

[X,Y] = meshgrid(SPL,freq);
[XX,YY] = meshgrid(min(SPL):1:max(SPL),freq);
ZZ = interp2(X,Y,phase_uw,XX,YY);

figure;
pcolor(XX,YY,ZZ);view(0,90);shading flat;hold on;
contour(XX,YY,ZZ,nrlines,'linewidth',2,'linecolor','k');
grid on;
set(gca,'layer','top','fontsize',14,'activepositionproperty','outerposition',...
    'ylim',[0 3],'xlim',[min(SPL)-5 max(SPL)+5]);
ylabel 'FREQUENCY (kHz)';
xlabel 'OVERALL LEVEL (dB SPL)';
axis square;
h = colorbar;
set(h,'fontsize',14);
ylabel(h,'PHASE (cycles)');

return;






































