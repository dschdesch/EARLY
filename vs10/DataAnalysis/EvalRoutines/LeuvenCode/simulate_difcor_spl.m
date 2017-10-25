function [CD,CP,difcor_df,BD,binauralSPLs] = simulate_difcor_spl(annum,unitnum,plotflag)

%A good example is:
%simulate_difcor_spl({'L16503','L16503'},{'19','17'})

datadir = 'C:\work\BinRev\MonRev\OLD\';
data = cell(1,2);
SPLs = cell(1,2);
indx = cell(1,2);
spec_crit = 3;
for i = 1:2
    data{i} = load(fullfile(datadir,[annum{i} '_' unitnum{i} '_MonRev']));
    SPLs{i} = unique(data{i}.ParamsOut.noiseSPLs);
end
[binauralSPLs,indx{1},indx{2}] = intersect(SPLs{1},SPLs{2});
nrSPLs = length(binauralSPLs);
cols = {[0 0 1],[0 1 0],[1 0 0],[0 1 1],[1 0 1],[0 0 0],[.8 .8 0],[0 0.5 0],[1 0.5 0],[0.5 0 0.5]};
time = data{1}.ParamsOut.Time;
h1freq = data{1}.ParamsOut.h1ffax;
dt = time(2)-time(1);

for i = 1:2
    h1{i} = data{i}.ParamsOut.h1filt(:,indx{i});
    h1df{i} = data{i}.ParamsOut.domfreq(indx{i});
    h1phase{i} = data{i}.ParamsOut.h1phase(:,indx{i});
    h1z{i} = data{i}.ParamsOut.h1zscore(:,indx{i});
    h1mag{i} = data{i}.ParamsOut.h1mag(:,indx{i});
end

%make the "ipsi" ear always the higher frequency of the two
if mean(h1df{1})<mean(h1df{2})
    temp = data{2};
    data{2} = data{1};
    data{1} = temp;clear temp;
    temp = SPLs{2};
    SPLs{2} = SPLs{1};
    SPLs{1} = temp;clear temp;
    temp = indx{2};
    indx{2} = indx{1};
    indx{1} = temp;clear temp;
    temp = h1{2};
    h1{2} = h1{1};
    h1{1} = temp;clear temp;
    temp = h1df{2};
    h1df{2} = h1df{1};
    h1df{1} = temp;clear temp;
    temp = h1phase{2};
    h1phase{2} = h1phase{1};
    h1phase{1} = temp;clear temp;
    temp = h1z{2};
    h1z{2} = h1z{1};
    h1z{1} = temp;clear temp;
    temp = h1mag{2};
    h1mag{2} = h1mag{1};
    h1mag{1} = temp;clear temp;
end

% %Get the pseudo-binaural CD and CP from the monaural phase-frequency
% %functions
for i = 1:2
    weights{i} = zeros(size(h1phase{i}));
    for j = 1:nrSPLs
        ind(i,j) = find(data{i}.ParamsOut.h1ffax==h1df{i}(j));
        try
        startind(i,j) = find(h1z{i}(1:ind(i,j)-1,j)>=spec_crit,1,'first')+1;
        catch
            foo;
        end
        stopind(i,j) = find(h1z{i}(ind(i,j)+1:end,j)>=spec_crit,1,'last')+ind(i,j)-1;
        weights{i}(startind(i,j):stopind(i,j),j) = 10.^(h1mag{i}(startind(i,j):stopind(i,j),j)/20);
        weights{i}(:,j) = weights{i}(:,j)/max(weights{i}(:,j));
        weights{i}(find(h1z{i}(startind(i,j):stopind(i,j),j)<spec_crit)+startind(i,j)-1,j)=0;
    end
end
for i = 1:nrSPLs
    startat = max(startind(:,i));
    stopat = min(stopind(:,i));
    Phase = [h1phase{1}(startat:stopat,i) h1phase{2}(startat:stopat,i)];
    Freq = h1freq(startat:stopat);
    BPhase_cycles = (Phase(:,2)-Phase(:,1))./(2*pi);
    Bweights = prod([weights{1}(startat:stopat,:) weights{2}(startat:stopat,:)],2);
    Bweights = Bweights/max(Bweights);
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
    [CP(i),CD(i),dummy1,dummy2] = fit_circ_lin(Freq,BPhase_cycles,Bweights');
end

%Get the xcorr at each SPL
for i = 1:nrSPLs
    difcor(:,i) = xcorr(h1{1}(:,i),h1{2}(:,i),'coeff');
    difcor_df(i) = local_difcor_spec(difcor(:,i));
end
itd = -[fliplr(-time) time(2:end)];
%scale the revcors for plotting
for i = 1:2
    for j = 1:nrSPLs
        h1{i}(:,j) = ((2*j)-2)+(h1{i}(:,j)/max(abs(h1{i}(:,j))));
    end
end
%scale the difcors for plotting
for j = 1:nrSPLs
    difcor(:,j) = ((2*j)-2)+difcor(:,j)./max(abs(difcor(:,j)));
    [bd(j).y,bd(j).ind] = max(difcor(:,j));
    bd(j).x = itd(bd(j).ind);
end
BD = [bd(:).x];


if plotflag
    figure;
    subplot(2,4,[1 5]);
    plot(time,h1{1},'b','linewidth',2);hold on;
    plot(time,h1{2},'r','linewidth',2);
    set(gca,'xlim',[0 10],'ylim',[-1.5 2*(nrSPLs-1)+1.5],'fontsize',16,'ytick',[0:2:2*(nrSPLs-1)+1],'yticklabel',binauralSPLs,'linewidth',1,'layer','top');
    box off;
    xlabel 'TIME (ms)';
    ylabel 'SOUND LEVEL (dB SPL)';
    
    subplot(2,4,[2 6]);
    plot([0 0],[-1.5 2*(nrSPLs-1)+1.5],'--','color',[0.6 0.6 0.6]);hold on;
    plot(itd,difcor,'k','linewidth',2);hold on;
    n = 1;
    for i = 1:nrSPLs
        plot(bd(i).x, bd(i).y,'marker','o','markersize',10,'markerfacecolor',cols{n},'markeredgecolor', cols{n});
        n = n+1;
        if n >10
            n = 1;
        end
    end
    set(gca,'xlim',[-4 4],'ylim',[-1.5 2*(nrSPLs-1)+1.5],'fontsize',16,'ytick',[0:2:2*(nrSPLs-1)+1],'yticklabel',binauralSPLs,'linewidth',1,'layer','top');
    box off;
    xlabel 'ITD (ms)';
    ylabel 'SOUND LEVEL (dB SPL)';
    
    %Dominant frequency of monaural revcors and peudo-binaural difcor as a
    %fuction of SPL
    subplot(2,4,3);
    plot(binauralSPLs,h1df{1},'b^','markerfacecolor','none','markersize',10,'linewidth',1);hold on;
    plot(binauralSPLs,h1df{2},'rv','markerfacecolor','none','markersize',10,'linewidth',1);hold on;
    n = 1;
    for i = 1:nrSPLs
        plot(binauralSPLs(i),difcor_df(i),'marker','o','markerfacecolor',cols{n},'markeredgecolor',cols{n},'markersize',10);
        n = n+1;
        if n >10
            n = 1;
        end
    end
    set(gca,'xlim',[min(binauralSPLs)-5 max(binauralSPLs)+5],'ylim',[0.7 0.9],'fontsize',16,'xtick',binauralSPLs(1:2:end),'xticklabel',binauralSPLs(1:2:end),'linewidth',1,'layer','top');
    box off;
    ylabel 'DOMINANT FREQUENCY (kHz)';
    xlabel 'SOUND LEVEL (dB SPL)';
    
    %BD and BD*DF as a function of SPL
    subplot(2,4,7);
    n = 1;
    for i = 1:nrSPLs
        plot(binauralSPLs(i),bd(i).x,'marker','o','markerfacecolor',cols{n},'markeredgecolor',cols{n},'markersize',10);hold on;
        n = n+1;
        if n >10
            n = 1;
        end
    end
    set(gca,'xlim',[min(binauralSPLs)-5 max(binauralSPLs)+5],'ylim',[0 1],'fontsize',16,'xtick',binauralSPLs(1:2:end),'xticklabel',binauralSPLs(1:2:end),'linewidth',1,'layer','top');
    box off;
    ylabel 'BEST DELAY (ms)';
    xlabel 'SOUND LEVEL (dB SPL)';
    
    %Characteristic Phase as a function of SPL
    subplot(2,4,4);
    n = 1;
    for i = 1:nrSPLs
        plot(binauralSPLs(i),CP(i),'marker','o','markerfacecolor',cols{n},'markeredgecolor',cols{n},'markersize',10);hold on;
        n = n+1;
        if n >10
            n = 1;
        end
    end
    set(gca,'xlim',[min(binauralSPLs)-5 max(binauralSPLs)+5],'ylim',[-0.5 0.5],'fontsize',16,'xtick',binauralSPLs(1:2:end),'xticklabel',binauralSPLs(1:2:end),'linewidth',1,'layer','top');
    box off;
    ylabel 'CHARACTERISTIC PHASE (cycles)';
    xlabel 'SOUND LEVEL (dB SPL)';
    
    %Characteristic Delay as a function of SPL
    subplot(2,4,8);
    n = 1;
    for i = 1:nrSPLs
        plot(binauralSPLs(i),CD(i),'marker','o','markerfacecolor',cols{n},'markeredgecolor',cols{n},'markersize',10);hold on;
        n = n+1;
        if n >10
            n = 1;
        end
    end
    set(gca,'xlim',[min(binauralSPLs)-5 max(binauralSPLs)+5],'ylim',[0 1],'fontsize',16,'xtick',binauralSPLs(1:2:end),'xticklabel',binauralSPLs(1:2:end),'linewidth',1,'layer','top');
    box off;
    ylabel 'CHARACTERISTIC DELAY (ms)';
    xlabel 'SOUND LEVEL (dB SPL)';
    
    
    %
    % subplot(3,3,3);
    for i = 1:2
        TC(i).x.data = data{i}.ParamsOut.TC.thr.freq;
        TC(i).y.data = data{i}.ParamsOut.TC.thr.thr;
        TC(i).x.fit = data{i}.ParamsOut.TC.fit.freq;
        TC(i).y.fit = data{i}.ParamsOut.TC.fit.thr;
    end
end
% semilogx(TC(1).x.data,TC(1).y.data,'b+','linewidth',1);hold on;
% semilogx(TC(1).x.fit,TC(1).y.fit,'b-','linewidth',2);hold on;
% semilogx(TC(2).x.data,TC(2).y.data,'r+','linewidth',1);hold on;
% semilogx(TC(2).x.fit,TC(2).y.fit,'r-','linewidth',2);hold on;

    function [DFdifcor] = local_difcor_spec(IN)
        diffcorsr = 1000/dt;
        NFFT = 2^21;
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
end

