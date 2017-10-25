function [] = PhilipFigure(AnNum,TCds,BFSds,SPLds)



dstc = dataset(AnNum,TCds);
dsbfs = struct(dataset(AnNum,BFSds));
for i = 1:length(SPLds)
    dsspl(i) = struct(dataset(AnNum,SPLds{i}));
end

TC = EvalTHR(dstc);

BFS.freq = dsbfs.Stimulus.IndepVar.Values;
BFS.freq = BFS.freq(~isnan(BFS.freq));
numvals = numel(BFS.freq);
ST = {dsbfs.Data.SpikeTimes{1:numvals}};
stimlen = dsbfs.Stimulus.Special.BurstDur;
for i = 1:numvals
    BFS.SpikeRate(i) = length(ST{i}(ST{i}<=stimlen))/(stimlen/1000);
end

for i = 1:length(SPLds)
    SPL(i).spl = dsspl(i).Stimulus.IndepVar.Values;
    SPL(i).freq = dsspl(i).Stimulus.Special.CarFreq;
    numvals = numel(SPL(i).spl);
    ST = {dsspl(i).Data.SpikeTimes{1:numvals}};
    stimlen = dsspl(i).Stimulus.Special.BurstDur;
    for j = 1:numvals
        SPL(i).SpikeRate(j) = length(ST{j}(ST{j}<=stimlen))/(stimlen/1000);
    end
end

cols = {'m','b','g','r'};

figure;set(gcf,'paperpositionmode','auto','color','none','units','normalized','position',[0.0650    0.1978    0.8400    0.5178]);

subplot(1,3,1);
semilogx(TC.fit.freq,TC.fit.thr,'k-','linewidth',3);hold on;
plot([TC.fit.cf TC.fit.cf],[-10 70],'--','color',[0.6 0.6 0.6],'linewidth',2);
for i = 1:length(SPLds)
    plot(SPL(i).freq,-5,'o','markersize',8,'color',cols{i},'markerfacecolor',cols{i},'linewidth',2);
end
set(gca,'layer','top','xlim',[100 4000],'ylim',[-10 70],'fontsize',18,'ytick',[0:20:60],'xtick',[100 1000],'xticklabel',[100 1000],'linewidth',2,'ticklength',[0.03 0.03],'tickdir','out');
ylabel 'THRESHOLD (dB SPL)';
xlabel 'FREQUENCY (Hz)';
box off;
axis square;

subplot(1,3,2);
plot([50 550],[50 550],'--','color',[0.6 0.6 0.6],'linewidth',2);hold on;
plot(BFS.freq,BFS.SpikeRate,'k-o','linewidth',2,'markersize',10);
for i = 1:length(SPLds)
    plot(SPL(i).freq,BFS.SpikeRate(BFS.freq==SPL(i).freq),'o','color',cols{i},'linewidth',2,'markersize',10,'markerfacecolor',cols{i});hold on;
end
set(gca,'fontsize',18,'xlim',[50 500],'xtick',[100:100:500],'ylim',[100 400],'ytick',[100:100:400],'linewidth',2,'ticklength',[0.03 0.03],'tickdir','out');
ylabel ('FIRING RATE (Spikes s^{-1})','interpreter','tex');
xlabel 'FREQUENCY (Hz)';
box off;
axis square;

subplot(1,3,3);
for i = 1:length(SPLds)
    plot([0 100],[SPL(i).freq SPL(i).freq],'--','color',cols{i},'linewidth',1);hold on;
    plot(SPL(i).spl,SPL(i).SpikeRate,'o-','color',cols{i},'linewidth',2,'markersize',10);hold on;
    plot(70,SPL(i).SpikeRate(SPL(i).spl==70),'o','markersize',10,'markerfacecolor',cols{i},'color',cols{i},'linewidth',2);hold on;
end
set(gca,'fontsize',18,'linewidth',2,'ticklength',[0.03 0.03],'tickdir','out');
box off;
axis square;
ylabel ('FIRING RATE (Spikes s^{-1})','interpreter','tex');
xlabel 'SOUND LEVEL (dB SPL)';