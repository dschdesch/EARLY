
MonRevData = load('C:\Users\Mark\Dropbox\Disparity_MATLAB\MonRevData.mat');
load('C:\Users\Mark\Dropbox\Disparity_MATLAB\TemchinQ10Data.mat');
datadir = 'C:\Users\Mark\Dropbox\MonRevBandwidth';
datafiles = dir(fullfile(datadir,'*Bandwidth*'));

figure;set(gcf,'paperpositionmode','auto','activepositionproperty','outerposition');

semilogx(MonRevData.CF(MonRevData.SR>=18)/1000,MonRevData.Q10(MonRevData.SR>=18),'k+','linewidth',2,'markersize',10);
hold on;
semilogx(MonRevData.CF(MonRevData.SR<18)/1000,MonRevData.Q10(MonRevData.SR<18),'r+','linewidth',2,'markersize',10);

set(gca,'yscale','log','fontsize',14,'linewidth',1,'xtick',[0.1 1 10],...
    'xticklabel',[0.1 1 10],'tickdir','out','ticklength',[0.02 0.02],...
'activepositionproperty','outerposition','xlim',[0.08 20],'ytick',[1 10],...
'yticklabel',[1 10],'ylim',[0.2 20]);
xlabel 'Characteristic Frequency (kHz)';
ylabel ('Q_{10dB}','interpreter','tex');
box off;
axis square;
plot(temchin.freq/1000,temchin.q10,'go-','markersize',10,'markerfacecolor','g');


for i = 1:length(datafiles)
    load(fullfile(datadir,datafiles(i).name));
    data.ID{i} = datafiles(i).name(1:end-14);
    data.CF(i) = CF;
    data.DF{i} = DF;
    data.ERBHz{i} = ERBHz;
    data.ERBSPL{i} = ERBSPL;
    data.LowerHz{i} = LowerHz;
    data.UpperHz{i} = UpperHz;
    data.M{i} = M;
    data.SPL{i} = SPL;
    data.SPLpHz{i} = SPLpHz;
    data.SR(i) = SR;
    data.Th(i) = Th;
    data.SymmetryRatio{i} = SymmetryRatio;
    data.lo_slope{i} = lo_slope;
    data.up_slope{i} = up_slope;
end
nrunits = length(data.CF);

%Figures
figure;set(gcf,'paperpositionmode','auto');
semilogx(data.CF(data.SR>=18),data.Th(data.SR>=18),'ko');hold on;
semilogx(data.CF(data.SR<18),data.Th(data.SR<18),'r*');hold on;


cmapy = jet;
cmapx = linspace(log10(100),log10(3e3),64)';
cmapx = repmat(cmapx,1,3);

figure;set(gcf,'paperpositionmode','auto');
for i = 1:nrunits
    rgbvals = [interp1(cmapx(:,1),cmapy(:,1),log10(data.CF(i)*1000))...
        interp1(cmapx(:,2),cmapy(:,2),log10(data.CF(i)*1000))...
        interp1(cmapx(:,3),cmapy(:,3),log10(data.CF(i)*1000))];
    plot(data.ERBSPL{i},data.DF{i},'-','color',rgbvals,'linewidth',2);hold on;
end
set(gca,'fontsize',14,'linewidth',1,'tickdir','out','ticklength',[0.02 0.02])
ylabel 'Dominant Frequency (kHz)';
xlabel 'Sound Level (dB SPL in ERB)';
box off;
axis square;



figure;set(gcf,'paperpositionmode','auto');
for i = 1:nrunits
    rgbvals = [interp1(cmapx(:,1),cmapy(:,1),log10(data.CF(i)*1000))...
         interp1(cmapx(:,2),cmapy(:,2),log10(data.CF(i)*1000))...
         interp1(cmapx(:,3),cmapy(:,3),log10(data.CF(i)*1000))];

        plot(data.ERBSPL{i},log2(data.DF{i}./data.CF(i)),'o-','color',rgbvals,'linewidth',2);hold on;

end
set(gca,'fontsize',14,'linewidth',1,'tickdir','out','ticklength',[0.02 0.02])
xlabel 'Sound Level (dB SPL in ERB)';
ylabel ({'Dominant Frequency', '(Octaves re CF)'});
box off;
axis square;


figure;set(gcf,'paperpositionmode','auto');
for i = 1:nrunits
    if data.SR(i)<18
        plot(data.ERBSPL{i}, (data.DF{i}*1000)./data.ERBHz{i},'r-','linewidth',2);hold on;
    else
        plot(data.ERBSPL{i}, (data.DF{i}*1000)./data.ERBHz{i},'k-','linewidth',2);hold on;
    end
end
set(gca,'fontsize',14,'linewidth',1,'tickdir','out','ticklength',[0.02 0.02])
xlabel 'Sound Level (dB SPL in ERB)';
ylabel ('Q_{ERB}','interpreter','tex');
box off;
axis square;

    figure;set(gcf,'paperpositionmode','auto');
SPLs = 30:10:90;
cols = {'g','c','m','k','r','y','b'};
for j = 1:length(SPLs)
    for i = 1:nrunits
        if ismember(SPLs(j),data.SPL{i})
            SPLind = find(data.SPL{i}==SPLs(j));
            [h(j)] = semilogx(data.CF(i), data.SymmetryRatio{i}(SPLind),'o','markersize',10,...
                'color',cols{j},'markerfacecolor',cols{j});hold on;
        end
    end
end
legend(h,{'30','40','50','60','70','80','90 dBSPL'},'location','northwest');
set(gca,'fontsize',14,'linewidth',1,'tickdir','out','ticklength',[0.02 0.02],...
    'xlim',[0.08 3],'xtick',[.1 1],'xticklabel',[0.1 1]);
xlabel 'CF (kHz)';
ylabel ('Symmetry Index');
box off;
axis square;
plot([0.08 3],[1 1],'--','color',[.6 .6 .6],'linewidth',2);
plot([0.7 0.7],[0 5],'--','color',[.6 .6 .6],'linewidth',2);

figure;set(gcf,'paperpositionmode','auto');
SPLs = 30:10:90;
for j = 1:length(SPLs)
    for i = 1:nrunits
        if ismember(SPLs(j),data.SPL{i})
            SPLind = find(data.SPL{i}==SPLs(j));
            [h(j)] = semilogx(data.CF(i), data.lo_slope{i}(SPLind),'o','markersize',10,...
                'color',cols{j},'markerfacecolor',cols{j});hold on;
        end
    end
end
legend(h,{'30','40','50','60','70','80','90 dBSPL'},'location','northwest');
set(gca,'fontsize',14,'linewidth',1,'tickdir','out','ticklength',[0.02 0.02],...
    'xlim',[0.08 3],'xtick',[.1 1],'xticklabel',[0.1 1]);
xlabel 'CF (kHz)';
ylabel ('Low-side slope (dB/oct)');
box off;
axis square;

figure;set(gcf,'paperpositionmode','auto');
SPLs = 30:10:90;
for j = 1:length(SPLs)
    for i = 1:nrunits
        if ismember(SPLs(j),data.SPL{i})
            SPLind = find(data.SPL{i}==SPLs(j));
            [h(j)] = semilogx(data.CF(i), data.up_slope{i}(SPLind),'o','markersize',10,...
                'color',cols{j},'markerfacecolor',cols{j});hold on;
        end
    end
end
legend(h,{'30','40','50','60','70','80','90 dBSPL'},'location','northwest');
set(gca,'fontsize',14,'linewidth',1,'tickdir','out','ticklength',[0.02 0.02],...
    'xlim',[0.08 3],'xtick',[.1 1],'xticklabel',[0.1 1]);
xlabel 'CF (kHz)';
ylabel ('Low-side slope (dB/oct)');
box off;
axis square;
