function [] = AN_CDCP_SPL

cd 'C:\work\BinRev\Analyses';

SPLs = 30:10:80;

figure;
for i = 1:length(SPLs)
    subplot(2,3,i);
    load(['AN_CD_CP_',num2str(SPLs(i))]);
    plot(CP(deltaF>0),CDcycles(deltaF>0),'k+')
    set(gca,'ylim',[-2 2],'xlim',[-0.5 0.5],'fontsize',16);
    box off;
    ylabel 'CD (cycles)';
    xlabel 'CP (cycles)';
    title ([num2str(SPLs(i)) 'dB SPL']);
    length(CDcycles(deltaF>0))
end