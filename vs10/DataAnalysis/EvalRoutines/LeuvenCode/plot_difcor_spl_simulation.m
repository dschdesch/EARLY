function [] = plot_difcor_spl_simulation

%First plot the AN data simulation
load('C:\work\BinRev\Analyses\AN_simulate_difcor_spl_pop.mat');
Xvals = [];
Yvals = [];
figure;
for i = 1:length(SPL)
    if mDF(i)<=1.5;
        x = SPL{i};if size(x,1)>size(x,2);x = x';end;
        y = BDDFnorm{i};if size(y,1)>size(y,2);y = y';end;
        ind = find(abs(y)<=0.5);
        y = y(ind);
        x = x(ind);
        if numel(x)>1
            plot(x,y,'color',[0.6 0.6 0.6],'linewidth',1);hold on;
            Xvals = [Xvals x];
            Yvals = [Yvals y];
        end
    end
end

spls = 30:10:80;
for i = 1:length(spls)
    yy(i) = mean(Yvals(Xvals==spls(i)));
    uu(i) = yy(i)+std(Yvals(Xvals==spls(i)));
    ll(i) = yy(i)-std(Yvals(Xvals==spls(i)));
end

plot(spls,yy,'ko-','linewidth',2,'markersize',10,'markerfacecolor','k');
plot([spls; spls],[uu;ll],'k-','linewidth',2);

%Now plot those MSO units where you have some SPL data for
%the revcors
MSOfiles = dir('C:\work\BinRev\Analyses\MSO_SPL*.mat');

for i = 1:length(MSOfiles)
    load(fullfile('C:\work\BinRev\Analyses',MSOfiles(i).name));
    yvals = BD.*difcor_df;
    yvals = yvals-yvals(uSPLs==70);
    plot(uSPLs,yvals,'r','linewidth',2);
end

set(gca,'xlim',[20 90],'fontsize',16,'linewidth',1,'layer','top');
box off;
xlabel 'SOUND LEVEL (dB SPL)';
ylabel ('\Delta(BD*DF) (cycles)','interpreter','tex');