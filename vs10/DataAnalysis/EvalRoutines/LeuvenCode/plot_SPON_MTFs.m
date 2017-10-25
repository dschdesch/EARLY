datapath = 'C:\Users\Mark\Dropbox\SPON\MTFs';

MTFfolders = dir(fullfile(datapath,'*_*'));

figure(1);
figure(2);
figure(3);
figure(4);

cols = {'k','r','b','g','m','c'};

for i = 1:length(MTFfolders)
    unitID = MTFfolders(i).name;
    
    MTFfiles = dir(fullfile(datapath,unitID,'*MTF*'));
    
    for j = 1:length(MTFfiles)
        
        
        load(fullfile(datapath,unitID,MTFfiles(j).name));
        
%         meanphase = wrapToPi((meanphase+0.5)*2*pi)/(2*pi);
        figure(1);set(gcf,'activepositionproperty','outerposition')
        semilogx(fm(Pval<0.001),VS(Pval<0.001),'ko','markersize',6,'linewidth',2);
        hold on;
        
        figure(2);set(gcf,'activepositionproperty','outerposition')
        semilogx(fm(Pval<0.001),spikerate(Pval<0.001),'k-','markersize',6,'linewidth',2);
        hold on;
        
        figure(3);set(gcf,'activepositionproperty','outerposition')
        semilogx(fm(Pval<0.001),syncrate(Pval<0.001),'k-','markersize',6,'linewidth',2);
        hold on;
        
        figure(4);set(gcf,'activepositionproperty','outerposition')
        semilogx(fm(Pval<0.001),meanphase(Pval<0.001),'ko','markersize',6,'linewidth',2);
        hold on;
        
        
        figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.0869    0.1100    0.3500    0.4978]);
        fm_Hz = round(fm*10)/10;
        [nrconds,nrreps]=size(raster_spikes);
        dY = 1/(nrconds*nrreps);
        dC = 1/nrconds;
        jj = 1;
        c_ix = 1;
        while jj<=nrconds
            if c_ix>6
                c_ix = 1;
            end
            for k = 1:nrreps
                Nsp = numel(raster_spikes{jj,k});
                plot([raster_spikes{jj,k};raster_spikes{jj,k}],repmat(((jj-1)*dC)+([k-1 k]*dY),Nsp,1)',cols{c_ix});
                hold on;
            end
            jj = jj+1;
            c_ix = c_ix+1;
        end
        set(gca,'fontsize',14,'ytick',[dC/2:dC:1],'yticklabel',fm_Hz,'tickdir','out');
        xlabel 'Time (ms)';
        ylabel 'Modulation Frequency (Hz)';
        box off;
        axis square;
    end
end
figure(1);set(gcf,'paperpositionmode','auto');
set(gca,'fontsize',14,'tickdir','out','activepositionproperty','outerposition',...
    'xlim',[1 2000],'xtick',[1 10 100 1000],'xticklabel',[1 10 100 1000],'ytick',[.98 .99 1],'ticklength',[0.02 0.02]);
axis square;box off;
xlabel 'Modulation Frequency (Hz)';
ylabel 'Vector Strength';

figure(2);set(gcf,'paperpositionmode','auto');
set(gca,'fontsize',14,'tickdir','out','activepositionproperty','outerposition',...
    'xlim',[1 2000],'xtick',[1 10 100 1000],'xticklabel',[1 10 100 1000],'ticklength',[0.02 0.02]);
axis square;box off;
xlabel 'Modulation Frequency (Hz)';
ylabel 'Firing Rate (spikes/s)';

figure(3);set(gcf,'paperpositionmode','auto');
set(gca,'fontsize',14,'tickdir','out','activepositionproperty','outerposition',...
    'xlim',[1 2000],'xtick',[1 10 100 1000],'xticklabel',[1 10 100 1000],'ticklength',[0.02 0.02]);
axis square;box off;
xlabel 'Modulation Frequency (Hz)';
ylabel 'Synchronized Rate (spikes/s)';

figure(4);set(gcf,'paperpositionmode','auto');
set(gca,'fontsize',14,'tickdir','out','activepositionproperty','outerposition',...
    'xlim',[1 2000],'xtick',[1 10 100 1000],'xticklabel',[1 10 100 1000],'ticklength',[0.02 0.02],'ylim',[-0.5 0.5],'ytick',-0.5:0.5:0.5);
axis square;box off;
xlabel 'Modulation Frequency (Hz)';
ylabel 'Phase (cycles)';






