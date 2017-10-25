datapath = 'C:\Users\Mark\Dropbox\SPON\RMs';

RMfolders = dir(fullfile(datapath,'*_*'));


for i = 5:length(RMfolders)
    unitID = RMfolders(i).name;
    
    RMfiles = dir(fullfile(datapath,unitID,'*RM*'));
    
    for j = 1:length(RMfiles)
        load(fullfile(datapath,unitID,RMfiles(j).name));
        figure(1);set(gcf,'paperpositionmode','auto');
        pcolor(freqs/1000,levels,OnRate);
        shading flat;
        
        set(gca,'fontsize',14,'tickdir','out','xscale','log','layer','top',...
            'linewidth',1,'xlim',[0.1 20],'xtick',[0.1 1 10],'xticklabel',[0.1 1 10],'ticklength',[0.02 0.02],'activepositionproperty','outerposition');
        xlabel 'Frequency (kHz)';
        ylabel 'Level (dB SPL)';
        box off;
        axis square;
        imagefilename = [fullfile(datapath,unitID) '_ON.eps'];
        print ('-f1','-depsc','-painters','-r600',imagefilename)
        
        
        figure(2);set(gcf,'paperpositionmode','auto');
        pcolor(freqs/1000,levels,OffRate);
        shading flat;
        
        set(gca,'fontsize',14,'tickdir','out','xscale','log','layer','top',...
            'linewidth',1,'xlim',[0.1 20],'xtick',[0.1 1 10],'xticklabel',[0.1 1 10],'ticklength',[0.02 0.02],'activepositionproperty','outerposition');
        xlabel 'Frequency (kHz)';
        ylabel 'Level (dB SPL)';
        box off;
        axis square;
        imagefilename = [fullfile(datapath,unitID) '_OFF.eps'];
        print ('-f2','-depsc','-painters','-r600',imagefilename)
        
        pause(1);
    end
end

