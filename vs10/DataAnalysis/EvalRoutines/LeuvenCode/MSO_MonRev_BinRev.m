function [] = MSO_MonRev_BinRev
datapath = 'C:\Users\Mark\Dropbox\BinRev';
anapath = 'C:\LeuvenDataAnalysis';
cd(datapath);
folderlist = dir('*_*');
count = 0;
for i  = 1:length(folderlist)
    foundit = 0;
    IPSI = [];CONTRA = [];kernels = [];
    unitname = folderlist(i).name;
    if exist(fullfile(datapath,unitname,[unitname '_MonRevCONTRA_70.00dB.mat']),'file')
        load(fullfile(datapath,unitname,[unitname '_MonRevCONTRA_70.00dB.mat']))
        foundit = 1;
    end
    if exist(fullfile(datapath,unitname,[unitname '_MonRevIPSI_70.00dB.mat']),'file')
        load(fullfile(datapath,unitname,[unitname '_MonRevIPSI_70.00dB.mat']))
        foundit = 1;
    end
    if foundit
        count = count+1;
        load(fullfile(datapath,unitname,[unitname '_BinRev_70.00dB.mat']));
        load(fullfile(datapath,unitname,[unitname '_NDF_70.00dB.mat']));
        if ~isempty(IPSI)
            if length(kernels.h1c)~=length(IPSI.h1c)
                if length(IPSI.h1c)>length(kernels.h1c)
                    monauralrev = IPSI.h1c;
                    binauralrev = interp1(kernels.Time,kernels.h1c(:,1),IPSI.Time)';
                    binauralrev(isnan(binauralrev)) = 0;
                else
                    binauralrev = kernels.h1c(:,1);
                    monauralrev = interp1(IPSI.Time,IPSI.h1c,kernels.Time)';
                    monauralrev(isnan(monauralrev)) = 0;
                end
            else
                binauralrev = kernels.h1c(:,1);
                monauralrev = IPSI.h1c;
            end
            ipsi.corr(count) = corr(binauralrev,monauralrev);
            ipsi.mon.h1{count} = IPSI.h1;
            ipsi.mon.h1c{count} = IPSI.h1c;
            ipsi.mon.h1wf{count} = IPSI.h1wf;
            ipsi.mon.time{count} = IPSI.Time;
            ipsi.mon.GD(count) = IPSI.h1coeffs.Monaural.Delay;
            ipsi.mon.DF(count) = IPSI.h1domfreq;
            ipsi.mon.mag{count} = IPSI.h1magmask;
            ipsi.mon.phase{count} = IPSI.h1phase_uw;
            ipsi.mon.freq{count} = IPSI.h1ffax;
        else
            ipsi.corr(count) = nan;
            ipsi.mon.h1{count} = [];
            ipsi.mon.h1c{count} = [];
            ipsi.mon.h1wf{count} = [];
            ipsi.mon.time{count} = [];
            ipsi.mon.GD(count) = nan;
            ipsi.mon.DF(count) = nan;
            ipsi.mon.mag{count} = [];
            ipsi.mon.phase{count} = [];
            ipsi.mon.freq{count} = [];
        end
        if ~isempty(CONTRA)
            if length(kernels.h1c)~=length(CONTRA.h1c)
                if length(CONTRA.h1c)>length(kernels.h1c)
                    monauralrev = CONTRA.h1c;
                    binauralrev = interp1(kernels.Time,kernels.h1c(:,2),CONTRA.Time)';
                    binauralrev(isnan(binauralrev)) = 0;
                else
                    binauralrev = kernels.h1c(:,2);
                    monauralrev = interp1(CONTRA.Time,CONTRA.h1c,kernels.Time)';
                    monauralrev(isnan(monauralrev)) = 0;
                end
            else
                binauralrev = kernels.h1c(:,2);
                monauralrev = CONTRA.h1c;
            end
            contra.corr(count) = corr(binauralrev,monauralrev);
            contra.mon.h1{count} = CONTRA.h1;
            contra.mon.h1c{count} = CONTRA.h1c;
            contra.mon.h1wf{count} = CONTRA.h1wf;
            contra.mon.time{count} = CONTRA.Time;
            contra.mon.GD(count) = CONTRA.h1coeffs.Monaural.Delay;
            contra.mon.DF(count) = CONTRA.h1domfreq;
            contra.mon.mag{count} = CONTRA.h1magmask;
            contra.mon.phase{count} = CONTRA.h1phase_uw;
            contra.mon.freq{count} = CONTRA.h1ffax;
        else
            contra.corr(count) = nan;
            contra.mon.h1{count} = [];
            contra.mon.h1c{count} = [];
            contra.mon.h1wf{count} = [];
            contra.mon.time{count} = [];
            contra.mon.GD(count) = nan;
            contra.mon.DF(count) = nan;
            contra.mon.mag{count} = [];
            contra.mon.phase{count} = [];
            contra.mon.freq{count} = [];
        end
        difcorDF(count) = ndf.difcor.peakhz/1000;
        ipsi.bin.h1{count} = kernels.h1(:,1);
        contra.bin.h1{count} = kernels.h1(:,2);
        ipsi.bin.h1c{count} = kernels.h1c(:,1);
        contra.bin.h1c{count} = kernels.h1c(:,2);
        ipsi.bin.h1wf{count} = kernels.h1wf(:,1);
        contra.bin.h1wf{count} = kernels.h1wf(:,2);
        ipsi.bin.time{count} = kernels.Time;
        contra.bin.time{count} = kernels.Time;
        ipsi.bin.GD(count) = kernels.h1coeffs.Monaural.Delay(1);
        contra.bin.GD(count) = kernels.h1coeffs.Monaural.Delay(2);
        ipsi.bin.DF(count) = kernels.h1domfreq(1);
        contra.bin.DF(count) = kernels.h1domfreq(2);
        ipsi.bin.mag{count} = kernels.h1magmask(:,1);
        contra.bin.mag{count} = kernels.h1magmask(:,2);
        ipsi.bin.phase{count} = kernels.h1phase_uw(:,1);
        contra.bin.phase{count} = kernels.h1phase_uw(:,2);
        ipsi.bin.freq{count} = kernels.h1ffax;
        contra.bin.freq{count} = kernels.h1ffax;
    end
end
cd (anapath);

A = [ipsi.bin.DF contra.bin.DF];
B = [ipsi.mon.DF contra.mon.DF];
COD_DF = Coeff_Determination(A(~isnan(B)),B(~isnan(B)));
[~,p_DF] = ttest(A,B);

A = [ipsi.bin.GD contra.bin.GD];
B = [ipsi.mon.GD contra.mon.GD];
COD_GD = Coeff_Determination(A(~isnan(B)),B(~isnan(B)));
[~,p_GD] = ttest(A,B);

for i = 5
    figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.1775    0.0522    0.5350    0.8567]);
    subplot(3,3,1);hold on;
    plot(ipsi.bin.time{i},ipsi.bin.h1wf{i},'b-','linewidth',3);
    plot(ipsi.mon.time{i},ipsi.mon.h1wf{i},'--','linewidth',3,'color',[.6 .6 .6]);
    set(gca,'xlim',[2 15],'tickdir','out','ticklength',[0.02 0.02],'fontsize',12,'linewidth',1);
    xlabel 'Time (ms)';
    ylabel 'Spikes/s/Pa';
    text(0.7,0.9,sprintf('\\rho: %2.2f',ipsi.corr(5)),'interpreter','tex','units','normalized','fontsize',12);
    
    subplot(3,3,2);hold on;
    plot(contra.bin.time{i},contra.bin.h1wf{i},'r-','linewidth',3);
    plot(contra.mon.time{i},contra.mon.h1wf{i},'--','linewidth',3,'color',[.6 .6 .6]);
    set(gca,'xlim',[2 15],'tickdir','out','ticklength',[0.02 0.02],'fontsize',12,'linewidth',1);
    xlabel 'Time (ms)';
    ylabel 'Spikes/s/Pa';
    text(0.7,0.9,sprintf('\\rho: %2.2f',contra.corr(5)),'interpreter','tex','units','normalized','fontsize',12);
    
    subplot(3,3,4);hold on;
    semilogx(ipsi.bin.freq{i},ipsi.bin.mag{i},'b-','linewidth',3);
    semilogx(ipsi.mon.freq{i},ipsi.mon.mag{i},'--','linewidth',3,'color',[.6 .6 .6]);
    set(gca,'xscale','log','xtick',[0.1 1],'xticklabel',[0.1 1],...
        'tickdir','out','ticklength',[0.02 0.02],'xlim',[0.05 3],'ylim',[-20 2],'fontsize',12,'linewidth',1);
    xlabel 'Frequency (kHz)';
    ylabel 'Gain (dB)';
    
    subplot(3,3,5);hold on;
    semilogx(contra.bin.freq{i},contra.bin.mag{i},'r-','linewidth',3);
    semilogx(contra.mon.freq{i},contra.mon.mag{i},'--','linewidth',3,'color',[.6 .6 .6]);
    set(gca,'xscale','log','xtick',[0.1 1],'xticklabel',[0.1 1],...
        'tickdir','out','ticklength',[0.02 0.02],'xlim',[0.05 3],'ylim',[-20 2],'fontsize',12,'linewidth',1);
    xlabel 'Frequency (kHz)';
    ylabel 'Gain (dB)';
    
    subplot(3,3,7);hold on;
    semilogx(ipsi.bin.freq{i},ipsi.bin.phase{i},'b-','linewidth',3);
    semilogx(ipsi.mon.freq{i},ipsi.mon.phase{i},'--','linewidth',3,'color',[.6 .6 .6]);
    set(gca,'xscale','linear','fontsize',12,'tickdir','out','ticklength',[0.02 0.02],'linewidth',1,'xlim',[0.1 1]);
    ylabel 'Phase (cycles)';
    xlabel 'Frequency (kHz)';
    
    subplot(3,3,8);hold on;
    semilogx(contra.bin.freq{i},contra.bin.phase{i},'r-','linewidth',3);
    semilogx(contra.mon.freq{i},contra.mon.phase{i},'--','linewidth',3,'color',[.6 .6 .6]);
    set(gca,'xscale','linear','fontsize',12,'linewidth',1,'tickdir','out','ticklength',[0.02 0.02],'xlim',[0.1 1]);
    ylabel 'Phase (cycles)';
    xlabel 'Frequency (kHz)';
    
    subplot(3,3,3);hold on;
    plot(ipsi.bin.DF,ipsi.mon.DF,'bo','markersize',8,'linewidth',1);
    plot(contra.bin.DF,contra.mon.DF,'ro','markersize',8,'linewidth',1);
    plot([0.1 1.1],[0.1 1.1],'k--','linewidth',1);
    set(gca,'xlim',[0.1 1.1],'ylim',[0.1 1.1],'fontsize',12,'tickdir','out','ticklength',[0.02 0.02],'tickdir','out','linewidth',1);
    ylabel 'Dom. Freq. {monaural} (kHz)';
    xlabel 'Dom. Freq. {binaural} (kHz)';
    text(0.6,0.3,{sprintf('R^2: %2.2f',COD_DF),sprintf('p: %2.2f',p_DF)},'fontsize',12,'interpreter','tex','units','normalized');
    
    subplot(3,3,6);hold on;
    plot(ipsi.bin.GD,ipsi.mon.GD,'bo','markersize',8,'linewidth',1);
    plot(contra.bin.GD,contra.mon.GD,'ro','markersize',8,'linewidth',1);
    plot([0 10],[0 10],'k--','linewidth',1);
    set(gca,'xlim',[4 9],'ylim',[4 9],'fontsize',12,'tickdir','out','ticklength',[0.02 0.02],'tickdir','out','linewidth',1);
    ylabel 'Group Delay {monaural} (ms)';
    xlabel 'Group Delay {binaural} (ms)';
    text(0.6,0.3,{sprintf('R^2: %2.2f',COD_GD),sprintf('p: %2.2f',p_GD)},'fontsize',12,'interpreter','tex','units','normalized');
    
    subplot(3,3,9);hold on;
    semilogy(ipsi.corr,difcorDF,'bo','markersize',8,'linewidth',1);
    semilogy(contra.corr,difcorDF,'ro','markersize',8,'linewidth',1);
    set(gca,'xlim',[0 1],'ylim',[0.1 1.1],'ytick',[0.1 1],'yticklabel',[0.1 1],'yscale','log','fontsize',12,'tickdir','out','ticklength',[0.02 0.02],'linewidth',1);
    ylabel 'Dom. Freq. {difcor} (kHz)';
    xlabel 'Correlation Coeff. {mon,bin}';
end













































