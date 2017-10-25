clear;clc;

load('C:\Users\Mark\Dropbox\Disparity_MATLAB\chinchilla_panorama.mat');
savedir = 'C:\Users\Mark\Dropbox\Disparity_MATLAB\';
delta_CF = -0.6:0.05:0.6;
N = length(delta_CF);
CONTRA_cf = 0.1.*(2.^(0:0.1:5));
CONTRA_cf = CONTRA_cf(CONTRA_cf<=2.5);
CONTRA_df = nan(size(CONTRA_cf));

IPSI_df = nan(length(CONTRA_cf),N);
CP = nan(length(CONTRA_cf),N);
CD = nan(length(CONTRA_cf),N);
difcorDF = nan(length(CONTRA_cf),N);
DF_tuningdiff_oct = nan(length(CONTRA_cf),N);

for j = 1:length(CONTRA_cf)
    IPSI_cf = CONTRA_cf(j).*(2.^delta_CF);
    
    CONTRA_phase = feval(phasemodel,CONTRA_cf(j)*ones(4096,1),FreqAx');
    CONTRA_weight = interp2(weightsX,weightsY,weightsZ,CONTRA_cf(j)*ones(4096,1),FreqAx');
    CONTRA_weight(isnan(CONTRA_weight)) = 0;
    [~,maxind] = max(CONTRA_weight);
    CONTRA_df(j) = FreqAx(maxind);
    
    for i = 1:N
        IPSI_phase = feval(phasemodel,IPSI_cf(i)*ones(4096,1),FreqAx');
        IPSI_weight = interp2(weightsX,weightsY,weightsZ,IPSI_cf(i)*ones(4096,1),FreqAx');
        IPSI_weight(isnan(IPSI_weight)) = 0;
        
        B_weight = CONTRA_weight.*IPSI_weight;
        B_weight = B_weight/max(B_weight);
        
        if ~all(isnan(B_weight))
            B_phase = CONTRA_phase-IPSI_phase;
%             [~,CP(j,i),CD(j,i)] = fit_circ_lin(FreqAx(B_weight>0),B_phase(B_weight>0),B_weight(B_weight>0));
            [~,maxind] = max(B_weight);
            difcorDF(j,i) = FreqAx(maxind);
            
            [~,maxind] = max(IPSI_weight);
            IPSI_df(j,i) = FreqAx(maxind);
            
            A.h1ffax = FreqAx;
            A.h1cmag = 20*log10([IPSI_weight CONTRA_weight]);
            B.difcor.power = 20*log10(B_weight);
            [~,maxind] = max(B_weight);
            B.difcor.peakhz = FreqAx(maxind)*1000;
            [tuning_diff] = get_tuning_difference(A,B);
            DF_tuningdiff_oct(j,i) = tuning_diff.tuning_diff.ipsivscontra_oct;
        end
    end
end

figure;pcolor(CONTRA_cf,delta_CF,CP');hold on;
set(gca,'xscale','log','clim',[-0.5 0.5],'fontsize',14,'layer','top','xtick',[0.1 1],...
    'xticklabel',[0.1 1],'xlim',[0.1 2]);
shading flat;
axis square;
box off;
ylabel ('\Delta CF (octaves)','interpreter','tex');
xlabel 'CONTRA Characteristic Frequency (kHz)';
cbh = colorbar;
set(cbh,'fontsize',12);
ylabel(cbh,'Characteristic Phase (cycles)');
colormap(flipud(jet));

figure;pcolor(CONTRA_cf,delta_CF,CD');hold on;
set(gca,'xscale','log','fontsize',14,'layer','top','xtick',[0.1 1],...
    'xticklabel',[0.1 1],'xlim',[0.1 2],'clim',[-1 1]);
shading flat;
axis square;
box off;
ylabel ('\Delta CF (octaves)','interpreter','tex');
xlabel 'CONTRA Characteristic Frequency (kHz)';
cbh = colorbar;
set(cbh,'fontsize',12);
ylabel(cbh,'Characteristic Delay (ms)');
colormap(flipud(jet));

figure;pcolor(CONTRA_cf,delta_CF,(CP./difcorDF)');hold on;
set(gca,'xscale','log','clim',[-1 1],'fontsize',14,'layer','top','xtick',[0.1 1],...
    'xticklabel',[0.1 1],'xlim',[0.1 2]);
shading flat;
axis square;
box off;
ylabel ('\Delta CF (octaves)','interpreter','tex');
xlabel 'CONTRA Characteristic Frequency (kHz)';
cbh = colorbar;
set(cbh,'fontsize',12);
ylabel(cbh,'Characteristic Phase (ms)');
colormap(flipud(jet));

figure;pcolor(CONTRA_cf,delta_CF,(CD.*difcorDF)');hold on;
set(gca,'xscale','log','fontsize',14,'layer','top','xtick',[0.1 1],...
    'xticklabel',[0.1 1],'xlim',[0.1 2],'clim',[-1.7 1.7]);
shading flat;
axis square;
box off;
ylabel ('\Delta CF (octaves)','interpreter','tex');
xlabel 'CONTRA Characteristic Frequency (kHz)';
cbh = colorbar;
set(cbh,'fontsize',12);
ylabel(cbh,'Characteristic Delay (cycles)');
colormap(flipud(jet));


save(fullfile(savedir,'Predicted_CP_CD_vs_deltaCF.mat'),'CONTRA_cf','IPSI_cf','delta_CF','CD','CP','difcorDF');










