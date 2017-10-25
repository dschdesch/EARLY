clc;clear;

MSOdatapath = 'C:\Users\Mark\Dropbox\BinRev';

MSOfolderlist = dir(fullfile(MSOdatapath,'*_*'));
nrMSOunits = length(MSOfolderlist);

AN = load('C:\Users\Mark\Dropbox\Disparity_MATLAB\Predicted_CP_CD_vs_deltaCF.mat');
count = 0;
for i = 1:nrMSOunits
    MSOunitID = MSOfolderlist(i).name;
    if ~isempty(dir(fullfile(MSOdatapath,MSOunitID,'*BinRev_70.00dB.mat')))
        count = count+1;
        load(fullfile(MSOdatapath,MSOunitID,[MSOunitID '_BinRev_70.00dB.mat']));
        load(fullfile(MSOdatapath,MSOunitID,[MSOunitID '_NDF_70.00dB.mat']));
        TD(count) = kernels.tuning_diff.ipsivscontra_oct;
        CF_contra(count) = kernels.h1cdomfreq(2);
        DifcorDF(count) = ndf.difcor.peakhz/1000;
        BD_ms(count) = ndf.bd.bd;
        BD_cycles(count) = BD_ms(count)*DifcorDF(count);
        totalCD_ms(count) = kernels.h1coeffs.Binaural.CD;
        totalCP_cycles(count) = kernels.h1coeffs.Binaural.CP;
        totalCD_cycles(count) = totalCD_ms(count)*DifcorDF(count);
        totalCP_ms(count) = totalCP_cycles(count)/DifcorDF(count);
    end
end

nrunits = count;

figure;pcolor(AN.CONTRA_cf,AN.delta_CF,AN.CP');hold on;
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
plot3(CF_contra,TD,ones(size(TD)),'ko');


figure;pcolor(AN.CONTRA_cf,AN.delta_CF,AN.CD');hold on;
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
plot3(CF_contra,TD,ones(size(TD)),'ko');


cochlearCP_cycles = interp2(AN.CONTRA_cf,AN.delta_CF,AN.CP',CF_contra,TD);
cochlearCP_ms = cochlearCP_cycles./DifcorDF;
cochlearCD_ms = interp2(AN.CONTRA_cf,AN.delta_CF,AN.CD',CF_contra,TD);
cochlearCD_cycles = cochlearCD_ms.*DifcorDF;

brainCD_ms = totalCD_ms-cochlearCD_ms;
brainCD_cycles = brainCD_ms.*DifcorDF;
brainCP_cycles = totalCP_cycles-cochlearCP_cycles;
brainCP_ms = brainCP_cycles./DifcorDF;


%Make a set of different models
X = [-1.5 1.5];
Y = X;
figure;
subplot(4,2,1);
plot(BD_cycles,DifcorDF.*(totalCD_ms+(totalCP_cycles./DifcorDF)),'o');
COD = Coeff_Determination(BD_cycles,DifcorDF.*(totalCD_ms+(totalCP_cycles./DifcorDF)));
hold on;
plot(X,Y,'r-');
title 'Total CD + Total CP';

subplot(4,2,2);
plot(BD_cycles,DifcorDF.*brainCD_ms,'o');
COD = Coeff_Determination(BD_cycles,DifcorDF.*brainCD_ms);
hold on;
plot(X,Y,'r-');
title 'Brain CD';

subplot(4,2,3);
plot(BD_cycles,brainCP_cycles,'o');
COD = Coeff_Determination(BD_cycles,brainCP_cycles);
hold on;
plot(X,Y,'r-');
title 'Brain CP';

subplot(4,2,4);
plot(BD_cycles,DifcorDF.*(brainCD_ms+(brainCP_cycles./DifcorDF)),'o');
COD = Coeff_Determination(BD_cycles,DifcorDF.*(brainCD_ms+(brainCP_cycles./DifcorDF)));
hold on;
plot(X,Y,'r-');
title 'Brain CD + Brain CP';

subplot(4,2,5);
plot(BD_cycles,DifcorDF.*cochlearCD_ms,'o');
COD = Coeff_Determination(BD_cycles,DifcorDF.*cochlearCD_ms);
hold on;
plot(X,Y,'r-');
title 'Cochlear CD';

subplot(4,2,6);
plot(BD_cycles,cochlearCP_cycles,'o');
COD = Coeff_Determination(BD_cycles,cochlearCP_cycles);
hold on;
plot(X,Y,'r-');
title 'Cochlear CP';

subplot(4,2,7);
plot(BD_cycles,DifcorDF.*(cochlearCD_ms+(cochlearCP_cycles./DifcorDF)),'o');
COD = Coeff_Determination(BD_cycles,DifcorDF.*(cochlearCD_ms+(cochlearCP_cycles./DifcorDF)));
hold on;
plot(X,Y,'r-');
title 'Cochlear CD + Cochlear CP';

figure;
plot(cochlearCD_cycles,brainCD_cycles,'o');

figure;
plot(wrapToPi(cochlearCP_cycles*2*pi)/(2*pi),wrapToPi(brainCP_cycles*2*pi)/(2*pi),'o');

%Fit a linear mixed-effects model by maximum-likelihood
y = BD_cycles';
X = [(DifcorDF.*brainCD_ms)' brainCP_cycles' (DifcorDF.*cochlearCD_ms)' cochlearCP_cycles'];

lm_full = fitlm(X,y);%Full 4-parameter model
% lm_bD_bP_cD = fitlm(X(:,[1 2 3]),y);%leave out cochlear phase
% lm_bD_bP_cP = fitlm(X(:,[1 2 4]),y);%leave out cochlear delay
% lm_bD_cD_cP = fitlm(X(:,[1 3 4]),y);%leave out brain phase
% lm_bP_cD_cP = fitlm(X(:,[2 3 4]),y);%leave out brain delay

Xbrain = sum(X(:,[1 2]),2);
Xcochlea = sum(X(:,[3 4]),2);
lm_brain_cochlea = fitlm([Xbrain Xcochlea],y);%Two-parameter (brain and cochlea) model
% lm_brain = fitlm(X(:,[1 2]),y);%leave out cochlea
% lm_cochlea = fitlm(X(:,[3 4]),y);%leave out brain

for i = 1:100
    indxA = randi([1 nrunits],nrunits,1);
    indxB = randi([1 nrunits],nrunits,1);
    indxC = randi([1 nrunits],nrunits,1);
    lm = fitlm([Xbrain(indxA) Xcochlea(indxB)],y(indxC));
end





