function [] = chinchilla_panorama(optionFlag)

datapath = 'C:\Users\Mark\Dropbox\MonRev';
savedir = 'C:\Users\Mark\Dropbox\Disparity_MATLAB';
ANdatalist = load('C:\Users\Mark\Dropbox\Disparity_MATLAB\MonRevData.mat');
CF = ANdatalist.CF/1000;
DF = ANdatalist.DF;

inds = find(CF<=3 & ~isnan(DF));
DF = DF(inds);
CF = CF(inds);
CF = CF+(0.002.*randn(size(CF)));
IDs = ANdatalist.unitID(inds);
for i = 1:length(IDs)
    animalNr(i) = str2num(IDs{i}(2:6));
end

uniqueAnimals = unique(animalNr);
for i = 1:length(uniqueAnimals)
    animalIndx{i} = find(animalNr==uniqueAnimals(i));
    nrunitsperanimal(i) = length(animalIndx{i});
end


switch optionFlag
    case 'population'
        Pinds = 1:length(CF);
        F = 0.09.*(2.^(0:0.15:5));
        F = F(F<=2.6);
    case 'single'
        Pinds = animalIndx{31};
        F = CF(Pinds);
end
nrF = length(F);


[MeanPhase,Weights,FreqAx,flag] = local_get_phase(F,CF(Pinds),IDs(Pinds),optionFlag);

F = F(logical(flag));
nrF = length(F);

if strcmp('single',optionFlag)
    [F,sortind] = sort(F,'ascend');
    MeanPhase = MeanPhase(sortind);
    Weights = Weights(sortind);
end

for i = 1:nrF-1
    CWeights = Weights{i}.*Weights{i+1};
    CWeights = CWeights/max(CWeights);
    CPhase = MeanPhase{i}-MeanPhase{i+1};
    [~,maxind] = max(CWeights);
    [CPw,CP_uw,CD] = fit_circ_lin(FreqAx(CWeights>0),CPhase(CWeights>0),CWeights(CWeights>0));
    targetphase = (FreqAx(maxind)*CD)+CP_uw;
    currentphase = MeanPhase{i+1}(maxind)-MeanPhase{i}(maxind);
    MeanPhase{i+1} = MeanPhase{i+1}-currentphase;
    MeanPhase{i+1} = MeanPhase{i+1}+targetphase;
end


figure;
for i = 1:nrF
    scatter3(FreqAx,MeanPhase{i},Weights{i},...
        30*ones(size(Weights{i})),Weights{i},'fill','o','markeredgecolor','none');
    hold on;
end


[X,Y] = meshgrid(F,FreqAx);
Z = [MeanPhase{:}];
figure;pcolor(X,Y,Z);
shading flat;
hold on;
contour(X,Y,Z,-10:0.5:2,'k','linewidth',2);
set(gca,'xlim',[0.1 3],'ylim',[0.05 5],'fontsize',14,'xtick',[.1 1],'ytick',[.1 1],...
    'xticklabel',[0.1 1],'yticklabel',[0.1 1],'ticklength',[0.02 0.02],...
    'tickdir','out','layer','top','xscale','log','yscale','log');
grid off;
box off;
ylabel 'Stimulus Frequency (kHz)';
xlabel 'Characteristic Frequency (kHz)';
axis square;

Zweights = [Weights{:}];
Zweights(Zweights==0)=nan;
figure;pcolor(X,Y,Zweights)
shading flat;
hold on;
plot([0.05 3],[0.05 3],'k-','linewidth',2);
set(gca,'xlim',[0.1 3],'ylim',[0.05 5],'fontsize',14,'xtick',[.1 1],'ytick',[.1 1],...
    'xticklabel',[0.1 1],'yticklabel',[0.1 1],'ticklength',[0.02 0.02],...
    'tickdir','out','layer','top','xscale','log','yscale','log');
grid off;
box off;
ylabel 'Stimulus Frequency (kHz)';
xlabel 'Characteristic Frequency (kHz)';
axis square;



Xhat = [X(:) Y(:)];
zz = Z(:);
ww = [Weights{:}];
ww = ww(:);
[phasemodel] = fitglm(Xhat,zz,'poly99','weights',ww);

zhat = feval(phasemodel,X,Y);
zhat(isnan(Z)) = nan;
figure;pcolor(X,Y,zhat);
shading flat;
hold on;
contour(X,Y,zhat,-10:0.5:2,'k','linewidth',2);
set(gca,'xlim',[0.1 3],'ylim',[0.05 5],'fontsize',14,'xtick',[.1 1],'ytick',[.1 1],...
    'xticklabel',[0.1 1],'yticklabel',[0.1 1],'ticklength',[0.02 0.02],...
    'tickdir','out','layer','top','xscale','log','yscale','log');
grid off;
box off;
ylabel 'Stimulus Frequency (kHz)';
xlabel 'Characteristic Frequency (kHz)';
axis square;

weightsX = 0.1:0.005:5;
weightsY = 0.05:0.005:5;
[weightsX,weightsY] = meshgrid(weightsX,weightsY);
weightsZ = interp2(X,Y,[Weights{:}],weightsX,weightsY);

figure;pcolor(weightsX,weightsY,20*log10(weightsZ));
shading flat;
hold on;
set(gca,'xlim',[0.1 3],'ylim',[0.05 5],'fontsize',14,'xtick',[.1 1],'ytick',[.1 1],...
    'xticklabel',[0.1 1],'yticklabel',[0.1 1],'ticklength',[0.02 0.02],...
    'tickdir','out','layer','top','xscale','log','yscale','log','clim',[-20 0]);
grid off;
box off;
ylabel 'Stimulus Frequency (kHz)';
xlabel 'Characteristic Frequency (kHz)';
axis square;

save(fullfile(savedir,'chinchilla_panorama.mat'),'X','Y','Z','zhat','phasemodel',...
    'weightsX','weightsY','weightsZ','FreqAx','F');




    function [P,W,Freq,flag] = local_get_phase(F,uCF,uIDs,optionflag)
        switch optionflag
            case 'population'
                totalunits = 0;
                count = 0;
                nF = length(F);
                for ii = 1:nF
                    loF = F(ii)*2^-0.075;
                    hiF = F(ii)*2^0.075;
                    idx = find(uCF>=loF & uCF<hiF);
                    nrunits = length(idx);
                    totalunits = totalunits+nrunits;
                    if nrunits>0
                        flag(ii) = 1;
                        count = count+1;
                        PhaseMatrix = nan(4096,nrunits);
                        WeightMatrix = zeros(4096,nrunits);
                        for jj = 1:nrunits
                            d = load(fullfile(datapath,uIDs{idx(jj)},[uIDs{idx(jj)} '_MonRev_70.00dB.mat']));
                            d = d.kernels;
                            PhaseMatrix(:,jj) = d.h1phase_uw;
                            WeightMatrix(:,jj) = d.h1coeffs.Monaural.weights;
                        end
                        W{count} = mean(WeightMatrix,2);
                        W{count} = W{count}/max(W{count});
                        P{count} = nanmean(PhaseMatrix,2);
                    else
                        flag(ii) = 0;
                    end
                end
                Freq = d.h1ffax;
            case 'single'
                nF = length(uCF);
                for ii = 1:nF
                    d = load(fullfile(datapath,uIDs{ii},[uIDs{ii} '_MonRev_70.00dB.mat']));
                    d = d.kernels;
                    W{ii} = d.h1coeffs.Monaural.weights;
                    W{ii} = W{ii}/max(W{ii});
                    P{ii} = d.h1phase_uw;
                end
                Freq = d.h1ffax;
                flag = ones(size(P));
        end
    end
end

















