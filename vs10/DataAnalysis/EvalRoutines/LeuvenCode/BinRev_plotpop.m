function [] = BinRev_plotpop
warning off all
anadir = 'C:\LeuvenDataAnalysis';
datadir = 'C:\work\BinRev\BinRev';
cd (datadir);
datafiles = dir('*_BinRev.mat');
nrdata = numel(datafiles);
ANdata = load('C:\LeuvenDataAnalysis\MonRevPopData.mat');

ANdfs = zeros(length(ANdata.AN),1);
for i = 1:length(ANdata.AN)
    ind = find(unique(ANdata.AN(i).noiseSPLs)==70,1);
    if ~isempty(ind)
        ANdfs(i) = ANdata.AN(i).df(ind);
    else
        ANdfs(i) = NaN;
    end
end
ANcfs = vertcat(ANdata.AN(:).cf);
ANfreqs = ANdfs;
ANfreqs(ANcfs>2000) = ANcfs(ANcfs>2000)/1000;
for i =1:nrdata
    cd (datadir);
    data = load (datafiles(i).name);
    ParamsOut = data.ParamsOut;clear data;
    if strcmp(ParamsOut.Peak_Trough,'P')
        PeakTrough(i) = 1;
    else
        PeakTrough(i) = 2;
    end
    if isfield(ParamsOut,'NDF')
        if isfield(ParamsOut.NDF,'NTDyPos_spline') && isfield(ParamsOut.NDF,'NTDyNeg_spline')
            NDFx{i} = ParamsOut.NDF.NTDx_spline;
            NDFy{i} = [ParamsOut.NDF.NTDyPos_spline;ParamsOut.NDF.NTDyNeg_spline];
            NDFy{i} = ((NDFy{i}./max(abs(NDFy{i}(:))))*0.5)-0.25;
            NDFz{i} = ParamsOut.NDF.Difcor_spline;
            
            NDFz{i} = ((NDFz{i}./max(abs(NDFz{i})))*0.5)-0.25;
            
            NDFxvals{i} = ParamsOut.NDF.NTDx;
            NDFyvals{i} = ParamsOut.NDF.Difcor_NTD;
            
            NDFtype(i) = 1;
            BD(i) = ParamsOut.NDF.BD_NTD;
            WD(i) = ParamsOut.NDF.WD_NTD;
            CD_NDF(i) = ParamsOut.NDF.BD_NTD_env;
            CP_NDF(i) = (BD(i)-CD_NDF(i))*ParamsOut.NDF.DifcorDomFreq_NTD;
            CP_NDF(i) = wrapToPi(CP_NDF(i)*2*pi)/(2*pi);
            
            DF(i) = ParamsOut.NDF.DifcorDomFreq_NTD;
            DifFreqAx{i} = ParamsOut.NDF.Difcor_Mag_Freq;
            DifS{i} = ParamsOut.NDF.Difcor_Mag_NTD;
            
        elseif isfield(ParamsOut,'TDF')
            
            NDFxvals{i} = ParamsOut.TDF.x.ITD;
            NDFyvals{i} = ParamsOut.TDF.Tone_Difcor;
            
            NDFx{i} = ParamsOut.TDF.x.ITD_spline;
            NDFy{i} = [ParamsOut.TDF.y.pITD_spline;ParamsOut.TDF.y.nITD_spline];
            NDFy{i} = ((NDFy{i}./max(abs(NDFy{i}(:))))*0.5)-0.25;
            NDFz{i} = ParamsOut.TDF.Tone_Difcor_Spline;
            
            NDFz{i} = ((NDFz{i}./max(abs(NDFz{i})))*0.5)-0.25;
            NDFtype(i) = 2;
            BD(i) = ParamsOut.TDF.BD;
            
            WD(i) = ParamsOut.TDF.WD;
            
            CD_NDF(i) = ParamsOut.TDF.BD_env;
            CP_NDF(i) = (BD(i)-CD_NDF(i))*ParamsOut.TDF.DFTonedifcor;
            CP_NDF(i) = wrapToPi(CP_NDF(i)*2*pi)/(2*pi);
            
            DF(i) = ParamsOut.TDF.DFTonedifcor;
            DifS{i} = ParamsOut.TDF.Tone_DifMagSpec;
            
        end
    elseif isfield(ParamsOut,'TDF')
        NDFxvals{i} = ParamsOut.TDF.x.ITD;
        NDFyvals{i} = ParamsOut.TDF.Tone_Difcor;
        
        NDFx{i} = ParamsOut.TDF.x.ITD_spline;
        NDFy{i} = [ParamsOut.TDF.y.pITD_spline;ParamsOut.TDF.y.nITD_spline];
        NDFy{i} = ((NDFy{i}./max(abs(NDFy{i}(:))))*0.5)-0.25;
        NDFz{i} = ParamsOut.TDF.Tone_Difcor_Spline;
       
        NDFz{i} = ((NDFz{i}./max(abs(NDFz{i})))*0.5)-0.25;
        NDFtype(i) = 2;
        BD(i) = ParamsOut.TDF.BD;
        WD(i) = ParamsOut.TDF.WD;
        CD_NDF(i) = ParamsOut.TDF.BD_env;
        CP_NDF(i) = (BD(i)-CD_NDF(i))*ParamsOut.TDF.DFTonedifcor;
        CP_NDF(i) = wrapToPi(CP_NDF(i)*2*pi)/(2*pi);
        
        DF(i) = ParamsOut.TDF.DFTonedifcor;
        DifS{i} = ParamsOut.TDF.Tone_DifMagSpec;
        
    end
    if isfield(ParamsOut,'TDF')
        CP(i) = ParamsOut.TDF.CP;
        CD(i) = ParamsOut.TDF.CD;
        CP_w(i) = ParamsOut.TDF.CP_w;
    else
        CP(i) = NaN;
        CD(i) = NaN;
        CP_w(i) = NaN;
    end
    pi_limit(i) = 0.5/DF(i);
    if isfield(ParamsOut,'Revcor')
        if isfield(ParamsOut,'GCFit')
            for j = 1:2
                F0(i,j) = abs(ParamsOut.GCFit.x{j}(3));
                c(i,j) = ParamsOut.GCFit.x{j}(4);
                tau(i,j) = ParamsOut.GCFit.x{j}(2);
            end
        else
            F0(i,:) = [NaN NaN];
            c(i,:) = [NaN NaN];
            tau(i,:) = [NaN NaN];
        end
        h1{i} = ParamsOut.Revcor.h1filt;
        Latency(i,:) = ParamsOut.Revcor.h1latency;
        Delta.y(i) = 2*ParamsOut.Revcor.h1magxcorpeakdoct;
        Delta.z(i) = ParamsOut.Revcor.domfreq(1)-ParamsOut.Revcor.domfreq(2);
        Delta.zoct(i) = log2(ParamsOut.Revcor.domfreq(1)/ParamsOut.Revcor.domfreq(2));
        Delta.x(i) = ParamsOut.Revcor.h1meanDomFreq;
        RevcorDF(i,:) = ParamsOut.Revcor.domfreq;
        DF_P(i) = ParamsOut.Revcor.DifcorDomFreq_Pred;
        if isfield(ParamsOut.Revcor,'ANFits')
            ANinds(i,:) = ParamsOut.Revcor.ANFits.ANinds;
            delaysamps(i,:) = ParamsOut.Revcor.ANFits.delaysamps;
            ANfilename{i} = ParamsOut.Revcor.ANFits.ANfilename;
            Q_3_AN(i,:) = ParamsOut.Revcor.ANFits.Q_3_AN;
            Q_6_AN(i,:) = ParamsOut.Revcor.ANFits.Q_6_AN;
            AddTime(i,:) = ParamsOut.Revcor.ANFits.AddTime;
            ANBD(i) = ParamsOut.Revcor.ANFits.ANBD;
            ANCP(i) = ParamsOut.Revcor.ANFits.ANCP;
            ANCD(i) = ParamsOut.Revcor.ANFits.ANCD;
            ANDF(i) = ParamsOut.Revcor.ANFits.ANDF;
        else
            %get the closest matching AN data and do the CD CP analyses on
            %those reponses
            cd (anadir);
            [ANinds(i,:),delaysamps(i,:),ANfilename{i}] = determine_AN_fits(h1{i},RevcorDF(i,:));
            [Q_3_AN(i,:),Q_6_AN(i,:),AddTime(i,:),ANBD(i),ANCP(i),ANCD(i),ANDF(i)] = predictCPCD(ANinds(i,:),delaysamps(i,:),ParamsOut);
            cd (datadir);
            ParamsOut.Revcor.ANFits.ANinds = ANinds(i,:);
            ParamsOut.Revcor.ANFits.delaysamps = delaysamps(i,:);
            ParamsOut.Revcor.ANFits.ANfilename = ANfilename{i};
            ParamsOut.Revcor.ANFits.Q_3_AN = Q_3_AN(i,:);
            ParamsOut.Revcor.ANFits.Q_6_AN = Q_6_AN(i,:);
            ParamsOut.Revcor.ANFits.AddTime = AddTime(i,:);
            ParamsOut.Revcor.ANFits.ANBD = ANBD(i);
            ParamsOut.Revcor.ANFits.ANCP = ANCP(i);
            ParamsOut.Revcor.ANFits.ANCD = ANCD(i);
            ParamsOut.Revcor.ANFits.ANDF = ANDF(i);
            %Update the saved data structure
            save(datafiles(i).name,'ParamsOut');
        end
        bbCP(i) = ParamsOut.Revcor.BPhaseFits.bbCP;
        bbCD(i) = ParamsOut.Revcor.BPhaseFits.bbCD;
%         bbCP_w(i) = ParamsOut.Revcor.bbCP_w;
        bP(i) = ParamsOut.Revcor.Best_Phase_ms_pred;
        BD_P(i) = ParamsOut.Revcor.BD_pred;
        BDest(i) = ParamsOut.Revcor.BDest;
        binauralBW(i) = max(ParamsOut.Revcor.h1ffax(~isnan(ParamsOut.Revcor.h1phasemask(:,1)) & ~isnan(ParamsOut.Revcor.h1phasemask(:,2))))-...
            min(ParamsOut.Revcor.h1ffax(~isnan(ParamsOut.Revcor.h1phasemask(:,1)) & ~isnan(ParamsOut.Revcor.h1phasemask(:,2))));
        if abs(BD(i))<=pi_limit(i)
            within_pi(i) = true;
        else
            within_pi(i) = false;
        end
        if abs(BD(i))<=pi_limit(i)/2
            within_pi_half(i) = true;
        else
            within_pi_half(i) = false;
        end
        Rho(i) = ParamsOut.Revcor.Difcor_xcor;
        Q_3(i,:) = ParamsOut.Revcor.qvals(1,:);
        Q_6(i,:) = ParamsOut.Revcor.qvals(2,:);
        
        FrontLatency(i,:) = ParamsOut.Revcor.h1frontlatency;
        GD(i,:) = -ParamsOut.Revcor.h1phasefit_coeffs(:,1)';
        
        minFreq(i) = min(ParamsOut.Revcor.h1ffax(~isnan(ParamsOut.Revcor.h1phasemask(:,1)) & ~isnan(ParamsOut.Revcor.h1phasemask(:,2))));
        maxFreq(i) = max(ParamsOut.Revcor.h1ffax(~isnan(ParamsOut.Revcor.h1phasemask(:,1)) & ~isnan(ParamsOut.Revcor.h1phasemask(:,2))));
        
        K(i) = ParamsOut.Revcor.k;
    else
        DF_P(i) = NaN;
        ANCP(i) = NaN;
        ANCD(i) = NaN;
        ANDF(i) = NaN;
        ANBD(i) = NaN;
        Q_3_AN(i,:) = [NaN NaN];
        Q_6_AN(i,:) = [NaN NaN];
        AddTime(i,:) = [NaN NaN];
        ANinds(i,:) = [NaN NaN];
        delaysamps(i,:) = [NaN NaN];
        ANfilename{i} = cell(1,2);
        F0(i,:) = [NaN NaN];
        c(i,:) = [NaN NaN];
        tau(i,:) = [NaN NaN];
        h1{i} = NaN;
        binauralBW(i) = NaN;
        Delta.x(i) = NaN;
        Delta.y(i) = NaN;
        Delta.z(i) = NaN;
        Delta.zoct(i) = NaN;
        RevcorDF(i,:) = [NaN NaN];
        bbCP(i) = NaN;
        bbCD(i) = NaN;
        bbCP_w(i) = NaN;
        bP(i) = NaN;
        BD_P(i) = NaN;
        BDest(i) = NaN;
        Rho(i) = NaN;
        Q_3(i,:) = NaN;
        Q_6(i,:) = NaN;
        within_pi(i) = false;
        within_pi_half(i) = false;
        Latency(i,:)=  [NaN NaN];
        FrontLatency(i,:) = [NaN NaN];
        GD(i,:) = [NaN NaN];

        minFreq(i) = NaN;
        maxFreq(i) = NaN;
        
        K(i) = NaN;
    end
    if isfield(ParamsOut,'NRHO')
        if isfield(ParamsOut.NRHO,'beta')
            nrhoP(i) = ParamsOut.NRHO.beta(3);
        else
            nrhoP(i) = NaN;
        end
        nrhoPpredicted(i) = ParamsOut.NRHO.predicted_beta(3);
    else
        nrhoP(i) = NaN;
        nrhoPpredicted(i) = NaN;
    end
end
cd (anadir);

%% Plot BD against difference in tuning
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
axes('ActivePositionProperty','outerposition');
xdata = Delta.y(within_pi & PeakTrough==1);%Difference in DF expressed as octaves (positive is ipsi-higher). Throw out LSO type and only look at cells in the pi limit.
ydata = BD(within_pi & PeakTrough==1);%BD - Throw out LSO type and only look at cells in the pi limit.

% xdata = Delta.y(~isnan(Delta.y) & PeakTrough==1);%Difference in DF expressed as octaves (positive is ipsi-higher). Throw out LSO type and only look at cells in the pi limit.
% ydata = BD(~isnan(Delta.y) & PeakTrough==1);%BD - Throw out LSO type and only look at cells in the pi limit.

% xdata2 = Delta.y(PeakTrough==2);
% ydata2 = WD(PeakTrough==2);
% 
% xdata = [xdata xdata2(~isnan(xdata2))];
% ydata = [ydata ydata2(~isnan(xdata2))];
% ind = find(xdata<0.45 | ydata>0);
% 
% xdata = xdata(ind);
% ydata = ydata(ind);

plot(xdata',ydata','k*','markersize',10);hold on;
rho = corr(xdata',ydata');
[beta,gof]=local_fit_linear(xdata',ydata');
[coeffpval,coefftval] = local_stats(xdata',ydata',beta,gof.dfe);
plot([-1 1],feval(beta,[-1 1]),'r--','linewidth',2);
set(gca,'fontsize',16,'ticklength',[0.02 0.02],'linewidth',1,'xlim',[-1.2 1.2]);
xlabel ({'\Delta DOMINANT', 'FREQUENCY (octaves)'},'interpreter','tex');
ylabel 'BEST DELAY (ms)';
box off;
axis square;
if coeffpval<0.001
    text(-0.8,1.2,sprintf('\\rho: %2.2f \n\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: <0.001',rho,beta.a,gof.dfe,coefftval),'interpreter','tex')
else
    text(-0.8,1.2,sprintf('\\rho: %2.2f \n\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',rho,beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex')
end


%% Plot FD, S-FD, CD, and CP against Delta-tuning
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.3244    0.4189    0.2819    0.3689]);
axes('ActivePositionProperty','outerposition');
DeltaFrontLatency = FrontLatency(:,2)-FrontLatency(:,1);
DeltaFrontLatency = DeltaFrontLatency';
xdata = Delta.y(~isnan(DeltaFrontLatency) & within_pi);
ydata = DeltaFrontLatency(~isnan(DeltaFrontLatency) & within_pi);
plot(xdata,ydata,'k*','linewidth',1,'markersize',10);hold on;
[beta,gof]=local_fit_linear(xdata',ydata');
[coeffpval,coefftval] = local_stats(xdata',ydata',beta,gof.dfe);
plot([min(xdata) max(xdata)],feval(beta,[min(xdata) max(xdata)]),'r--','linewidth',2);
%add some text describing the regression coefficient and stats
if coeffpval<0.001
    text(-0.8,4,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: <0.001',beta.a,gof.dfe,coefftval),'interpreter','tex')
else
    text(-0.8,4,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex')
end
set(gca,'xlim',[-1.2 1.2],'linewidth',1,'fontsize',16,'ticklength',[0.02 0.02]);
xlabel ({'\Delta DOMINANT', 'FREQUENCY (octaves)'},'interpreter','tex');
ylabel ({'\Delta SIGNAL-FRONT' ,'DELAY (ms)'},'interpreter','tex');
axis square;
box off;



figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.3244    0.4189    0.2819    0.3689]);
axes('ActivePositionProperty','outerposition');
DeltaGDC = (GD(:,2)-FrontLatency(:,2))-(GD(:,1)-FrontLatency(:,1));
ydata = DeltaGDC(~isnan(DeltaGDC) & within_pi');
xdata =  Delta.y(~isnan(DeltaGDC) & within_pi');
plot(xdata,ydata,'k*','linewidth',1,'markersize',12);hold on;
[beta,gof]=local_fit_linear(xdata',ydata);
[coeffpval,coefftval] = local_stats(xdata',ydata,beta,gof.dfe);
plot([min(xdata) max(xdata)],feval(beta,[min(xdata) max(xdata)]),'r--','linewidth',2);
%add some text describing the regression coefficient and stats
if coeffpval<0.001
    text(-0.8,1.5,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: <0.001',beta.a,gof.dfe,coefftval),'interpreter','tex')
else
    text(-0.8,1.5,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex')
end
set(gca,'xlim',[-1 1],'linewidth',1,'fontsize',16,'ticklength',[0.02 0.02]);
xlabel ({'\Delta DOMINANT', 'FREQUENCY (octaves)'},'interpreter','tex');
ylabel ('\Delta FILTER DELAY (ms)','interpreter','tex');
axis square;
box off;

figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.3244    0.4189    0.2819    0.3689]);
axes('ActivePositionProperty','outerposition');
DeltaGD = GD(:,2)-GD(:,1);
ydata = DeltaGD(~isnan(DeltaGD) & within_pi');
xdata =  Delta.y(~isnan(DeltaGD) & within_pi');
plot(xdata,ydata,'k*','linewidth',1,'markersize',12);hold on;
[beta,gof]=local_fit_linear(xdata',ydata);
[coeffpval,coefftval] = local_stats(xdata',ydata,beta,gof.dfe);
plot([min(xdata) max(xdata)],feval(beta,[min(xdata) max(xdata)]),'r--','linewidth',2);
%add some text describing the regression coefficient and stats
if coeffpval<0.001
    text(-0.8,1.5,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: <0.001',beta.a,gof.dfe,coefftval),'interpreter','tex')
else
    text(-0.8,1.5,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex')
end
set(gca,'xlim',[-1 1],'linewidth',1,'fontsize',16,'ticklength',[0.02 0.02]);
xlabel ({'\Delta DOMINANT', 'FREQUENCY (octaves)'},'interpreter','tex');
ylabel ({'CHARACTERISTIC' ,'DELAY (ms)'},'interpreter','tex');
axis square;
box off;

figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.3244    0.4189    0.2819    0.3689]);
axes('ActivePositionProperty','outerposition');
ydata = bbCP(~isnan(bbCP) & within_pi)./DF(~isnan(bbCP) & within_pi);
xdata =  Delta.y(~isnan(bbCP) & within_pi);
plot(xdata,ydata,'k*','linewidth',1,'markersize',12);hold on;
[beta,gof]=local_fit_linear(xdata',ydata');
[coeffpval,coefftval] = local_stats(xdata',ydata',beta,gof.dfe);
plot([min(xdata) max(xdata)],feval(beta,[min(xdata) max(xdata)]),'r--','linewidth',2);
%add some text describing the regression coefficient and stats
if coeffpval<0.001
    text(-0.8,1.5,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: <0.001',beta.a,gof.dfe,coefftval),'interpreter','tex')
else
    text(-0.8,1.5,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex')
end
set(gca,'xlim',[-1 1],'linewidth',1,'fontsize',16,'ticklength',[0.02 0.02]);
xlabel ({'\Delta DOMINANT', 'FREQUENCY (octaves)'},'interpreter','tex');
ylabel ({'CHARACTERISTIC' ,'PHASE (ms)'},'interpreter','tex');
axis square;
box off;

%% Plot signal-front delay versus filter delay for the same units
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.3244    0.4189    0.2819    0.3689]);
axes('ActivePositionProperty','outerposition');
DeltaFrontLatency = FrontLatency(:,2)-FrontLatency(:,1);
DeltaGDC = (GD(:,2)-FrontLatency(:,2))-(GD(:,1)-FrontLatency(:,1));
xdata = DeltaGDC(~isnan(DeltaGDC) & within_pi');
ydata = DeltaFrontLatency(~isnan(DeltaFrontLatency) & within_pi');
plot(xdata,ydata,'k*','linewidth',1,'markersize',10);hold on;
[beta,gof]=local_fit_linear(xdata,ydata);
[coeffpval,coefftval] = local_stats(xdata,ydata,beta,gof.dfe);
plot([min(xdata) max(xdata)],feval(beta,[min(xdata) max(xdata)]),'r--','linewidth',2);
%add some text describing the regression coefficient and stats
if coeffpval<0.001
    text(0,5,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: <0.001',beta.a,gof.dfe,coefftval),'interpreter','tex')
else
    text(0,5,sprintf('\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex')
end
set(gca,'xlim',[-3 2],'linewidth',1,'fontsize',16,'ticklength',[0.02 0.02]);
xlabel ('\Delta FILTER DELAY (ms)','interpreter','tex');
ylabel ({'\Delta SIGNAL-FRONT' ,'DELAY (ms)'},'interpreter','tex');
axis square;
box off;


%% Plot BD (measured with NTD) against the model prediction based on the CD
%% and CP from the binaural revcors
figure;set(gcf,'paperpositionmode','auto');
axes('ActivePositionProperty','outerposition');
plot([0 0],[-3 3],'--',[-3 3],[0 0],'--',[-3 3],[-3 3],'color',[.6 .6 .6],'linewidth',1);hold on;

%all factors in the model: FD, S-FD, and CP;
plot(BD,BD_P,'ko','markersize',10,'linewidth',2);hold on;
rho = corr(BD(~isnan(BD_P))',BD_P(~isnan(BD_P))');
[beta,gof]=local_fit_linear(BD(~isnan(BD_P))',BD_P(~isnan(BD_P))');
[coeffpval,coefftval] = local_stats(BD(~isnan(BD_P))',BD_P(~isnan(BD_P))',beta,gof.dfe);
plot([min(BD) max(BD)],feval(beta,[min(BD) max(BD)]),'k--','linewidth',2);
if coeffpval<0.001
    text(-2.5,2.2,sprintf('\\rho: %2.2f \n\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',rho,beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex','color','k');
else
    text(-2.5,2.2,sprintf('\\rho: %2.2f \n\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',rho,beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex','color','k');
end


% %only FD and CP included;
% FD_CP_pred = ((GD(:,2)-FrontLatency(:,2))-(GD(:,1)-FrontLatency(:,1)))+(bbCP./DF)';
% plot(BD,FD_CP_pred,'r^','markersize',10,'linewidth',2);hold on;
% rho = corr(BD(~isnan(FD_CP_pred))',FD_CP_pred(~isnan(FD_CP_pred)));
% [beta,gof]=local_fit_linear(BD(~isnan(FD_CP_pred))',FD_CP_pred(~isnan(FD_CP_pred)));
% [coeffpval,coefftval] = local_stats(BD(~isnan(FD_CP_pred))',FD_CP_pred(~isnan(FD_CP_pred)),beta,gof.dfe);
% plot([min(BD) max(BD)],feval(beta,[min(BD) max(BD)]),'r--','linewidth',2);
% if coeffpval<0.001
%     text(1.8,-0.7,sprintf('\\rho: %2.2f \n\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: <0.001',rho,beta.a,gof.dfe,coefftval),'interpreter','tex','color','r');
% else
%     text(1.8,-0.7,sprintf('\\rho: %2.2f \n\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',rho,beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex','color','r');
% end
% 
% %only S-FD and CP included;
% SFD_CP_pred = (FrontLatency(:,2)-FrontLatency(:,1))+(bbCP./DF)';
% plot(BD,SFD_CP_pred,'bv','markersize',10,'linewidth',2);hold on;
% rho = corr(BD(~isnan(SFD_CP_pred))',SFD_CP_pred(~isnan(FD_CP_pred)));
% [beta,gof]=local_fit_linear(BD(~isnan(SFD_CP_pred)),SFD_CP_pred(~isnan(SFD_CP_pred)));
% [coeffpval,coefftval] = local_stats(BD(~isnan(SFD_CP_pred)),SFD_CP_pred(~isnan(SFD_CP_pred)),beta,gof.dfe);
% plot([min(BD) max(BD)],feval(beta,[min(BD) max(BD)]),'b--','linewidth',2);
% if coeffpval<0.001
%     text(1.8,-2.1,sprintf('\\rho: %2.2f \n\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: <0.001',rho,beta.a,gof.dfe,coefftval),'interpreter','tex','color','b');
% else
%     text(1.8,-2.1,sprintf('\\rho: %2.2f \n\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',rho,beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex','color','b');
% end
% 
% %only CP included;
% CP_pred = (bbCP./DF)';
% plot(BD,CP_pred,'gsq','markersize',10,'linewidth',2);hold on;
% rho = corr(BD(~isnan(CP_pred))',CP_pred(~isnan(CP_pred)));
% [beta,gof]=local_fit_linear(BD(~isnan(CP_pred)),CP_pred(~isnan(CP_pred)));
% [coeffpval,coefftval] = local_stats(BD(~isnan(CP_pred)),CP_pred(~isnan(CP_pred)),beta,gof.dfe);
% plot([min(BD) max(BD)],feval(beta,[min(BD) max(BD)]),'g--','linewidth',2);
% if coeffpval<0.001
%     text(1.8,-3.5,sprintf('\\rho: %2.2f \n\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: <0.001',rho,beta.a,gof.dfe,coefftval),'interpreter','tex','color','g');
% else
%     text(1.8,-3.5,sprintf('\\rho: %2.2f \n\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',rho,beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex','color','g');
% end


set(gca,'xlim',[-3 3],'ylim',[-3 3],'ticklength',[0.02 0.02],'fontsize',16);
xlabel ('BEST DELAY_{\itNDF}\rm (ms)','interpreter','tex');
ylabel ('BEST DELAY_{\itpredicted}\rm (ms)','interpreter','tex');
box off;
axis square;


%% 
% xdata = (FrontLatency(:,2)-FrontLatency(:,1));
% ydata = ((GD(:,2)-FrontLatency(:,2))-(GD(:,1)-FrontLatency(:,1)));
% [beta,gof]=local_fit_linear(xdata(~isnan(xdata) & ~isnan(ydata)),ydata(~isnan(xdata) & ~isnan(ydata)));
% [coeffpval,coefftval] = local_stats(xdata(~isnan(xdata) & ~isnan(ydata)),ydata(~isnan(xdata) & ~isnan(ydata)),beta,gof.dfe);
% figure;
% 
% plot(xdata,ydata,'ko');

%% Plot Filter delay contra vs ipsi
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.3244    0.4189    0.2819    0.3689]);

xdata = GD(:,1)-FrontLatency(:,1);
ydata = GD(:,2)-FrontLatency(:,2);

mycmap = [[zeros(32,1);(linspace(0,1,16))';ones(16,1)] [(linspace(0,1,16))';ones(32,1);(linspace(1,0,16))'] [ones(16,1); linspace(1,0,16)'; zeros(32,1)]];
colormap(mycmap);
colormapyvals = mycmap;
minCF = min(DF(~isnan(xdata) & ~isnan(ydata)));
maxCF = max(DF(~isnan(xdata) & ~isnan(ydata)));
colormapxvals = logspace(log10(minCF),log10(maxCF),64);

for i = 1:length(xdata)
    if ~isnan(xdata(i))
        currentrgb = [interp1(colormapxvals,colormapyvals(:,1),DF(i)) interp1(colormapxvals,colormapyvals(:,2),DF(i)) interp1(colormapxvals,colormapyvals(:,3),DF(i))];
        plot(xdata(i),ydata(i),'color',currentrgb,'markersize',10,'marker','*');hold on;
    end
end

plot([0 ceil(max([xdata;ydata]))],[0 ceil(max([xdata;ydata]))],'color',[.6 .6 .6]);
set(gca,'xlim',[0 ceil(max([xdata;ydata]))],'ylim',[0 ceil(max([xdata;ydata]))],'linewidth',1,'ticklength',[0.02 0.02],'fontsize',16);
box off;
axis square;
xlabel 'IPSI. FILTER DELAY (ms)';
ylabel 'CONTRA. FILTER DELAY (ms)';
%Run a two-tailed paired t-test on the ipsi vs contra filter delay
%estimates
[H,P,C,S] = ttest2(xdata(~isnan(xdata) & ~isnan(ydata)),ydata(~isnan(xdata) & ~isnan(ydata)),'vartype','unequal');
%Display some text descriibg the stats
text(1, ceil(max([xdata;ydata]))-1, sprintf('\\itt\\rm_{(%2.1f)}: %2.2f \n\\itp\\rm: %2.2f',S.df,S.tstat,P),'interpreter','tex','fontsize',12);
%Add a colorbar
cba = colorbar;
set(cba,'fontsize',12,'linewidth',1,'ytick',[1 64],'yticklabel',{sprintf('%2.1f',minCF), sprintf('%2.1f',maxCF)},'ylim',[1 64]);
ylabel (cba,'DIFCOR DOMINANT FREQUENCY (kHz)');

%% Plot a histogram of differnce in filter delays (contra-ipsi)
xvals = [-4:0.1:4];
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.3244    0.4189    0.2819    0.3689]);
yvals = hist(ydata-xdata,xvals);
plot([0 0],[0 ceil(max(yvals))+1],'--','color',[.6 .6 .6]);hold on;
stairs(xvals,yvals,'k-','linewidth',2);
set(gca,'linewidth',1,'fontsize',16,'xlim',[-4 4],'ylim',[0 ceil(max(yvals))+1]);
xlabel ('\Delta FILTER DELAY (ms)','interpreter','tex');
ylabel 'N';
box off;
axis square;
%run a t-test on the difference in filter delay with the null hypothesis
%that the mean difference is zero
[H,P,C,S] = ttest(ydata-xdata,0);
text(-3, ceil(max(yvals))-1, sprintf('\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.2f',S.df,S.tstat,P),'interpreter','tex','fontsize',12);


%% Plot Signal-Front delay contra vs ipsi
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.3244    0.4189    0.2819    0.3689]);

xdata = FrontLatency(:,1);
ydata = FrontLatency(:,2);

mycmap = [[zeros(32,1);(linspace(0,1,16))';ones(16,1)] [(linspace(0,1,16))';ones(32,1);(linspace(1,0,16))'] [ones(16,1); linspace(1,0,16)'; zeros(32,1)]];
colormap(mycmap);
colormapyvals = mycmap;
minCF = min(DF(~isnan(xdata) & ~isnan(ydata)));
maxCF = max(DF(~isnan(xdata) & ~isnan(ydata)));
colormapxvals = logspace(log10(minCF),log10(maxCF),64);

for i = 1:length(xdata)
    if ~isnan(xdata(i))
        currentrgb = [interp1(colormapxvals,colormapyvals(:,1),DF(i)) interp1(colormapxvals,colormapyvals(:,2),DF(i)) interp1(colormapxvals,colormapyvals(:,3),DF(i))];
        plot(xdata(i),ydata(i),'color',currentrgb,'markersize',10,'marker','*');hold on;
    end
end

plot([2 ceil(max([xdata;ydata]))],[2 ceil(max([xdata;ydata]))],'color',[.6 .6 .6]);
set(gca,'xlim',[2 ceil(max([xdata;ydata]))],'ylim',[2 ceil(max([xdata;ydata]))],'linewidth',1,'ticklength',[0.02 0.02],'fontsize',16);
box off;
axis square;
xlabel 'IPSI. SIGNAL-FRONT DELAY (ms)';
ylabel 'CONTRA. SIGNAL-FRONT DELAY (ms)';
%Run a two-tailed paired t-test on the ipsi vs contra filter delay
%estimates
[H,P,C,S] = ttest2(xdata(~isnan(xdata) & ~isnan(ydata)),ydata(~isnan(xdata) & ~isnan(ydata)),'vartype','unequal');
%Display some text descriibg the stats
text(3, ceil(max([xdata;ydata]))-1, sprintf('\\itt\\rm_{(%2.1f)}: %2.2f \n\\itp\\rm: %2.2f',S.df,S.tstat,P),'interpreter','tex','fontsize',12);
%Add a colorbar
cba = colorbar;
set(cba,'fontsize',12,'linewidth',1,'ytick',[1 64],'yticklabel',{sprintf('%2.1f',minCF), sprintf('%2.1f',maxCF)},'ylim',[1 64]);
ylabel (cba,'DIFCOR DOMINANT FREQUENCY (kHz)');

%% Plot a histogram of differnce in signal-front delays (contra-ipsi)
xvals = [-4:0.1:4];
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.3244    0.4189    0.2819    0.3689]);
yvals = hist(ydata-xdata,xvals);
plot([0 0],[0 ceil(max(yvals))+1],'--','color',[.6 .6 .6]);hold on;
stairs(xvals,yvals,'k-','linewidth',2);
set(gca,'linewidth',1,'fontsize',16,'xlim',[-4 4],'ylim',[0 ceil(max(yvals))+1]);
xlabel ('\Delta SIGNAL-FRONT DELAY (ms)','interpreter','tex');
ylabel 'N';
box off;
axis square;
%run a t-test on the difference in filter delay with the null hypothesis
%that the mean difference is zero
[H,P,C,S] = ttest(ydata-xdata,0);
text(-3, ceil(max(yvals))-1, sprintf('\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.2f',S.df,S.tstat,P),'interpreter','tex','fontsize',12);

%% PLot a histogram of BDs

xvals = [-2:0.05:2];
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.3244    0.4189    0.2819    0.3689]);
yvals = hist(BD.*DF,xvals);
plot([0 0],[0 ceil(max(yvals))+1],'--','color',[.6 .6 .6]);hold on;
plot([0.5 0.5],[0 ceil(max(yvals))+1],'-','color','r');hold on;
plot([0.25 0.25],[0 ceil(max(yvals))+1],'-','color','g');hold on;
plot([0.125 0.125],[0 ceil(max(yvals))+1],'-','color','b');hold on;
stairs(xvals,yvals,'k-','linewidth',2);
set(gca,'linewidth',1,'fontsize',16,'xlim',[-2 2],'ylim',[0 ceil(max(yvals))+1]);
xlabel ('BD*DF (cycles)','interpreter','tex');
ylabel 'N';
box off;
axis square;
%run a t-test on the BDs with the null hypothesis
%that the mean is zero
[H,P,C,S] = ttest(BD.*DF,0);
text(-1.5, ceil(max(yvals))-1, sprintf('\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.2f',S.df,S.tstat,P),'interpreter','tex','fontsize',12);

%% Plot Q_3dB for the contra ear against the ipsi ear
Q_3_new = [Q_3(~isnan(Q_3(:,1)) & ~isnan(Q_3(:,2)),1) Q_3(~isnan(Q_3(:,1)) & ~isnan(Q_3(:,2)),2)];
ColV = Rho(~isnan(Q_3(:,1)) & ~isnan(Q_3(:,2)));
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
axes('ActivePositionProperty','outerposition');
scatter(Q_3_new(:,1),Q_3_new(:,2),100,ColV,'+');hold on;
colormap (flipud(hot));
set(gca,'clim',[0 1],'linewidth',1,'fontsize',16,'layer','top','ticklength',[0.02 0.02],'xlim',[0 6],'ylim',[0 6]);
xlabel ('IPSI \itQ_{3dB}','interpreter','tex');
ylabel ('CONTRA \itQ_{3dB}','interpreter','tex');
axis square;box off;
plot([0 6],[0 6],'-','linewidth',1,'color',[0.6 0.6 0.6]);
cbh = colorbar;
set(cbh,'fontsize',14);
ylabel (cbh,'\rho','interpreter','tex');
% [h,p,ci,stats] = ttest2(Q_3_new(:,1),Q_3_new(:,2)',0.05);
set(gca,'position',[0.1822    0.1100    0.5608    0.8150]);

%% Plot Q_6dB for the contra ear against the ipsi ear
Q_6_new = [Q_6(~isnan(Q_6(:,1)) & ~isnan(Q_6(:,2)),1) Q_6(~isnan(Q_6(:,1)) & ~isnan(Q_6(:,2)),2)];
ColV = Rho(~isnan(Q_6(:,1)) & ~isnan(Q_6(:,2)));
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
axes('ActivePositionProperty','outerposition');
scatter(Q_6_new(:,1),Q_6_new(:,2),100,ColV,'+');hold on;
colormap (flipud(hot));
set(gca,'clim',[0 1],'linewidth',1,'fontsize',16,'layer','top','ticklength',[0.02 0.02],'xlim',[0 4],'ylim',[0 4]);
xlabel ('IPSI \itQ_{6dB}','interpreter','tex');
ylabel ('CONTRA \itQ_{6dB}','interpreter','tex');
axis square;box off;
plot([0 4],[0 4],'-','linewidth',1,'color',[0.6 0.6 0.6]);
cbh = colorbar;
set(cbh,'fontsize',14);
ylabel (cbh,'\rho','interpreter','tex');
% [h,p,ci,stats] = ttest2(Q_6_new(:,1),Q_6_new(:,2)',0.05);
set(gca,'position',[0.1822    0.1100    0.5608    0.8150]);


%% Delta DF vs mean DF
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
axes('ActivePositionProperty','outerposition');
xdata = DF_P(~isnan(DF_P))';
% xdata = Delta.x(~isnan(Delta.x) & ~isnan(Delta.y))';
ydata = Delta.y(~isnan(Delta.x) & ~isnan(Delta.y))';
[xdata,ind] = sort(xdata);
ydata = ydata(ind);
plot(xdata,ydata,'k*','markersize',10);hold on;
set(gca,'linewidth',1,'fontsize',16,'layer','top','ticklength',[0.02 0.02],'ylim',[-1.0 1.0],'xlim',[0.1 1.2]);
xlabel ({'MEAN', 'DOMINANT FREQUENCY (kHz)'},'interpreter','tex');
ylabel ({'\Delta DOMINANT', 'FREQUENCY (octs)'},'interpreter','tex');
axis normal;box off;
plot([0 1.2],[0 0],'-','linewidth',1,'color',[0.6 0.6 0.6]);
% set(gca,'position',[0.1822    0.15    0.5608    0.815]);
set(gca,'ylim',[-1.2 1.2]);

% histax = axes('ActivePositionProperty','outerposition','position',[0.8 0.15 0.1 0.815]);
% yyn = hist(Delta.y,-1.2:0.1:1.2);
% plot(histax,yyn,-1.2:0.1:1.2,'k-','linewidth',2);hold on;
% box off;
% set(gca,'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'ylim',[-1.0 1.0],'ytick',[-1:0.5:1],'xdir','reverse','yaxislocation','right');
% xlabel ('\itn','interpreter','tex');
% plot([0 max(get(histax,'xlim'))],[0 0],'-','color',[0.6 0.6 0.6]);hold on;

[ydatamean,ydatastd] = mylowessbootstrap(xdata,ydata,0.3,'lowess',500);
plot(xdata,ydatamean,'r-','linewidth',2);
plot(xdata,ydatamean+ydatastd,'r--','linewidth',1);
plot(xdata,ydatamean-ydatastd,'r--','linewidth',1);

%% CD vs CP for narrowband estimates
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
axes('ActivePositionProperty','outerposition');
plot(CP_w,CD,'k*','markersize',10);hold on;
set(gca,'linewidth',1,'fontsize',16,'layer','top','ticklength',[0.02 0.02],'ylim',[-1.5 1.5],'xlim',[-0.5 0.5]);
xlabel ('CP_{\ittone\rm} (cycles)','interpreter','tex');
ylabel ('CD_{\ittone\rm} (ms)','interpreter','tex');
axis square;box off;
plot([-0.5 0.5],[0 0],'-','linewidth',1,'color',[0.6 0.6 0.6]);
plot([0 0],[-1.5 1.5],'-','linewidth',1,'color',[0.6 0.6 0.6]);
set(gca,'position',[0.1822    0.1100    0.5608    0.8150]);

%% CD vs CP for broadband estimates
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
axes('ActivePositionProperty','outerposition');
plot(bbCP,bbCD,'k*','markersize',10);hold on;
set(gca,'linewidth',1,'fontsize',16,'layer','top','ticklength',[0.02 0.02],'ylim',[-1.5 1.5],'xlim',[-0.5 0.5]);
xlabel ('CP_{\itnoise\rm} (cycles)','interpreter','tex');
ylabel ('CD_{\itnoise\rm} (ms)','interpreter','tex');
axis square;box off;
plot([-0.5 0.5],[0 0],'-','linewidth',1,'color',[0.6 0.6 0.6]);
plot([0 0],[-1.5 1.5],'-','linewidth',1,'color',[0.6 0.6 0.6]);
set(gca,'position',[0.1822    0.1100    0.5608    0.8150]);

%% Plot the predicted CP vs the noise CP - a torus - and make a movie (needs MATLAB 2013b)

% X = rand(1000,1);Y = rand(1000,1);
% [rho,pval,tscore]=circ_corrcc(X*2*pi,X*2*pi);

[rho,pval,tscore]=circ_corrcc(bbCP(~isnan(bbCP))*2*pi,CP_NDF(~isnan(bbCP))*2*pi);

figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.1956    0.1367    0.4219    0.7222],'color','w');
ax2 = axes;
set(ax2,'visible','off');
pos = get(gca,'position');
%Display some text on the figure indicating the correlation
if pval > 0.001
    text(0.1, 0.9, sprintf('\\rho: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',rho,length(bbCP_w(~isnan(bbCP_w)))-1,tscore,pval),'interpreter','tex','fontsize',15);
else
    text(0.1, 0.9, sprintf('\\rho: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: <0.001',rho,length(bbCP_w(~isnan(bbCP_w)))-1,tscore),'interpreter','tex','fontsize',15);
end

torus_a = 5;
torus_c = 10;
[u,v] = meshgrid(0:10:360);
x = (torus_c+torus_a*cosd(v)).*cosd(u);
y = (torus_c+torus_a*cosd(v)).*sind(u);
z = torus_a*sind(v);

ax = axes('position',pos);
mesh_h = mesh(x,y,z,'facelighting','gouraud');
set(mesh_h,'facealpha',0.7);
hold on;
props = {'CameraViewAngle','DataAspectRatio','PlotBoxAspectRatio'};
set(ax,props,get(ax,props));

%Draw identity line;
Uvals = 0:1:360;
Vvals = 0:1:360;
Xvals = (torus_c+torus_a*cosd(Vvals)).*cosd(Uvals);
Yvals = (torus_c+torus_a*cosd(Vvals)).*sind(Uvals);
Zvals = torus_a*sind(Vvals);
plot3(Xvals,Yvals,Zvals,'k-','linewidth',3);

Uvals = bbCP(~isnan(bbCP))*360;
Vvals = CP_NDF(~isnan(bbCP))*360;
% Uvals = X*360;
% Vvals = X*360;
Xvals = (torus_c+torus_a*cosd(Vvals)).*cosd(Uvals);
Yvals = (torus_c+torus_a*cosd(Vvals)).*sind(Uvals);
Zvals = torus_a*sind(Vvals);
scatter3(Xvals,Yvals,Zvals,200,Zvals,'filled');
set(ax,'visible','off');
axis equal;


set(gcf,'renderer','opengl','currentaxes',ax);
opengl software;
VV = VideoWriter('Torus.mp4','MPEG-4');
VV.FrameRate = 15;
VV.Quality = 100;
open(VV);

azimuth  = 1:2:360;
T = linspace(0,1,length(azimuth));
elevation  = 40*cos(2*pi*2.*T);
for i = 1:length(azimuth)
    set(ax,'view',[azimuth(i) elevation(i)]);
    frame = getframe(gcf);
    writeVideo(VV,frame);
end
close(VV);


%% PLot predicted CD vs noise CD
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
axes('ActivePositionProperty','outerposition');
plot(CD_NDF,bbCD,'k*','markersize',10);hold on;
rho = corr(CD_NDF(~isnan(bbCD))',bbCD(~isnan(bbCD))');
[beta,gof]=local_fit_linear(CD_NDF(~isnan(bbCD))',bbCD(~isnan(bbCD))');
[coeffpval,coefftval] = local_stats(CD_NDF(~isnan(bbCD))',bbCD(~isnan(bbCD))',beta,gof.dfe);
plot([-2 2],[-2 2],'color',[.6 .6 .6],'linewidth',1)


plot([min(CD_NDF(~isnan(bbCD))) max(CD_NDF(~isnan(bbCD)))],feval(beta,[min(CD_NDF(~isnan(bbCD))) max(CD_NDF(~isnan(bbCD)))]),'r--','linewidth',2);
text(1,-1,sprintf('\\rho: %2.2f \n\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',rho,beta.a,gof.dfe,coefftval,coeffpval),'interpreter','tex','color','r');
set(gca,'linewidth',1,'ticklength',[0.02 0.02],'fontsize',16,'layer','top','xlim',[-2 2],'ylim',[-2 2]);
axis square;
box off;
xlabel ('CD_{\itNDF\rm} (ms)','interpreter','tex');
ylabel ('CD_{\itpredicted\rm} (ms)','interpreter','tex');

%% BD vs dominant frequency
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
axes('ActivePositionProperty','outerposition');
plot([0 1.2],[0 0],'-','color',[0.6 0.6 0.6]);hold on;
scatter(DF(PeakTrough==1),BD(PeakTrough==1),100,'k+'); hold on;
scatter(DF(PeakTrough==2),BD(PeakTrough==2),100,'b+'); hold on;
pi_lim.x = [0:0.01:1.5];
pi_lim.y = [-0.5./(pi_lim.x);0.5./(pi_lim.x)];
qpi_lim.y = [-0.125./(pi_lim.x);0.125./(pi_lim.x)];
plot(pi_lim.x,pi_lim.y,'r-');
plot(pi_lim.x,qpi_lim.y,'g-');
set(gca,'xlim',[0.1 1.4],'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'ylim',[-5 5],'ytick',[-4:2:4]);
xlabel ({'Difcor','DOMINANT FREQUENCY (kHz)'});
ylabel 'BEST DELAY (ms)';
axis normal;box off;
% set(gca,'position',[0.1822    0.15    0.5608    0.815]);
% histax = axes('ActivePositionProperty','outerposition','position',[0.8 0.15 0.1 0.815]);
% yyn = hist(BD,-5:0.4:5);
% plot([0 max(yyn)*1.2],[0 0],'-','color',[0.6 0.6 0.6]);hold on;
% plot(histax,yyn,-5:0.4:5,'k-','linewidth',2);hold on;
% box off;
% set(gca,'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'ylim',[-5 5],'ytick',[-4:2:4],'xdir','reverse','yaxislocation','right');
% xlabel ('\itn','interpreter','tex');

%% BD*DF(cycles) vs dominant frequency
figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
axes('ActivePositionProperty','outerposition');
plot([0 0],[0 1.5],'-','color',[0.6 0.6 0.6]);hold on;
scatter(BD(PeakTrough==1).*DF(PeakTrough==1),DF(PeakTrough==1),100,'ko'); hold on;
scatter(BD(PeakTrough==2).*DF(PeakTrough==2),DF(PeakTrough==2),100,'ko'); hold on;
plot([0.5 0.5],[0 1.5],'r-');hold on;
plot([0.25 0.25],[0 1.5],'b-');hold on;
plot([0.125 0.125],[0 1.5],'g-');hold on;

set(gca,'ylim',[0.1 1.5],'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'xlim',[-1 1]);
ylabel ({'Difcor','DOMINANT FREQUENCY (kHz)'});
xlabel 'BD*DF (cycles)';
axis normal;box off;
% set(gca,'position',[0.1822    0.15    0.5608    0.815]);
% histax = axes('ActivePositionProperty','outerposition','position',[0.8 0.15 0.1 0.815]);
% yyn = hist(BD,-5:0.4:5);
% plot([0 max(yyn)*1.2],[0 0],'-','color',[0.6 0.6 0.6]);hold on;
% plot(histax,yyn,-5:0.4:5,'k-','linewidth',2);hold on;
% box off;
% set(gca,'fontsize',16,'linewidth',1,'ticklength',[0.02 0.02],'ylim',[-5 5],'ytick',[-4:2:4],'xdir','reverse','yaxislocation','right');
% xlabel ('\itn','interpreter','tex');




%% NRHO exponent vs. Dominant frequency
% 
% figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2600    0.2644    0.2675    0.3878]);
% axes('ActivePositionProperty','outerposition');
% scatter(DF,nrhoP,100,Rho,'+','linewidth',1);hold on;
% colormap (flipud(hot));
% set(gca,'clim',[0 1],'linewidth',1,'fontsize',16,'layer','top','ticklength',[0.02 0.02]);
% ylabel ('NRHO EXPONENT','interpreter','tex');
% xlabel ('DOMINANT FREQUENCY (kHz)','interpreter','tex');
% axis square;box off;
% plot([0 1.4],[1 1],'--','linewidth',1,'color',[0.6 0.6 0.6]);
% cbh = colorbar;
% set(cbh,'fontsize',14);
% ylabel (cbh,'\rho','interpreter','tex');
% set(gca,'position',[0.1822    0.1100    0.5608    0.8150]);

%% Plot a series of noise-delay curves

NDFlist = find(NDFtype==1);
NDFlist = NDFlist(BD(NDFlist)>0);%restrict to positive BDs
NDFlist = NDFlist(PeakTrough(NDFlist)==1);
% NDFlist = NDFlist(NDFtype(NDFlist)==1);%restrict to real NDFs
%Get all DFs for those NDFs selected
DFlist = DF(NDFlist);
%Put them in order of DF
[DFlist,indx] = sort(DFlist,'descend');
NDFlist = NDFlist(indx);

%Define a list of DFs to search for
DFsearch = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
for i = 1:length(DFsearch)
    %Find a fiber matching this DF
    dif = abs(DFlist-DFsearch(i));
    ind = find(dif==min(dif));
    if length(ind)>1
        ind = ind(2);
    end
    newNDFlist(i) = NDFlist(ind);
    newDFlist(i) = DFlist(ind);
end
DFlist = newDFlist;clear newDFlist;
NDFlist = newNDFlist;clear newNDFlist;

figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.0869    0.1100    0.2800    0.7544]);
ax1 = axes('ActivePositionProperty','outerposition');
count = 0;
for i = 1:length(DFlist)
    plot(ax1,NDFx{NDFlist(i)},(NDFy{NDFlist(i)}(1,:)-mean(mean(NDFy{NDFlist(i)})))+(0.5*i),'b-','linewidth',3);hold on;
    plot(ax1,NDFx{NDFlist(i)},(NDFy{NDFlist(i)}(2,:)-mean(mean(NDFy{NDFlist(i)})))+(0.5*i),'m-','linewidth',3);hold on;
end
plot([0 0],[0 5.5],'k--','linewidth',1);
set(gca,'linewidth',1,'fontsize',20,'ticklength',[0.02 0.02],'layer','top','xlim',[-5 5],'ylim',[0 5.5],'ytick',[1:5],'yticklabel',[0.2:0.2:1]);
box off;
xlabel 'ITD (ms)';
ylabel 'DOMINANT FREQUENCY (kHz)';

figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.0869    0.1100    0.2800    0.7544]);
ax2 = axes('ActivePositionProperty','outerposition');
colorord = get(gca,'colororder');
colorord = [colorord;colorord];
for i = length(DFlist):-1:1
    plot(ax2,NDFx{NDFlist(i)},(NDFz{NDFlist(i)}*0.5)+(0.5*i),'k-','linewidth',3);hold on;
    plot(ax2,BD(NDFlist(i)),0.125+(i*0.5),'o','color',colorord(i,:),'markersize',15,'linewidth',3,'markerfacecolor',colorord(i,:));
end
plot([0 0],[0 5.5],'k--','linewidth',1);
set(gca,'linewidth',1,'fontsize',20,'ticklength',[0.02 0.02],'layer','top','xlim',[-5 5],'ylim',[0 5.5],'ytick',[1:5],'yticklabel',[0.2:0.2:1]);
box off;
xlabel 'ITD (ms)';
ylabel 'DOMINANT FREQUENCY (kHz)';

    function [beta_out,gof_out]=local_fit_linear(varargin)
        if nargin ==3
            x_in = varargin{1};
            y_in = varargin{2};
            weight_in = varargin{3};
            weights = weight_in;
            weights = weights./max(weights);
            s = fitoptions('Method','LinearLeastSquares','Lower',[-Inf -Inf],'Upper',[Inf Inf],'Weights',weights,'robust','on');
        elseif nargin ==2
            x_in = varargin{1};
            y_in = varargin{2};
            s = fitoptions('Method','LinearLeastSquares','Lower',[-Inf -Inf],'Upper',[Inf Inf],'robust','on');
        end
        f = fittype({'x','1'},'coefficients',{'a','b'},'options',s);
        [beta_out,gof_out] = fit(x_in(:),y_in,f);
    end
    function [coeffpval,coefftval] = local_stats(xdata,ydata,beta,dfe)
        r=ydata-feval(beta,xdata);%residuals
        SSres=sum(r.^2);%sum square of residuals
        coeffSE = sqrt(SSres/dfe)/sqrt(sum((xdata-mean(xdata)).^2));%standard error of the regression coefficient
        coefftval = beta.a/coeffSE;
        coeffpval=2*(1-tcdf(abs(coefftval),dfe));
    end
end
































