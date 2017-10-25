function [] = CPanalysis

datadir = 'C:\work\BinRev\BinRev';
anadir = 'C:\LeuvenDataAnalysis';

cd (datadir)
datafiles = dir('*_BinRev.mat');
nrdata = numel(datafiles);

for i = 1:nrdata
    cd (datadir);
    data = load (datafiles(i).name);
    ParamsOut = data.ParamsOut;clear data;
    if isfield(ParamsOut,'Revcor')
        h1DF = ParamsOut.Revcor.domfreq;
        freqax = ParamsOut.Revcor.h1ffax;
        h1Phase = ParamsOut.Revcor.h1phasemask_uw;
        BD(i) = ParamsOut.NDF.BD_NTD;
        DifCorDF(i) = ParamsOut.NDF.DifcorDomFreq_NTD;
        GD = -ParamsOut.Revcor.h1phasefit_coeffs(:,1)';
        if strcmp(ParamsOut.Peak_Trough,'P')
            PT(i) = 1;
        else
            PT(i) = 2;
        end
        kvals = [-10:1:10];
        %Find the phase for each ear at its CF
        for j = 1:2
            dif = (abs(freqax-h1DF(j)));
            [~,ind(j)] = min(dif);
            Phase(j) = h1Phase(ind(j),j);
        end
        %predicted CP without considering k
        predicted_CP(i) = (2*pi)-(Phase(2)-Phase(1));%radians
        
        predicted_CP(i) = predicted_CP(i)/(2*pi);
        
%         predicted_CP(i) = 1-wrapToPi(Phase(2)- Phase(1))/(2*pi);
        diffDF(i) = h1DF(2)-h1DF(1);
        
        CD = GD(2) - GD(1);%contra - ipsi (i.e., contra leading should be a positive CD).
        CP(i) = (BD(i)-CD)*DifCorDF(i);
        CP(i) = CP(i)*2*pi;
%         CP(i) = wrapToPi(CP(i)*2*pi);
        
%         predicted_CP(i) = wrapToPi(predicted_CP(i));
        
        
%         We know that BD = CD+(CP/DF);
        for j = 1:length(kvals)
            CPest(j) = (predicted_CP(i)+kvals(j));
            BDest(j) = CD+(CPest(j)/DifCorDF(i));
            BDest_error(j) = abs(BD(i)-BDest(j));
        end
        [minerror,minind] = min(BDest_error);
        predicted_CP(i) = CPest(minind);
        predicted_BD(i) = BDest(minind);
    else
        diffDF(i) = NaN;
        predicted_CP(i) = NaN;
        CP(i) = NaN;
        BD(i) = NaN;
        PT(i) = NaN;
        DirCorDF(i) = NaN;
        predicted_BD(i) = NaN;
    end
end
cd(anadir);

figure;set(gcf,'paperpositionmode','auto');
plot(diffDF,diffIF,'k*','markersize',10);hold on;
set(gca,'xlim',[-0.3 0.7],'ylim',[-0.3 0.7],'fontsize',16);
axis square;
box off;
xlabel ('\Delta DOMINANT FREQ. (kHz)','interpreter','tex');
ylabel ('\Delta INSTANTANEOUS FREQ. (kHz)','interpreter','tex');
%add the identity line;
plot([-0.3 0.7],[-0.3 0.7],'--','color',[0.6 0.6 0.6]);
%Calculate the correlation coefficient between delta DF and delta IF
rho = corr(diffDF(~isnan(diffDF))',diffIF(~isnan(diffIF))');
%Linear regression
[beta,gof] = local_fit_linear(diffDF(~isnan(diffDF))',diffIF(~isnan(diffIF))');
[pval,tscore] = local_stats(diffDF(~isnan(diffDF))',diffIF(~isnan(diffIF))',beta,gof.dfe);
plot([-0.3 0.7],feval(beta,[-0.3 0.7]),'r--','linewidth',1);
%display some text
if pval<0.001
    text(-0.2,0.5,sprintf('\\rho: %2.2f \n\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: <0.001',...
        rho,beta.a,gof.dfe,tscore),'interpreter','tex','color','r','fontsize',12);
else
    text(-0.2,0.5,sprintf('\\rho: %2.2f \n\\beta_{1}: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',...
        rho,beta.a,gof.dfe,tscore,pval),'interpreter','tex','color','r','fontsize',12);
end


%Now the CP analysis

[rho,pval,tscore]=circ_corrcc(CP_NDF(~isnan(CP_NDF))*2*pi,predicted_CP(~isnan(CP_NDF))*2*pi);%circular-circular correlation coefficient

figure;set(gcf,'paperpositionmode','auto');
plot(CP_NDF,predicted_CP,'k*','markersize',12);hold on;
set(gca,'xlim',[-0.5 0.5],'ylim',[-0.5 0.5],'ticklength',[0.02 0.02],'fontsize',16,'linewidth',1);
xlabel ('CHARACTERISTIC PHASE_{(\itNDF\rm)} \rm(cycles)','interpreter','tex');
ylabel ({'PREDICTED','CHARACTERISTIC PHASE (cycles)'},'interpreter','tex');
box off;
axis square;
plot([-0.5 0.5],[0 0],'--','color',[.6 .6 .6]);
plot([0 0],[-0.5 0.5],'--','color',[.6 .6 .6]);
plot([-0.5 0.5],[-0.5 0.5],'-','color',[.6 .6 .6],'linewidth',1);



figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.1956    0.1367    0.4219    0.7222],'color','w');
ax2 = axes;
set(ax2,'visible','off');
pos = get(gca,'position');
%Display some text on the figure indicating the correlation
if pval > 0.001
    text(0.1, 0.9, sprintf('\\rho: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: %2.3f',rho,length(CP_NDF(~isnan(CP_NDF)))-1,tscore,pval),'interpreter','tex','fontsize',15);
else
    text(0.1, 0.9, sprintf('\\rho: %2.2f \n\\itt\\rm_{(%2.0f)}: %2.2f \n\\itp\\rm: <0.001',rho,length(CP_NDF(~isnan(CP_NDF)))-1,tscore),'interpreter','tex','fontsize',15);
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

Uvals = CP_NDF(~isnan(CP_NDF))*360;
Vvals = predicted_CP(~isnan(CP_NDF))*360;
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


    function [beta_out,gof_out]=local_fit_linear(varargin)
        if nargin ==3
            x_in = varargin{1};
            y_in = varargin{2};
            weight_in = varargin{3};
            wt = weight_in;
            wt = wt./max(wt);
            s = fitoptions('Method','LinearLeastSquares','Lower',[-Inf -Inf],'Upper',[Inf Inf],'Weights',wt,'robust','Bisquare');
        elseif nargin ==2
            x_in = varargin{1};
            y_in = varargin{2};
            s = fitoptions('Method','LinearLeastSquares','Lower',[-Inf -Inf],'Upper',[Inf Inf],'robust','Bisquare');
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

























