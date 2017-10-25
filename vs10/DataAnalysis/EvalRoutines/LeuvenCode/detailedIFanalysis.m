function [] = detailedIFanalysis

datadir = 'C:\work\BinRev\BinRev';
anadir = 'C:\LeuvenDataAnalysis';

cd (datadir)
datafiles = dir('*_BinRev.mat');
nrdata = numel(datafiles);
figure;set(gcf,'units','normalized','position',[0.3244 0.1333 0.3500 0.7522]);
axh = axes;

for i = 1:nrdata
    cd (datadir);
    data = load (datafiles(i).name);
    ParamsOut = data.ParamsOut;clear data;
    if isfield(ParamsOut,'Revcor')
        h1DF = ParamsOut.Revcor.domfreq;
        Time  = ParamsOut.Revcor.IFTime;
        dt = Time(2)-Time(1);
        BD = ParamsOut.NDF.BD_NTD;
        CD_NDF(i) = ParamsOut.NDF.BD_NTD_env;
        CP_NDF(i) = (BD-CD_NDF(i))*ParamsOut.NDF.DifcorDomFreq_NTD;
        CP_NDF(i) = wrapToPi(CP_NDF(i)*2*pi)/(2*pi);
        h1 = ParamsOut.Revcor.h1filt;
        IP = unwrap(angle(hilbert(h1)));
        IF = [0 0; abs((1/dt)/(2*pi)*diff(IP))];%instantaneous frequency (Hz)
        
        envz = ParamsOut.Revcor.h1filtenvzscore;
        GD = -ParamsOut.Revcor.h1phasefit_coeffs(:,1)';
        %         GD = ParamsOut.Revcor.h1frontlatency;
        for j = 1:2
            dif = (abs(Time-GD(j)));
            [minval,ind(j)] = min(dif);
        end
        predicted_CP(i) = IP(ind(2),2)-IP(ind(1),1);
%         predicted_CP(i) = predicted_CP(i)/(2*pi);
        %This compensates everything for group delay
%         if GD(1)>GD(2);
%             indx = round((GD(1)-GD(2))/dt);
%             sIP(:,1) = IP(indx+1:end,1);
%             sIP(:,2) = IP(1:end-indx,2);
%             senvz(:,1) = envz(indx+1:end,1);
%             senvz(:,2) = envz(1:end-indx,2);
%             sIF(:,1) = IF(indx+1:end,1);
%             sIF(:,2) = IF(1:end-indx,2);
%             %             sTime = Time(1:length(sIF));
%         else
%             indx = round((GD(2)-GD(1))/dt);
%             sIP(:,2) = IP(indx+1:end,2);
%             sIP(:,1) = IP(1:end-indx,1);
%             senvz(:,2) = envz(indx+1:end,2);
%             senvz(:,1) = envz(1:end-indx,1);
%             sIF(:,2) = IF(indx+1:end,2);
%             sIF(:,1) = IF(1:end-indx,1);
% %             sTime = Time(1:length(sIF));
%         end
%         
%         for j = 1:2
%             [maxval,ind] = max(senvz(:,j));
%             startval = find(senvz(1:ind,j)<3,1,'last')+1;
%             stopval = find(senvz(ind:end,j)<3,1,'first')+ind-2;
%             senvz(1:startval-1,j) = 0;
%             senvz(stopval+1:end,j) = 0;
%             sIP(1:startval-1,j) = NaN;
%             sIP(stopval+1:end,j) = NaN;
%             sIF(1:startval-1,j) = NaN;
%             sIF(stopval+1:end,j) = NaN;
%             senvz(:,j) = senvz(:,j)./max(senvz(:,j));
%         end
%         
%         subplot(3,1,1);
%         plot(senvz,'linewidth',1);
%         subplot(3,1,2);
%         plot(sIP,'linewidth',1);
%         subplot(3,1,3);
%         plot(sIF,'linewidth',1);
%         
        diffDF(i) = h1DF(2)-h1DF(1);
%         dIF = sIF(:,2)-sIF(:,1);
%         
%         
%         weights = prod(senvz,2);
%         
%         diffIF(i) = sum(dIF(~isnan(dIF)).*weights(~isnan(dIF)))/sum(weights(~isnan(dIF)));
%         diffIP = sIP(:,1)-sIP(:,2);
        
%         cd(anadir);
%         predicted_CP(i) = circ_mean(diffIP(~isnan(diffIP)),weights(~isnan(diffIP)));
% predicted_CP(i) = circ_mean(diffIP(~isnan(diffIP)));
        predicted_CP(i) = wrapToPi(predicted_CP(i));
        predicted_CP(i) = predicted_CP(i)./(2*pi);
        predicted_CP(i) = predicted_CP(i)+ParamsOut.Revcor.k;
        
        %         savedfilename = datafiles(i).name;
        %         save (savedfilename,'ParamsOut');
        %         disp (['Saved: ' savedfilename]);
        clear sIF;
        clear sIP;
        clear senvz;
        clear sTime;
    else
%         diffIF(i) = NaN;
        diffDF(i) = NaN;
        predicted_CP(i) = NaN;
        CP_NDF(i) = NaN;
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

























