function [] = BFS_plotpop

currdir = cd;
datadir = 'C:\work\BinRev\BFS';
cd (datadir);
datafiles = dir('*_BFS.mat');
CDdata = []; CPdata = []; NDFdata = []; TDFdata = [];
for i =1:length(datafiles)
    load (datafiles(i).name);
    CDdata = [CDdata dataout.CD];
    CPdata = [CPdata dataout.CP];
    NDFdata = [NDFdata dataout.NDF];
    TDFdata = [TDFdata dataout.TDF];
end

currdir = cd;
datadir = 'C:\work\BinRev\BinRev';
cd (datadir);

datafiles = dir('*_BinRev.mat');
meanF = []; diffF = []; BD = [];diffO = [];Difcorxcor = [];
Latency_I =  [];Latency_C = [];
GD_I = [];GD_C = [];
for i =1:length(datafiles)
    load (datafiles(i).name);
    meanF = [meanF ParamsOut.DifcorDomFreq_NTD];
    diffF = [diffF ParamsOut.h1magxcorpeakdfreq];
    diffO = [diffO ParamsOut.h1magxcorpeakdoct];
    BD = [BD ParamsOut.BD_NTD];
    Difcorxcor = [Difcorxcor ParamsOut.Difcor_xcor];
    Latency_I = [Latency_I ParamsOut.h1latency(1)];
    Latency_C = [Latency_C ParamsOut.h1latency(2)];
    GD_I = [GD_I ParamsOut.h1phasefit_coeffs(1,1)];
    GD_C = [GD_C ParamsOut.h1phasefit_coeffs(2,1)];
end

currdir = cd;

figure;
set(gcf,'paperpositionmode','auto');
scatter3(CPdata,CDdata,Difcorxcor,100,Difcorxcor,'filled');
colormap (flipud(hot));view (0,90);
set(gca,'fontsize',16);
xlabel 'CP (Cycles)';
ylabel 'CD (ms)';
axis square;hcol = colorbar;
set(hcol,'fontsize',16);

figure;
set(gcf,'paperpositionmode','auto');
scatter3(Latency_I,Latency_C,BD,100,BD,'filled');
colormap (jet);view (0,90);
set(gca,'clim',[-1 1],'xlim',[4 10],'ylim',[4 10]);
set(gca,'fontsize',16);
xlabel 'Ipsi. Latency (ms)';
ylabel 'Contra. Latency (ms)';
axis square;hcol = colorbar;
set(hcol,'fontsize',16);
ylabel (hcol,'BD (ms)');
hold on
plot([4 10],[4 10],'k--','linewidth',1);

figure;
set(gcf,'paperpositionmode','auto');
scatter3(NDFdata,TDFdata,Difcorxcor,100,Difcorxcor,'filled');
colormap (flipud(hot));view (0,90);
set(gca,'fontsize',16);
xlabel 'Noise DF (kHz)';
ylabel 'Tone DF (kHz)';
axis square;
cd (currdir);hcol = colorbar;
set(hcol,'fontsize',16);

figure;
set(gcf,'paperpositionmode','auto');
scatter3(meanF,diffF,Difcorxcor,100,Difcorxcor,'filled');
colormap (flipud(hot));view (0,90);
set(gca,'fontsize',16);
xlabel 'Mean Freq (kHz)';
ylabel 'Freq Difference (kHz)';
axis square;hcol = colorbar;
set(hcol,'fontsize',16);

figure;
set(gcf,'paperpositionmode','auto');
scatter3(meanF,diffO,Difcorxcor,100,Difcorxcor,'filled');
colormap (flipud(hot));view (0,90);
set(gca,'fontsize',16);
xlabel 'Mean Freq (kHz)';
ylabel 'Freq Difference (Octaves)';
axis square;hcol = colorbar;
set(hcol,'fontsize',16);

figure;
set(gcf,'paperpositionmode','auto');
scatter3(meanF,diffF,Difcorxcor,100,Difcorxcor,'filled');
colormap (flipud(hot));view (0,90);
set(gca,'fontsize',16);
xlabel 'BD (ms)';
ylabel 'Freq Difference (kHz)';
axis square;hcol = colorbar;
set(hcol,'fontsize',16);

figure;
set(gcf,'paperpositionmode','auto');
scatter3(meanF,BD,Difcorxcor,100,Difcorxcor,'filled');
colormap (flipud(hot));view (0,90);
set(gca,'fontsize',16);
ylabel 'BD (ms)';
xlabel 'DF (kHz)';
axis square;hcol = colorbar;
set(hcol,'fontsize',16);
