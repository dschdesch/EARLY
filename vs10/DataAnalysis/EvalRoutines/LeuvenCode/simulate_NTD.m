function [] = simulate_NTD(CD,CP)

%simulates a noise delay function as a SAM tone.
%Inputs
%CD = envelope delay in ms (+ = contra leading, - = ipsi leading)
%CP = fine structure delay in cycles

sr = 50e3; %sample rate in Hz;
nf = sr/2;
dt = 1/sr; %sample period in s;

Fcar = 1000;
Fmod = Fcar/6;

Nsamps = (1/Fmod)*sr;

t = 0:dt:dt*(Nsamps-1);

env = (1+(sin(2*pi*Fmod.*t - (pi/2))));
s  = cos(2*pi*Fcar.*t - CP*(2*pi)).*env;

totalzeros = Nsamps;

CDsamps = (CD/1000)*sr;
if CD>0
    negzeros = (totalzeros/2)+CDsamps;
    poszeros = (totalzeros/2)-CDsamps;
elseif CD<0
    poszeros = (totalzeros/2)+abs(CDsamps);
    negzeros = (totalzeros/2)-abs(CDsamps);
else
    poszeros = totalzeros/2;
    negzeros = totalzeros/2;
end

s = [zeros(1,negzeros) s zeros(1,poszeros)];
env = [zeros(1,negzeros) env zeros(1,poszeros)];

newNsamps = length(s);

iat = [-dt*((newNsamps/2)):dt:-dt 0:dt:dt*((newNsamps/2)-1)];

[BDy,indBD] = max(s);
BDx = iat(indBD);
BD = CD+(CP/(Fcar/1000));


figure;set(gcf,'paperpositionmode','auto','units','normalized','position',[0.2838    0.3300    0.3544    0.3156]);
plot([0 0],[-2.5 2.5],'--','color',[.6 .6 .6],'linewidth',1);hold on;
plot([-4 4],[0 0],'--','color',[.6 .6 .6],'linewidth',1);hold on;
plot([CD CD],[-2.5 2.5],'b-','linewidth',2);
plot([CD+(CP/(Fcar/1000)) CD+(CP/(Fcar/1000))],[-2.5 2.5],'r-','linewidth',2);
plot(iat*1000,s,'k-',iat*1000,env,'b--','linewidth',2);
plot(BDx*1000,BDy,'o','markerfacecolor','r','markeredgecolor','r','markersize',12);

set(gca,'xlim',[-3 3],'ylim',[-2.5 2.5],'fontsize',16,'ytick',[-2.5 0 2.5],'yticklabel',{'-','0','+'},'ticklength',[0.01 0.01],'linewidth',1);
xlabel 'INTER-AURAL TIME DIFFERENCE (ms)';
ylabel ('\Delta FIRING RATE','interpreter','tex');
box off;
axis normal;

text(1,-2,sprintf('CD: %2.2f ms \nCP: %2.2f cycles \nBD: %2.2f ms',CD,CP,BD),'fontsize',12);


%take the FFT and try to re-compute the CP
NFFT = 2^13;
spec = fft(s,NFFT);
mag = 20*log10(spec);
phase = angle(spec);
specfreq = sr*linspace(0,0.5,NFFT/2)/1000; % freq in kHz



















