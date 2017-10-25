function []= plot_spike_wave(dataIN,fs);

dt = 1/fs;
T = 0:dt:(numel(dataIN)*dt)-dt;


%The whole spike trace
figure;set(gcf,'units','normalized','paperpositionmode','auto','position',[0.0119    0.4033    0.9775    0.2633]);
plot(T,dataIN,'k-','linewidth',1);hold on;
plot(T(T>=2.5 & T<=3.5),dataIN(T>=2.5 & T<=3.5),'r-','linewidth',1);
box off;
set(gca,'fontsize',20,'ticklength',[0.01 0.01],'ytick',[],'xtick',[0:1:6],'ylim',[-0.01 0.1]);
xlabel 'Time (s)';

%Get a splatter plot and mean spike

[dum,ind] = findpeaks(dataIN,'minpeakheight',0.02);

for i = 1:length(ind)
    spike{i} = dataIN(ind(i)-37:ind(i)+37);
end

T2 = [-37*dt:dt:0 dt:dt:dt*37];

spikes = vertcat(spike{:});
figure;set(gcf,'units','normalized','paperpositionmode','auto','position',[0.2    0.4033    0.235    0.2633]);
plot(T2*1000,spikes,'color',[.6 .6 .6],'linewidth',0.5);hold on;
plot(T2*1000,mean(spikes,1),'k-','linewidth',2);
set(gca,'fontsize',20,'ticklength',[0.01 0.01],'ytick',[],'xtick',[-1 0 1],'ylim',[-0.01 0.1],'xlim',[-1.5 1.5]);
xlabel 'Time (ms)';
box off;


figure;set(gcf,'units','normalized','paperpositionmode','auto','position',[0.0119    0.4033    0.9775    0.2633]);
plot(T,dataIN,'r-','linewidth',1);
box off; 
set(gca,'fontsize',20,'ticklength',[0.01 0.01],'ytick',[],'xtick',[2.5 3.5],'xlim',[2.5 3.5],'ylim',[-0.01 0.1]);
xlabel 'Time (s)';



figure;set(gcf,'units','normalized','paperpositionmode','auto','position',[0.0119    0.4033    0.9775    0.2633]);
plot(T, sin(2*pi*20.*T),'b-','linewidth',1);hold on;
plot(T, sin(2*pi*21.*T),'g-','linewidth',1);hold on;
set(gca,'fontsize',20,'ticklength',[0.01 0.01],'ytick',[],'xtick',[2.5 3 3.5],'xlim',[2.5 3.5],'ylim',[-1.1 1.1]);
xlabel 'Time (s)';