clear;clc;

load ('C:\ExpData\L16553\L16553_177.mat');

nconds = size(spikes,1);

t = 0:0.2:110; %Time vector (millseconds)

mu = cell(1,nconds);
sigma = cell(1,nconds);
CV = cell(1,nconds);

for i = 1 : nconds
    %concatenate, nan-padd, and transpose the spike times matrix - one column per rep
    s = catpad(1,spikes{i,:})';
    
    [mu{i},sigma{i},CV{i}] = VCN_regularity(0:0.2:110,s,'exact',5);
end

figure;

for i = 1 : nconds
    
    dBval = stim_param.SPL(i);
    
    subplot(3,3,i);
    
    plot(t,CV{i},'r-','linewidth',1);
    
    set(gca,'xlim',[0 60],'fontsize',12,'tickdir','out','ylim',[0 1]);
    
    xlabel 'Time (ms)';
    
    ylabel ' Coeffcient of Variation';
    
    title (sprintf('%2.0f dB SPL',dBval));
    
    box off;
    
end