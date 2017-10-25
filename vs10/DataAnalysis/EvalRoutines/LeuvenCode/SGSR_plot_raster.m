function [] = SGSR_plot_raster(AnNum,SeqID)
cols = {'k','r','b','g','m','c'};

dataobj = dataset(AnNum,SeqID);
data = struct(dataobj);
RepDur = data.Stimulus.Special.RepDur;

switch data.ID.StimType
    case 'FSLOG'
        IndepVar = data.Stimulus.Special.CarFreq;
        IndepVar_unit = 'Hz';
    case 'SPL'
        IndepVar = data.Stimulus.IndepVar.Values;
        IndepVar = IndepVar(~isnan(IndepVar));
        IndepVar_unit = 'dB SPL';
end

figure;set(gcf,'paperpositionmode','auto','units','normalized',...
    'position',[0.0869    0.0533    0.3500    0.8556],...
    'activepositionproperty','outerposition');
[nrconds,nrreps]=size(data.Data.SpikeTimes);
dY = 1/(nrconds*nrreps);
dC = 1/nrconds;
jj = 1;
c_ix = 1;
while jj<=nrconds
    if c_ix>6
        c_ix = 1;
    end
    for k = 1:nrreps
        Nsp = numel(data.Data.SpikeTimes{jj,k});
        plot([data.Data.SpikeTimes{jj,k};data.Data.SpikeTimes{jj,k}],repmat(((jj-1)*dC)+([k-1 k]*dY),Nsp,1)',cols{c_ix});
        hold on;
    end
    jj = jj+1;
    c_ix = c_ix+1;
end
set(gca,'fontsize',14,'ytick',[dC/2:dC:1],'yticklabel',round(IndepVar),'tickdir','out','xlim',[0 RepDur],...
    'activepositionproperty','outerposition');
xlabel 'Time (ms)';
switch data.ID.StimType
    case 'FSLOG'
        ylabel 'Carrier Frequency (Hz)';
    case 'SPL'
        ylabel 'Sound Level (dB SPL)';
end
box off;
