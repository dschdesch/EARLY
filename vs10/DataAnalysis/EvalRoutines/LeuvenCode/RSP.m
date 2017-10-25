function [] = RSP(varargin)

annum = varargin{1};
for i = 2:nargin
    ds{i-1} = dataset(annum,varargin{i});
end

for i = 1:length(ds)
    xvals(:,i) = ds{i}.Stimulus.IndepVar.Values;
    nrconds = length(ds{i}.Stimulus.IndepVar.Values);
    nrreps = ds{i}.Stimulus.StimParam.stimcntrl.repcount;
    stimdur = ds{i}.Stimulus.StimParam.indiv.stim{1}.duration;
    yvals(:,i) = repmat(ds{i}.Stimulus.StimParam.indiv.stim{1}.spl,length(ds{i}.Stimulus.IndepVar.Values),1);
    S = ds{i}.SPT;
    for j = 1:nrconds
        spiketimes = [S{j,:}];
        spiketimes = spiketimes(spiketimes<=stimdur);
        zvals(j,i) = length(spiketimes)/nrreps/(stimdur/1000);
    end
end
figure;surf(xvals,yvals,zvals)
shading interp
view (0,90)
set(gca,'xscale','log')