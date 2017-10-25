function [] = AddGCfits

datadir = 'C:\work\BinRev\BinRev';
anadir = 'C:\LeuvenDataAnalysis';

cd (datadir)
datafiles = dir('*_BinRev.mat');
nrdata = numel(datafiles);
for i = 53:nrdata
    data = load (datafiles(i).name);
    ParamsOut = data.ParamsOut;clear data;
    if isfield(ParamsOut,'Revcor')
        h1{i} = ParamsOut.Revcor.h1filt; %#ok<AGROW>
        time  = ParamsOut.Revcor.Time;
        envz = ParamsOut.Revcor.h1filtenvzscore;
        Latency(i,:) = ParamsOut.Revcor.h1latency;
        cd (anadir);
        for j = 1:2
            latencyind = find(time==Latency(i,j));
            startind = find(envz(1:latencyind,j)<3,1,'last')+1;
            stopind = find(envz(latencyind:end,j)<3,1,'first')-1+latencyind;
            tt = time(startind:stopind)-time(startind);
            hh = h1{i}(startind:stopind,j);
            [x{j},resnorm{j},residual{j},exitflag{j},output{j}] = fitGCtoIR(tt,hh);
        end
        cd(datadir);
        ParamsOut.GCFit.x = x;
        ParamsOut.GCFit.residual = residual;
        ParamsOut.GCFit.exitflag = exitflag;
        ParamsOut.GCFit.resnorm = resnorm;
        ParamsOut.GCFit.output = output;
        savedfilename = datafiles(i).name;
        save (savedfilename,'ParamsOut');
        disp (['Saved: ' savedfilename]);
    end
end