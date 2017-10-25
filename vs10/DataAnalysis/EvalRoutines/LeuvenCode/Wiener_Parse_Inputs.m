function [datastruct,params] = Wiener_Parse_Inputs(MaxOrder,Nchan,IN)

params.MaxOrder = MaxOrder;
params.Nchan = Nchan;
params.AnNum = IN{1};
params.UnNum = IN{2};
for i = 3:length(IN)
    if strcmp(IN{i-1},'ndp')
        if iscell(IN{i})
            ndpdata = IN{i};
        else
            ndpdata = IN(i);
        end
    elseif strcmp(IN{i-1},'ndn')
        if iscell(IN{i})
            ndndata = IN{i};
        else
            ndndata = IN(i);
        end
    elseif strcmp(IN{i-1},'bfs')
        if iscell(IN{i})
            bfsdata = IN{i};
        else
            bfsdata = IN(i);
        end
    elseif strcmp(IN{i-1},'ndz')
        if iscell(IN{i})
            ndzdata = IN{i};
        else
            ndzdata = IN(i);
        end
    elseif strcmp(IN{i-1},'thr')
        if iscell(IN{i})
            thrdata = IN{i};
        else
            thrdata = IN(i);
        end
    elseif strcmp(IN{i-1},'nrho')
        if iscell(IN{i})
            nrhodata = IN{i};
        else
            nrhodata = IN(i);
        end
    elseif strcmp(IN{i-1},'ndipsi')
        if iscell(IN{i})
            ipsidata = IN{i};
        else
            ipsidata = IN(i);
        end
    elseif strcmp(IN{i-1},'ndcontra')
        if iscell(IN{i})
            contradata = IN{i};
        else
            contradata = IN(i);
        end
    elseif strcmp(IN{i-1},'spl')
        if iscell(IN{i})
            spldata = IN{i};
        else
            spldata = IN(i);
        end
    end
end
if exist('ndzdata','var')
    for i = 1:length(ndzdata)
        try
            dataobj = dataset(params.AnNum,[ndzdata{i} '-NTD']);
        catch
            dataobj = dataset(params.AnNum,[ndzdata{i} '-NSPL']);
        end
        datastruct.dsz{i} = struct(dataobj);
        datastruct.dsz{i}.dataobj = dataobj;
    end
end
if exist('ndndata','var')
    for i = 1:length(ndndata)
        dataobj = dataset(params.AnNum,ndndata{i});
        datastruct.dsn{i} = struct(dataobj);
        datastruct.dsn{i}.dataobj = dataobj;
    end
end
if exist('ndpdata','var')
    for i = 1:length(ndpdata)
        dataobj = dataset(params.AnNum,ndpdata{i});
        datastruct.dsp{i} = struct(dataobj);
        datastruct.dsp{i}.dataobj = dataobj;
    end
end
if exist('bfsdata','var')
    for i = 1:length(bfsdata)
        dataobj = dataset(params.AnNum,bfsdata{i});
        datastruct.dsbfs{i} = struct(dataobj);
        datastruct.dsbfs{i}.dataobj = dataobj;
    end
end
if exist('thrdata','var')
    for i = 1:length(thrdata)
        dataobj = dataset(params.AnNum,[thrdata{i} '-THR']);
        datastruct.dsthr{i} = struct(dataobj);
        datastruct.dsthr{i}.dataobj = dataobj;
    end
end
if exist('nrhodata','var')
    for i = 1:length(nrhodata)
        dataobj = dataset(params.AnNum,[nrhodata{i} '-NRHO']);
        datastruct.dsnrho{i} = struct(dataobj);
        datastruct.dsnrho{i}.dataobj = dataobj;
    end
end
if exist('ipsidata','var')
    for i = 1:length(ipsidata)
        try
            dataobj = dataset(params.AnNum,[ipsidata{i} '-NTD']);
        catch
            dataobj = dataset(params.AnNum,[ipsidata{i} '-NSPL']);
        end
        datastruct.dsipsi{i} = struct(dataobj);
        datastruct.dsipsi{i}.dataobj = dataobj;
    end
end
if exist('contradata','var')
    for i = 1:length(contradata)
        try
            dataobj = dataset(params.AnNum,[contradata{i} '-NTD']);
        catch
            dataobj = dataset(params.AnNum,[contradata{i} '-NSPL']);
        end
        datastruct.dscontra{i} = struct(dataobj);
        datastruct.dscontra{i}.dataobj = dataobj;
    end
end
if exist('spldata','var')
    for i = 1:length(spldata)
        dataobj = dataset(params.AnNum,spldata{i});
        datastruct.dsspl{i} = struct(dataobj);
        datastruct.dsspl{i}.dataobj = dataobj;
    end
end

%% Make sure you keep track of ITD conventions here.
%Usually, since we're typically recording on the right - the contra channel
%is 1.
if exist('ndzdata','var')
    params.contrachan = datastruct.dsz{1}.Stimulus.StimParam.stimcntrl.contrachan;
    [~, params.dt] = StimSam(datastruct.dsz{1}.dataobj,1);
    params.dt = 1e-3*params.dt; %us -> ms
    params.gatems = datastruct.dsz{1}.Stimulus.StimParam.indiv.stim{1}.fall;
    params.maxST = datastruct.dsz{1}.Stimulus.Special.BurstDur(1)-params.gatems; %Maximum spike time in ms
elseif exist('ndndata','var')
    params.contrachan = datastruct.dsn{1}.Stimulus.StimParam.stimcntrl.contrachan;
    [~, params.dt] = StimSam(datastruct.dsn{1}.dataobj,1);
    params.dt = 1e-3*params.dt; %us -> ms
    params.gatems = datastruct.dsn{1}.Stimulus.StimParam.indiv.stim{1}.fall;
    params.maxST = datastruct.dsn{1}.Stimulus.Special.BurstDur(1)-params.gatems; %Maximum spike time in ms
elseif exist('ndpdata','var')
    params.contrachan = datastruct.dsp{1}.Stimulus.StimParam.stimcntrl.contrachan;
    [~, params.dt] = StimSam(datastruct.dsp{1}.dataobj,1);
    params.dt = 1e-3*params.dt; %us -> ms
    params.gatems = datastruct.dsp{1}.Stimulus.StimParam.indiv.stim{1}.fall;
    params.maxST = datastruct.dsp{1}.Stimulus.Special.BurstDur(1)-params.gatems; %Maximum spike time in ms
elseif exist('bfsdata','var')
    params.contrachan = datastruct.dsbfs{1}.Stimulus.StimParam.stimcntrl.contrachan;
    [~, params.dt] = StimSam(datastruct.dsbfs{1}.dataobj,1);
    params.dt = 1e-3*params.dt; %us -> ms
    params.gatems = datastruct.dsbfs{1}.Stimulus.StimParam.indiv.stim{1}.fall;
    params.maxST = datastruct.dsbfs{1}.Stimulus.Special.BurstDur(1)-params.gatems; %Maximum spike time in ms
end
params.correctITD = 1;%Flag to apply correction to ITD conventions - this will always be the case in SGSR data since they are always wrong - i.e., SGSR takes the recording side into account for ITDs but applies the wrong convention.
params.minST = 50;%Minimum spike time in ms

if exist('thrdata','var')
    TC = EvalTHR(datastruct.dsthr{1}.dataobj);
end
%% Set some parameters specific for the revcors
if exist('ndzdata','var')%only if those data exist
    params.spec_crit = 2; %Z-score(f) criterion for "significance" in the frequency domain
    params.maxTime = params.maxST;
    params.minTime = params.gatems;
    params.Nboot = 10;
    params.maxSample = ceil(params.maxTime/params.dt);
    params.minSample = ceil(params.minTime/params.dt);
    if exist('TC','var') && TC.fit.cf>1000
        params.Nwin = 1.5*2^10;
    else
        params.Nwin = 2^11;
    end
    params.RCwin = params.Nwin*params.dt;%Length of reverse-correlation window in ms
end
%%
mfilelist = dbstack;
if strcmp(mfilelist(end).name,'Binaural_Wiener') && Nchan==1
    params.saveYN = 0;
    params.dirname = 'C:\Users\Mark\Dropbox\BinRev';
else
    params.saveYN = 1;
    if strcmp(mfilelist(end).name,'Binaural_Wiener')
        params.dirname = 'C:\Users\Mark\Dropbox\BinRev';
    elseif strcmp(mfilelist(end).name,'Monaural_Wiener')
        params.dirname = 'C:\Users\Mark\Dropbox\MonRev';
    elseif strcmp(mfilelist(end).name,'update_AN_Revcors')
        params.dirname = 'C:\Users\Mark\Dropbox\MonRev';
    end
end

