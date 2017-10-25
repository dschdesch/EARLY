function [] = get_CP_CD_from_NTD

datapath = 'C:\Users\Mark\Dropbox\BinRev';
folderlist = dir(fullfile(datapath,'*_*'));
nrunits = length(folderlist);
%main loop over all units
%find all BFS and NDF data files and add estimate of CP and CD
for i = 1:nrunits
    BFSlist = [];NDFlist = [];
    unitID = folderlist(i).name;
    if ~isempty(dir(fullfile(datapath,unitID,[unitID '_BFS*.mat'])))
        BFSlist = dir(fullfile(datapath,unitID,[unitID '_BFS*.mat']));
    end
    if ~isempty(dir(fullfile(datapath,unitID,[unitID '_NDF*.mat'])))
        NDFlist = dir(fullfile(datapath,unitID,[unitID '_NDF*.mat']));
    end
    if ~isempty(NDFlist)
        nrNDFfiles = length(NDFlist);
        for j = 1:nrNDFfiles
            ndf = load(fullfile(datapath,unitID,NDFlist(j).name));
            ndf = ndf.ndf;
            if isfield(ndf,'difcor')
                [CP,CD] = local_CP_CD(ndf,1);
            elseif isfield(ndf,'positive')
                [CP,CD,pseudodifcor] = local_CP_CD(ndf,3);
                ndf.pseudodifcor.rate = pseudodifcor;
            elseif isfield(ndf,'negative')
                [CP,CD,pseudodifcor] = local_CP_CD(ndf,4);
                ndf.pseudodifcor.rate = pseudodifcor;
            end
            ndf.bd.CD = CD;
            ndf.bd.CP = CP;
            save(fullfile(datapath,unitID,NDFlist(j).name),'ndf');
        end
    end
    if ~isempty(BFSlist)
        nrBFSfiles = length(BFSlist);
        for j = 1:nrBFSfiles
            bfs = load(fullfile(datapath,unitID,BFSlist(j).name));
            bfs = bfs.bfs;
            [CP,CD] = local_CP_CD(bfs,2);
            bfs.difcor.CD = CD;
            bfs.difcor.CP = CP;
            save(fullfile(datapath,unitID,BFSlist(j).name),'bfs');
        end
    end
end
    function [CP,CD,difcor] = local_CP_CD(IN,flag)
        switch flag
            case 1
                difcor = IN.difcor.rate.spline.mean;
                spec = fft(difcor,2^13);
                freq = IN.difcor.freq/1000;
            case 2
                difcor = IN.difcor.ratespline;
                spec = fft(difcor,2^13);
                freq = IN.difcor.freq/1000;
            case 3
                Pos = IN.positive.rate.spline.mean;
                Neg = mean(Pos)-(Pos-mean(Pos));
                Neg = max(Neg,0);
                difcor = Pos-Neg;
                spec = fft(difcor,2^13);
                freq = IN.positive.freq/1000;
            case 4
                Neg = IN.negative.rate.spline.mean;
                Pos = mean(Neg)-(Neg-mean(Neg));
                Pos = max(Pos,0);
                difcor = Pos-Neg;
                spec = fft(difcor,2^13);
                freq = IN.negative.freq/1000;
        end
        spec = spec(1:4096);
        mag = 20*log10(abs(spec));
        mag = mag-max(mag);
        mag(mag<-10)=nan;
        weight = 10.^(mag/20);
        phase = angle(spec)/(2*pi);
        weight(isnan(weight)) = 0;
        weight(freq<0.08) = 0;
        [CP,CD] = fit_circ_lin(freq(weight>0),phase(weight>0)',weight(weight>0));
        correctionFx = abs(min(IN.xvals.ITDxspline));
        CD = CD-correctionFx;
    end
end