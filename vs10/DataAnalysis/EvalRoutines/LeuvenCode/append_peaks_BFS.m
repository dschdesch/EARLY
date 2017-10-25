function [] = append_peaks_BFS

datapath = 'C:\Users\Mark\Dropbox\BinRev\';

cd(datapath);

unitfolders = dir;
unitfolders = unitfolders(4:end-1);

for i = 1:length(unitfolders)
    unitname = unitfolders(i).name;
    cd(unitname);
    files = dir([datapath unitname '\*BFS*.mat']);
    for j = 1:length(files)
        bfs = [];
        load(files(j).name);
        if ~isfield(bfs.difcor,'peaks') || ~isfield(bfs.difcor,'troughs')
            [bfs] = add_extra_peaks(bfs);
            save(files(j).name,'bfs');
            fprintf(['Updated ' files(j).name '\n']);
        end
    end
    cd(datapath);
end

    function [data] = add_extra_peaks(data)
        BD = data.difcor.bd;
        Xspline = data.xvals.ITDxspline;
        vals = data.difcor.ratespline-min(data.difcor.ratespline);
        [peakval,peakind] = findpeaks(vals);
        [troughval,troughind] = findpeaks(-vals);
        troughval = -troughval;
        if min(peakind)<min(troughind)
            troughind = [1 troughind];
            troughval = [vals(1) troughval];
        end
        if max(peakind)>max(troughind)
            troughind = [troughind length(vals)];
            troughval = [troughval vals(end)];
        end
        Dpeaks = Xspline(peakind);
        Dtroughs = Xspline(troughind);
        keeppeak = zeros(size(peakval));
        keeptrough = zeros(size(troughval));
        for ii = 1:length(peakval)
            clo = (peakval(ii)-troughval(find(troughind<peakind(ii),1,'last')))/...
                (peakval(ii)+troughval(find(troughind<peakind(ii),1,'last')));
            chi = (peakval(ii)-troughval(find(troughind>peakind(ii),1,'first')))/...
                (peakval(ii)+troughval(find(troughind>peakind(ii),1,'first')));
            if clo>0.5
                keeppeak(ii) = 1;
                keeptrough(find(troughind<peakind(ii),1,'last'))=1;
            elseif chi>0.5
                keeppeak(ii) = 1;
                keeptrough(find(troughind>peakind(ii),1,'first'))=1;
            end
        end
        Dpeaks = Dpeaks(logical(keeppeak));
        data.difcor.troughs = Dtroughs(logical(keeptrough));
        data.difcor.peaks = Dpeaks(Dpeaks~=BD);
    end
end