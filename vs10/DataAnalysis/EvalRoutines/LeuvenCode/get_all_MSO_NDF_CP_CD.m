datapath = 'C:\Users\Mark\Dropbox\BinRev';
folderlist = dir(fullfile(datapath,'*_*'));
nrunits = length(folderlist);
%main loop over all units
%find all BFS and NDF data files and add estimate of CP and CD
count = 0;
for i = 1:nrunits
    BFSlist = [];NDFlist = [];
    ndf = []; bfs = [];
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
            count = count+1;
            ndf = load(fullfile(datapath,unitID,NDFlist(j).name));
            ndf = ndf.ndf;
            CD(count) = ndf.bd.CD;
            CP(count) = ndf.bd.CP;
            Type(count) = ndf.bd.bdtype;
            unitNum{count} = unitID;
            SPL(count) = ndf.SPL;
            switch Type(count)
                case 1
                    DF(count) = ndf.difcor.peakhz/1000;
                case 2
                    DF(count) = ndf.positive.peakhz/1000;
                case 3
                    DF(count) = ndf.negative.peakhz/1000;
            end
        end
    end
    if ~isempty(BFSlist)
        nrBFSfiles = length(BFSlist);
        for j = 1:nrBFSfiles
            count = count+1;
            bfs = load(fullfile(datapath,unitID,BFSlist(j).name));
            bfs = bfs.bfs;
            CD(count) = bfs.difcor.CD;
            CP(count) = bfs.difcor.CP;
            Type(count) = 4;
            unitNum{count} = unitID;
            SPL(count) = bfs.SPL(1);
            DF(count) = bfs.difcor.peakhz/1000;
        end
    end
end
save('C:\Users\Mark\Dropbox\Disparity_MATLAB\MSO_NDF_CP_CD.mat','CD','CP','Type','unitNum','SPL','DF');