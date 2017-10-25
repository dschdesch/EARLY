function [] = early_updateBFSdatasets
%Applies an update to BFS data structures
%Add monaural estimates of phase locking
%23/6/16

MSOdatapath = 'C:\Users\Mark\Dropbox\BinRev';
MSOfolderlist = dir(fullfile(MSOdatapath,'*_*'));
params.dirname = 'C:\Users\Mark\Dropbox\BinRev';
istart = 250;
for i = istart:length(MSOfolderlist)
    MSOunitID = MSOfolderlist(i).name;
    if ~isempty(dir(fullfile(MSOdatapath,MSOunitID,'*BFS*')))
        clear inputstruct bfs datastruct ITDx;
        filenames = dir(fullfile(MSOdatapath,MSOunitID,'*BFS*'));
        inputfile = dir(fullfile(MSOdatapath,MSOunitID,'*inputs*'));
        load(fullfile(MSOdatapath,MSOunitID,inputfile.name));
        load(fullfile(MSOdatapath,MSOunitID,filenames(1).name));
        if strcmp(inputstruct{1},'EARLY')
            AnNum = inputstruct{2};
            params.AnNum = AnNum;
            UnNum = inputstruct{3};
            params.UnNum = UnNum;
            earlydatadir = 'C:\ExpData';
            if exist(fullfile(earlydatadir,AnNum),'dir')
                thisandir = fullfile(earlydatadir,AnNum);
            elseif exist(fullfile(earlydatadir,'Mark',AnNum),'dir')
                thisandir = fullfile(earlydatadir,'Mark',AnNum);
            elseif exist(fullfile(earlydatadir,'Philip',AnNum),'dir')
                thisandir = fullfile(earlydatadir,'Philip',AnNum);
            else
                error ('Can''t find the data');
            end
            for j = 1:length(inputstruct{5})
                data{j} = load(fullfile(thisandir,[AnNum '_' num2str(inputstruct{5}(j)) '.mat']));
            end
            params.dt = 1e3/data{1}.stim_param.Fsam;
            ITDx = bfs.xvals.ITDx;
            [BFS] = analyse_BBFC(data,params,ITDx);
            for nn = 1:length(BFS)
                savethedata('BFS',BFS{nn},BFS{nn}.SPL);
            end
        else
            disp('SGSR dataset... skipping...');
        end
    end
    clear data;
    pause(.1);
    i
end

    function [] = savethedata(datatype,varargin)
        % Save the data structure
        UnNum = params.UnNum;
        AnNum = params.AnNum;
        savedirname = fullfile(params.dirname,[AnNum '_' UnNum]);
        if ~exist(savedirname,'dir')
            mkdir(savedirname);
        end
        
        savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('BFS_%2.2fdB.mat',varargin{2})]);
        bfs = varargin{1};
        save (savedfilename,'bfs');
        disp (['Saved: ' savedfilename]);
        clear bfs;
    end
end