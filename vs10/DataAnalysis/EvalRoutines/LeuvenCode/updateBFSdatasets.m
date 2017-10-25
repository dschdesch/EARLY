function [] = updateBFSdatasets
%Applies an update to BFS data structures
%Add monaural estimates of phase locking
%23/6/16

MSOdatapath = 'C:\Users\Mark\Dropbox\BinRev';
MSOfolderlist = dir(fullfile(MSOdatapath,'*_*'));
params.dirname = 'C:\Users\Mark\Dropbox\BinRev';
for i = 1:length(MSOfolderlist)
    MSOunitID = MSOfolderlist(i).name;
    if ~isempty(dir(fullfile(MSOdatapath,MSOunitID,'*BFS*')))
        clear inputstruct bfs datastruct ITDx;
        filenames = dir(fullfile(MSOdatapath,MSOunitID,'*BFS*'));
        inputfile = dir(fullfile(MSOdatapath,MSOunitID,'*inputs*'));
        load(fullfile(MSOdatapath,MSOunitID,inputfile.name));
        load(fullfile(MSOdatapath,MSOunitID,filenames(1).name));
        if ~strcmp(inputstruct{1},'EARLY')
            AnNum = inputstruct{1};
            params.AnNum = AnNum;
            UnNum = inputstruct{2};
            params.UnNum = UnNum;
%             if strcmp('L15007',AnNum) && strcmp('7',UnNum)
                
                for j = 3:length(inputstruct)
                    if strcmp(inputstruct{j-1},'bfs')
                        if iscell(inputstruct{j})
                            bfsdata = inputstruct{j};
                        else
                            bfsdata = inputstruct(j);
                        end
                    end
                end
                if exist('bfsdata','var')
                    for j = 1:length(bfsdata)
                        dataobj = dataset(params.AnNum,bfsdata{j});
                        datastruct.dsbfs{j} = struct(dataobj);
                        datastruct.dsbfs{j}.dataobj = dataobj;
                    end
                end
                params.contrachan = datastruct.dsbfs{1}.Stimulus.StimParam.stimcntrl.contrachan;
                [~, params.dt] = StimSam(datastruct.dsbfs{1}.dataobj,1);
                params.dt = 1e-3*params.dt; %us -> ms
                params.gatems = datastruct.dsbfs{1}.Stimulus.StimParam.indiv.stim{1}.fall;
                params.maxST = datastruct.dsbfs{1}.Stimulus.Special.BurstDur(1)-params.gatems; %Maximum spike time in ms
                params.correctITD = 1;%Flag to apply correction to ITD conventions - this will always be the case in SGSR data since they are always wrong - i.e., SGSR takes the recording side into account for ITDs but applies the wrong convention.
                params.minST = 50;%Minimum spike time in ms
                
                ITDx = bfs.xvals.ITDx;
                [BFS] = bfs_analysis(datastruct.dsbfs,params,ITDx);
                for nn = 1:length(BFS)
                    savethedata('BFS',BFS{nn},BFS{nn}.SPL);
                end
%             end
        else
            disp('EARLY dataset... skipping...');
        end
    end
    clear data;clear BFS;
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