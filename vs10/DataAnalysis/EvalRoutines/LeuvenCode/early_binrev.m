function [] = early_binrev(earlyFLAG,AnNum,UnNum,varargin)
%This is all a big hack job to get the EARLY data processed in the same way as the
%SGSR data. Some not so pretty code here.

%convert data
% convert_data(AnNum);

% varargin is [MOVN],[BFS],[NITD]
inputstruct = cell(1,6);
inputstruct{1} = 'EARLY';
inputstruct{2} = AnNum;
inputstruct{3} = UnNum;
inputstruct{4} = varargin{1};
inputstruct{5} = varargin{2};
inputstruct{6} = varargin{3};

NDF = [];

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
params.MaxOrder = 1;
params.minST = 50;
params.Nboot = 10;
params.spec_crit = 5;
params.AnNum=AnNum;
params.UnNum=UnNum;
params.contrachan = 1;
params.dirname = ['C:\ExpData\' AnNum];
Nin = length(varargin);
savethedata('inputs');
for i = 1:Nin
    if ~isempty(varargin{i})
        for j = 1:length(varargin{i})
            data{j} = load(fullfile(thisandir,[AnNum '_' num2str(varargin{i}(j)) '.mat']));
        end
        datatype = data{1}.stim_param.StimType;
        if strcmp(datatype,'NITD')
            if abs(data{1}.stim_param.Corr)
                rhoflag = 1; %This is a correlated noise NITD dataset
            else
                rhoflag = 0; %This is an uncorrelated noise NITD dataset for revcors
            end
        end
        params.dt = 1e3/data{1}.stim_param.Fsam;
        switch datatype
            case 'MOVN'
                [NDF] = analyse_MOVN(data,params);
                for nn = 1:length(NDF)
                    savethedata('NDF',NDF{nn},NDF{nn}.SPL);
                end
            case 'NITD'
                if rhoflag
                    [NDF] = analyse_MOVN(data,params);
                    savethedata('NDF',NDF,NDF.SPL);
                else
                    [KERNELS,bayesian_model,NRHO,IAP,NDF] = early_get_binrev(data,params,NDF);
                    %quick hack here to update BD estimate if it's ~2pi from
                    %the prediction
                    for ii = 1:length(KERNELS)
                        kernelspls(ii) = KERNELS{ii}.SPL;
                    end
                    for ii = 1:length(NDF)
                        ndfspls(ii) = NDF{ii}.SPL;
                    end
                    for ii = 1:length(kernelspls)
                        ndfind = find(ndfspls==kernelspls(ii));
                        if ~isempty(ndfind)
                            if isfield(NDF{ndfind},'difcor');
                                [NDF{ndfind}] = early_updateBD(NDF{ndfind},KERNELS{ii});
                            end
                        end
                    end
                    for nn = 1:length(KERNELS)
                        savethedata('kernels',KERNELS{nn},KERNELS{nn}.SPL);
                        if ~isempty(bayesian_model{nn})
                            savethedata('BM',bayesian_model{nn},KERNELS{nn}.SPL);
                        end
                        if ~isempty(NRHO{nn})
                            savethedata('NRHO',NRHO{nn},KERNELS{nn}.SPL);
                        end
                        if ~isempty(IAP{nn})
                            savethedata('IAP',IAP{nn},KERNELS{nn}.SPL);
                        end
                    end
                    for nn = 1:length(NDF)
                        savethedata('NDF',NDF{nn},NDF{nn}.SPL);
                    end
                end
            case 'BBFC'
                if ~isempty(NDF)
                    [BFS] = analyse_BBFC(data,params,NDF{1}.xvals.ITDx);
                else
                    [BFS] = analyse_BBFC(data,params,-5:0.1:5);
                end
                for nn = 1:length(BFS)
                    savethedata('BFS',BFS{nn},BFS{nn}.SPL);
                end
        end
        clear data;
    end
end
%% Save
    function [] = savethedata(datatype,varargin)
        % Save the data structure
        UnNum = params.UnNum;
        AnNum = params.AnNum;
        savedirname = fullfile(params.dirname,[AnNum '_' UnNum]);
        if ~exist(savedirname,'dir')
            mkdir(savedirname);
        end
        switch datatype
            case 'inputs'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_inputs.mat']);
                save (savedfilename,'inputstruct');
                disp (['Saved: ' savedfilename]);
            case 'IAP'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('IAP_%2.2fdB.mat',varargin{2})]);
                iap = varargin{1};
                save (savedfilename,'iap');
                disp (['Saved: ' savedfilename]);
            case 'BM'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('BM_%2.2fdB.mat',varargin{2})]);
                bm = varargin{1};
                save (savedfilename,'bm');
                disp (['Saved: ' savedfilename]);
            case 'kernels'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('BinRev_%2.2fdB.mat',varargin{2})]);
                kernels = varargin{1};
                save (savedfilename,'kernels');
                disp (['Saved: ' savedfilename]);
            case 'NDF'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('NDF_%2.2fdB.mat',varargin{2})]);
                ndf = varargin{1};
                save (savedfilename,'ndf');
                disp (['Saved: ' savedfilename]);
                clear ndf;
            case 'BFS'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('BFS_%2.2fdB.mat',varargin{2})]);
                bfs = varargin{1};
                save (savedfilename,'bfs');
                disp (['Saved: ' savedfilename]);
                clear bfs;
            case 'NRHO'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('NRHO_%2.2fdB.mat',varargin{2})]);
                nrho = varargin{1};
                save (savedfilename,'nrho');
                disp (['Saved: ' savedfilename]);
                clear nrho;
        end
    end
end