anapath = 'C:\LeuvenDataAnalysis';
datapath = 'C:\Users\Mark\Dropbox\BinRev';
cd (anapath);
folderlist = dir(fullfile(datapath,'*_*'));
%main loop over all unit folders
for i = 1:length(folderlist)
    unitID = folderlist(i).name;
    if ~isempty(dir(fullfile(datapath,unitID,[unitID, '_*BinRev*'])));
        datafiles = dir(fullfile(datapath,unitID,[unitID, '_*BinRev*']));
        for j = 1:length(datafiles)
            load(fullfile(datapath,unitID,datafiles(j).name));
            if ~isfield(kernels,'h1cdomfreq');
                [~,maxind] = max(kernels.h1cmag);
                h1cdomfreq = kernels.h1ffax(maxind);
                kernels.h1cdomfreq = h1cdomfreq;
                save(fullfile(datapath,unitID,datafiles(j).name),'kernels');
                clear kernels;
            end
        end
    end
    
    
%     if ~isempty(dir(fullfile(datapath,unitID,[unitID, '_*CONTRA*'])));
%         datafiles = dir(fullfile(datapath,unitID,[unitID, '_*CONTRA*']));
%         for j = 1:length(datafiles)
%             load(fullfile(datapath,unitID,datafiles(j).name));
%             [~,maxind] = max(CONTRA.h1cmag);
%             h1cdomfreq = CONTRA.h1ffax(maxind);
%             CONTRA.h1cdomfreq = h1cdomfreq;
%             save(fullfile(datapath,unitID,datafiles(j).name),'CONTRA');
%             clear CONTRA;
%         end
%         clear datafiles;
%     end
%     if ~isempty(dir(fullfile(datapath,unitID,[unitID, '_*IPSI*'])));
%         datafiles = dir(fullfile(datapath,unitID,[unitID, '_*IPSI*']));
%         for j = 1:length(datafiles)
%             load(fullfile(datapath,unitID,datafiles(j).name));
%             [~,maxind] = max(IPSI.h1cmag);
%             h1cdomfreq = IPSI.h1ffax(maxind);
%             IPSI.h1cdomfreq = h1cdomfreq;
%             save(fullfile(datapath,unitID,datafiles(j).name),'IPSI');
%             clear IPSI;
%         end
%         clear datafiles;
%     end
end