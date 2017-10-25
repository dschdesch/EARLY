anapath = 'C:\LeuvenDataAnalysis';
datapath = 'C:\Users\Mark\Dropbox\MonRev';
cd (anapath);
folderlist = dir(fullfile(datapath,'*_*'));
%main loop over all unit folders
for i = 1:length(folderlist)
    unitID = folderlist(i).name;
    clear inputstruct kernels;
    if ~isempty(dir(fullfile(datapath,unitID,[unitID, '_*MonRev*'])));
        datafiles = dir(fullfile(datapath,unitID,[unitID, '_*MonRev*']));
        for j = 1:length(datafiles)
            load(fullfile(datapath,unitID,datafiles(j).name));
            [~,maxind] = max(kernels.h1cmag);
            h1cdomfreq = kernels.h1ffax(maxind);
            kernels.h1cdomfreq = h1cdomfreq;
            save(fullfile(datapath,unitID,datafiles(j).name),'kernels');
            clear kernels;
        end
        clear datafiles;
%         load(fullfile(datapath,unitID,[unitID '_inputs.mat']));
%         Monaural_Wiener(1,0,inputstruct{:});
    end
end