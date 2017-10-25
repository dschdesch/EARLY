datadir = 'C:\Users\Mark\Dropbox\MonRev';

folderlist = dir(fullfile(datadir,'L*'));

for i = 1:length(folderlist)
    ID = folderlist(i).name;
    load(fullfile(datadir,ID,[ID '_inputs.mat']));
    AnNum = inputstruct(1);
    UnNum = inputstruct(2);
    if ~isempty(dir(fullfile(datadir,ID,[ID '_MonRev*'])))
        quick_WK_level(AnNum{1},UnNum{1});
    end
end