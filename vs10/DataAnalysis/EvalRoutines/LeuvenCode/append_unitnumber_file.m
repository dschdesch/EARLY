
datapath = 'C:\Users\Mark\Dropbox\BinRev\';

cd(datapath);

unitfolders = dir;
unitfolders = unitfolders(4:end-1);

for i = 1:length(unitfolders)
    unitname = unitfolders(i).name;
    cd(unitname);
    files = dir([datapath unitname '\*.mat']);
    
    for j = 1:length(files)
        oldname = [datapath unitname '\' files(j).name];
        if ~strcmp(files(j).name(1),'L') && ~strcmp(files(j).name(1),'F') 
            newname = [datapath unitname '\' unitname '_' files(j).name];
            movefile(oldname,newname);
        end
    end
    cd(datapath);
end