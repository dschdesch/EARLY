function [] = converth1h2data

load('MonRevPopData.mat');
nrfiles = length(AN);
for i = 594:nrfiles
    DatasetIDs = AN(i).DatasetIDs;
    ind = regexp(AN(i).filename,'_');
    unitnum = AN(i).filename(ind(1)+1:ind(2)-1);
    dirtosave = [DatasetIDs.Animal '_' unitnum];
    if isfield(DatasetIDs,'noise') && isfield(DatasetIDs,'spl')
        Monaural_Wiener(0,2,DatasetIDs.Animal,DatasetIDs.noise,'thr',DatasetIDs.thr,'spl',DatasetIDs.spl,'savedirname',dirtosave);
    elseif isfield(DatasetIDs,'noise')
        Monaural_Wiener(0,2,DatasetIDs.Animal,DatasetIDs.noise,'thr',DatasetIDs.thr,'savedirname',dirtosave);
    else
        Monaural_Wiener(0,2,DatasetIDs.Animal,'thr',DatasetIDs.thr,'savedirname',dirtosave);
    end
    fprintf('Completed %1.0f of %1.0f files...\n',i,nrfiles);
end