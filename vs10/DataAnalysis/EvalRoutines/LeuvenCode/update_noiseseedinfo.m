function [] = update_noiseseedinfo

ANdata = load('C:\LeuvenDataAnalysis\MonRevPopData.mat');
nrdata = length(ANdata.AN);
for i = 1:nrdata
    load(fullfile('C:\work\BinRev\MonRev',ANdata.AN(i).filename));
    if isfield(ParamsOut,'StimulusInfo')
        ParamsOut = rmfield(ParamsOut,{'noiseseed','SPLDataIDs'});
        
%         count = 0;
%         for j = 1:length(ParamsOut.StimulusInfo)
%             nseed(j) = ParamsOut.StimulusInfo{j}.StimParam.RandomSeed;
%             if isfield(ParamsOut.StimulusInfo{j}.StimParam,'SPL')
%                 count = count+1;
%                 spl(j) = ParamsOut.StimulusInfo{j}.StimParam.SPL(1);
%                 dataID{j} = ParamsOut.DatasetIDs.noise{count};
%             else
%                 if count==0
%                     count=1;
%                 end
%                 allspls = ParamsOut.StimulusInfo{j}.StimParam.startSPL(1):ParamsOut.StimulusInfo{j}.StimParam.stepSPL(1):ParamsOut.StimulusInfo{j}.StimParam.endSPL(1);
%                 spl(j) = allspls(j);
%                 dataID{j} = ParamsOut.DatasetIDs.noise{count};
%             end
%         end
%         uSPLs = ParamsOut.noiseSPLs;
%         for j = 1:length(uSPLs)
%             ParamsOut.noiseseed{j} = nseed(spl==uSPLs(j));
%             ParamsOut.SPLDataIDs{j} = {dataID{spl==uSPLs(j)}};
%         end
        
        save(fullfile('C:\work\BinRev\MonRev',ANdata.AN(i).filename),'ParamsOut');
    end
%     clear count;
%     clear allspls;
%     clear nseed;
%     clear spl;
    clear ParamsOut;
end