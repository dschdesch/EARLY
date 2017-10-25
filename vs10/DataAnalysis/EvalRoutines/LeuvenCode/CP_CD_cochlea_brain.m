clc;clear;

MSOdatapath = 'C:\Users\Mark\Dropbox\BinRev';

MSOfolderlist = dir(fullfile(MSOdatapath,'*_*'));
nrMSOunits = length(MSOfolderlist);
count = 0;
for i = 1:nrMSOunits
    unitID = MSOfolderlist(i).name;
    if ~isempty(dir(fullfile(MSOdatapath,unitID,[unitID '_AN_MSO_match_70.00dB.mat'])));
        count = count+1;
        load(fullfile(MSOdatapath,unitID,[unitID '_AN_MSO_match_70.00dB.mat']));
        load(fullfile(MSOdatapath,unitID,[unitID '_BinRev_70.00dB.mat']));
        load(fullfile(MSOdatapath,unitID,[unitID '_NDF_70.00dB.mat']));
        
        data.tuning_diff_oct(count) = kernels.tuning_diff.ipsivscontra_oct;
        data.difcorDF(count) = ndf.difcor.peakhz/1000;
        data.BD(count) = ndf.bd.bd;
        data.CP(count) = ndf.bd.CP;
        data.CD(count) = ndf.bd.CD;
        data.BD_cycles(count) = data.BD(count)*data.difcorDF(count);
        
        data.CD_cochlea.mean(count) = CD_cochlea.mean;
        data.CD_brain.mean(count) = CD_brain.mean;
        data.CD_cochlea.cycles.mean(count) = CD_cochlea.mean*data.difcorDF(count);
        data.CD_brain.cycles.mean(count) = CD_brain.mean*data.difcorDF(count);
        data.CP_cochlea.mean(count) = CP_cochlea.mean;
        data.CP_brain.mean(count) = CP_brain.mean;
        
        data.CD_cochlea.sem(count) = CD_cochlea.sem;
        data.CD_brain.sem(count) = CD_brain.sem;
        data.CD_cochlea.cycles.sem(count) = CD_cochlea.sem*data.difcorDF(count);
        data.CD_brain.cycles.sem(count) = CD_brain.sem*data.difcorDF(count);
        data.CP_cochlea.sem(count) = CP_cochlea.sem;
        data.CP_brain.sem(count) = CP_brain.sem;
        
        data.CD_cochlea.std(count) = CD_cochlea.std;
        data.CD_brain.std(count) = CD_brain.std;
        data.CD_cochlea.cycles.std(count) = CD_cochlea.std*data.difcorDF(count);
        data.CD_brain.cycles.std(count) = CD_brain.std*data.difcorDF(count);
        data.CP_cochlea.std(count) = CP_cochlea.std;
        data.CP_brain.std(count) = CP_brain.std;
    end
end


figure;
axes;hold on;
plot(data.CD_cochlea.cycles.mean, data.CD_brain.cycles.mean,'ko','markersize',8,'linewidth',2);

% plot([data.CD_cochlea.cycles.mean-(2*data.CD_cochlea.cycles.sem) ;data.CD_cochlea.cycles.mean+(2*data.CD_cochlea.cycles.sem)],...
%     [data.CD_brain.cycles.mean ;data.CD_brain.cycles.mean],'r-','linewidth',1);
% plot([data.CD_cochlea.cycles.mean ;data.CD_cochlea.cycles.mean],...
%     [data.CD_brain.cycles.mean-(2*data.CD_brain.cycles.sem) ;data.CD_brain.cycles.mean+(2*data.CD_brain.cycles.sem)],...
%     'r-','linewidth',1);

figure;
axes;hold on;
plot(data.CP_cochlea.mean, data.CP_brain.mean,'ko','markersize',8,'linewidth',2);
% plot([data.CP_cochlea.mean-(1*data.CP_cochlea.sem) ;data.CP_cochlea.mean+(1*data.CP_cochlea.sem)],...
%     [data.CP_brain.mean ;data.CP_brain.mean],'r-','linewidth',1);
% plot([data.CP_cochlea.mean ;data.CP_cochlea.mean],...
%     [data.CP_brain.mean-(1*data.CP_brain.sem) ;data.CP_brain.mean+(1*data.CP_brain.sem)],...
%     'r-','linewidth',1);

































