clear;clc;

MSOdatapath = 'C:\Users\Mark\Dropbox\BinRev';

MSOfolderlist = dir(fullfile(MSOdatapath,'*_*'));
nrMSOunits = length(MSOfolderlist);

for i = 1:nrMSOunits
    MSOunitID = MSOfolderlist(i).name;
    if ~isempty(dir(fullfile(MSOdatapath,MSOunitID,'*NDF*')))
        filenames = dir(fullfile(MSOdatapath,MSOunitID,'*NDF*'));
        for j = 1:length(filenames)
            load(fullfile(MSOdatapath,MSOunitID,filenames(j).name));
            %             switch ndf.bd.bdtype
            %                 case 1
            %                     DifcorDF = ndf.difcor.peakhz/1000;
            %                 case 2
            %                     DifcorDF = ndf.positive.peakhz/1000;
            %                 case 3
            %                     DifcorDF = ndf.negative.peakhz/1000;
            %             end
            %             BD_ms = ndf.bd.bd;
            %             BD_cycles = BD_ms*DifcorDF;
            %             totalCD_ms = ndf.bd.CD;
            totalCP_cycles = ndf.bd.CP;
            %             BDest = DifcorDF*(totalCD_ms+(totalCP_cycles/DifcorDF));
            %             BDerr = BD_cycles-BDest;
            %             if abs(BDerr)>0.5
            %                 ndf.bd.CP = ndf.bd.CP+round(BDerr);
            %                 save(fullfile(MSOdatapath,MSOunitID,filenames(j).name),'ndf');
            %             end
            if abs(totalCP_cycles)>0.5
                ndf.bd.CP = wrapToPi(ndf.bd.CP*2*pi)/(2*pi);
                save(fullfile(MSOdatapath,MSOunitID,filenames(j).name),'ndf');
            end
        end
    end
end