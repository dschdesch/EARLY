function [] = convert_data(animal_ID);

% clear
% 
% animal_ID = 'L16544';
pathtofile = 'C:\ExpData\';
list_dir = dir([pathtofile,animal_ID]);

pattern = 'EarlyDS';

for idir=1:length(list_dir)
    if (length(strfind(list_dir(idir).name, pattern))~=0)
        name = list_dir(idir).name;
        idx = name(end-12:end-8);
        
        
        
        %        if ~exist(newfolder)
        %            mkdir(newfolder)
        %        end
        
        
        
        idx=str2num(idx);
        %         done=0;
        %         while done==0
        %
        %             %because THR sometime screws everything,
        %             %catch error, increment idx until it finds a normal stuff
        %             try
        %                 D = read(dataset,animal_ID,idx);
        %                 stim_param =  D.stimparam;
        %                 spikes = D.spiketimes;
        %                 done=1
        %             catch err
        %                 idx=idx+1
        %                 idir
        %             end
        %         end
        
        new_name = [animal_ID,'_',num2str(idx)];
        newfolder = [pathtofile,animal_ID,'\'];
        newname = [newfolder,new_name,'.mat'];
        
        D = read(dataset,animal_ID,idx);
        stim_param =  D.stimparam;
        
        if (length(strfind(stim_param.StimType, 'THR'))~=0) & (~isobject(D.Data))
            thr = D.Data.Thr;
            freq = stim_param.Presentation.X.PlotVal;
            save(newname,'thr','freq','stim_param');
        elseif (length(strfind(stim_param.StimType, 'THR'))==0)
            spikes = D.spiketimes;
            save(newname,'spikes','stim_param');
        end
        
        
        
        
        
    end
end
