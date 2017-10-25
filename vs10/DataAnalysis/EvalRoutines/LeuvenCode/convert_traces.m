

clear

animal_ID = 'L16552';
pathtofile = 'C:\ExpData\';
list_dir = dir([pathtofile,animal_ID]);

pattern = 'EarlyDS';

for idir=1:length(list_dir)
    
    if (length(strfind(list_dir(idir).name, pattern))~=0)
        name = list_dir(idir).name;
        idx = name(end-12:end-8)
        if str2num(idx)>=129
            new_name = [animal_ID,'_',num2str(str2num(idx)),'_trace'];
            newfolder = [pathtofile,animal_ID,'\',new_name];
            if ~exist(newfolder)
                mkdir(newfolder)
            end
            
            
            
            [newfolder,'\',new_name,'.mat'],idir
            D = read(dataset,animal_ID,str2num(idx));
            
            stim_param = D.stimparam;
            %         Nconds = prod(stim_param.Ncond_XY);
            Nconds =stim_param.Presentation.Ncond;
            Nreps = stim_param.Nrep;
            
            
            if isfield(D.data(1),'RX6_analog_1');
                try
                    trace = anadata(D, 1, 1, 1);
                catch err
                    continue
                end
                nsamples = length(trace);
                for icond=1:Nconds
                    
                    traces=zeros(Nreps,nsamples);
                    for irep=1:Nreps
                        try
                            tmp = anadata(D, 1, icond, irep);
                            if min(size(tmp))<1
                                continue;
                            end
                            if size(tmp,2)>1
                                'there is a problem' %it seems that a single trial has all repetation in matrix
                            end
                            traces(irep,:)=tmp(:,1);
                        catch err %throws an error if interupeted in middle
                            'interrupted'
                        end
                        
                        clear tmp;
                    end
                    save([newfolder,'\',new_name,num2str(icond),'.mat'],'traces','stim_param')
                end
                
                
                
            end
        end
    end
end
