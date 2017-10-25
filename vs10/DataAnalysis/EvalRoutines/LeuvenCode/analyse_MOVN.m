function [NDF] = analyse_MOVN(data,params)


NFFT = 2^13;

datatype = data{1}.stim_param.StimType;

switch datatype %This is a hack to deal with some old early data with static NITD curves
    case 'MOVN'
        %First find out which SPLs we have
        for i = 1:length(data)
            allSPLs(i) = data{i}.stim_param.SPL;
        end
        %find unique SPLs
        [uSPLs] = unique(allSPLs);
        nrSPLs = length(uSPLs);
        
        for SPLind = 1:nrSPLs %main loop over SPL conditions
            
            thesedata = data(allSPLs == uSPLs(SPLind));
            
            minITD = thesedata{1}.stim_param.ITD1;
            maxITD = thesedata{1}.stim_param.ITD2;
            
            ITDspeed = thesedata{1}.stim_param.startSpeed;
            
            dt = params.dt;
            burstdur = (maxITD-minITD)/ITDspeed*1e3;
            repdur = thesedata{1}.stim_param.ISI;
            SPL = thesedata{1}.stim_param.SPL;
            
            for idx = 1:length(thesedata)
                rho(idx) = thesedata{idx}.stim_param.Corr;
            end
            
            for rhodx = 1:length(thesedata)
                [nrcond,nrreps] = size(thesedata{rhodx}.spikes);
                nrcond=1;
                count(rhodx) = 0;
                for idx = 1:nrcond
                    for jdx = 1:nrreps
                        if ~isempty(thesedata{rhodx}.spikes{idx,jdx})
                            count(rhodx) = count(rhodx)+1;
                        end
                    end
                end
                totaldur(rhodx) = burstdur*count(rhodx);
            end
            
            for rhodx = 1:length(thesedata)
                ST = cat(2,thesedata{rhodx}.spikes{:});
                ST = ST(ST<=burstdur);
                
                %make spike times into units of ITD
                SITD = 1e-3*((ST/burstdur)*(maxITD-minITD)-maxITD);
                
                %correction is needed here: ITDs in MOVN are the wrong convention
                SITD = SITD*-1;
                
                if strcmp(params.AnNum,'F15075')%This animal used the recording side in set-up, so we need to re-correct here.
                    SITD = -SITD;
                end
                
                %xvals - to match typical resolution used in SGSR datasets
                ITDx = linspace(minITD,maxITD,41)/1000;
                nrbins = length(ITDx);
                %yvals
                yvals{rhodx} = hist(SITD,ITDx);
                %correct for edge effects
                yvals{rhodx}(1) = yvals{rhodx}(1)*2;
                yvals{rhodx}(end) = yvals{rhodx}(end)*2;
                
                spikerate{rhodx} = yvals{rhodx}/((totaldur(rhodx)/1000)/nrbins);
            end
            if ~isempty(find(rho==1,1))
                if numel(find(rho==1))>1
                    NDF{SPLind}.positive.rate.mean = mean(vertcat(spikerate{rho==1}))';
                else
                    NDF{SPLind}.positive.rate.mean = spikerate{rho==1}';
                end
            end
            if ~isempty(find(rho==-1,1))
                if numel(find(rho==-1))>1
                    NDF{SPLind}.negative.rate.mean = mean(vertcat(spikerate{rho==-1}))';
                else
                    NDF{SPLind}.negative.rate.mean = spikerate{rho==-1}';
                end
            end
        end
    case 'NITD'
        SPL = data{1}.stim_param.SPL;
        for rhodx = 1:length(data)
            rho(rhodx) = data{rhodx}.stim_param.Corr;
            burstdur = data{rhodx}.stim_param.BurstDur;
            ITDx = data{rhodx}.stim_param.ITD;
            [nrcond,nrreps] = size(data{rhodx}.spikes);
            for i = 1:nrcond
                ST = [data{rhodx}.spikes{i,:}];
                ST = ST<=burstdur;
                spikerate{rhodx}(i) = length(ST)/(nrreps*(burstdur/1000));
            end
        end
        
        if ~isempty(find(rho==1,1))
            if numel(find(rho==1))>1
                NDF.positive.rate.mean = mean(vertcat(spikerate{rho==1}))';
            else
                NDF.positive.rate.mean = spikerate{rho==1}';
            end
        end
        if ~isempty(find(rho==-1,1))
            if numel(find(rho==-1))>1
                NDF.negative.rate.mean = mean(vertcat(spikerate{rho==-1}))';
            else
                NDF.negative.rate.mean = spikerate{rho==-1}';
            end
        end
end
for SPLind = 1:nrSPLs
    if isfield(NDF{SPLind},'positive') && isfield(NDF{SPLind},'negative')
        NDF{SPLind}.difcor.rate.mean = NDF{SPLind}.positive.rate.mean-NDF{SPLind}.negative.rate.mean;
    end
    NDF{SPLind}.xvals.ITDx = ITDx;
end

switch datatype
    case 'MOVN'
        for SPLind = 1:nrSPLs
            NDF{SPLind}.xvals.ITDxspline = minITD/1000:dt:maxITD/1000;
        end
    case 'NITD'
        NDF.xvals.ITDxspline = linspace(min(ITDx),max(ITDx),1000);
end

for SPLind = 1:nrSPLs
    if isfield(NDF{SPLind},'positive')
        NDF{SPLind}.positive.rate.spline.mean = spline(NDF{SPLind}.xvals.ITDx,NDF{SPLind}.positive.rate.mean,NDF{SPLind}.xvals.ITDxspline);
    end
    if isfield(NDF{SPLind},'negative')
        NDF{SPLind}.negative.rate.spline.mean = spline(NDF{SPLind}.xvals.ITDx,NDF{SPLind}.negative.rate.mean,NDF{SPLind}.xvals.ITDxspline);
    end
    if isfield(NDF{SPLind},'positive') && isfield(NDF{SPLind},'negative')
        NDF{SPLind}.difcor.rate.spline.mean = spline(NDF{SPLind}.xvals.ITDx,NDF{SPLind}.difcor.rate.mean,NDF{SPLind}.xvals.ITDxspline);
    end
    
    %Get the BD
    if isfield(NDF{SPLind},'positive') && isfield(NDF{SPLind},'negative')
        %You have the difcor - therefore BD is just the peak
        BDtype = 1;
        [dummy,ind] = max(NDF{SPLind}.difcor.rate.spline.mean);
    elseif isfield(NDF{SPLind},'positive')
        %You only have the positive correlation NDF - take the peak as BD
        BDtype = 2;
        [dummy,ind] = max(NDF{SPLind}.positive.rate.spline.mean);
    elseif isfield(NDF{SPLind},'negative')
        %You only have the negative correlation NDF - take the minimum as "BD"
        BDtype = 3;
        [dummy,ind] = min(NDF{SPLind}.negative.rate.spline.mean);%Probably this will never be used... but just here for completeness
    end
    BD = NDF{SPLind}.xvals.ITDxspline(ind);
    NDF{SPLind}.bd.bdtype = BDtype;
    NDF{SPLind}.bd.bd = BD;
    NDF{SPLind}.SPL = uSPLs(SPLind);
    
    %Find other peaks and valleys
    if isfield(NDF{SPLind},'positive') && isfield(NDF{SPLind},'negative')
        vals = NDF{SPLind}.difcor.rate.spline.mean-min(NDF{SPLind}.difcor.rate.spline.mean);
        [peakval,peakind] = findpeaks(vals);
        
        [troughval,troughind] = findpeaks(-vals);
        troughval = -troughval;
        
        if min(peakind)<min(troughind)
            troughind = [1 troughind];
            troughval = [vals(1) troughval];
        end
        if max(peakind)>max(troughind)
            troughind = [troughind length(NDF{SPLind}.difcor.rate.spline.mean)];
            troughval = [troughval vals(length(NDF{SPLind}.difcor.rate.spline.mean))];
        end
        Dpeaks = NDF{SPLind}.xvals.ITDxspline(peakind);
        Dtroughs = NDF{SPLind}.xvals.ITDxspline(troughind);
        keeppeak = zeros(size(peakval));
        keeptrough = zeros(size(troughval));
        for i = 1:length(peakval)
            clo = (peakval(i)-troughval(find(troughind<peakind(i),1,'last')))/...
                (peakval(i)+troughval(find(troughind<peakind(i),1,'last')));
            chi = (peakval(i)-troughval(find(troughind>peakind(i),1,'first')))/...
                (peakval(i)+troughval(find(troughind>peakind(i),1,'first')));
            if clo>0.5
                keeppeak(i) = 1;
                keeptrough(find(troughind<peakind(i),1,'last'))=1;
            elseif chi>0.5
                keeppeak(i) = 1;
                keeptrough(find(troughind>peakind(i),1,'first'))=1;
            end
        end
        Dpeaks = Dpeaks(logical(keeppeak));
        Dtroughs = Dtroughs(logical(keeptrough));
        Dpeaks = Dpeaks(Dpeaks~=BD);
    elseif isfield(NDF{SPLind},'positive')
        vals = NDF{SPLind}.positive.rate.spline.mean;
        [peakval,peakind] = findpeaks(vals);
        Ppeaks = NDF{SPLind}.xvals.ITDxspline(peakind);
        [troughval,troughind] = findpeaks(-vals);
        troughval = -troughval;
        Ptroughs = NDF{SPLind}.xvals.ITDxspline(troughind);
        if min(peakind)<min(troughind)
            troughind = [1 troughind];
            troughval = [vals(1) troughval];
        end
        if max(peakind)>max(troughind)
            troughind = [troughind length(NDF{SPLind}.positive.rate.spline.mean)];
            troughval = [troughval vals(length(NDF{SPLind}.positive.rate.spline.mean))];
        end
        keeppeak = zeros(size(peakval));
        keeptrough = zeros(size(troughval));
        for i = 1:length(peakval)
            clo = (peakval(i)-troughval(find(troughind<peakind(i),1,'last')))/...
                (peakval(i)+troughval(find(troughind<peakind(i),1,'last')));
            chi = (peakval(i)-troughval(find(troughind>peakind(i),1,'first')))/...
                (peakval(i)+troughval(find(troughind>peakind(i),1,'first')));
            if clo>0.5
                keeppeak(i) = 1;
                keeptrough(find(troughind<peakind(i),1,'last'))=1;
            elseif chi>0.5
                keeppeak(i) = 1;
                keeptrough(find(troughind>peakind(i),1,'first'))=1;
            end
        end
        Ppeaks = Ppeaks(logical(keeppeak));
        Ptroughs = Ptroughs(logical(keeptrough));
        Ppeaks = Ppeaks(Ppeaks~=BD);
    elseif isfield(NDF{SPLind},'negative')
        vals = NDF{SPLind}.negative.rate.spline.mean;
        [peakval,peakind] = findpeaks(vals);
        Npeaks = NDF{SPLind}.xvals.ITDxspline(peakind);
        [troughval,troughind] = findpeaks(-vals);
        troughval = -troughval;
        Ntroughs = NDF{SPLind}.xvals.ITDxspline(troughind);
        if min(peakind)<min(troughind)
            troughind = [1 troughind];
            troughval = [vals(1) troughval];
        end
        if max(peakind)>max(troughind)
            troughind = [troughind length(NDF{SPLind}.negative.rate.spline.mean)];
            troughval = [troughval vals(length(NDF{SPLind}.negative.rate.spline.mean))];
        end
        keeppeak = zeros(size(peakval));
        keeptrough = zeros(size(troughval));
        for i = 1:length(peakval)
            clo = (peakval(i)-troughval(find(troughind<peakind(i),1,'last')))/...
                (peakval(i)+troughval(find(troughind<peakind(i),1,'last')));
            chi = (peakval(i)-troughval(find(troughind>peakind(i),1,'first')))/...
                (peakval(i)+troughval(find(troughind>peakind(i),1,'first')));
            if clo>0.5
                keeppeak(i) = 1;
                keeptrough(find(troughind<peakind(i),1,'last'))=1;
            elseif chi>0.5
                keeppeak(i) = 1;
                keeptrough(find(troughind>peakind(i),1,'first'))=1;
            end
        end
        Npeaks = Npeaks(logical(keeppeak));
        Ntroughs = Ntroughs(logical(keeptrough));
        Ntroughs = Ntroughs(Ntroughs~=BD);
    end
    
    
    
    
    
    
    if isfield(NDF{SPLind},'positive') && isfield(NDF{SPLind},'negative')
        NDF{SPLind}.bd.dpeaks = Dpeaks;
        NDF{SPLind}.bd.dtroughs = Dtroughs;
    elseif isfield(NDF{SPLind},'positive')
        NDF{SPLind}.bd.ppeaks = Ppeaks;
        NDF{SPLind}.bd.ptroughs = Ptroughs;
    elseif isfield(NDF{SPLind},'negative')
        NDF{SPLind}.bd.npeaks = Npeaks;
        NDF{SPLind}.bd.ntroughs = Ntroughs;
    end
    
    if strcmp(datatype,'NITD')
        dt = diff(NDF{SPLind}.xvals.ITDxspline([1 2]));
    end
    
    %Get the frequency domain signals
    nf = 0.5/(dt/1000);
    freq = linspace(0,1,NFFT/2)*nf;
    if isfield(NDF{SPLind},'positive') && isfield(NDF{SPLind},'negative')
        vals = NDF{SPLind}.difcor.rate.spline.mean;
        hanwin = hann(length(vals));
        spec = fft((vals'-mean(vals)).*hanwin,NFFT);
        Dphase = angle(spec(1:NFFT/2));
        Dspec = 20*log10(abs(spec(1:NFFT/2)));
        Dspec = Dspec-max(Dspec);
        Dphase(Dspec<-20) = nan;
        [dummy,ind] = max(Dspec);
        DspecPeakHz = freq(ind);
    end
    if isfield(NDF{SPLind},'negative')
        vals = NDF{SPLind}.negative.rate.spline.mean;
        hanwin = hann(length(vals));
        spec = fft((vals'-mean(vals)).*hanwin,NFFT);
        Nphase = angle(spec(1:NFFT/2));
        Nspec = 20*log10(abs(spec(1:NFFT/2)));
        Nspec = Nspec-max(Nspec);
        Nphase(Nspec<-20) = nan;
        [dummy,ind] = max(Nspec);
        NspecPeakHz = freq(ind);
    end
    if isfield(NDF{SPLind},'positive')
        vals = NDF{SPLind}.positive.rate.spline.mean;
        hanwin = hann(length(vals));
        spec = fft((vals'-mean(vals)).*hanwin,NFFT);
        Pphase = angle(spec(1:NFFT/2));
        Pspec = 20*log10(abs(spec(1:NFFT/2)));
        Pspec = Pspec-max(Pspec);
        Pphase(Pspec<-20) = nan;
        [dummy,ind] = max(Pspec);
        PspecPeakHz = freq(ind);
    end
    
    
    
    if isfield(NDF{SPLind},'positive')
        NDF{SPLind}.positive.power = Pspec;
        NDF{SPLind}.positive.phase = Pphase;
        NDF{SPLind}.positive.peakhz = PspecPeakHz;
        NDF{SPLind}.positive.freq = freq;
    end
    if isfield(NDF{SPLind},'negative')
        NDF{SPLind}.negative.power = Nspec;
        NDF{SPLind}.negative.phase = Nphase;
        NDF{SPLind}.negative.peakhz = NspecPeakHz;
        NDF{SPLind}.negative.freq = freq;
    end
    if isfield(NDF{SPLind},'positive') && isfield(NDF{SPLind},'negative')
        NDF{SPLind}.difcor.power = Dspec;
        NDF{SPLind}.difcor.phase = Dphase;
        NDF{SPLind}.difcor.peakhz = DspecPeakHz;
        NDF{SPLind}.difcor.freq = freq;
    end
end
return;