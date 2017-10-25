function [NDF] = NDF_analysis(dsp,dsn,params)

%returns the noise-delay functions, difcor, and their spline fits

%some useful parameters
NFFT = 2^13;
dt = params.dt;
maxST = params.maxST;
minST = params.minST;
stepST = (maxST-minST)/10;

%If there is more than one dataset - are they the same SPL or different?
if ~isempty(dsp)
    nrdsp = length(dsp);
    SPLp = zeros(1,nrdsp);
    for i = 1:nrdsp
        SPLp(i) = dsp{i}.Stimulus.StimParam.SPL(1);
    end
else
    SPLp = [];
end
if ~isempty(dsn)
    nrdsn = length(dsn);
    SPLn = zeros(1,nrdsn);
    for i = 1:nrdsn
        SPLn(i) = dsn{i}.Stimulus.StimParam.SPL(1);
    end
else
    SPLn = [];
end
%All unique SPLs presented
uSPLs = unique([SPLn SPLp]);
nrSPLs = length(uSPLs);

%Main loop over all unique SPL conditions
for j = 1:nrSPLs
    
    %indices for positive and negative at this SPL
    if ~isempty(SPLp)
        possplind = find(SPLp==uSPLs(j));
    else
        possplind =[];
    end
    if ~isempty(SPLn)
        negsplind = find(SPLn==uSPLs(j));
    else
        negsplind = [];
    end
    
    %Which noise delays were played?
    %Get list of all delay values
    if ~isempty(possplind)
        Pxvals = cell(1,length(possplind));
        for i = 1:length(possplind)
            Pxvals{i} = dsp{possplind(i)}.Stimulus.IndepVar.Values;
        end
    else
        Pxvals = [];
    end
    if ~isempty(negsplind)
        Nxvals = cell(1,length(negsplind));
        for i = 1:length(negsplind)
            Nxvals{i} = dsn{negsplind(i)}.Stimulus.IndepVar.Values;
        end
    else
        Nxvals = [];
    end
    %Find out if any of them were missed. If so, discard the last one recorded
    %(likely these data are compromised).
    if ~isempty(Pxvals)
        for i = 1:length(Pxvals)
            if any(isnan(Pxvals{i}))
                ind = find(isnan(Pxvals{i}),1,'first');
                Pxvals{i} = Pxvals{i}(1:ind-2);
            end
        end
    end
    if ~isempty(Nxvals)
        for i = 1:length(Nxvals)
            if any(isnan(Nxvals{i}))
                ind = find(isnan(Nxvals{i}),1,'first');
                Nxvals{i} = Nxvals{i}(1:ind-2);
            end
        end
    end
    
    
    %Make one vector of ITD values, based on all available data
    if ~isempty(Pxvals) && ~isempty(Nxvals)
        [theseNxvals,Nindsort] = sort(Nxvals{1},'ascend');
        [thesePxvals,Pindsort] = sort(Pxvals{1},'ascend');
        [Xvals,Pind,Nind] = intersect(round(thesePxvals),round(theseNxvals));
        Pind = Pindsort(Pind);
        Nind = Nindsort(Nind);
    elseif ~isempty(Pxvals)
        Xvals = Pxvals{1};
        Pind = 1:length(Pxvals{1});
    elseif ~isempty(Nxvals)
        Xvals = Nxvals{1};
        Nind = 1:length(Nxvals{1});
    end
    if params.correctITD%Correct the ITD convention for SGSR
        Xvals = -Xvals;
    end
    Xvals = Xvals/1000; %us --> ms
    %Get the spike data for each value of X
    N = length(Xvals);
    Xspline = min(Xvals):dt:max(Xvals);
    if ~isempty(Nxvals)
        for i = 1:N
            ST = dsn{negsplind}.Data.SpikeTimes{Nind(i)}; %new method
            NRate.reps{i} = get_bs_rate(ST); %new method
        end
        [NRate.mean,NRate.std,NRate.spline.mean,NRate.spline.std] = bs_rate(vertcat(NRate.reps{:}));
    end
    if ~isempty(Pxvals)
        for i = 1:N
            ST = dsp{possplind}.Data.SpikeTimes{Pind(i)}; %new method
            PRate.reps{i} = get_bs_rate(ST); %new method
        end
        [PRate.mean,PRate.std,PRate.spline.mean,PRate.spline.std] = bs_rate(vertcat(PRate.reps{:}));
    end
    if ~isempty(Pxvals) && ~isempty(Nxvals)
        [DRate.mean,DRate.std,DRate.spline.mean,DRate.spline.std] = bs_rate(PRate.reps(1,:),NRate.reps(1,:));
    end
    
    %Get the BD
    if ~isempty(Pxvals) && ~isempty(Nxvals)
        %You have the difcor - therefore BD is just the peak
        BDtype(j) = 1;
        [~,ind] = max(DRate.spline.mean);
    elseif ~isempty(Pxvals)
        %You only have the positive correlation NDF - take the peak as BD
        BDtype(j) = 2;
        [~,ind] = max(PRate.spline.mean);
    elseif ~isempty(Nxvals)
        %You only have the negative correlation NDF - take the minimum as "BD"
        BDtype(j) = 3;
        [~,ind] = min(NRate.spline.mean);%Probably this will never be used... but just here for completeness
    end
    BD(j) = Xspline(ind);
    
    %Find other peaks and valleys
    if ~isempty(Pxvals) && ~isempty(Nxvals)
        vals = DRate.spline.mean-min(DRate.spline.mean);
        [peakval,peakind] = findpeaks(vals);
        [troughval,troughind] = findpeaks(-vals);
        troughval = -troughval;
        if min(peakind)<min(troughind)
            troughind = [1;troughind];
            troughval = [vals(1);troughval];
        end
        if max(peakind)>max(troughind)
            troughind = [troughind;length(DRate.spline.mean)];
            troughval = [troughval;vals(length(DRate.spline.mean))];
        end
        Dpeaks{j} = Xspline(peakind);
        Dtroughs{j} = Xspline(troughind);
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
        Dpeaks{j} = Dpeaks{j}(logical(keeppeak));
        Dtroughs{j} = Dtroughs{j}(logical(keeptrough));
        Dpeaks{j} = Dpeaks{j}(Dpeaks{j}~=BD(j));
    elseif ~isempty(Pxvals)
        vals = PRate.spline.mean;
        [peakval,peakind] = findpeaks(vals);
        [troughval,troughind] = findpeaks(-vals);
        troughval = -troughval;
        if min(peakind)<min(troughind)
            troughind = [1;troughind];
            troughval = [vals(1);troughval];
        end
        if max(peakind)>max(troughind)
            troughind = [troughind;length(PRate.spline.mean)];
            troughval = [troughval;vals(length(PRate.spline.mean))];
        end
        Ppeaks{j} = Xspline(peakind);
        Ptroughs{j} = Xspline(troughind);
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
        Ppeaks{j} = Ppeaks{j}(logical(keeppeak));
        Ptroughs{j} = Ptroughs{j}(logical(keeptrough));
        Ppeaks{j} = Ppeaks{j}(Ppeaks{j}~=BD(j));
    elseif ~isempty(Nxvals)
        vals = NRate.spline.mean;
        [peakval,peakind] = findpeaks(vals);
        [troughval,troughind] = findpeaks(-vals);
        troughval = -troughval;
        if min(peakind)<min(troughind)
            troughind = [1;troughind];
            troughval = [vals(1);troughval];
        end
        if max(peakind)>max(troughind)
            troughind = [troughind;length(NRate.spline.mean)];
            troughval = [troughval;vals(length(NRate.spline.mean))];
        end
        Npeaks{j} = Xspline(peakind);
        Ntroughs{j} = Xspline(troughind);
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
        Npeaks{j} = Npeaks{j}(logical(keeppeak));
        Ntroughs{j} = Ntroughs{j}(logical(keeptrough));
        Ntroughs{j} = Ntroughs(Ntroughs{j}~=BD(j));
    end
    
    %Get the frequency domain signals
    nf = 0.5/(dt/1000);
    freq = linspace(0,1,NFFT/2)*nf;
    if ~isempty(Pxvals) && ~isempty(Nxvals)
        vals = DRate.spline.mean;
        hanwin = hann(length(vals));
        spec = fft((vals-mean(vals)).*hanwin,NFFT);
        Dphase{j} = angle(spec(1:NFFT/2));
        Dspec{j} = 20*log10(abs(spec(1:NFFT/2)));
        Dspec{j} = Dspec{j}-max(Dspec{j});
        Dphase{j}(Dspec{j}<-20) = nan;
        [~,ind] = max(Dspec{j});
        DspecPeakHz(j) = freq(ind);
    end
    if ~isempty(Nxvals)
        vals = NRate.spline.mean;
        hanwin = hann(length(vals));
        spec = fft((vals-mean(vals)).*hanwin,NFFT);
        Nphase{j} = angle(spec(1:NFFT/2));
        Nspec{j} = 20*log10(abs(spec(1:NFFT/2)));
        Nspec{j} = Nspec{j}-max(Nspec{j});
        Nphase{j}(Nspec{j}<-20) = nan;
        [~,ind] = max(Nspec{j});
        NspecPeakHz(j) = freq(ind);
    end
    if ~isempty(Pxvals)
        vals = PRate.spline.mean;
        hanwin = hann(length(vals));
        spec = fft((vals-mean(vals)).*hanwin,NFFT);
        Pphase{j} = angle(spec(1:NFFT/2));
        Pspec{j} = 20*log10(abs(spec(1:NFFT/2)));
        Pspec{j} = Pspec{j}-max(Pspec{j});
        Pphase{j}(Pspec{j}<-20) = nan;
        [~,ind] = max(Pspec{j});
        PspecPeakHz(j) = freq(ind);
    end
    
    %Build the output structure
    
    if ~isempty(Pxvals)
        NDF{j}.positive.rate = PRate;
        NDF{j}.positive.power = Pspec{j};
        NDF{j}.positive.phase = Pphase{j};
        NDF{j}.positive.peakhz = PspecPeakHz(j);
        NDF{j}.positive.freq = freq;
    end
    if ~isempty(Nxvals)
        NDF{j}.negative.rate = NRate;
        NDF{j}.negative.power = Nspec{j};
        NDF{j}.negative.phase = Nphase{j};
        NDF{j}.negative.peakhz = NspecPeakHz(j);
        NDF{j}.negative.freq = freq;
    end
    if ~isempty(Pxvals) && ~isempty(Nxvals)
        NDF{j}.difcor.rate = DRate;
        NDF{j}.difcor.power = Dspec{j};
        NDF{j}.difcor.phase = Dphase{j};
        NDF{j}.difcor.peakhz = DspecPeakHz(j);
        NDF{j}.difcor.freq = freq;
    end
    NDF{j}.xvals.ITDx = Xvals;
    NDF{j}.xvals.ITDxspline = Xspline;
    NDF{j}.bd.bdtype = BDtype(j);
    NDF{j}.bd.bd = BD(j);
    if ~isempty(Pxvals) && ~isempty(Nxvals)
        NDF{j}.bd.dpeaks = Dpeaks{j};
        NDF{j}.bd.dtroughs = Dtroughs{j};
    elseif ~isempty(Pxvals)
        NDF{j}.bd.ppeaks = Ppeaks{j};
        NDF{j}.bd.ptroughs = Ptroughs{j};
    elseif ~isempty(Nxvals)
        NDF{j}.bd.npeaks = Npeaks{j};
        NDF{j}.bd.ntroughs = Ntroughs{j};
    end
    NDF{j}.SPL = uSPLs(j);
    clear PRate NRate DRate;
end

%% Local functions
    function [rateout] = get_bs_rate(spikes)
        w = minST;
        count = 1;
        while w<maxST
            numsp = length(spikes(spikes>w & spikes<w+stepST));
            rateout(count) = numsp/(stepST/1000);
            w = w+stepST;
            count = count+1;
        end
    end

    function [varargout] = bs_rate(varargin)
        if nargin==2
            rP = varargin{1};
            rN = varargin{2};
            for ii = 1:length(rP)
                for jj = 1:length(rP{ii})
                    for kk = 1:length(rN{ii})
                        rD(ii,jj,kk) = rP{ii}(jj)-rN{ii}(kk);
                    end
                end
            end
            nrest = length(rP{1})*length(rN{1});
            rD = reshape(rD,length(rP),nrest);
            for ii = 1:nrest
                rDspline(:,ii) = spline(Xvals,rD(:,ii),Xspline);
            end
            varargout{1} = mean(rD,2);
            varargout{2} = std(rD,0,2);
            varargout{3} = mean(rDspline,2);
            varargout{4} = std(rDspline,0,2);
        elseif nargin==1
            varargout{1} = mean(varargin{1},2);
            varargout{2} = std(varargin{1},0,2);
            for ii = 1:size(varargin{1},2)
                rspline(:,ii) = spline(Xvals,varargin{1}(:,ii),Xspline);
            end
            varargout{3} = mean(rspline,2);
            varargout{4} = std(rspline,0,2);
        end
    end
end