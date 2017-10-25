function [] = Binaural_Wiener(MaxOrder,varargin)
%-------------------------------------------------------------------------
%Binaural revcor analyses
%-------------------------------------------------------------------------
%% 0. Parse input arguments
inputstruct = varargin;
[data,params] = Wiener_Parse_Inputs(MaxOrder,2,varargin);
% savethedata('inputs');
%% 1. Noise-Delay Function analyses
if isfield(data,'dsp') && isfield(data,'dsn')
    [NDF] = NDF_analysis(data.dsp,data.dsn,params);
elseif isfield(data,'dsp')
    [NDF] = NDF_analysis(data.dsp,[],params);
elseif isfield(data,'dsn')
    [NDF] = NDF_analysis([],data.dsn,params);
else
    NDF = [];
end
NDFspls = nan(size(NDF));
for i = 1:length(NDF)
    NDFspls(i) = NDF{i}.SPL;
    savethedata('NDF',NDF{i},NDFspls(i));
end
%% 2. Binaural-Beat analyses
if isfield(data,'dsbfs')
    if ~isempty(NDF)
        [BFS] = bfs_analysis(data.dsbfs,params,NDF{1}.xvals.ITDx);
    else
        [BFS] = bfs_analysis(data.dsbfs,params,-5:0.2:5);
    end
    BFSspls = nan(size(BFS));
    for i = 1:length(BFS)
        BFSspls(i) = BFS{i}.SPL;
        savethedata('BFS',BFS{i},BFSspls(i));
    end
else
    BFS = [];
end
%% 3. NRHO Function analyses
if isfield(data,'dsnrho')
    [NRHO] = nrho_analysis(params.contrachan,data.dsnrho);
else
    NRHO = [];
end
NRHOspls = nan(size(NRHO));
for i = 1:length(NRHO)
    NRHOspls(i) = NRHO(i).SPL;
    savethedata('NRHO',NRHO(i),NRHOspls(i));
end
%% 4. THR analysis
if isfield(data,'dsthr')
    TC.data = cell(1,length(data.dsthr));
    for i=1:length(data.dsthr)
        TC.data{i} = EvalTHR(data.dsthr{i}.dataobj);
        TC.data{i}.channel = data.dsthr{i}.Stimulus.Special.ActiveChan;%1 = L, 2 = R, 0 = D;
    end
    savethedata('TC');
else
    TC = [];
end
%% 5. Monaural Weiner analyses (h0,h1,h2)
if isfield(data,'dsipsi')
    nripsidata = length(data.dsipsi);
    ipsidsid = cell(1,nripsidata);
    for i = 1:nripsidata
        ipsidsid{i} = data.dsipsi{i}.ID.SeqID;
        switch data.dsipsi{i}.ID.StimType
            case 'NTD'
                ipsiind = regexp(ipsidsid{i},'-NTD');
            case 'NSPL'
                ipsiind = regexp(ipsidsid{i},'-NSPL');
        end
        ipsistring{i} = ipsidsid{i}(1:ipsiind-1);
    end
    [IPSI] = Monaural_Wiener(MaxOrder,0,params.AnNum,params.UnNum,'ndz',ipsistring);
    savethedata('IPSI',IPSI.SPL);
else
    IPSI = [];
end
if isfield(data,'dscontra')
    nrcontradata = length(data.dscontra);
    contradsid = cell(1,nrcontradata);
    for i = 1:nrcontradata
        contradsid{i} = data.dscontra{i}.ID.SeqID;
        switch data.dscontra{i}.ID.StimType
            case 'NTD'
                contraind = regexp(contradsid{i},'-NTD');
            case 'NSPL'
                contraind = regexp(contradsid{i},'-NSPL');
        end
        contrastring{i} = contradsid{i}(1:contraind-1);
    end
    [CONTRA] = Monaural_Wiener(MaxOrder,0,params.AnNum,params.UnNum,'ndz',contrastring);
    savethedata('CONTRA',CONTRA.SPL);
else
    CONTRA = [];
end
%% 5. Binaural Wiener analyses (h0,h1,h2) for both ears and the cross-second-order kernel(h2X)
if isfield(data,'dsz')%only if those data exist
    nrzdata = length(data.dsz);
    dszspl = cell(1,nrzdata);
    for i = 1:nrzdata
        switch data.dsz{i}.ID.StimType
            case 'NTD'
                %Its an NTD dataset
                dszspl{i} = data.dsz{i}.Stimulus.StimParam.SPL(1);
            case 'NSPL'
                %Its not an NTD dataset - try NSPL
                dszspl{i} = data.dsz{i}.Stimulus.IndepVar.Values;
        end
    end
    allspls = cat(1,dszspl{:});
    allspls = allspls(~isnan(allspls));
    uniquespls = unique(allspls);
    nrresponses = length(allspls);
    nrspls = length(uniquespls);
    whichdata = [];
    for i = 1:nrzdata
        whichdata = [whichdata; ones(length(dszspl{i}(~isnan(dszspl{i}))),1)*i];
    end
    fprintf('Found %2.0f binaural revcor datasets at %2.0f SPLs\n',nrresponses,nrspls);
    for kk = 1:nrspls
        thisspl = uniquespls(kk);
        thesedata = find(allspls==thisspl);
        nrthesedata = length(thesedata);
        NrReps = zeros(1,nrthesedata);
        thiswv = cell(1,nrthesedata);
        thesespikes = cell(1,nrthesedata);
        for i = 1:nrthesedata
            %Get the stimulus waveform and sampling period
            thiswv{i} = StimSam(data.dsz{whichdata(thesedata(i))}.dataobj,1);
            %Truncate to the relevant part
            thiswv{i} = thiswv{i}(params.minSample:params.maxSample,:);
            %Re-scale stimulus waveform in units of Pascals
            [thiswv{i},signalPOWER] = rescalewv(thiswv{i},thisspl,'Pa');
            %get spike times
            thesespikes{i} = data.dsz{whichdata(thesedata(i))}.Data.SpikeTimes;
            thesespikes{i} = thesespikes{i}((dszspl{whichdata(thesedata(i))}==thisspl),:);
            %How many reps are added from this dataset?
            NrReps(i) = size(thesespikes{i},2);
            %concatenate all spikes from all reps from this dataset
            thesespikes{i} = cat(2,thesespikes{i}{:});
            %limit spike times to only include driven spikes
            thesespikes{i} = thesespikes{i}(thesespikes{i}>=params.minST & thesespikes{i}<=params.maxST);
            %Correct for gate time
            thesespikes{i} = thesespikes{i}-params.minTime;
            %Correct for concatenation of noise tokens
            thesespikes{i} = thesespikes{i}+((params.maxTime-params.minTime)*(i-1));
        end
        %Total input noise duration - all reps, all noise tokens
        totaldur = sum(NrReps)*(params.maxST-params.minST);
        %All spike times, relative to the concatenated input noise
        spikes = cat(2,thesespikes{:})';clear thesespikes;
        %Concatenated input noise
        wv = cat(1,thiswv{:});clear thiswv;
        %Make sure the order is [ipsi contra]
        [~,~,wv] = flipchannels(params.AnNum,params.contrachan,wv);
        %-----h1 estimation------------
        %Inform the user something is happening
        fprintf('Calculating h1...\n');
        %Get the zeroth and first-order kernels
        [h0,h1,Time,TotalSpikes] = get_h0_h1_kernel(spikes,wv,params,signalPOWER,totaldur,'SpSpPa');
        %Inform the user something is happening
        fprintf('Completed h1 calculations!\n');
        %-----h2 and h2X estimation----
        %If you want the second order kernels
        if params.MaxOrder==2
            %Inform the user something is happening
            fprintf('Calculating h2...\n');
            %Get the second-order kernel
            [h2,h2X] = get_h2_kernel(spikes,wv,params,signalPOWER,totaldur,'SpSpPa2');
            %Inform the user something is happening
            fprintf('Completed h2 calculations!\n');
        else
            %Inform the user something didn't happen
            fprintf('h2 not requested\n');
        end
        %Bootstrapping
        if  params.MaxOrder==2
            %Bootstrap h1 h2 and h2X
            [nullh1,nullh2,nullh2X] = bootstrap_h1_h2(spikes,wv,params,totaldur,thisspl);
        else
            %Bootstrap h1 only
            [nullh1] = bootstrap_h1_h2(spikes,wv,params,totaldur,thisspl);
        end
        %Get the spectra and the Wiener filtered versions
        [h1mag,nullh1mag,ffax,h1phase,h1zscore,h1phasemask,h1magmask,h1wf,nullh1wf] = Phase_Freq_analysis(h1,nullh1,params);
        %Get the dominant frequency for each revcor and the bandwidth measures
        [h1bw,qvals,domfreq,domind] = Bandwidth_analysis(h1mag,h1zscore,params.spec_crit,ffax);
        %Gather all the information on the kernels
        kernels = struct('h0',h0,'h1',h1,'nullh1',nullh1,'h1wf',h1wf,'nullh1wf',nullh1wf,...
            'h1magmask',h1magmask,'h1mag',h1mag,'h1phasemask',h1phasemask,'h1phase',h1phase,...
            'nullh1mag',nullh1mag,'h1ffax',ffax,'h1zscore',h1zscore,'Time',Time,...
            'TotalSpikes',TotalSpikes,'h1bw',h1bw,'h1qvals',qvals,'h1domfreq',domfreq,...
            'h1domind',domind,'SPL',thisspl);
        %clear some memory
        clear h0 h1 nullh1 h1wf nullh1wf h1magmask h1mag h1phasemask h1phase nullh1mag...
            ffax h1zscore Time TotalSpikes h1bw qvals domfreq domind coeffs tuning_diff SPL;
        if MaxOrder==2
            kernels.h2 = h2; kernels.nullh2 = nullh2; kernels.h2X = h2X; kernels.nullh2X = nullh2X;
            clear h2 nullh2 h2X nullh2X;
        end
        %Low pass filter the kernels
        kernels = low_pass_kernels(kernels,params,[]);
        %Use singular value decomposition to clean up h2 if it exists
        if MaxOrder==2
            [kernels.H2,kernels.H2X] = decompose_h2(kernels,params);%Do the singular value decomposition and return denoised kernels
        end
        %Get the envelope and instantaneous frequency info, and clean up the h1 and h2 kernels
        [kernels] = envelope_if_analysis(kernels,signalPOWER,params);
        %Do the circular-linear regression on the monaural phase-frequency data and on the binaural phase frequency data
        [kernels] = Phase_Freq_Fit(kernels);
        %Get the tuning difference between ipsi and contra, using the cross correlation method
        if ismember(thisspl,NDFspls)
            [~,ndfind] = ismember(thisspl,NDFspls);
            [kernels] = get_tuning_difference(kernels,NDF{ndfind});
        else
            kernels.tuning_diff = [];
        end
        %------Predict the binaural responses (if you have them at this SPL)---------------
        if ismember(thisspl,NDFspls)
            %Make the Bayesian model
            [bayesian_model] = Make_Bayesian_Model(spikes,wv,kernels,params.dt/1000);
            %Make predicted noise-delay curve
            [~,ndfind] = ismember(thisspl,NDFspls);
            [NDF{ndfind},bayesian_model] = NDF_prediction(params.AnNum,params.contrachan,NDF{ndfind},params.dt/1000,kernels,bayesian_model);
            savethedata('NDF',NDF{ndfind},thisspl);
            savethedata('BM',thisspl);
            [updatedBD,NDF{ndfind}] = correctNDF_BD(kernels,NDF{ndfind});
            if updatedBD
                savethedata('NDF',NDF{ndfind},thisspl);
            end
            if ismember(thisspl,NRHOspls)
                [~,nrhoind] = ismember(thisspl,NRHOspls);
                %Make predicted NRHO function
                if isfield(data,'dsnrho')
                    %If we have the real data, use the same parameters
                    rho = data.dsnrho{1}.Stimulus.IndepVar.Values;
                    rho = rho(~isnan(rho));%make sure you only take values of rho for which we actually have data
                    nrhoitd = data.dsnrho{1}.Stimulus.StimParam.delay;%delay was probably always applied to the left ear, but check at next stage
                else
                    %Otherwise use some defaults
                    nrhoitd = NDF{ndfind}.bd.bd;
                    rho = (-1:0.1:1)';
                end
                if ~isempty(NRHO)
                    temp = NRHO_prediction(params.AnNum,params.contrachan,rho,params.dt/1000,kernels,bayesian_model,nrhoitd,NRHO(nrhoind));
                else
                    temp = NRHO_prediction(params.AnNum,params.contrachan,rho,params.dt/1000,kernels,bayesian_model,nrhoitd,NRHO);
                end
                NRHO(nrhoind).predicted = temp.predicted;
                NRHO(nrhoind).COD = temp.COD;
                clear temp;
                savethedata('NRHO',NRHO(nrhoind),thisspl);
            end
            interauralphase = -0.375:0.125:0.5;
            [IAP] = IAP_prediction(params.AnNum,params.contrachan,params.dt/1000,kernels,bayesian_model,NDF{ndfind}.xvals.ITDx,interauralphase);
            savethedata('IAP',thisspl);
        else
            bayesian_model = [];
            IAP = [];
        end
        if params.MaxOrder==2
            %Remove some of the bootstrapped data prior to saving (to avoid huge data files)
            kernels = rmfield(kernels,{'nullh2','nullh2X'});
            %Save the second-order kernels as single-precision floating-point
            %numbers for storage efficiency
            kernels.h2 = single(kernels.h2);kernels.h2X = single(kernels.h2X);...
                kernels.H2 = single(kernels.H2);kernels.H2X = single(kernels.H2X);...
                kernels.h2c = single(kernels.h2c);kernels.h2Xc = single(kernels.h2Xc);
        end
        savethedata('kernels',thisspl);
        clear kernels;
    end
else
    bayesian_model = [];
    kernels = [];
    IAP = [];
end
% plotstruct = struct('kernels',kernels,'BFS',BFS,'NDF',NDF,'NRHO',NRHO,'BM',bayesian_model,'IAP',IAP,'CONTRA',CONTRA,'IPSI',IPSI,'TC',TC);
% plot_Binaural_Wiener(plotstruct);
% clear plotstruct;
%% Save the data structure
    function [] = savethedata(datatype,varargin)
        % Save the data structure
        UnNum = params.UnNum;
        AnNum = params.AnNum;
        savedirname = fullfile(params.dirname,[AnNum '_' UnNum]);
        if ~exist(savedirname,'dir')
            mkdir(savedirname);
        end
        switch datatype
            case 'inputs'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_inputs.mat']);
                save (savedfilename,'inputstruct');
                disp (['Saved: ' savedfilename]);
            case 'IAP'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('IAP_%2.2fdB.mat',varargin{1})]);
                save (savedfilename,'IAP');
                disp (['Saved: ' savedfilename]);
            case 'BM'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('BM_%2.2fdB.mat',varargin{1})]);
                save (savedfilename,'bayesian_model');
                disp (['Saved: ' savedfilename]);
            case 'kernels'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('BinRev_%2.2fdB.mat',varargin{1})]);
                save (savedfilename,'kernels');
                disp (['Saved: ' savedfilename]);
            case 'NDF'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('NDF_%2.2fdB.mat',varargin{2})]);
                ndf = varargin{1};
                save (savedfilename,'ndf');
                disp (['Saved: ' savedfilename]);
                clear ndf;
            case 'NRHO'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('NRHO_%2.2fdB.mat',varargin{2})]);
                nrho = varargin{1};
                save (savedfilename,'nrho');
                disp (['Saved: ' savedfilename]);
                clear nrho;
            case 'BFS'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('BFS_%2.2fdB.mat',varargin{2})]);
                bfs = varargin{1};
                save (savedfilename,'bfs');
                disp (['Saved: ' savedfilename]);
                clear bfs;
            case 'TC'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_TC.mat']);
                save (savedfilename,'TC');
                disp (['Saved: ' savedfilename]);
            case 'IPSI'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('MonRevIPSI_%2.2fdB.mat',varargin{1})]);
                save (savedfilename,'IPSI');
                disp (['Saved: ' savedfilename]);
            case 'CONTRA'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('MonRevCONTRA_%2.2fdB.mat',varargin{1})]);
                save (savedfilename,'CONTRA');
                disp (['Saved: ' savedfilename]);
        end
    end
end