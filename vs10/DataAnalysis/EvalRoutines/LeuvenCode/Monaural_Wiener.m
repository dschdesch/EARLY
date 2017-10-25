function [OUT] = Monaural_Wiener(MaxOrder,plotYN,varargin)
%Calculates 1st and 2nd order Wiener kernels from monaural spike data

%M. Sayles 2015. Leuven.
%-------------------------------------------------------------------------
[data,params] = Wiener_Parse_Inputs(MaxOrder,1,varargin);
inputstruct = varargin;
if params.saveYN
    savethedata('inputs');
end
%% Tuning curve analysis
if isfield(data,'dsthr')
    TC = EvalTHR(data.dsthr{1}.dataobj);
    if params.saveYN
    savethedata('TC');
    end
else
    TC = [];
end
%% Rate-Level function analysis
if isfield(data,'dsspl')
    RLV = get_rate_level(data.dsspl);
    if params.saveYN
    savethedata('RLV');
    end
else
    RLV = [];
end

%% 1. Get revcors (h1: "First-order" Wiener kernel) for both ears
if isfield(data,'dsz')
    nrzdata = length(data.dsz);
    dszspl = cell(1,nrzdata);
    for i = 1:nrzdata
        switch data.dsz{i}.ID.StimType
            case {'NTD','NRHO'}
                %Its an NTD or an NRHO dataset
                dszspl{i} = data.dsz{i}.Stimulus.StimParam.SPL(1);
            case 'NSPL'
                %Its an NSPL dataset
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
    fprintf('Found %2.0f monaural revcor datasets at %2.0f SPLs\n',nrresponses,nrspls);
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
        %Inform the user something is happening
        fprintf('Calculating h1...\n');
        %Get the zeroth and first-order kernels
        [h0,h1,Time,TotalSpikes] = get_h0_h1_kernel(spikes,wv,params,signalPOWER,totaldur,'SpSpPa');
        %Inform the user something is happening
        fprintf('Completed h1 calculations!\n');
        if TotalSpikes<1500
            fprintf('Not enough spikes!\n');
            continue;
        end
        if params.MaxOrder==2
            %Inform the user something is happening
            fprintf('Calculating h2...\n');
            %Get the second-order kernel
            [h2] = get_h2_kernel(spikes,wv,params,signalPOWER,totaldur,'SpSpPa2');
            %Inform the user something is happening
            fprintf('Completed h2 calculations!\n');
        else
            %Inform the user something didn't happen
            fprintf('h2 not requested\n');
        end
        %Bootstrapping
        if  params.MaxOrder==2
            %Bootstrap h1 and h2
            [nullh1,nullh2] = bootstrap_h1_h2(spikes,wv,params,totaldur,thisspl);
        else
            %Bootstrap h1 only
            [nullh1] = bootstrap_h1_h2(spikes,wv,params,totaldur,thisspl);
        end
        %Get the spectra and the Wiener filtered versions
        [h1mag,nullh1mag,ffax,h1phase,h1zscore,h1phasemask,h1magmask,h1wf,nullh1wf] = Phase_Freq_analysis(h1,nullh1,params);
        %Get the dominant frequency for each revcor and the bandwidth measures
        [h1bw,qvals,domfreq,domind] = Bandwidth_analysis(h1mag,h1zscore,params.spec_crit,ffax);
        %Do the circular-linear regression on the monaural phase-frequency data
        kernels = struct('h0',h0,'h1',h1,'nullh1',nullh1,'h1wf',h1wf,'nullh1wf',nullh1wf,...
            'h1magmask',h1magmask,'h1mag',h1mag,'h1phasemask',h1phasemask,'h1phase',h1phase,...
            'nullh1mag',nullh1mag,'h1ffax',ffax,'h1zscore',h1zscore,'Time',Time,...
            'TotalSpikes',TotalSpikes,'h1bw',h1bw,'h1qvals',qvals,'h1domfreq',domfreq,...
            'h1domind',domind,'SPL',thisspl);
        %clear some memory
        clear h0 h1 nullh1 h1wf nullh1wf h1magmask h1mag h1phasemask h1phase nullh1mag...
            ffax h1zscore Time TotalSpikes h1bw qvals domfreq domind coeffs;
        if params.MaxOrder==2
            kernels.h2 = h2; kernels.nullh2 = nullh2;
            clear h2 nullh2;
        end
        %Low pass filter the kernels
        kernels = low_pass_kernels(kernels,params,TC);
        %Use singular value decomposition to clean up h2 if it exists
        if MaxOrder==2
            [kernels.H2] = decompose_h2(kernels,params);%Do the singular value decomposition and return denoised kernels
        end
        %Get the envelope and instantaneous frequency info, and clean up the h1 and h2 kernels
        [kernels] = envelope_if_analysis(kernels,signalPOWER,params);
        if params.MaxOrder==2
            %Remove some of the bootstrapped data prior to saving (to avoid huge data files)
            kernels = rmfield(kernels,'nullh2');
            %Save the second-order kernels as single-precision floating-point
            %numbers
            kernels.h2 = single(kernels.h2);kernels.H2 = single(kernels.H2); kernels.h2c = single(kernels.h2c);
        end
        %Do the circular-linear regression on the monaural phase-frequency data
        [kernels] = Phase_Freq_Fit(kernels);
        %Figures if required
        if plotYN
            plot_Monaural_Wiener(kernels);
        end
        db = dbstack;
        if strcmp(db(end).name,'Binaural_Wiener')
            OUT = kernels;
            return;
        else
            if params.saveYN
                savethedata('kernels');
            end
            clear kernels;
        end
    end
end
%% Local functions
    function [] = savethedata(datatype)
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
            case 'kernels'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_' sprintf('MonRev_%2.2fdB.mat',thisspl)]);
                save (savedfilename,'kernels');
                disp (['Saved: ' savedfilename]);
            case 'RLV'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_RLV.mat']);
                save (savedfilename,'RLV');
                disp (['Saved: ' savedfilename]);
            case 'TC'
                savedfilename = fullfile(savedirname, [AnNum '_' UnNum '_TC.mat']);
                save (savedfilename,'TC');
                disp (['Saved: ' savedfilename]);
        end
    end
end