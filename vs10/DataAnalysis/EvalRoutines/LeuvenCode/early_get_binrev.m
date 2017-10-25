function [kernels,BM,NRHO,IAP,NDF] = early_get_binrev(data,params,NDF)
anadir = 'C:\EARLY\vs10\DataAnalysis';
cd ('C:\Early_StimDefLeuven\helpers');
for i = 1:length(data)
    %re-make the stimulus waveform without the ear calibration
    data{i}.stim_param = noiseDelayStim_nocalib(data{i}.stim_param);
end
cd (anadir);

nchan = 2;
fs = data{1}.stim_param.Fsam;
dt = 1e3/fs;
burstdur = data{1}.stim_param.BurstDur;
Nsamps = ceil(burstdur/dt);
params.dt = dt;
params.Nchan = nchan;
params.Nwin = 2^nextpow2(ceil(25/dt));

for idx = 1:length(data)
    SPL(idx) = data{idx}.stim_param.SPL;
end

uSPLs = unique(SPL);
nrSPLs = length(uSPLs);


for i = 1:length(NDF)
   ndfspl(i) = NDF{i}.SPL; 
end

for SPLidx = 1:nrSPLs
    
    ndfindx = find(ndfspl==uSPLs(SPLidx));
    thesedata = data(SPL==uSPLs(SPLidx));
    
    for idx = 1:length(thesedata)
        wv{idx}(:,1) = thesedata{idx}.stim_param.Waveform(1,1).Samples;
        wv{idx}(:,2) = thesedata{idx}.stim_param.Waveform(1,2).Samples;
        wv{idx} = wv{idx}(1,:);
        [wv{idx},signalPOWER] = rescalewv(wv{idx},uSPLs(SPLidx),'Pa');
        
        %we're always recording on the right, and early is organized ([Left,
        %Right])... and we want ([ipsi,contra])... therefore...
        wv{idx} = wv{idx}(:,[2 1]);
        
        Nreps = size(thesedata{idx}.spikes,2);
        Nreps_presented = Nreps;
        for i = 1:Nreps
            if isempty(thesedata{idx}.spikes{i})
                Nreps_presented = i-1;
                break;
            end
        end
        totaldur(idx) = burstdur*Nreps_presented;
        ST{idx} = cat(2,thesedata{idx}.spikes{:});
        ST{idx} = ST{idx}(ST{idx}<=burstdur);
        ST{idx} = ST{idx}(ST{idx}>params.minST);
        if size(ST{idx},1)==1
            ST{idx} = ST{idx}';
        end
        ST{idx} = ST{idx}+((idx-1)*burstdur);
    end
    totaldur = sum(totaldur);
    ST = cat(1,ST{:});
    wv = cat(1,wv{:});
    [h0,h1,Time,TotalSpikes] = get_h0_h1_kernel(ST,wv,params,signalPOWER,totaldur,'SpSpPa');
    
    if params.MaxOrder==2
        %Inform the user something is happening
        fprintf('Calculating h2...\n');
        %Get the second-order kernel
        [h2,h2X] = get_h2_kernel(ST,wv,params,signalPOWER,totaldur,'SpSpPa2');
        %Inform the user something is happening
        fprintf('Completed h2 calculations!\n');
    else
        %Inform the user something didn't happen
        fprintf('h2 not requested\n');
    end
    %Bootstrapping
    if  params.MaxOrder==2
        %Bootstrap h1 h2 and h2X
        [nullh1,nullh2,nullh2X] = bootstrap_h1_h2(ST,wv,params,totaldur,uSPLs(SPLidx));
    else
        %Bootstrap h1 only
        [nullh1] = bootstrap_h1_h2(ST,wv,params,totaldur,uSPLs(SPLidx));
    end
    %Get the spectra and the Wiener filtered versions
    [h1mag,nullh1mag,ffax,h1phase,h1zscore,h1phasemask,h1magmask,h1wf,nullh1wf] = Phase_Freq_analysis(h1,nullh1,params);
    %Get the dominant frequency for each revcor and the bandwidth measures
    [h1bw,qvals,domfreq,domind] = Bandwidth_analysis(h1mag,h1zscore,params.spec_crit,ffax);
    %Gather all the information on the kernels
    kernels{SPLidx} = struct('h0',h0,'h1',h1,'nullh1',nullh1,'h1wf',h1wf,'nullh1wf',nullh1wf,...
        'h1magmask',h1magmask,'h1mag',h1mag,'h1phasemask',h1phasemask,'h1phase',h1phase,...
        'nullh1mag',nullh1mag,'h1ffax',ffax,'h1zscore',h1zscore,'Time',Time,...
        'TotalSpikes',TotalSpikes,'h1bw',h1bw,'h1qvals',qvals,'h1domfreq',domfreq,...
        'h1domind',domind,'SPL',uSPLs(SPLidx));
    %clear some memory
    clear h0 h1 nullh1 h1wf nullh1wf h1magmask h1mag h1phasemask h1phase nullh1mag...
        ffax h1zscore Time TotalSpikes h1bw qvals domfreq domind coeffs tuning_diff;
    if params.MaxOrder==2
        kernels{SPLidx}.h2 = h2; kernels{SPLidx}.nullh2 = nullh2; kernels{SPLidx}.h2X = h2X; kernels{SPLidx}.nullh2X = nullh2X;
        clear h2 nullh2 h2X nullh2X;
    end
    %Low pass filter the kernels
    kernels{SPLidx} = low_pass_kernels(kernels{SPLidx},params,[]);
    %Use singular value decomposition to clean up h2 if it exists
    if params.MaxOrder==2
        [kernels{SPLidx}.H2,kernels{SPLidx}.H2X] = decompose_h2(kernels{SPLidx},params);%Do the singular value decomposition and return denoised kernels
    end
    %Get the envelope and instantaneous frequency info, and clean up the h1 and h2 kernels
    [kernels{SPLidx}] = envelope_if_analysis(kernels{SPLidx},signalPOWER,params);
    %Do the circular-linear regression on the monaural phase-frequency data and on the binaural phase frequency data
    [kernels{SPLidx}] = Phase_Freq_Fit(kernels{SPLidx});
    %Get the tuning difference between ipsi and contra, using the cross correlation method
    
    if ~isempty(ndfindx)
        [kernels{SPLidx}] = get_tuning_difference(kernels{SPLidx},NDF{ndfindx});
    else
        kernels{SPLidx}.tuning_diff = [];
    end
    
    if ~isempty(ndfindx)
        %------Predict the binaural responses---------------
        %Make the Bayesian model
        [BM{SPLidx}] = Make_Bayesian_Model(ST,wv,kernels{SPLidx},params.dt/1000);
        %Make predicted noise-delay curve
        [NDF{ndfindx},BM{SPLidx}] = NDF_prediction(params.AnNum,params.contrachan,NDF{ndfindx},params.dt/1000,kernels{SPLidx},BM{SPLidx});
        %use some defaults
        nrhoitd = NDF{ndfindx}.bd.bd;
        rho = (-1:0.1:1)';
        NRHO{SPLidx} = [];
        [NRHO{SPLidx}] = NRHO_prediction(params.AnNum,params.contrachan,rho,params.dt/1000,kernels{SPLidx},BM{SPLidx},nrhoitd,NRHO{SPLidx});
        interauralphase = -0.375:0.125:0.5;
        [IAP{SPLidx}] = IAP_prediction(params.AnNum,params.contrachan,params.dt/1000,kernels{SPLidx},BM{SPLidx},NDF{ndfindx}.xvals.ITDx,interauralphase);
    else
        NRHO{SPLidx} = [];
        IAP{SPLidx} = [];
        BM{SPLidx} = [];
    end
    
    
    if params.MaxOrder==2
        %Remove some of the bootstrapped data prior to saving (to avoid huge data files)
        kernels{SPLidx} = rmfield(kernels{SPLidx},{'nullh2','nullh2X'});
        %Save the second-order kernels as single-precision floating-point
        %numbers for storage efficiency
        kernels{SPLidx}.h2 = single(kernels{SPLidx}.h2);kernels{SPLidx}.h2X = single(kernels{SPLidx}.h2X);...
            kernels{SPLidx}.H2 = single(kernels{SPLidx}.H2);kernels{SPLidx}.H2X = single(kernels{SPLidx}.H2X);...
            kernels{SPLidx}.h2c = single(kernels{SPLidx}.h2c);kernels{SPLidx}.h2Xc = single(kernels{SPLidx}.h2Xc);
    end
    clear wv totaldur ST;
end









