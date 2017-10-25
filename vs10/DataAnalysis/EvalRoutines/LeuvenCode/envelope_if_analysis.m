function [kernels] = envelope_if_analysis(kernels,signalPOWER,params)
maxorder = params.MaxOrder;
Nchan = params.Nchan;
Nwin = params.Nwin;
dt = params.dt;
env_crit = 2;
gate_time = 0.5;%ms
gate_samps = round(gate_time/dt);
x = linspace(0.25,0,gate_samps)';
gate = cos(x*2*pi).^2;
gate(1) = 0;
Nsamps = length(kernels.h1wf);
kernels.h1c = zeros(size(kernels.h1wf));
if maxorder==2
    kernels.h2c = zeros(size(kernels.h2));
    if Nchan==2
        kernels.h2Xc = zeros(size(kernels.h2X));
    end
end
h1envif.zcriterion = env_crit;
startind = zeros(Nchan,1);
stopind = zeros(Nchan,1);
window = zeros(Nwin,Nchan);
for i = 1:Nchan
    z = hilbert(kernels.h1wf(:,i));
    h1envif.IF(:,i) = (1000/dt)/(2*pi)*diff(unwrap(angle(z)));%instantaneous frequency (Hz)
    h1envif.ENV(:,i) = abs(z);%Hilbert envelope of the revcors
    znull = hilbert(permute(kernels.nullh1wf(:,i,:),[1 3 2]));
    nullenv = abs(znull);
    %convert to z scores
    h1envif.z(:,i) = (h1envif.ENV(:,i)-mean(nullenv,2))./std(nullenv,1,2);
    %Apply smoothing (Triangular smoothing function)
    h1envif.z(:,i) = Trifilter(h1envif.z(:,i)',15)';
    [h1envif.ENVmax_sps(i),peakind] = max(h1envif.ENV(:,i));
    %express the peak as a modulation factor of the mean firing rate (Henry
    %& Heinz, 2016) - facilitates comparison across units with different
    %mean firing rates
    h1envif.ENVmax_mod(i) = h1envif.ENVmax_sps(i)*(sqrt(signalPOWER)/kernels.h0);
    if ~isempty(find(h1envif.z(1:peakind,i)<env_crit,1,'last'))
        startind(i) = min(length(h1envif.z(:,i))-(2*gate_samps)-1, max(gate_samps+1,find(h1envif.z(1:peakind,i)<env_crit,1,'last')));
    else
        startind(i) = gate_samps+1;
    end
    if ~isempty(find(h1envif.z(peakind+1:end,i)<env_crit,1,'first'))
        stopind(i) = max(min(length(h1envif.z(:,i))-gate_samps,find(h1envif.z(peakind+1:end,i)<env_crit,1,'first')+peakind),startind(i)+1);
    else
        stopind(i) = max(length(h1envif.z(:,i))-gate_samps,startind(i)+1);
    end
    h1envif.frontlatency(i) = startind(i)*dt;
    h1envif.endlatency(i) = stopind(i)*dt;
    h1envif.z([1:startind(i)-1-gate_samps stopind(i)+1+gate_samps:end],i) = NaN;
    h1envif.IF([1:startind(i)-1-gate_samps stopind(i)+1+gate_samps:end],i) = NaN;
    window(:,i) = [zeros(startind(i)-1-gate_samps,1); gate; ones(stopind(i)-startind(i)+1,1); flipud(gate); zeros(Nsamps-stopind(i)-gate_samps,1)];
    kernels.h1c(:,i) = kernels.h1wf(:,i).*window(:,i);
    kernels.h1c_mod(:,i) = kernels.h1c(:,i)*(sqrt(signalPOWER)/kernels.h0);
    if maxorder ==2
        window2D = window(:,i)*window(:,i)';
        kernels.h2c(:,:,i) = kernels.h2(:,:,i).*window2D;
    end
end
maxsize = max(stopind)+gate_samps+1;
if maxorder==2 && Nchan==2
    window2D = window(:,2)*window(:,1)';
    kernels.h2Xc = kernels.h2X.*window2D;
    kernels.h2Xc = kernels.h2Xc(1:maxsize,1:maxsize);
end
if maxorder==2
    kernels.h2c = kernels.h2c(1:maxsize,1:maxsize,:);
end
kernels.h1envif = h1envif;
return;