function [nrho] = NRHO_prediction(AnNum,contrachan,rho,dt_s,kernels,BM,ITD,nrho)
SPL = kernels.SPL;
N = length(rho);
spikerate = zeros(N,1);
noisedur = 1000;
if length(ITD)==2
    %in this case, "ITD" is really the delays (left, right) from the NRHO
    %dataset
    ITDear = find(ITD~=0);
    if isempty(ITDear)
        ITD = 0;%easy
    else
        ITD = ITD(ITDear);
        if (contrachan==1 && ITDear==1) || (contrachan==2 && ITDear==2)
            ITD = -1*ITD;
        end
    end
end

for i = 1:N
    wv = MakeGaussNoiseBand(10, 20000, noisedur, 1/dt_s, 2, rho(i), ITD, 0, 999,contrachan);
    wv = rescalewv(wv,SPL,'Pa');%Scale it as pressure
    [dummy,dummy,wv] = flipchannels(AnNum,contrachan,wv);%Make sure the order is [ipsi contra]
    s = zeros(size(wv));%allocate memory
    for j = 1:2
        s(:,j) = conv(wv(:,j),kernels.h1c(:,j),'same');%Do the convolution for each channel
        s(:,j) = s(:,j)/std(s(:,j));%Scale the output relative to the rms for each channel
    end
    s = min(s,4);s = max(s,-4);%Clip anything beyond 4 sigma
    spikeprob = interp2(BM.X,BM.Y,BM.IOprob,s(:,1),s(:,2));%Get the instantaneous firing probability
    spikerate(i) = sum(spikeprob)/(noisedur/1000);%Get the mean firing rate
end

nrho.predicted.rate = spikerate;
nrho.predicted.itd = ITD;
nrho.predicted.rho = rho;
nrho.predicted.powerfit = fit_NRHO(rho,spikerate);

if isfield(nrho,'data')
    if nrho.data.itd==nrho.predicted.itd
        [nrho.COD] = Coeff_Determination(nrho.data.rate,nrho.predicted.rate);
    else
        disp 'COD between measured and predicted NRHO not computed: different ITDs were used!';
    end
else
    nrho.COD = nan;
end

end