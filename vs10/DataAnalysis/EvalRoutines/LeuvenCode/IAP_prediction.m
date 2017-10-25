function [out] = IAP_prediction(AnNum,contrachan,dt_s,kernels,BM,ITD,IAP)
SPL = kernels.SPL;
Nphase = length(IAP);
Nitd = length(ITD);
spikerate = zeros(Nphase,Nitd);
noisedur = 1000;

for i = 1:Nphase
    for k = 1:Nitd
        wv = MakeGaussNoiseBand(10, 20000, noisedur, 1/dt_s, 2, 1, ITD(k), IAP(i), 999,contrachan);
        wv = rescalewv(wv,SPL,'Pa');%Scale it as pressure
        [dummy,dummy,wv] = flipchannels(AnNum,contrachan,wv);%Make sure the order is [ipsi contra]
        s = zeros(size(wv));%allocate memory
        for j = 1:2
            s(:,j) = conv(wv(:,j),kernels.h1c(:,j),'same');%Do the convolution for each channel
            s(:,j) = s(:,j)/std(s(:,j));%Scale the output relative to the rms for each channel
        end
        s = min(s,4);s = max(s,-4);%Clip anything beyond 4 sigma
        spikeprob = interp2(BM.X,BM.Y,BM.IOprob,s(:,1),s(:,2));%Get the instantaneous firing probability
        spikerate(i,k) = sum(spikeprob)/(noisedur/1000);%Get the mean firing rate
    end
end
out.spikerate = spikerate;
out.IAP = IAP;
out.ITD = ITD;
end