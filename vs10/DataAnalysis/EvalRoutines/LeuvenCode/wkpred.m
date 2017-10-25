% function psth = wkpred(stimulus,kernels,fs) - Ken
function psth = wkpred(stimulus,kernels)

psth.h0 = ones(length(stimulus),1)*kernels.h0;
for i = 1:2
    psth.h1(:,i) = conv(stimulus(:,i),kernels.h1(end:-1:1,i),'same');
end
[psth.h2,psth.h2X] = conv_2Dkernel(stimulus,kernels.h2,kernels.h2X);