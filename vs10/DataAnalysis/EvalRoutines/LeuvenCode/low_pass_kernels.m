function [kernels] = low_pass_kernels(kernels,params,TC)

%Apply some low pass filtering to the kernels
domfreq = kernels.h1domfreq;
maxorder = params.MaxOrder;
dt = params.dt;
Nboot = size(kernels.nullh1,3);
Nchan = size(kernels.h1,2);
if isempty(TC)
    fcut(1) = max(domfreq)*1000*(2^1.5);
else
    fcut(1) = TC.fit.cf*(2^1.5);
end
fcut(2) = fcut(1)+1000;
if maxorder==2
    if Nchan==2
        [kernels.h2X] = my_low_pass(kernels.h2X,dt,2,fcut(1),fcut(2));
        for j = 1:Nboot
            kernels.nullh2X(:,:,j) = my_low_pass(kernels.nullh2X(:,:,j),dt,2,fcut(1),fcut(2));
        end
    end
    for i = 1:Nchan
        if isempty(TC)
            fcut(1) = min(500/dt-2000,domfreq(i)*1000*(2^1.5));
        else
            fcut(1) = min(500/dt-2000,TC.fit.cf*(2^1.5));
        end
        fcut(2) = min(500/dt-1000,fcut(1)+1000);
        kernels.h1wf(:,i) = my_low_pass(kernels.h1wf(:,i),dt,1,fcut(1),fcut(2));
        [kernels.h2(:,:,i)] = my_low_pass(kernels.h2(:,:,i),dt,2,fcut(1),fcut(2));
        for j = 1:Nboot
            kernels.nullh1wf(:,i,j) = my_low_pass(kernels.nullh1wf(:,i,j),dt,1,fcut(1),fcut(2));
            [kernels.nullh2(:,:,i,j)] = my_low_pass(kernels.nullh2(:,:,i,j),dt,2,fcut(1),fcut(2));
        end
    end
else
    for i = 1:Nchan
        if isempty(TC)
            fcut(1) = min(500/dt-2000,domfreq(i)*1000*(2^1.5));
        else
            fcut(1) = min(500/dt-2000,TC.fit.cf*(2^1.5));
        end
        fcut(2) = min(500/dt-1000,fcut(1)+1000);
        kernels.h1wf(:,i) = my_low_pass(kernels.h1wf(:,i),dt,1,fcut(1),fcut(2));
        for j = 1:Nboot
            kernels.nullh1wf(:,i,j) = my_low_pass(kernels.nullh1wf(:,i,j),dt,1,fcut(1),fcut(2));
        end
    end
end

return;