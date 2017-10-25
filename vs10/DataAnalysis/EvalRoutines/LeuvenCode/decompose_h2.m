function [H2,H2X] = decompose_h2(kernels,params)
Nchan = params.Nchan;
H2 = zeros(size(kernels.h2));
if Nchan==2
    H2X = zeros(size(kernels.h2X));
end
Npoints = size(kernels.h2,1);
Nboot = size(kernels.nullh2,4);
for i = 1:Nchan%Loop over the two channels
    progress_bar = waitbar(0,'Denoising h2...');
    knull = zeros(Npoints,Nboot);
    for j = 1:Nboot%Loop over number of bootstrap reps
        [Unull,Snull,Vnull] = svd(kernels.nullh2(:,:,i,j));
        Snull = diag(Snull);
        knull(:,j) = sign(diag(Unull).*diag(Vnull)).*Snull;
        waitbar(j/Nboot);
    end
    knull = sort(knull,'descend');
    knull = knull([1 Npoints],:);
    poscrit = mean(knull(1,:))+3*std(knull(1,:));
    negcrit = mean(knull(2,:))+3*std(knull(2,:));
    clear knull Snull Unull Vnull;
    
    [U,S,V] = svd(kernels.h2(:,:,i));
    S = diag(S);
    k = sign(diag(U).*diag(V)).*S;
    U = repmat(abs(k'),Npoints,1).*U;%The weighted vectors
    posind = find(k>poscrit);
    negind = find(k<negcrit);
    clear k S;
    if or(~isempty(posind),~isempty(negind))
        H2(:,:,i) = U(:,[posind;negind])*V(:,[posind;negind])';
    end
    clear U V posind negind;
    close(progress_bar);
end
%For the cross kernel
if Nchan==2
    progress_bar = waitbar(0,'Denoising cross-h2...');
    for j = 1:Nboot%Loop over number of bootstrap reps
        [Unull,Snull,Vnull] = svd(kernels.nullh2X(:,:,j));
        Snull = diag(Snull);
        knull(:,j) = sign(diag(Unull).*diag(Vnull)).*Snull;
        waitbar(j/Nboot);
    end
    knull=sort(knull,'descend');
    knull = knull([1 Npoints],:);
    poscrit = mean(knull(1,:))+3*std(knull(1,:));
    negcrit = mean(knull(2,:))+3*std(knull(2,:));
    clear knull Snull Unull Vnull;
    
    [U,S,V] = svd(kernels.h2X);
    S = diag(S);
    k = sign(diag(U).*diag(V)).*S;
    U = repmat(abs(k'),Npoints,1).*U;%The weighted vectors
    posind = find(k>poscrit);
    negind = find(k<negcrit);
    clear k S;
    if or(~isempty(posind),~isempty(negind))
        H2X = U(:,[posind;negind])*V(:,[posind;negind])';
    end
    clear U V posind negind;
    close(progress_bar);
end
return;
