function [otvect]=Trifilter(invect,nfw)
nfwi = 2*floor(nfw/2) + 1;
filt = zeros(1,nfwi);
summ = 0;
nfw2 = floor(nfwi/2);
for jj=1:nfw2
    filt(jj) = jj;
    filt(nfwi+1-jj) = jj;
    summ = summ + 2*jj;
end
nfw3 = nfw2 + 1;
filt(nfw3) = nfw3;
summ = summ + nfw3;
filt = filt./summ;
otvect = zeros(size(invect));
for i = 1:size(invect,1)
    svect = size(invect,2) + 2*nfw2;
    vect1 = zeros(1,svect);
    vect1(1:nfw2) = invect(i,1)*ones(1,nfw2);
    vect1(nfw3:svect-nfw2) = invect(i,:);
    vect1(svect-nfw2+1:svect) = invect(i,size(invect,2))*ones(1,nfw2);
    vect2 = conv(vect1, filt);
    otvect(i,:) = vect2(2*nfw2+1:svect);
end
end