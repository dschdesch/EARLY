load 'cirevcordata';

X = meanh1domfreq;
Y = domfreqdiff_oct;

ind = find(abs(Y)<=1 & X<=1.4);

X = X(ind);
Y = Y(ind);
[X,ind] = sort(X);
Y = Y(ind);

[Ylowessfit,YlowessSTD] = mylowessbootstrap(X',Y',0.3,'lowess',200);


figure;
plot(X,Y,'k+','markersize',10,'linewidth',1);

