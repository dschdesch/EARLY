function density = mykde2d(data1,data2,n,minx,miny,maxx,maxy)
MIN_XY = [minx miny];
MAX_XY = [maxx maxy];
[bandwidth,density,X,Y] = kde2d([data1 data2],n,MIN_XY,MAX_XY);
end