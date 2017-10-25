function [h1,h2,h2X,h1null,h2null,h2Xnull] = normalize_kernels(ndata,nrspikes,h1,h2,h2X,h1null,h2null,h2Xnull)
h1 = h1/ScaleFactor;
h2 = h2/ScaleFactor;
h2X = h2X/ScaleFactor;
h1null = h1null/ScaleFactor;
h2null = h2null/ScaleFactor;
h2Xnull = h2Xnull/ScaleFactor;
end