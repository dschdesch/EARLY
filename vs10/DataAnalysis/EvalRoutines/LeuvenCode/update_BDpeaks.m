function [] = update_BDpeaks

BinRevDir = 'C:\work\BinRev\BinRev';

cd (BinRevDir);

datafiles = dir('*.mat');
nrdata = length(datafiles);

for i = 1:nrdata
    cd (BinRevDir);
    S = load(datafiles(i).name);
    if isfield(S.ParamsOut,'NDF')%If we have the noise data
        figure(1);
        yvals = S.ParamsOut.NDF.Difcor_spline;
        yvals = yvals/max(abs(yvals));
        xvals = S.ParamsOut.NDF.NTDx_spline;
        BD = S.ParamsOut.NDF.BD_NTD;
        plot(xvals,yvals,'k-','linewidth',1);hold on;
        ind = find(xvals==BD);
        plot(BD,yvals(ind),'rp','markersize',15);
        [dummy,lowind] = findpeaks(smooth(yvals(1:ind),20),'sortstr','descend','minpeakheight',0.2);
        lowval = yvals(lowind);
        if ~isempty(lowval)
            lowval = lowval(1);
            lowind = lowind(1);
            plot(xvals(lowind),lowval,'g^','markersize',15);
        end
        [dummy,hiind] = findpeaks(smooth(yvals(ind:end),20),'sortstr','descend','minpeakheight',0.2);
        hival = yvals(hiind+ind-1);
        if ~isempty(hival)
            hival = hival(1);
            hiind = hiind(1);
            plot(xvals(hiind+ind-1),hival,'bv','markersize',15);
        end
    else%Only tone data
        figure(1);
        yvals = S.ParamsOut.TDF.Tone_Difcor_Spline;
        yvals = yvals/max(abs(yvals));
        xvals = S.ParamsOut.TDF.x.ITD_spline;
        BD = S.ParamsOut.TDF.BD;
        plot(xvals,yvals,'k-','linewidth',1);hold on;
        ind = find(xvals==BD);
        plot(BD,yvals(ind),'rp','markersize',15);
        [dummy,lowind] = findpeaks(smooth(yvals(1:ind),20),'sortstr','descend','minpeakheight',0.2);
        lowval = yvals(lowind);
        if ~isempty(lowval)
            lowval = lowval(1);
            lowind = lowind(1);
            plot(xvals(lowind),lowval,'g^','markersize',15);
        end
        [dummy,hiind] = findpeaks(smooth(yvals(ind:end),20),'sortstr','descend','minpeakheight',0.2);
        hival = yvals(hiind+ind-1);
        if ~isempty(hival)
            hival = hival(1);
            hiind = hiind(1);
            plot(xvals(hiind+ind-1),hival,'bv','markersize',15);
        end
    end
    clf;
    ParamsOut = S.ParamsOut;
    if isfield(ParamsOut,'NDF')
        ParamsOut.NDF.LSP = [];
        ParamsOut.NDF.USP = [];
    elseif isfield(ParamsOut,'TDF')
        ParamsOut.TDF.LSP = [];
        ParamsOut.TDF.USP = [];
    end
    %     cd (BinRevDir);
    %     save (datafiles(i).name,'ParamsOut');
end