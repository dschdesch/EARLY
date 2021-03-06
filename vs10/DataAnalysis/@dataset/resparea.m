function [H,G]=resparea(D, figh, P)
% dataset/resparea - plot Response Area of a dataset
%    resparea(D) displays a Response Area of the spike times in dataset D.
%
%    resparea(D,figh) uses figure handle figh for plotting
%    (default = [] -> gcf). 
%
%    resparea(D, figh, P) uses parameters P for displaying the resparea.
%    P is typically a dataviewparam object or a valid 2nd input argument to
%    the dataviewparam constructor method, such as a parameter filename.
%
%    resparea is a standard "dataviewer", meaning that it may serve as
%    viewer for online data analysis during data collection. In addition,
%    the plot generated by all dataviewers allow an interactive change of
%    analysis parameter view the Params|Edit pullodwn menu (Ctr-Q).
%    For details on dataviewers, see dataviewparam.
%
%    See also dataviewparam, dataset/enableparamedit.

% recursive call using default plot params
if numel(D) > 1 
    clf;
    for ii=1:numel(D),
        resparea(D(ii));
    end
    return;
end

%===========single D from here=============

% handle the special case of parameter queries. Do this immediately to 
% avoid endless recursion with dataviewparam.
if isvoid(D) && isequal('params', figh),
    [H,G] = local_ParamGUI;
    return;
end

% open a new figure or use existing one?
if nargin<2 || isempty(figh),
    open_new = isempty(get(0,'CurrentFigure'));
    figh=gcf; 
else,
    open_new = isSingleHandle(figh);
end

% parameters
if nargin<3, P = []; end
if isempty(P), % use default paremeter set for this dataviewer
    P = dataviewparam(mfilename); 
end

% delegate the real work to local fcn
H = local_resparea(D, figh, open_new, P);

% enable parameter editing when viewing offline
if isSingleHandle(figh, 'figure'), enableparamedit(D, P, figh); end;


%============================================================
%============================================================
function data_struct = local_resparea(D, figh, open_new, P);
% the real work for the Response Area

if ~has2varparams(D)
    if isSingleHandle(figh, 'figure') % offline
        errordlg('No Response Area for datasets with no 2 independent parameters.','Not applicable'); 
        close(figh)
    else  % online
        warning('No Response Area for datasets with no 2 independent parameters.','Not applicable'); 
    end
    return;
end

% prepare plot
if isSingleHandle(figh, 'figure')
    figure(figh); clf; axes1 = gca;
    if open_new, placefig(figh, mfilename, D.Stim.GUIname); end % restore previous size 
else
    axes1 = axes('parent', figh);
end

% Check varied stimulus Params
Pres = D.Stim.Presentation;
P = struct(P); P = P.Param;
isortPlot = P.iCond(P.iCond<=Pres.Ncond); % limit to actual Ncond
if isortPlot==0, isortPlot = 1:Pres.Ncond; end;
Ncond = numel(isortPlot);
AW = P.Anwin;

Chan = 1; % digital channel
Nrep = NrepRec(D);
SPT = spiketimes(D, Chan, 'no-unwarp');
BurstDur = GenericStimparams(D,'BurstDur');
for icond=1:Ncond,
    if isequal('burstdur', AW)
        bdur = max(BurstDur(icond,:)); % burst dur in ms
        aw = [0 bdur]; 
    else
        aw = AW;
    end
    spt = AnWin(SPT(icond, :),aw);
    Nsp(icond) = numel([spt{:}]); 
    Rate(icond) = 1e3*Nsp(icond)./Nrep(icond)./bdur;
    
    data_struct.spt{icond} = spt;
    data_struct.Nsp(icond) = Nsp(icond);
    data_struct.Rate(Rate) = Rate(icond);
end; 

data_struct.Chan = 1;
data_struct.BurstDur = BurstDur;
data_struct.aw = aw;



X = D.Stim.Presentation.X;
Y = D.Stim.Presentation.Y;

[Nx, Ny] = DealElements(D.Stim.Ncond_XY);
ra = zeros(Ny,Nx);
for iy=1:Ny, % plot Rate against X for single Y value
    icond = (1:Nx)+Nx*(iy-1); % select conditions having current Y value
    ra(iy,:) = Rate(icond); %fill upside down
    %xplot(X.PlotVal(icond), Rate(icond), lico(iy));
    LegStr{iy} = sprintf(Y.FormatString, Y.PlotVal(icond(1)));
    xLabel=X.PlotVal(icond);
    yLabel(iy)=Y.PlotVal(icond(1));
end

data_struct.X = X;
data_struct.Y = Y;

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[749.9 4252]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[45 115]);
box(axes1,'on');
hold(axes1,'all');

image(xLabel,yLabel,ra,'Parent',axes1,'CDataMapping','scaled');
set(gca,'YDir','normal');
% Create xlabel
xlabel([lower(X.ParName) ' ('  X.ParUnit ')'],'fontsize',10);
% Create ylabel
ylabel([lower(Y.ParName) ' ('  Y.ParUnit ')'],'fontsize',10);
colorbar('peer',axes1);
H = Rate;

data_struct.xLabel = xLabel;
data_struct.yLabel = yLabel;
data_struct.ra = ra;
data_struct.xlabel = [lower(X.ParName) ' ('  X.ParUnit ')'];
data_struct.ylabel = [lower(Y.ParName) ' ('  Y.ParUnit ')'];

if nargout<1, clear H; end; % suppress unwanted echoing
title(IDstring(D, 'full'), 'fontsize', 12, 'fontweight', 'bold', 'interpreter', 'none');

function [T,G] = local_ParamGUI
% Returns the GUI for specifying the analysis parameters.
P = GUIpanel('resparea','');
iCond = ParamQuery('iCond', 'iCond:', '0', '', 'integer',...
    'Condition indices for which to calculate the Response Area. 0 means: all conditions.', 20);
Anwin = ParamQuery('Anwin', 'analysis window:', 'burstdur', '', 'anwin',...
    'Analysis window (in ms) [t0 t1] re the stimulus onset. The string "burstdur" means [0 t], in which t is the burst duration of the stimulus.');
P = add(P, iCond);
P = add(P, Anwin, below(iCond));
P = marginalize(P,[4 4]);
G = GUIpiece([mfilename '_parameters'],[],[0 0],[10 10]);
G = add(G,P);
G = marginalize(G,[10 10]);
% list all parameters in a struct
T = VoidStruct('iCond/Anwin');