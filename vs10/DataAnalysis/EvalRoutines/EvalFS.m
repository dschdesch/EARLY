function ArgOut = EvalFS(varargin)
%EVALFS   plot composite curve, interaural phase versus stimulus frequency, 
%         monaural phase versus stimulus frequency and ratecurve for two
%         cells with different CF. 
%         Datasets for both cells must contain responses of these cells to
%         monospectral tones. 
%   Data = EVALFS(ds1, ds2, ... )
%
%   Optional properties and their values can be given as a comma-separated
%   list. To view list of all possible properties and their default values,
%   use 'list' as only property. 
%
%   The result is clipped to 2 cycles by default. To change this, use the
%   property 'display'. Use Inf as value to display all cycles.

%B. Van de Sande 19-08-2004
%K. Spiritus     02-10-2007

%% ---------------- CHANGELOG -----------------------
%  Wed Jul 13 2011  Abel   
%  - Added support for new userdata syntax


%Parameter Lijst ...
DefParam.plot     = 'YES';
DefParam.plotmode = 'ALL';
DefParam.corbinwidth = 0.05;
DefParam.cormaxlag   = 30;
DefParam.coranwin    = [0 -1];
DefParam.corrunav    = 0.5;
DefParam.corcutoff   = 0.10;
DefParam.cyclenbin   = 64;
DefParam.fftrunav    = 25;
DefParam.linregtype  = 'W'; %Weighted or Non-weigthed ...
DefParam.display       = 2; % limit display to a certain amount of periods periods: yes or no

%By Abel: Specific options for THR using the userdata DB
DefParam.ignoreuserdata = false;    %Ignore userdata when looking for THR information (true/false)
DefParam.thrsequencenr = [];        %Force sequence NR for the THR curve (instead of last or NR from userdata DB)

%by Abel: force no diagonals (for identical datasets)
DefParam.nodiag = false;


%----------------------------------template--------------------------------------
Template.ds1.filename         = '';          
Template.ds1.icell            = NaN;         
Template.ds1.iseq             = NaN;         
Template.ds1.seqid            = '';          
Template.ds2.filename         = '';          
Template.ds2.icell            = NaN;         
Template.ds2.iseq             = NaN;         
Template.ds2.seqid            = '';          
Template.tag                  = 0;           %General purpose tag field
Template.createdby            = mfilename;   %Name of MATLAB function that generated the data
Template.thr1.cf               = NaN;         %Characteristic frequency retrieved from threshold curve
Template.thr1.sr               = NaN;         %Spontaneous rate retrieved from threshold curve
Template.thr1.thr              = NaN;         %Threshold at characteristic frequency
Template.thr1.q10              = NaN;         %Q10 retrieved from threshold curve
Template.thr1.bw               = NaN;         %Width 10dB above threshold (Hz)
Template.thr2.cf               = NaN;         %Characteristic frequency retrieved from threshold curve
Template.thr2.sr               = NaN;         %Spontaneous rate retrieved from threshold curve
Template.thr2.thr              = NaN;         %Threshold at characteristic frequency
Template.thr2.q10              = NaN;         %Q10 retrieved from threshold curve
Template.thr2.bw               = NaN;         %Width 10dB above threshold (Hz)
Template.DF                    = NaN;
Template.BestITD               = NaN;
Template.DeltaCF               = NaN;
Template.DeltaCoDist           = NaN;
Template.BinFreq               = NaN;
Template.R                     = NaN;
Template.Phase                 = NaN;
Template.RaySig                = NaN;
Template.N                     = NaN;


%Parameters evalueren ...
[ExpFolder, Cell1, Cell2, ds1, ds2, CalcParam, PlotParam, Param] = EvalParam(varargin, DefParam);

%By Abel: Calculate CF (from userdata) and format THR title string
exp_name = ds1.ID.Experiment.ID.Name;
Thr1 = getThr4Cell(ds1.ID.Experiment, Cell1.CellNr, Param.ignoreuserdata, Param.thrsequencenr);
CellParam(1).CF = Thr1.cf;
CellParam(1).SA = Thr1.sr;
CellParam(1).CD = greenwood(CellParam(1).CF);

Thr2 = getThr4Cell(ds2.ID.Experiment, Cell2.CellNr, Param.ignoreuserdata, Param.thrsequencenr);
CellParam(2).CF = Thr2.cf;
CellParam(2).SA = Thr2.sr;
CellParam(2).CD = greenwood(CellParam(2).CF);

%Gemeenschappelijke onafhankelijke variabelen opsporen ...
[StimFreq, ComSub] = findComSub(ds1, ds2); NSub = size(ComSub,2);
if isempty(ComSub)
    error('No common subsequences.'); 
elseif NSub == 1
    error('Only one common subsequence.');
else
    fprintf('Common subsequences are : '); 
    fprintf('%dHz,', ds1.Stim.Fcar(ComSub(1,:))); 
    fprintf('\b.\n'); 
end

%Voor elke gemeenschappelijke toon het CrossCorrelogram berekenen en daarop cyclehistogram nagaan.
%Bovendien ook op monaurale gegevens het cyclehistogram berekenen ...
DisSub = [];
fprintf('Calculating Correlogram and Cyclehistogram for : ');


[spt1,count1] = spiketimes(ds1);
%Remove all cells that contain no spikes
remove = find(~count1);
[spt2,count2] = spiketimes(ds2);
remove = [remove find(~count2)]; 
remove = unique(remove);
StimFreq(remove) = [];
ComSub(:,remove) = [];
NSub = size(ComSub,2);

for n = 1:NSub
    fprintf([int2str(StimFreq(n)) 'Hz,']);
    SpkTr1 = spt1(ComSub(1, n), :);
	SpkTr2 = spt2(ComSub(2, n), :);
    N = length(SpkTr1);
	M = length(SpkTr2);
    AnDur = abs(diff(CalcParam.coranwin));
	
	if isequal(SpkTr1, SpkTr2)
		[Y, X] = SPTCORR(...
			anwin(SpkTr1, CalcParam.coranwin),...
			'nodiag', CalcParam.cormaxlag, CalcParam.corbinwidth);
		NRep = N*(N-1);
		Y = (1000*Y)/(AnDur*NRep);
	elseif (Param.nodiag)
		%By Abel: Force no diagonals on 2 different datasets
		%see code from SPTCORR()
		%- Get spikes in analysis window
		spt1 = anwin(SpkTr1, CalcParam.coranwin);
		spt2 = anwin(SpkTr2, CalcParam.coranwin);
		maxlag = CalcParam.cormaxlag;
		binwidth = CalcParam.corbinwidth;
		%- SPTCORR without normalization
		%  Y = histogram of the spike-time differences
		%      between the spike pairs from SPT1 and SPT2
		%  X = bin centers
		[Y, X] = SPTCORR([spt1{:}], [spt2{:}],...
			maxlag, binwidth);
		%- Substract diagonal terms
		for irep=1:length(spt1)
			Yd = SPTCORR(spt1{irep}, spt2{irep},...
				maxlag, binwidth);
			Y = Y - Yd;
		end
		NRep = N*M;
		Y = (1000*Y)/(AnDur*NRep);
    else
        if n==17
            disp(n)
        end
		[Y, X] = SPTCORR(...
			anwin(SpkTr1, CalcParam.coranwin),...
			anwin(SpkTr2, CalcParam.coranwin),...
			CalcParam.cormaxlag, CalcParam.corbinwidth);
		NRep = N*M;
		Y = (1000*Y)/(AnDur*NRep);
	end
	CrossCor(n) = CollectInStruct(X, Y);
    
    %N = size(ds1.spt, 2); M = size(ds2.spt, 2); if isequal(ds1.spt, ds2.spt), NRep = N*(N-1); else NRep = N*M; end
    BinCycleHist(n) = SGSR_cyclehist(X, Y, CalcParam.coranwin(2)-CalcParam.coranwin(1), NRep, StimFreq(n), CalcParam.cyclenbin);
    
    if all(CrossCor(n).Y == 0)
        fprintf('\b discarded,'); 
        DisSub = [DisSub n]; 
    end
end

fprintf('\b.\n');


if ~isempty(DisSub) %Indien er subsequenties wegvallen ... 
    StimFreq(DisSub) = [];
    CrossCor(DisSub) = [];
    BinCycleHist(DisSub) = [];
    ComSub(:, DisSub) = [];
end

if length(StimFreq) == 1
    error('Only one common subsequence left.'); 
end

for n = 1:length(ds1.Stim.Fcar)
    MonCycleHist1(n) = SGSR_cyclehist(ds1, n, -1, CalcParam.cyclenbin, CalcParam.coranwin(1), CalcParam.coranwin(2));
    if sum(MonCycleHist1(n).Y) == 0
       fprintf('Subssequence with frequency %dHz has no spikes in dataset 1\n',ds1.Stim.Fcar(n));
       MonCycleHist1(n).Y = [];
    end
end

for n = 1:length(ds2.Stim.Fcar)
    MonCycleHist2(n) = SGSR_cyclehist(ds2, n, -1, CalcParam.cyclenbin, CalcParam.coranwin(1), CalcParam.coranwin(2));
    if sum(MonCycleHist1(n).Y) == 0
       fprintf('Subssequence with frequency %dHz has no spikes in dataset 2\n',ds1.Stim.Fcar(n));
       MonCycleHist2(n).Y = [];
    end
end



%Composite Curve berekenen ...
X = CrossCor(1).X;
Y = mean(cat(1, CrossCor.Y));

%Indien MaxLag groter is dan 10ms dan zoeken naar hoofdpiek en niet naar die van de
%aanwezige harmonieken op 10ms (100Hz) en 20ms(50Hz)...
if CalcParam.cormaxlag > 10
    BestITD = getmaxloc(X, Y, CalcParam.corrunav, [-7, +7]);
else
    BestITD = getmaxloc(X, Y, CalcParam.corrunav); 
end

CFdiff  = log2(CellParam(2).CF/CellParam(1).CF);
CDdiff  = CellParam(2).CD - CellParam(1).CD;

CC = CollectInStruct(X, Y, BestITD, CFdiff, CDdiff);

%Gegevens uit de composite curve halen ...
%CAVE: Fourier.m houdt geen rekening met het feit dat er harmonische pieken aanwezig zijn
%in de Composite Curve ...
[XFft, YFft, YFftAv, DF] = spectana(X, Y, CalcParam.fftrunav);
FFT.X = XFft; 
FFT.Y = [YFft; YFftAv]; 
FFT.DF = DF;

%SuperImpose Curve berekenen ... van elk crosscorrelogram 2 cycli overhouden en elk normaliseren naar
%zijn maximum ... Bovendien worden slechts die correlogrammen geplot die boven bepaalde cutoff rate komen ...
Y = cat(1, CrossCor.Y);
Ymax = max(Y, [], 2); 
iExc = find(Ymax < CalcParam.corcutoff);
Y(iExc, :) = []; 

Freq = StimFreq; 
Freq(iExc) = [];
Period = repmat(1000./Freq, 1, size(Y, 2));
Xtemp = repmat(X, size(Y, 1), 1);
display = CalcParam.display * Period / 2;
iExc = find((Xtemp <= -display) | (Xtemp >= display));
Y(iExc)= NaN;

%No normalization ... 16-02-2004
%Ymax = max(Y, [], 2);
%Y = Y ./ repmat(Ymax, 1, size(Y, 2));

if isempty(Y)
    warning('No curves left for super impose plot.'); 
    Y = repmat(NaN, 1, length(X)); 
end
SI = CollectInStruct(X,Y);

%Interaurale VS magnitude curve berekenen ...
X = StimFreq';
Y = cat(2, BinCycleHist.R);
pRayleigh = cat(2, BinCycleHist.pRaySig);
MaxR = max(Y);

IR = CollectInStruct(X, Y, pRayleigh, MaxR);

%Monaurale VS magnitude curven berekenen ...
X = ds1.Stim.Fcar(1:ds1.Stim.Presentation.Ncond)';
Y = cat(2, MonCycleHist1.R);
pRayleigh = cat(2, MonCycleHist1.pRaySig);
MaxR = max(Y);

MR(1) = CollectInStruct(X, Y, pRayleigh, MaxR);

X = ds2.Stim.Fcar(1:ds2.Stim.Presentation.Ncond)';
Y = cat(2, MonCycleHist2.R);
pRayleigh = cat(2, MonCycleHist2.pRaySig);
MaxR = max(Y);

MR(2) = CollectInStruct(X, Y, pRayleigh, MaxR);

%Response area curves berekenen ... 
MRA(1).X = ds1.Stim.Fcar(1:ds1.Stim.Presentation.Ncond)';
MRA(1).Y = getrate(ds1, 1:ds1.Stim.Presentation.Ncond, CalcParam.coranwin(1), CalcParam.coranwin(2));

MRA(2).X = ds2.Stim.Fcar(1:ds2.Stim.Presentation.Ncond)';
MRA(2).Y = getrate(ds2, 1:ds2.Stim.Presentation.Ncond, CalcParam.coranwin(1), CalcParam.coranwin(2));

IRA.X = StimFreq';
IRA.Y = mean(cat(1, BinCycleHist.Y), 2);
IRA.Max = max(IRA.Y(:));

%Interaurale Synchronicity Rate berekenen ...
X = StimFreq';
Y = IRA.Y' .* IR.Y;
pRayleigh = IR.pRayleigh;
Max = max(Y); 
SRleft = Y(1); SRright = Y(end);
NormSRleft = SRleft/Max; NormSRright = SRright/Max;

ISR = CollectInStruct(X, Y, pRayleigh, Max, SRleft, SRright, NormSRleft, NormSRright);

%Interaurale phase curve berekenen ...
X = StimFreq';
Y = cat(2, BinCycleHist.Ph);
Y = unwrap(Y * (2*pi))/2/pi;

%Gegevens uit interaurale phase curve halen ...
pRayleigh = cat(2, BinCycleHist.pRaySig);
iSign = find(pRayleigh <= 0.001);
if length(X(iSign)) > 1
    if strncmpi(CalcParam.linregtype, 'W', 1)
        %Gewogen gemiddelde naar syncrate ...
        P = linregfit(X(iSign), Y(iSign), ISR.Y(iSign));
        [pLinReg, MSerror, DF] = signlinreg(P, X(iSign), Y(iSign), ISR.Y(iSign));
    else
        P = polyfit(X(iSign), Y(iSign), 1); 
        [pLinReg, MSerror, DF] = signlinreg(P, X(iSign), Y(iSign));
    end
    CD = P(1)*1000; CP = P(2);
else
    [CD, CP, pLinReg, MSerror, DF] = deal(NaN); 
end

IPC = CollectInStruct(X, Y, pRayleigh, CD, CP, pLinReg, MSerror, DF);

%Monaurale Synchronicity Rate berekenen ...
X = ds1.Stim.Fcar(1:ds1.Stim.Presentation.Ncond)';
idx = find(ismember(MRA(1).X, intersect(MRA(1).X, MR(1).X)));
Y = MRA(1).Y(idx) .* MR(1).Y;
pRayleigh = MR(1).pRayleigh;
Max = max(Y);

MSR(1) = CollectInStruct(X, Y, pRayleigh, Max);

X = ds2.Stim.Fcar(1:ds2.Stim.Presentation.Ncond)';
idx = find(ismember(MRA(2).X, intersect(MRA(2).X, MR(2).X)));
Y = MRA(2).Y(idx) .* MR(2).Y;
pRayleigh = MR(2).pRayleigh;
Max = max(Y);

MSR(2) = CollectInStruct(X, Y, pRayleigh, Max);

%Monaurale phase curven berekenen ...
X = ds1.Stim.Fcar(1:ds1.Stim.Presentation.Ncond)';
Y = cat(2, MonCycleHist1.Ph);
Y = unwrap(Y * (2*pi))/2/pi;
pRayleigh = cat(2, MonCycleHist1.pRaySig);

iSign = find(pRayleigh <= 0.001);
if length(X(iSign)) > 1,
    if strncmpi(CalcParam.linregtype, 'W', 1)
        %Gewogen gemiddelde naar syncrate ...
        P = linregfit(X(iSign), Y(iSign), MSR(1).Y(iSign));
        [pLinReg, MSerror, DF] = signlinreg(P, X(iSign), Y(iSign), MSR(1).Y(iSign));
    else, 
        P = polyfit(X(iSign), Y(iSign), 1); 
        [pLinReg, MSerror, DF] = signlinreg(P, X(iSign), Y(iSign));
    end
    Slope = P(1)*1000; YInterSect = P(2);
else, [Slope, YInterSect, pLinReg, MSerror, DF] = deal(NaN); end

MPC(1) = CollectInStruct(X, Y, pRayleigh, Slope, YInterSect, pLinReg, MSerror, DF);

X = ds2.Stim.Fcar(1:ds2.Stim.Presentation.Ncond)';
Y = cat(2, MonCycleHist2.Ph);
Y = unwrap(Y * (2*pi))/2/pi;
pRayleigh = cat(2, MonCycleHist2.pRaySig);

iSign = find(pRayleigh <= 0.001);
if length(X(iSign)) > 1,
    if strncmpi(CalcParam.linregtype, 'W', 1)
        %Gewogen gemiddelde naar syncrate ...
        P = linregfit(X(iSign), Y(iSign), MSR(2).Y(iSign));
        [pLinReg, MSerror, DF] = signlinreg(P, X(iSign), Y(iSign), MSR(2).Y(iSign));
    else, 
        P = polyfit(X(iSign), Y(iSign), 1);
        [pLinReg, MSerror, DF] = signlinreg(P, X(iSign), Y(iSign));
    end
    Slope = P(1)*1000; YInterSect = P(2);
else [Slope, YInterSect, pLinReg, MSerror, DF] = deal(NaN); end   

MPC(2) = CollectInStruct(X, Y, pRayleigh, Slope, YInterSect, pLinReg, MSerror, DF);

%Monaurale phaseverschil ...
X = intersect(MPC(1).X, MPC(2).X);
idx1 = find(ismember(round(MPC(1).X), round(X)));
idx2 = find(ismember(round(MPC(2).X), round(X)));
Y = MPC(1).Y(idx1) - MPC(2).Y(idx2);
idx = find(~isnan(Y)); X = X(idx); Y = Y(idx)';

iRayleigh = find(all([MPC(1).pRayleigh(idx1); MPC(2).pRayleigh(idx2)] <= 0.001));

if size(X) == size(Y')
    Y=Y';
end
P = polyfit(X, Y, 1);
Slope = P(1)*1000; YInterSect = P(2);
idx = find(ismember(X, IPC.X));
[pCorrCoef, CorrCoef] = signcorr(IPC.Y, Y); 

PD = CollectInStruct(X, Y, iRayleigh, Slope, YInterSect, pCorrCoef, CorrCoef);

%Gegevens herorganiseren ...
BurstDur  = ds1.Stim.BurstDur;
RepDur    = ds1.Stim.ISI;
StimParam = CollectInStruct(BurstDur, RepDur, StimFreq);

Param = CollectInStruct(CellParam, StimParam, CalcParam);
CellInfo = CollectInStruct(exp_name, Cell1, Cell2);
CalcData = CollectInStruct(CrossCor, BinCycleHist, MonCycleHist1, MonCycleHist2, CC, SI, FFT, IPC, IR, MR, MPC, MRA, IRA, ISR, MSR, PD);

% Extract the comments if available to add next to the plots
comment1 = '';
comment2 = '';
if isfield(ds1.ID,'comment')
    comment1 = ds1.ID.comment;
end
if isfield(ds2.ID,'comment')
    comment2 = ds2.ID.comment;
end

if strcmpi(PlotParam.plot, 'YES')
    if any(strcmpi(PlotParam.plotmode, {'ALL', 'CC'})) 
        %Overzicht van alle CrossCorrelogrammen met gecentralizeerd de CompositeCurve ...
        PlotCC(CC, SI, CrossCor, FFT, CellInfo, Param, comment1,comment2);
    end
    if any(strcmpi(PlotParam.plotmode, {'ALL', 'IPC'}))
        %Overzicht van alle CycleHistogrammen met gecentraliseerd de Interaurale phase curve ...
        PlotIPC(IPC, IR, BinCycleHist, CellInfo, Param, comment1,comment2);
    end
    if any(strcmpi(PlotParam.plotmode, {'ALL', 'MPC'}))
        %Overzicht van alle CycleHistogrammen met gecentraliseerd de monaurale phase curve voor beide cellen ...
        PlotMPC(1, MPC, MonCycleHist1, CellInfo, Param, comment1);
        PlotMPC(2, MPC, MonCycleHist2, CellInfo, Param, comment2);
    end
    if any(strcmpi(PlotParam.plotmode, {'ALL', 'PRA'}))
        %Monaurale phaseverschil vergelijken met interaurale phase ...
        PlotPRA(MRA, IRA, MPC, IPC, IR, MR, MSR, ISR, PD, CellInfo, Param, comment1,comment2);
    end
end

if nargout >= 1
    [CalcData.ds1, CalcData.ds2, CalcData.CalcParam, CalcData.PlotParam, CalcData.stim, CalcData.thr1, ...
        CalcData.thr2, CalcData.DF, CalcData.BestITD, CalcData.DeltaCF, CalcData.DeltaCoDist, ...
        CalcData.BinFreq, CalcData.R, CalcData.Phase, CalcData.RaySig, CalcData.N] = deal(ds1, ds2,...
        CalcParam, PlotParam, StimParam, Thr1, Thr2, FFT.DF, BestITD, CFdiff, CDdiff, cat(2, Param.StimParam.StimFreq), ...
        cat(2, BinCycleHist.R), cat(2, BinCycleHist.Ph), cat(2, BinCycleHist.pRaySig), cat(2, BinCycleHist.N));

    ArgOut = structtemplate(CalcData, Template, 'reduction', 'off');
end

%-------------------------------------------------------------------------------------------------------------------%
%                                             LOCALE FUNCTIES                                                       % 
%-------------------------------------------------------------------------------------------------------------------%
%-----------%
% EVALPARAM %
%-----------%
function [DataFolder, Cell1, Cell2, ds1, ds2, CalcParam, PlotParam, Param] = EvalParam(ParamList, DefParam)
%EVALFS(ds1, ds2, ...)

if length(ParamList) < 2, error('Wrong number of input parameters.'); end

ds1 = ParamList{1};
ds2 = ParamList{2};

%Nagaan of datasets wel van zelfde experiment afkomstig zijn ...
exp_name1 = ds1.ID.Experiment.ID.Name;
exp_name2 = ds2.ID.Experiment.ID.Name;
if ~strcmp(exp_name1, exp_name2)
	error('Datasets are from different experiment.'); 
end
DataFolder = folder(ds1.ID.Experiment);

Cell1 = CreateCellStruct(ds1);
Cell2 = CreateCellStruct(ds2);

ParamList(1:2) = [];

%Optionele lijst van property/values evalueren ...
Param = checkproplist(DefParam, ParamList{:});
    
%ParameterLijst evalueren ...
if ~isnumeric(Param.display)
    error('Property display must be numeric.'); 
end

Param.plot = upper(Param.plot);
if ~any(strcmpi(Param.plot, {'YES', 'NO'}))
    error('Property Plot must be YES or NO.'); 
end

Param.plotmode = upper(Param.plotmode);
if ~any(strcmpi(Param.plotmode, {'ALL', 'IPC', 'CC', 'MPC', 'PRA'})), error('Wrong mode. Mode should be ''ALL'', ''IPC'', ''CC'', ''MPC'' or ''PRA''.'); end

Param.linregtype = upper(Param.linregtype);
if ~any(strncmpi(Param.linregtype, {'W', 'N'}, 1)), error('Wrong value for property linregtype.'); end

if (length(Param.corcutoff) ~= 1) | (Param.corcutoff < 0), error('Wrong value for property corcutoff.'); end

%Nagaan of datasets wel dezelfde stimulustype hebben en of ze wel voor dit project bedoeld zijn ...
if ~any(strcmp('FS', { ds1.Stim.StimType, ds2.Stim.StimType }))
    error('Wrong datasets.');
end

%Nagaan of datasets dezelfde stimulusparameters hebben ...
if CheckStimParam(ds1, ds2), error('Stimulus Parameters are different for both datasets.'); end

%Binaurale stimuli zijn niet toegelaten als ze verschillende frequenties
%voor elk oor hebben
if strcmp(ds1.Stim.DAC,'Both') || strcmp(ds2.Stim.DAC,'Both')
    if sum(ds1.Stim.Fcar(:,1)==ds1.Stim.Fcar(:,2))==0  && ...
            sum(ds2.Stim.Fcar(:,1)==ds2.Stim.Fcar(:,2))==0
        error('SGSR:ERROR','EvalFS can only be used for monaural stimuli!');
    end
end

%Einde analysewindow kan -1 zijn, wat staat voor de stimulusduur ...
if Param.coranwin(2) == -1, Param.coranwin(2) = ds1.Stim.BurstDur; end
    
%PlotParam-structuur opstellen ...
CalcParam = getfields(Param, {'corbinwidth', 'cormaxlag', 'coranwin', 'corrunav', 'corcutoff', 'cyclenbin', 'fftrunav', 'linregtype', 'display'});
PlotParam = getfields(Param, {'plot', 'plotmode'});

%----------------%
% CHECKSTIMPARAM %
%----------------%
function Err = CheckStimParam(ds1, ds2)

Err = 0;

%Kijken of stimulusparameters gelijk zijn tussen de cellen en of beide cellen zelfde SPL hebben ...
if ~isequal(ds1.Stim.BurstDur, ds2.Stim.BurstDur) | ...
   ~isequal(ds1.Stim.ISI, ds2.Stim.ISI)     | ...
   ~isequal(ds1.Stim.SPL, ds2.Stim.SPL)
    Err = 1;
    return
end

%------------%
% FINDCOMSUB %
%------------%
% Finds the Common Independant variables in this case freq
function [ComIndepVar, ComSub] = findComSub(varargin)

ds1freq = varargin{1}.Stim.Fcar;
ds2freq = varargin{2}.Stim.Fcar;
[common_indexes indices_loc] = ismember(round(ds1freq),round(ds2freq));
ComIndepVar = ds1freq(common_indexes(:,1));

if isempty(ComIndepVar), ComSub = []; return; end

ComSub(1,:) = find(common_indexes(:,1)==1);
ComSub(2,:) = indices_loc(ComSub(1,:));

%--------%
% PLOTCC %
%--------%
function PlotCC(ComCurve, SuperImpose, CrossCor, FFT, CellInfo, Param, comment1, comment2)

DataID  = [ CellInfo.exp_name ' <' CellInfo.Cell1.dsID '> & <' CellInfo.Cell2.dsID '>' ];

Interface = figure('Name', ['EvalFS: Composite Curve for ' DataID], ...
    'NumberTitle', 'off', ...
    'PaperType', 'a4letter', ...
    'PaperUnits', 'normalized', ...
    'PaperPosition', [0.05 0.05 0.90 0.90], ...
    'PaperOrientation', 'landscape');

Ax_CC = axes('Position', [0.195 0.59 0.3025 0.37]); 
line(ComCurve.X, ComCurve.Y, 'LineStyle', '-', 'Color', 'r', 'Marker', 'none');
title(['Composite Curve for ' DataID], 'FontSize', 12);
xlabel('Delay(ms)'); ylabel('Coincidence Rate(spk/sec)');

MinX = -5; MaxX = +5;
MinY = 0; MaxY = 1.5; 
axis([MinX MaxX MinY MaxY]); 

text(MinX, MaxY, {sprintf('DF : %.0fHz', FFT.DF); ...
                  sprintf('BestITD : %.3fms', ComCurve.BestITD)}, ...
                  'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

line([0 0], [MinY MaxY], 'LineStyle', ':', 'Color', 'k');
line([ComCurve.BestITD ComCurve.BestITD], [MinY MaxY], 'LineStyle', ':', 'Color', 'k');

Ax_SI = axes('Position', [0.5425 0.59 0.3025 0.37]);
line(SuperImpose.X, SuperImpose.Y, 'LineStyle', '-', 'Marker', '.');
title(['Superimpose Curve for ' DataID], 'FontSize', 12);
xlabel('Delay(ms)'); ylabel('Coincidence Rate (spk/sec)');

if size(SuperImpose.Y, 1) == 1,
    iMinX = min(find(~isnan(SuperImpose.Y)));
    iMaxX = max(find(~isnan(SuperImpose.Y)));
else
    iMinX = min(find(any(~isnan(SuperImpose.Y))));
    iMaxX = max(find(any(~isnan(SuperImpose.Y))));
end    
MinX = SuperImpose.X(iMinX); MaxX = SuperImpose.X(iMaxX);
MinY = 0; MaxY = max(SuperImpose.Y(:)); 
if isempty(MinX), MinX = 0; end
if isempty(MaxX), MaxX = 1; end
if MinY == MaxY
    MaxY = 1;
end
axis([MinX MaxX MinY MaxY]); 

%Weergeven van plotparameters en cellparameters naast Composite Curve ...
printinfo([0.02 0.59 0.125 0.37], {['<' CellInfo.Cell1.dsID '>']; ...
        sprintf('CF = %.0fHz', Param.CellParam(1).CF); ...
        sprintf('SR = %.0fspk/sec', Param.CellParam(1).SA); ...
        sprintf('Comment: %s',comment1); ...
        ''; ...
        ['<' CellInfo.Cell2.dsID '>']; ...
        sprintf('CF = %.0fHz', Param.CellParam(2).CF); ...
        sprintf('SR = %.0fspk/sec', Param.CellParam(2).SA); ...
        sprintf('Comment: %s',comment2); ...
        ''; ...
        ['\DeltaCF = ' sprintf('%.3f', log2(Param.CellParam(2).CF/Param.CellParam(1).CF)) 'oct']; ...
        ['\DeltaCoDist = ' sprintf('%.2f', (Param.CellParam(2).CD - Param.CellParam(1).CD)) 'mm']});
printinfo([0.855  0.59 0.125 0.37], {'Stimulus Parameters:'; ...
        sprintf('StimDur = %dms', Param.StimParam.BurstDur); ...
        sprintf('RepDur  = %dms', Param.StimParam.RepDur); ...
        ''; ...
        'CrossCorrelation:'; ...
        sprintf('BinWidth = %.3fms', Param.CalcParam.corbinwidth); ...
        sprintf('MaxLag = %dms', Param.CalcParam.cormaxlag); ...
        sprintf('ANWIN = [%dms-%dms]', Param.CalcParam.coranwin(1), Param.CalcParam.coranwin(2))});


%Alle crosscorrelogrammen plotten ...
PlotInfo.XLabel = 'Delay(ms)'; PlotInfo.XLim = [-5 +5];
PlotInfo.YLabel = 'Rate(spk/sec)';
PlotInfo.Text = '[''StimFreq = '' num2str(PlotInfo.StimFreq(n)) ''Hz'']';
PlotInfo.StimFreq = Param.StimParam.StimFreq;
PlotInfo.Type = 'PLOT';
PlotSmallCurves(CrossCor, PlotInfo);

%---------%
% PLOTIPC %
%---------%
function PlotIPC(IPC, IR, BinCycleHist, CellInfo, Param, comment1, comment2)

DataID  = [ CellInfo.exp_name ' <' CellInfo.Cell1.dsID '> & <' CellInfo.Cell2.dsID '>' ];

%Weergeven van de Interaurale Phase Curve ...
Interface = figure('Name', ['EvalFS: Interaural Phase Curve for ' DataID], ...
    'NumberTitle', 'off', ...
    'PaperType', 'a4letter', ...
    'PaperUnits', 'normalized', ...
    'PaperPosition', [0.05 0.05 0.90 0.90], ...
    'PaperOrientation', 'landscape');

Ax_IPC = axes('Position', [0.195 0.59 0.3025 0.37]); 
line(IPC.X, IPC.Y, 'LineStyle', '-', 'Color', 'r', 'Marker', 'o');
title(['Interaural Phase Curve for ' DataID], 'FontSize', 12);
xlabel('Stimulus Frequency(Hz)'); ylabel('Interaural Phase(cycles)');

MinX = 0; MaxX = IPC.X(end);
MinY = floor(min(IPC.Y)); MaxY = ceil(max(IPC.Y));
axis([MinX MaxX MinY MaxY]);

text(MinX, MaxY, {sprintf('CD : %.3fms', IPC.CD); ...
                  sprintf('CP : %.3f', IPC.CP); ...
                  sprintf('MSerror : %.3f', IPC.MSerror); ...
                  sprintf('pLinReg : %.3f', IPC.pLinReg)}, ...
                  'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

line([0, MaxX], [IPC.CP, polyval([IPC.CD/1000 IPC.CP], MaxX)], 'Color', 'k', 'LineStyle', ':', 'Marker', 'none');

%Niet significante waarden donker plotten
idx = find(IPC.pRayleigh > 0.001); line(IPC.X(idx), IPC.Y(idx), 'Color', [1 0 0], 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 24);

Ax_IR = axes('Position', [0.5425 0.59 0.3025 0.37]); 
line(IR.X, IR.Y, 'LineStyle', '-', 'Color', 'r', 'Marker', 'o');
title(['Interaural R for ' DataID], 'FontSize', 12);
xlabel('Stimulus Frequency(Hz)'); ylabel('R');

MinX = IR.X(1); MaxX = IR.X(end);
MinY = 0; MaxY = 1;
axis([MinX MaxX MinY MaxY]);

text(MinX, MaxY, sprintf('Rmax = %.3f', IR.MaxR), ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
         
%Niet significante waarden donker plotten
idx = find(IR.pRayleigh > 0.001); line(IR.X(idx), IR.Y(idx), 'Color', [1 0 0], 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 24);

%Weergeven van plotparameters en cellparameters naast Interaurale Phase Curve ...
printinfo([0.02 0.59 0.125 0.37], {['<' CellInfo.Cell1.dsID '>']; ...
        sprintf('CF = %.0fHz', Param.CellParam(1).CF); ...
        sprintf('SR = %.0fspk/sec', Param.CellParam(1).SA); ...
        sprintf('Comment: %s',comment1); ...
        ''; ...
        ['<' CellInfo.Cell2.dsID '>']; ...
        sprintf('CF = %.0fHz', Param.CellParam(2).CF); ...
        sprintf('SR = %.0fspk/sec', Param.CellParam(2).SA); ...
        sprintf('Comment: %s',comment2); ...
        ''; ...
        ['\DeltaCF = ' sprintf('%.3f', log2(Param.CellParam(2).CF/Param.CellParam(1).CF)) 'oct']; ...
        ['\DeltaCoDist = ' sprintf('%.2f', (Param.CellParam(2).CD - Param.CellParam(1).CD)) 'mm']});
printinfo([0.855  0.59 0.125 0.37], {'Stimulus Parameters:'; ...
        sprintf('StimDur = %dms', Param.StimParam.BurstDur); ...
        sprintf('RepDur  = %dms', Param.StimParam.RepDur); ...
        ''; ...
        'CrossCorrelation:'; ...
        sprintf('BinWidth = %.3fms', Param.CalcParam.corbinwidth); ...
        sprintf('MaxLag = %dms', Param.CalcParam.cormaxlag); ...
        sprintf('ANWIN = [%dms-%dms]', Param.CalcParam.coranwin(1), Param.CalcParam.coranwin(2)); ...
        ''; ...
        'CycleHistogram:'; ...
        sprintf('Nr. of Bins = %d', Param.CalcParam.cyclenbin)});

%Plotten van alle cyclehistogrammen
PlotInfo.XLabel = 'Cycle'; PlotInfo.XLim = [0 1];
PlotInfo.YLabel = 'Rate(spk/sec)';
PlotInfo.Text = ['{[''BF = '' int2str(PlotInfo.BinFreq(n)) ''Hz'' ];' ...
                 ' [''R = '' sprintf(''%.3f'', PlotInfo.R(n))];' ...
                 ' [''\phi = '' sprintf(''%.2f'', PlotInfo.Phase(n))];' ...
                 ' [''RaySig = '' num2str(PlotInfo.RaySig(n))];' ...
                 ' [''N = '' int2str(PlotInfo.N(n))]}' ];
PlotInfo.BinFreq = Param.StimParam.StimFreq;
PlotInfo.R       = cat(1, BinCycleHist.R);
PlotInfo.Phase   = cat(1, BinCycleHist.Ph);
PlotInfo.RaySig  = cat(1, BinCycleHist.pRaySig);
PlotInfo.N       = cat(1, BinCycleHist.N);
PlotInfo.Type = 'BAR';
PlotSmallCurves(BinCycleHist, PlotInfo);

%---------%
% PLOTMPC %
%---------%
function PlotMPC(CellNr, MPC, MonCycleHist, CellInfo, Param, comment)

DataID = eval(sprintf('[CellInfo.exp_name '' <'' CellInfo.Cell%d.dsID ''>''];', CellNr));

%Weergeven van de monaurale phase curve ...
MinX = 0; MaxX = MPC(CellNr).X(end);
MinY = floor(min(MPC(CellNr).Y)); MaxY = ceil(max(MPC(CellNr).Y));

PlotInfo.Title  = ['Monaural Phase Curve for ' DataID];
PlotInfo.XLabel = 'Stimulus Frequency(Hz)'; PlotInfo.XLim = [MinX MaxX];
PlotInfo.YLabel = 'Phase(cycles)'; PlotInfo.YLim = [MinY MaxY];
PlotInfo.Text   = {sprintf('Slope : %.3fms', MPC(CellNr).Slope); ...
                   sprintf('YInterSect : %.3f', MPC(CellNr).YInterSect); ...
                   sprintf('MSerror : %.3f', MPC(CellNr).MSerror); ...
                   sprintf('pLinReg : %.3f', MPC(CellNr).pLinReg)};
PlotInfo.Style  = 'ro-';
Interface = PlotMainCurve(struct('X', MPC(CellNr).X, 'Y', MPC(CellNr).Y), PlotInfo);

%Niet significante waarden donker plotten
idx = find(MPC(CellNr).pRayleigh > 0.001); 
line(MPC(CellNr).X(idx), MPC(CellNr).Y(idx), 'Color', [1 0 0], 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 24);

line([0, MaxX], [MPC(CellNr).YInterSect, polyval([MPC(CellNr).Slope/1000 MPC(CellNr).YInterSect], MaxX)], 'Color', 'k', 'LineStyle', ':', 'Marker', 'none');

MinX = 0; MaxX = MPC(CellNr).X(end);
MinY = floor(min(MPC(CellNr).YInterSect)); MaxY = ceil(max(MPC(CellNr).Y));
if isnan(MinY)
    MinY = 0;
end
axis([MinX MaxX MinY MaxY]);

%Weergeven van plotparameters en cellparameters  ...
s = eval(sprintf('[ ''<'' CellInfo.Cell%d.dsID ''>'']', CellNr));
printinfo([0.02 0.59 0.125 0.37], {s; ...
        sprintf('CF = %.0fHz', Param.CellParam(CellNr).CF); ...
        sprintf('SR = %.0fspk/sec', Param.CellParam(CellNr).SA);...
        sprintf('Comment: %s',comment)}); ...
printinfo([0.855  0.59 0.125 0.37], {'Stimulus Parameters:'; ...
        sprintf('StimDur = %dms', Param.StimParam.BurstDur); ...
        sprintf('RepDur  = %dms', Param.StimParam.RepDur); ...
        ''; ...
        'CycleHistogram:'; ...
        sprintf('ANWIN = [%dms-%dms]', Param.CalcParam.coranwin); ...
        sprintf('Nr. of Bins = %d', Param.CalcParam.cyclenbin)});

%Plotten van alle cyclehistogrammen
PlotInfo.XLabel = 'Cycle'; PlotInfo.XLim = [0 1];
PlotInfo.YLabel = 'Rate(spk/sec)';
PlotInfo.Text = ['{[''BF = '' int2str(PlotInfo.BinFreq(n)) ''Hz'' ];' ...
                 ' [''R = '' sprintf(''%.3f'', PlotInfo.R(n))];' ...
                 ' [''\phi = '' sprintf(''%.3f'', PlotInfo.Phase(n))];' ...
                 ' [''RaySig = '' num2str(PlotInfo.RaySig(n))];' ...
                 ' [''N = '' int2str(PlotInfo.N(n))]}' ];
PlotInfo.BinFreq = MPC(CellNr).X;
PlotInfo.R       = cat(1, MonCycleHist.R);
PlotInfo.Phase   = cat(1, MonCycleHist.Ph);
PlotInfo.RaySig  = cat(1, MonCycleHist.pRaySig);
PlotInfo.N       = cat(1, MonCycleHist.N);
PlotInfo.Type = 'BAR';
PlotSmallCurves(MonCycleHist, PlotInfo);

%---------%
% PLOTPRA %
%---------%
function PlotPRA(MRA, IRA, MPC, IPC, IR, MR, MSR, ISR, PD, CellInfo, Param, comment1, comment2);

DataID  = [ CellInfo.exp_name ' <' CellInfo.Cell1.dsID '> & <' CellInfo.Cell2.dsID '>' ];
DataID1 = [ CellInfo.exp_name ' <' CellInfo.Cell1.dsID '>' ];
DataID2 = [ CellInfo.exp_name ' <' CellInfo.Cell2.dsID '>' ];

Interface = figure('Name', ['EvalFS: Response Area, Phase, R and SyncRate for ' DataID ], ...
    'NumberTitle', 'off', ...
    'PaperType', 'a4letter', ...
    'PaperUnits', 'normalized', ...
    'PaperPosition', [0.05 0.05 0.90 0.90], ...
    'PaperOrientation', 'landscape');

StimFreqs = unique(cat(2, MRA.X));
MinFreq = min(StimFreqs);
MaxFreq = max(StimFreqs);

%Eerst ratecurve weergeven, dit voor beide monaurale gegevens, maar ook voor de gesimuleerde
%coincidence-detector ...
Ax_Mon = axes('Position', [0.05 0.55 0.30 0.40]); set(Ax_Mon, 'XLim', [MinFreq MaxFreq]);
Hdl_LineMon1 = line(MRA(1).X, MRA(1).Y, 'Color', 'b', 'LineStyle', '-', 'Marker', '^');
Hdl_LineMon2 = line(MRA(2).X, MRA(2).Y, 'Color', 'b', 'LineStyle', '-', 'Marker', 'v');

xlabel('Stimulus Frequency(Hz)'); ylabel('Monaural Rate(spk/sec)');

line([MinFreq MaxFreq], [Param.CellParam(1).SA Param.CellParam(1).SA], 'Color', 'k', 'LineStyle', ':', 'Marker', 'none');
line(MinFreq, Param.CellParam(1).SA, 'Color', 'b', 'LineStyle', 'none', 'Marker', '^');
line([MinFreq MaxFreq], [Param.CellParam(2).SA Param.CellParam(2).SA], 'Color', 'k', 'LineStyle', ':', 'Marker', 'none');
line(MinFreq, Param.CellParam(2).SA, 'Color', 'b', 'LineStyle', 'none', 'Marker', 'v');

Ax_Bin = axes('Position', get(Ax_Mon, 'Position')); set(Ax_Bin, 'XLim', [MinFreq MaxFreq], 'Color', 'none', 'YAxisLocation', 'right');
Hdl_LineBin = line(IRA.X, IRA.Y, 'Color', 'g', 'LineStyle', '-', 'Marker', 'o');

ylabel('Interaural Rate(spk/sec)');

%Daarna phasecurven weergeven ...
subplot('Position', [0.05 0.05 0.30 0.40]); hold on;
idx1 = find(MPC(1).pRayleigh <= 0.001);
plot(MPC(1).X(idx1), MPC(1).Y(idx1), 'b^-');
idx2 = find(MPC(2).pRayleigh <= 0.001);
plot(MPC(2).X(idx2), MPC(2).Y(idx2), 'bv-');

plot(PD.X(PD.iRayleigh), PD.Y(PD.iRayleigh), 'rx-');
idx = find(IPC.pRayleigh <= 0.001);
plot(IPC.X(idx), IPC.Y(idx), 'go-');

xlabel('Stimulus Frequency(Hz)'); ylabel('Phase(cycles)');

MinX = 0; MaxX = max([MPC(1).X(end) MPC(2).X(end)]);
line([MinX, MaxX], [MPC(1).YInterSect, polyval([MPC(1).Slope/1000 MPC(1).YInterSect], MaxX)], 'Color', 'k', 'LineStyle', ':', 'Marker', 'none');
line([MinX, MaxX], [MPC(2).YInterSect, polyval([MPC(2).Slope/1000 MPC(2).YInterSect], MaxX)], 'Color', 'k', 'LineStyle', ':', 'Marker', 'none');
line([MinX, MaxX], [IPC.CP, polyval([IPC.CD/1000 IPC.CP], MaxX)], 'Color', 'k', 'LineStyle', ':', 'Marker', 'none');
line([MinX, MaxX], [PD.YInterSect, polyval([PD.Slope/1000 PD.YInterSect], MaxX)], 'Color', 'k', 'LineStyle', ':', 'Marker', 'none');
xlim([MinX MaxX]);

legend({[DataID1]; ...
        [DataID2]; ...
        ['\Delta'];...
        ['CC']},2);

%R curven plotten ...
Ax_IR = axes('Position', [0.45 0.55 0.30 0.40]); 
idx = find(MR(1).pRayleigh <= 0.001);
line(MR(1).X(idx), MR(1).Y(idx), 'LineStyle', '-', 'Color', 'b', 'Marker', '^');
idx = find(MR(2).pRayleigh <= 0.001);
line(MR(2).X(idx), MR(2).Y(idx), 'LineStyle', '-', 'Color', 'b', 'Marker', 'v');
idx = find(IR.pRayleigh <= 0.001);
line(IR.X(idx), IR.Y(idx), 'LineStyle', '-', 'Color', 'g', 'Marker', 'o');
xlabel('Stimulus Frequency(Hz)'); ylabel('R');

MinX = min([MR(1).X(1) MR(2).X(1) IR.X(1)]); MaxX = max([MR(1).X(end) MR(2).X(end) IR.X(end)]);
MinY = 0; MaxY = 1;
axis([MinX MaxX MinY MaxY]);

%SyncRate curven plotten ...
Ax_MSR = axes('Position', [0.45 0.05 0.30 0.40]); 
idx = find(MSR(1).pRayleigh <= 0.001);
Hdl_MSR1 = line(MSR(1).X(idx), MSR(1).Y(idx), 'LineStyle', '-', 'Color', 'b', 'Marker', '^');
idx = find(MSR(2).pRayleigh <= 0.001);
Hdl_MSR2 = line(MSR(2).X(idx), MSR(2).Y(idx), 'LineStyle', '-', 'Color', 'b', 'Marker', 'v');

xlabel('Stimulus Frequency(Hz)'); ylabel('Monaurale SyncRate(R*Rate)');

MinX = min([MSR(1).X MSR(2).X]); MaxX = max([MSR(1).X MSR(2).X]);
xlim([MinX MaxX]);

Ax_ISR = axes('Position', get(Ax_MSR, 'Position')); set(Ax_ISR, 'XLim', [MinX MaxX], 'Color', 'none', 'YAxisLocation', 'right');
idx = find(ISR.pRayleigh <= 0.001);
Hdl_ISR = line(ISR.X(idx), ISR.Y(idx), 'LineStyle', '-', 'Color', 'g', 'Marker', 'o');

ylabel('Interaurale SyncRate(R*Rate)');

%Extra gegevens weergeven ...
printinfo([0.80 0.55 0.15 0.40], {sprintf('%s', DataID); ...
         ''; ...
         ['<' CellInfo.Cell1.dsID '>']; ...
        sprintf('CF = %.0fHz', Param.CellParam(1).CF); ...
        sprintf('SR = %.0fspk/sec', Param.CellParam(1).SA); ...
        sprintf('Comment: %s',comment1); ...
        ''; ...
        ['<' CellInfo.Cell2.dsID '>']; ...
        sprintf('CF = %.0fHz', Param.CellParam(2).CF); ...
        sprintf('SR = %.0fspk/sec', Param.CellParam(2).SA); ...
        sprintf('Comment: %s',comment2); ...
        ''; ...
        ['\DeltaCF = ' sprintf('%.3f', log2(Param.CellParam(2).CF/Param.CellParam(1).CF)) 'oct']; ...
        ['\DeltaCoDist = ' sprintf('%.2f', (Param.CellParam(2).CD - Param.CellParam(1).CD)) 'mm']});
printinfo([0.80 0.05 0.15 0.40], {'Rate:'; ...
        sprintf('ANWIN = [%dms-%dms]', Param.CalcParam.coranwin(1), Param.CalcParam.coranwin(2)); ...
        ''; ...
        'CycleHistogram:'; ...
        sprintf('ANWIN = [%dms-%dms]', Param.CalcParam.coranwin(1), Param.CalcParam.coranwin(2)); ...
        sprintf('Nr. of Bins = %d', Param.CalcParam.cyclenbin); ...
        ''; ...
        'IPC-\DeltaMonPhase Corr:'; ...
        sprintf('CorrCoeff : %.3f', PD.CorrCoef); ...
        sprintf('pCorrCoef : %.3f', PD.pCorrCoef)});

%---------------%
% PLOTMAINCURVE %
%---------------%
function Hdl = PlotMainCurve(PlotData,PlotInfo)

Hdl = figure('Name', ['EvalFS: ' PlotInfo.Title ], ...
    'NumberTitle', 'off', ...
    'PaperType', 'a4letter', ...
    'PaperUnits', 'normalized', ...
    'PaperPosition', [0.05 0.05 0.90 0.90], ...
    'PaperOrientation', 'landscape');

subplot('Position', [0.25 0.59 0.5 0.37]); plot(PlotData.X, PlotData.Y, PlotInfo.Style);
txt_title = title(PlotInfo.Title); set(txt_title, 'FontSize', 12);
xlabel(PlotInfo.XLabel); ylabel(PlotInfo.YLabel);
axis([ PlotInfo.XLim, PlotInfo.YLim ]);
hold on;

txt_info = text(PlotInfo.XLim(1), PlotInfo.YLim(2), PlotInfo.Text);
set(txt_info, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');


%-----------------%
% PLOTSMALLCURVES %
%-----------------%
function PlotSmallCurves(PlotData, PlotInfo)

XStartSpace = 0.055; YStartSpace = 0.055;
XSpacing = 0.015; YSpacing = 0.020; Width = 0.50; Length = 1;

N       = length(PlotData);   %Aantal correlogrammen te plotten
[nL,nW] = fact2(N);       %Beste organisatie van plots nagaan 

PlotWidth = (Width - (nW*YSpacing) - YStartSpace)/nW; 
PlotLength = (Length - (nL*XSpacing) - XStartSpace)/nL;

for n = 1:N
    Nrow = nW - floor((n-1)/nL); 
    Ncol = rem(n,nL); if Ncol == 0, Ncol = nL; end
    X = (Ncol-1)*PlotLength + (Ncol-1) * XSpacing + XStartSpace;
    Y = (Nrow-1)*PlotWidth + (Nrow-1) * YSpacing + YStartSpace;
    Pos  = [X Y PlotLength PlotWidth];
    
    axhdl(n) = axes('Parent', gcf, 'Position', Pos); 
    switch PlotInfo.Type
    case 'PLOT', plot(PlotData(n).X, PlotData(n).Y, 'b-');
    case 'BAR' 
        barhdl = bar(PlotData(n).X, PlotData(n).Y, 1);
        set(barhdl, 'EdgeColor', [0.7 0.7 0.7], 'FaceColor', [0.7 0.7 0.7]);
    end    
    set(axhdl(n), 'Box', 'on', 'XtickLabel', [], 'YTickLabel', []);
    
    MaxY(n) = max(PlotData(n).Y); if isnan(MaxY(n)), MaxY(n) = 1; end
    if MaxY(n) == 0
        MaxY(n) = 1;
    end
    axis([PlotInfo.XLim, 0, MaxY(n)]);
    
    if Nrow == 1
        set(axhdl(n), 'XTickLabel', [PlotInfo.XLim(1) mean(PlotInfo.XLim) PlotInfo.XLim(2)]);
        xlabel(PlotInfo.XLabel);
    end
    
    txt_info(n) = text(PlotInfo.XLim(1), MaxY(n), eval(PlotInfo.Text));
    set(txt_info(n), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 8);
end
MaxY = max(MaxY);
for n = 1:N
    Nrow = nW - floor((n-1)/nL); 
    Ncol = rem(n,nL); if Ncol == 0, Ncol = nL; end
    
    set(axhdl(n),'YLim', [0, MaxY] );
    set(txt_info(n), 'Position', [PlotInfo.XLim(1), MaxY]);
    
    if Ncol == 1
        YTickNum = (0:2)' *(MaxY/2);
        YTickLabel = []; for i = 1:length(YTickNum), YTickLabel = strvcat(YTickLabel, sprintf('%.1f', YTickNum(i))); end
        set(axhdl(n), 'YTick', YTickNum, 'YTickLabel', YTickLabel);
        axes(axhdl(n));
    end
    if (Ncol == 1) & (Nrow == floor(nW/2)), ylabel(PlotInfo.YLabel); end
end


function fileName = getFileName(ds)
fileName = ds.fileName;
if (length(fileName) >= 6)
	fileName = fileName(1:6);
end



