function Freq=EvalFrequencyHARHAR(figh, Prefix, P);
% EvalFrequencyStepper - compute frequency series from Frequency stepper GUI
%   Freq=EvalFrequencyStepper(figh) reads frequency-stepper specs
%   from paramqueries StartFreq, StepFrequency, EndFrequency, AdjustFreq,
%   in the GUI figure with handle figh (see FrequencyStepper), and converts
%   them to the individual frequencies of the frequency-stepping series.
%   Any errors in the user-specified values results in an empty return 
%   value Freq, while an error message is displayed by GUImessage.
%
%   EvalFrequencyStepper(figh, 'Foo') uses prefix Foo for the query names,
%   i.e., FooStartFreq, etc. The prefix defaults to ''.
%
%   EvalFrequencyStepper(figh, Prefix, P) does not read the queries, but
%   extracts them from struct P which was previously returned by GUIval.
%   This is the preferred use of EvalFrequencyStepper, because it leaves
%   the task of reading the parameters to the generic GUIval function. The
%   first input arg figh is still needed for error reporting.
%
%   See StimGUI, FrequencyStepper, GUIval, GUImessage.

if nargin<2, Prefix=''; end
if nargin<3, P = []; end
Freq = []; % allow premature return

if isempty(P), % obtain info from GUI. Non-preferred method; see help text.
    EXP = getGUIdata(figh,'Experiment');
    Q = getGUIdata(figh,'Query');
    StartNHN = read(Q([Prefix 'StartNHN']));
    [StepFreq, StepFreqUnit] = read(Q([Prefix 'StepNHN']));
    EndFreq = read(Q([Prefix 'EndNHN']));
    AdjustFreq = read(Q([Prefix 'AdjustFreq']));
    P = CollectInStruct(StartFreq, StepFreq, StepFreqUnit, EndFreq, AdjustFreq);
else,
    P = dePrefix(P, Prefix);
    EXP = P.Experiment;
end

% check if stimulus frequencies are within range [CT.minFreq CT.maxFreq]
somethingwrong=1;
if any((P.CF/P.StartNHN)<EXP.minStimFreq),
    GUImessage(figh, {'Start frequency violates min stim frequency'...
        ['of ' num2str(EXP.minStimFreq) ' Hz']},'error', [Prefix 'StartNHN']);
elseif any((P.CF/P.EndNHN)<EXP.minStimFreq),
    GUImessage(figh, {'End frequency violates min stim frequency'...
        ['of ' num2str(EXP.minStimFreq) ' Hz']},'error', [Prefix 'EndNHN']);
elseif any((P.CF/P.StartNHN)>EXP.maxStimFreq),
    GUImessage(figh, {'Start frequency exceeds max stim frequency'...
        ['of ' num2str(EXP.maxStimFreq) ' Hz']},'error', [Prefix 'StartNHN']);
elseif any((P.CF/P.EndNHN)>EXP.maxStimFreq),
    GUImessage(figh, {'End frequency exceeds max stim frequency'...
        ['of ' num2str(EXP.maxStimFreq) ' Hz']},'error', [Prefix 'EndNHN']);
else, % passed all the tests..
    somethingwrong=0;
end
if somethingwrong, return; end

% delegate the computation to generic EvalStepper
StepMode = 'Linear';
P.StepFreqUnit = 'Hz';

[Freq, Mess]=EvalStepper((P.CF/P.StartNHN), (P.CF/P.StepNHN), (P.CF/P.EndNHN), StepMode, ...
    P.AdjustFreq, [EXP.minStimFreq EXP.maxStimFreq], EXP.maxNcond);
if isequal('nofit', Mess),
    Mess = {'Stepsize does not exactly fit Frequency bounds', ...
        'Adjust Frequency parameters or toggle Adjust button.'};
elseif isequal('largestep', Mess)
    Mess = 'Step size exceeds frequency range';
elseif isequal('cripple', Mess)
    Mess = 'Different # frequency steps in the two DA channels';
elseif isequal('toomany', Mess)
    Mess = {['Too many (>' num2str(EXP.maxNcond) ') frequency steps.'] 'Increase stepsize or decrease range'};
end
GUImessage(figh,Mess, 'error',{[Prefix 'StartFreq'] [Prefix 'StepFreq'] [Prefix 'EndFreq'] });




