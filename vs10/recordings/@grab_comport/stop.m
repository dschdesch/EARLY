function stop(G, MessMode);
% grab_comport/stop - stop event grabbing
%    stop(G) stops data grabbing by grab_comport object G.
% 
%    Two steps:
%       1. clock the DA progress and store it in G
%       2. call generic action/stop
%
%    See also sortConditions, GUImessage, stimpresent/DAstatus.

eval(IamAt('indent')); % indent action method below

G = download(G);
if isequal('finished', status(G)), % just too late; cancel the stop action
    return;
end
SPS = seqplaystatus;
G.datagrabber = stopstatus(G.datagrabber, SPS.isam_abs);
upload(G);
stop(G.action); % generic action/stop - includes uploading



