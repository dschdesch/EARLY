function alpha = RaySign(r,n,roundflag)
% RaySign - significance of circular distr. according to Raleigh's z
%
% alpha = RaySign(r,n) returns the min. level of significance z=n*r^2.
%       These values are interpolated from a table (see reference below).
%       r are (complex) vector strengths, n is the number observations used
%       to calculate r.
%
% r may be an array, in which case alpha is calculated for each r(j) at the
%  provided n. If both r and n are arrays (then, same size required!),
%  alpha(j) is for r(j) and n(j).
%
% This function is an alternative implementation from the function
% RayleighSign.m, somewhat faster. It uses TriScatteredInterp.m to
% calculate alpha for a given r and n from the table (see ref. below).
% With this, alpha is continuous between 0.0001 and 0.5. To get at the
% same alpha's as from RayleighSign.m, these (interpolated) alpha's are
% disretized. To avoid this discretization (and make the function even
% faster), use syntax:
%           alpha = RaySign(r,n,false)
%
% T = RaySign('table') returns the used table & TriScatteredInterp-function
% in struct T.
% 
% --- Reference ---:
%   Jerrold Zar, Biostatistical Analysis, 
%   table 571, who refers to Greenwood and Durand (1955), Ann Math Stat 26
%
% See also RayleighSign, TriScatteredInterp

if nargin<3 | isempty(roundflag), roundflag = true; end

%store the table + the TriScatteredInterp function in a persistent
persistent SigniTable
if isempty(SigniTable),
   SigniTable = local_SigniTable;
end

if isequal(r,'table'), % just return the table
   alpha = SigniTable;
   return;
end

%check dimensions of n and r
if any(size(n)>1) & ~isequal(size(n), size(r)),
    error('The dimensions of r and n are not compatible.');
end
[r,n]=SameSize(r,n);

r = abs(r); % discard phase in r

%restrict a lot of observations to within table
n(n>max(SigniTable.n)) = max(SigniTable.n);

z = n.*r.^2;
alpha = nan(size(r)); %init of alpha

inumb = n<min(SigniTable.n); %too few observations --> alpha = max(SigniTable.alpha)
alpha(inumb) = max(SigniTable.alpha);

alpha(~inumb) = SigniTable.TSI(n(~inumb),z(~inumb)); %calculates alpha's

%correct for z-values outside the table
indx = isnan(alpha);
if any(indx),
    imax = z>=min(SigniTable.z_a_n(:,end));
    alpha(indx & imax) = min(SigniTable.alpha);

    imin = z<=max(SigniTable.z_a_n(:,1));
    alpha(indx & imin) = max(SigniTable.alpha);
end

if roundflag==1,
    %round to next-above alpha, the for-loop is faster than arrayfun
    for ii = 1:numel(alpha),
        indx = find(alpha(ii)-SigniTable.alpha<=0,1,'last');
        alpha(ii) = SigniTable.alpha(indx);
    end
end

%---------------local--------------------
function s = local_SigniTable
% table adapted from TAP. 
% See Jerrold Zar, Biostatistical Analysis, 
% table 571, who referes to Greenwood and Durand (1955), Ann Math Stat 26
%
%NOTE: table is different from RayleighSign.m in that the n=inf is missing
%       Inf is not allowed in TriScatteredInterp.m

alpha = [0.5	0.2	0.1	0.05	0.02	0.01	0.005	0.002	0.001];
n = [6:30 32:2:50 55:5:75 80 90 100:20:200 300 500 6e6];

z_a_n = [0.734	1.639	2.274	2.865	3.576	4.058	4.491	4.985	5.297
0.727	1.634	2.278	2.885	3.627	4.143	4.617	5.181	5.556
0.723	1.631	2.281	2.899	3.665	4.205	4.71	5.322	5.743
0.719	1.628	2.283	2.91	3.694	4.252	4.78	5.43	5.885
0.717	1.626	2.285	2.919	3.716	4.229	4.835	5.514	5.996
0.715	1.625	2.287	2.926	3.735	4.319	4.879	5.582	6.085
0.713	1.623	2.288	2.932	3.75	4.344	4.916	5.638	6.158
0.711	1.622	2.289	2.937	3.763	4.365	4.947	5.685	6.219
0.71	1.621	2.29	2.941	3.774	4.383	4.973	5.725	6.271
0.709	1.62	2.291	2.945	3.784	4.398	4.996	5.759	6.316
0.708	1.62	2.292	2.948	3.792	4.412	5.015	5.789	6.354
0.707	1.619	2.292	2.951	3.799	4.423	5.033	5.815	6.388
0.706	1.619	2.293	2.954	3.806	4.434	5.048	5.838	6.418
0.705	1.618	2.293	2.956	3.811	4.443	5.061	5.858	6.445
0.705	1.618	2.294	2.958	3.816	4.451	5.074	5.877	6.469
0.704	1.617	2.294	2.96	3.821	4.459	5.085	5.893	6.491
0.704	1.617	2.295	2.961	3.825	4.466	5.095	5.908	6.51
0.703	1.616	2.295	2.963	3.829	4.472	5.104	5.922	6.528
0.703	1.616	2.295	2.964	3.833	4.478	5.112	5.935	6.544
0.702	1.616	2.296	2.966	3.836	4.483	5.12	5.946	6.559
0.702	1.616	2.296	2.967	3.839	4.488	5.127	5.957	6.573
0.702	1.615	2.296	2.968	3.842	4.492	5.133	5.966	6.586
0.701	1.615	2.296	2.969	3.844	4.496	5.139	5.975	6.598
0.701	1.615	2.297	2.97	3.847	4.5	5.145	5.984	6.609
0.701	1.615	2.297	2.971	3.849	4.504	5.15	5.992	6.619
0.7	1.614	2.297	2.972	3.853	4.51	5.159	6.006	6.637
0.7	1.614	2.297	2.974	3.856	4.516	5.168	6.018	6.654
0.7	1.614	2.298	2.975	3.859	4.521	5.175	6.03	6.668
0.699	1.614	2.298	2.976	3.862	4.525	5.182	6.039	6.681
0.699	1.613	2.298	2.977	3.865	4.529	5.188	6.048	6.692
0.699	1.613	2.298	2.978	3.867	4.533	5.193	6.056	6.703
0.698	1.613	2.299	2.979	3.869	4.536	5.198	6.064	6.712
0.698	1.613	2.299	2.979	3.871	4.539	5.202	6.07	6.721
0.698	1.613	2.299	2.98	3.873	4.542	5.206	6.076	6.729
0.698	1.613	2.299	2.981	3.874	4.545	5.21	6.082	6.736
0.697	1.612	2.299	2.982	3.878	4.55	5.218	6.094	6.752
0.697	1.612	2.3	2.983	3.881	4.555	5.225	6.104	6.765
0.697	1.612	2.3	2.984	3.883	4.559	5.231	6.113	6.776
0.696	1.612	2.3	2.985	3.885	4.562	5.235	6.12	6.786
0.696	1.612	2.3	2.986	3.887	4.565	5.24	6.127	6.794
0.696	1.611	2.3	2.986	3.889	4.567	5.243	6.132	6.801
0.696	1.611	2.301	2.987	3.891	4.572	5.249	6.141	6.813
0.695	1.611	2.301	2.988	3.893	4.575	5.254	6.149	6.822
0.695	1.611	2.301	2.99	3.896	4.58	5.262	6.16	6.837
0.695	1.611	2.301	2.99	3.899	4.584	5.267	6.168	6.847
0.695	1.61	2.301	2.991	3.9	4.586	5.271	6.174	6.855
0.694	1.61	2.302	2.992	3.902	4.588	5.274	6.178	6.861
0.694	1.61	2.302	2.992	3.903	4.59	5.276	6.182	6.865
0.694	1.61	2.302	2.993	3.906	4.595	5.284	6.193	6.879
0.694	1.61	2.302	2.994	3.908	4.599	5.29	6.201	6.891
0.6931	1.6094	2.3026	2.9957	3.912	4.6052	5.2983	6.2146	6.9078];

A = sort(repmat(alpha(:),numel(n),1),'descend');
N = repmat(n(:), numel(alpha),1);
Z = reshape(z_a_n,numel(n)*numel(alpha),1);
TSI = TriScatteredInterp([N Z] ,A,'natural'); %return alpha, given n and z

s = CollectInStruct(alpha, n, z_a_n, TSI);










